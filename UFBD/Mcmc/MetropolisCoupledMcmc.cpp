#include "MetropolisCoupledMcmc.hpp"
#include "Msg.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Serialize.hpp"
#include "ThreadPool.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"
#include "WriteTSV.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <chrono>
#include <future>
#include <iomanip>
#include <iostream>
#include <sstream>

MetropolisCoupledMcmc::MetropolisCoupledMcmc(unsigned long ng, int thin, std::vector<PhylogeneticModel*> m, unsigned int masterSeed) : numCycles(ng), thinning(thin),
    models(m),
    numModels(models.size()),
    coldModelIdx(-1),
    deltaT(0.2),
    threadPool(std::max(1, std::min((int)m.size(), UserSettings::userSettings().getNumCores()))){
    swapRng.setSeed(masterSeed + numModels);
    currLnL.reserve(numModels); //needs to be reserve for pushback
    newLnL.resize(numModels);
    currLnP.reserve(numModels); //needs to be reserve for pushback
    newLnP.resize(numModels);
    indices.reserve(numModels); //needs to be reserve for pushback
    lnProposalRatio.resize(numModels);
    lnLikelihoodRatio.resize(numModels);
    lnPriorRatio.resize(numModels);
    lnAcceptanceProbabilities.resize(numModels);
    
    UserSettings& settings = UserSettings::userSettings();
    treeOut = settings.getTreeOutput();
    paramOut = settings.getParamOutput();
    writeTrees = (treeOut.empty() == false && models[0]->getTree()->getNumBackbone() > 0);
    gen = 0;

    int N = settings.getNumCores();
    for(int i = 0; i < numModels; i++){
        int per = N / numModels + (i < N % numModels ? 1 : 0);
        if(per < 1) per = 1;
        chainPools.push_back(new ThreadPool(per));
    }
    ThreadPool::shared().setChainCap(std::max(1, N / std::max(1, numModels)));
    chainCalibGen = 0;
    chainSeqT = 0.0;
    chainParT = 0.0;
    chainDecision = -1;
}

MetropolisCoupledMcmc::~MetropolisCoupledMcmc(void){
    for(ThreadPool* p : chainPools)
        delete p;
}

void MetropolisCoupledMcmc::printMoveDiagnostics(int rep){
    for(size_t i = 0; i < models.size(); i++){
        std::cout << "  run " << rep << " chain " << i
                  << ((int)i == coldModelIdx ? " [cold]  " : " [heated]  ");
        models[i]->print();
    }
}

Tree* MetropolisCoupledMcmc::getTree(void){
    int i = (coldModelIdx >= 0) ? coldModelIdx : 0;
    return models[i]->getTree();
}

bool MetropolisCoupledMcmc::treeHasFossils(void){
    int i = (coldModelIdx >= 0) ? coldModelIdx : 0;
    return models[i]->treeIncludesFossils();
}

double MetropolisCoupledMcmc::calcHeating(int idx){
    if(idx == 0)
        return 1.0;
    return (1 / (1 + deltaT * idx));
}

void MetropolisCoupledMcmc::init(void) {
    if(numCycles >= std::numeric_limits<unsigned long>::max())
        Msg::error("number of cycles requested is greater than the largest possible value");
    RandomVariable::setActiveInstance(&swapRng);
    int idx = 0;
    for(PhylogeneticModel* m : models){
        currLnL.push_back(m->lnLikelihood());
        currLnP.push_back(m->lnPriorProbability());
        indices.push_back(idx);
        idx++;
    }
    gen = 0;

    UserSettings& settings = UserSettings::userSettings();
    if(settings.getSigma2Param() == Sigma2Param::PNCP && settings.getPncpTuningGens() > 0){
        for(PhylogeneticModel* m : models)
            m->setChainLabel(runLabel);
        tuning = true;
        advance(settings.getPncpTuningGens());
        for(PhylogeneticModel* m : models)
            m->freezePncpTuning();
        tuning = false;
        gen = 0;
    }
}

void MetropolisCoupledMcmc::finalize(void) {
    params.closeTSV();
    if(writeTrees)
        trees.appendDataTSV("END;");
    trees.closeTSV();
    if(writeLatent)
        latent.closeTSV();
}

void MetropolisCoupledMcmc::advance(unsigned long nGens) {
    RandomVariable::setActiveInstance(&swapRng);
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    unsigned long target = gen + nGens;
    while (gen < target){
        gen++;
        unsigned long n = gen;
        if(n < 10000 && n % 50 == 0)
            updateDeltaT();
        
        bool canParallel = (numModels > 1 && threadPool.size() > 1);
        bool parallelChains = canParallel && (chainDecision < 0 ? chainCalibGen >= 100 : chainDecision == 1);
        std::chrono::steady_clock::time_point cwT0 = std::chrono::steady_clock::now();
        if(parallelChains){
            std::vector<std::function<void()>> tasks;
            tasks.reserve(numModels);
            for(int i = 0; i < numModels; i++){
                tasks.push_back([this, i](){
                    RandomVariable::setActiveInstance(models[i]->getRng());
                    ThreadPool::setCurrent(chainPools[i]);
                    lnProposalRatio[i] = models[i]->update();
                    if(lnProposalRatio[i] == -INFINITY){
                        newLnL[i] = currLnL[i];
                        newLnP[i] = currLnP[i];
                        lnAcceptanceProbabilities[i] = -INFINITY;
                    }else{
                        newLnL[i] = models[i]->lnLikelihood();
                        newLnP[i] = models[i]->lnPriorProbability();
                        lnLikelihoodRatio[i] = newLnL[i] - currLnL[i];
                        lnPriorRatio[i]      = newLnP[i] - currLnP[i];
                        int    chainIdx = indices[i];
                        double heat     = (chainIdx == 0) ? 1.0 : calcHeating(chainIdx);
                        lnAcceptanceProbabilities[i] = heat * (lnLikelihoodRatio[i] + lnPriorRatio[i])
                                                       + lnProposalRatio[i];
                    }
                    ThreadPool::setCurrent(nullptr);
                });
            }
            threadPool.parallelTasks(tasks);
        }else{
            for(int i = 0; i < numModels; i++){
                RandomVariable::setActiveInstance(models[i]->getRng());
                lnProposalRatio[i] = models[i]->update();

                if(lnProposalRatio[i] == -INFINITY){
                    newLnL[i] = currLnL[i];
                    newLnP[i] = currLnP[i];
                    lnAcceptanceProbabilities[i] = -INFINITY;
                    continue;
                }

                newLnL[i]          = models[i]->lnLikelihood();
                newLnP[i]          = models[i]->lnPriorProbability();

                lnLikelihoodRatio[i] = newLnL[i] - currLnL[i];
                lnPriorRatio[i]      = newLnP[i] - currLnP[i];

                int    chainIdx = indices[i];
                double heat     = (chainIdx == 0) ? 1.0 : calcHeating(chainIdx);
                lnAcceptanceProbabilities[i] = heat * (lnLikelihoodRatio[i] + lnPriorRatio[i])
                                               + lnProposalRatio[i];
            }
        }
        if(canParallel && chainDecision < 0){
            double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - cwT0).count();
            if(parallelChains) chainParT += dt; else chainSeqT += dt;
            if(++chainCalibGen >= 200)
                chainDecision = (chainParT < chainSeqT) ? 1 : 0;
        }


        // accept or reject the proposed state
        coldModelIdx = -1;
        for (int i = 0; i < numModels; i++) {
            if (indices[i] == 0) coldModelIdx = i;
            if (log(rng.uniformRv()) < lnAcceptanceProbabilities[i]) {
                currLnL[i] = newLnL[i];
                currLnP[i] = newLnP[i];
                models[i]->updateForAcceptance();
            } else {
                models[i]->updateForRejection();
            }
        }
        
        if (verbose && n % thinning == 0){
            const size_t numAccepted = std::count(recentAcceptRej.begin(), recentAcceptRej.end(), true);
            const double acceptanceRate = recentAcceptRej.empty() ? 0.0 : static_cast<double>(numAccepted) / recentAcceptRej.size();
            std::ostringstream os;
            os << std::fixed << std::setprecision(2) << "chain " << runLabel << "  " << n << " -- posterior "
               << (currLnL[coldModelIdx] + currLnP[coldModelIdx]) << " likelihood " << currLnL[coldModelIdx]
               << " | swap accept " << acceptanceRate << "\n";
            ChainRunner::logLine(os.str());
        }
            
        //choose two chains and swap
        int chain1 = (int)(rng.uniformRv() * numModels);
        int chain2 = -1;
        do{
            chain2 = (int)(rng.uniformRv() * numModels);
        }while(chain1 == chain2);
    
        int idx1 = indices[chain1];
        int idx2 = indices[chain2];
        
        double chain1LnPost = currLnL[chain1] + currLnP[chain1];
        double chain2LnPost = currLnL[chain2] + currLnP[chain2];
        
        double chain1Heating = calcHeating(idx1);
        double chain2Heating = calcHeating(idx2);
        
        double lnProbSwap = chain2Heating * chain1LnPost + chain1Heating * chain2LnPost;
        double lnProbStay = chain1Heating * chain1LnPost + chain2Heating * chain2LnPost;
        double lnAcceptanceSwap =lnProbSwap - lnProbStay;
        
        //accept or reject swap
        if(log(rng.uniformRv()) < lnAcceptanceSwap){
            indices[chain1] = idx2;
            indices[chain2] = idx1;
            recentAcceptRej.push_back(true);
        }else{
            recentAcceptRej.push_back(false);
        }
        
        if (recentAcceptRej.size() > 10000)
            recentAcceptRej.pop_front();
        
        if (tuning == false && (n == 1 || n == numCycles || n % thinning == 0))
            sample(n);
    }
}

void MetropolisCoupledMcmc::sample(unsigned long n) {
    if(n == 1){
        params.addFilepath(paramOut, true);
        std::vector<std::string> cn = {"n", "posterior", "likelihood", "prior"};
        std::vector<std::string> headStr = models[coldModelIdx]->getParameterNames();
        cn.insert( cn.end(), headStr.begin(), headStr.end() );
        params.addColumnNamesTSV(cn);

        std::vector<std::string> latNames = models[coldModelIdx]->getLatentNames();
        if(latNames.empty() == false){
            latentOut = paramOut;
            size_t dp = latentOut.rfind(".log");
            latentOut.insert(dp != std::string::npos ? dp : latentOut.size(), "_latent");
            latent.addFilepath(latentOut, true);
            std::vector<std::string> lcn = {"n"};
            lcn.insert(lcn.end(), latNames.begin(), latNames.end());
            latent.addColumnNamesTSV(lcn);
            writeLatent = true;
            latentNms = latNames;
            latentCols.assign(latNames.size(), std::vector<double>());
        }

        if(writeTrees){
            trees.addFilepath(treeOut, true);
            trees.appendDataTSV("#NEXUS");
            trees.appendDataTSV("BEGIN TREES;");
        }

        traceNms.clear();
        traceNms.push_back("posterior");
        traceNms.push_back("likelihood");
        traceNms.push_back("prior");
        for(const std::string& s : headStr)
            traceNms.push_back(s);
        traceCols.assign(traceNms.size(), std::vector<double>());
    }

    double cl = currLnL[coldModelIdx];
    double cp = currLnP[coldModelIdx];
    std::vector<double> dat = {(double)n, cl + cp, cl, cp};
    std::vector<double> parmStr = models[coldModelIdx]->getParameterString();
    dat.insert( dat.end(), parmStr.begin(), parmStr.end() );
    params.appendDataTSV(dat);

    if(writeLatent){
        std::vector<double> ls = models[coldModelIdx]->getLatentString();
        std::vector<double> ldat = {(double)n};
        ldat.insert(ldat.end(), ls.begin(), ls.end());
        latent.appendDataTSV(ldat);
        for(size_t j = 0; j < latentCols.size() && j < ls.size(); j++)
            latentCols[j].push_back(ls[j]);
    }

    if(writeTrees)
        trees.appendDataTSV("\tTREE STATE_" + std::to_string(n) + " = [&R] " + models[coldModelIdx]->getTree()->getBackboneNewickString(models[coldModelIdx]->treeIncludesFossils()));

    std::vector<double> tv = {cl + cp, cl, cp};
    tv.insert(tv.end(), parmStr.begin(), parmStr.end());
    for(size_t j = 0; j < traceCols.size() && j < tv.size(); j++)
        traceCols[j].push_back(tv[j]);

    writeCheckpoint();
}

void MetropolisCoupledMcmc::writeCheckpoint(void) {
    std::string path = paramOut + ".ckp";
    std::string tmp = path + ".tmp";
    std::ofstream os(tmp);
    os << std::setprecision(17);
    os << gen << ' ' << thinning << ' ' << deltaT << ' ' << numModels << '\n';
    for(int i = 0; i < numModels; i++)
        os << indices[i] << ' ';
    os << '\n';
    swapRng.writeState(os);
    Serialize::writeBoolDeque(os, recentAcceptRej);
    for(int i = 0; i < numModels; i++){
        os << currLnL[i] << ' ' << currLnP[i] << '\n';
        models[i]->getRng()->writeState(os);
        models[i]->writeState(os);
    }
    os.flush();
    os.close();
    std::rename(tmp.c_str(), path.c_str());
}

bool MetropolisCoupledMcmc::loadCheckpoint(void) {
    std::string path = paramOut + ".ckp";
    std::ifstream is = openCheckpoint(path);
    int nm, storedSf;
    is >> gen >> storedSf >> deltaT >> nm;
    reconcileThinning(storedSf, thinning);
    indices.assign(nm, 0);
    coldModelIdx = -1;
    for(int i = 0; i < nm; i++){
        is >> indices[i];
        if(indices[i] == 0)
            coldModelIdx = i;
    }
    swapRng.readState(is);
    Serialize::readBoolDeque(is, recentAcceptRej);
    currLnL.assign(nm, 0.0);
    currLnP.assign(nm, 0.0);
    for(int i = 0; i < nm; i++){
        is >> currLnL[i] >> currLnP[i];
        models[i]->getRng()->readState(is);
        models[i]->readState(is);
    }
    requireCheckpointIntact(is, path);
    return true;
}

RandomVariable* MetropolisCoupledMcmc::resumeRng(void) {
    return &swapRng;
}

std::vector<std::string> MetropolisCoupledMcmc::resumeParameterNames(void) {
    return models[coldModelIdx]->getParameterNames();
}

std::vector<std::string> MetropolisCoupledMcmc::resumeLatentNames(void) {
    return models[coldModelIdx]->getLatentNames();
}

void MetropolisCoupledMcmc::updateDeltaT(void) {
    if(recentAcceptRej.size() > 50){
        const size_t numAccepted = std::count(recentAcceptRej.begin(), recentAcceptRej.end(), true);
        const double acceptanceRate = static_cast<double>(numAccepted) / recentAcceptRej.size();
        
        constexpr double targetAcceptance = 0.23;
        constexpr double lowerAcceptance = targetAcceptance - 0.1;
        constexpr double upperAcceptance = targetAcceptance + 0.1;
        if(acceptanceRate < lowerAcceptance){
            deltaT *= 0.99;
        }else if (acceptanceRate > upperAcceptance){
            deltaT *= 1.01;
        }
    }
}
