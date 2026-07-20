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
#include <cmath>
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
    deltaT(UserSettings::userSettings().getDeltaTemperature()),
    resampleEvery(UserSettings::userSettings().getResampleEvery()),
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

    chainCalibGen = 0;
    chainSeqT = 0.0;
    chainParT = 0.0;
    chainDecision = -1;
    swapAdaptCount = 0;
    swapAdaptAcc = 0;
    swapAdaptAtt = 0;
    numSwapSweeps = 0;
    roundTrips = 0;
    lastEnd.assign(numModels, -1);
    emaAcc.assign(std::max(0, numModels - 1), 0.5);
    betas.assign(numModels, 1.0);
    for(int k = 1; k < numModels; k++)
        betas[k] = 1.0 / (1.0 + deltaT * k);
}

MetropolisCoupledMcmc::~MetropolisCoupledMcmc(void){
    for(ThreadPool* p : chainPools)
        delete p;
}

void MetropolisCoupledMcmc::printMoveDiagnostics(int rep){
    std::ios_base::fmtflags fmt = std::cout.flags();
    std::streamsize prec = std::cout.precision();
    std::cout << std::defaultfloat << std::setprecision(6);
    for(size_t i = 0; i < models.size(); i++){
        std::cout << "  run " << rep << " chain " << i
                  << ((int)i == coldModelIdx ? " [cold]  " : " [heated]  ");
        models[i]->print();
    }
    double Ehat = 0.0, Lhat = 0.0;
    for(size_t k = 0; k < emaAcc.size(); k++){
        double r = 1.0 - emaAcc[k];
        Lhat += r;
        Ehat += r / std::max(emaAcc[k], 1e-9);
    }
    std::cout << "  run " << rep << " PT: roundTrips " << roundTrips
              << " Lambda " << Lhat << " tauDEO " << (1.0 / (2.0 + 2.0 * Ehat))
              << " | beta";
    for(int k = 0; k < numModels; k++)
        std::cout << " " << betas[k];
    std::cout << " | adjAccept";
    for(size_t k = 0; k < emaAcc.size(); k++)
        std::cout << " " << emaAcc[k];
    std::cout << "\n";
    std::cout.flags(fmt);
    std::cout.precision(prec);
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
    return betas[idx];
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
    if(chainPools.empty()){
        UserSettings& s = UserSettings::userSettings();
        int nRuns = std::max(1, s.getNumRuns());
        int budget = s.getNumCores() / nRuns + ((runLabel < s.getNumCores() % nRuns) ? 1 : 0);
        for(int i = 0; i < numModels; i++){
            int per = budget / numModels + ((i < budget % numModels) ? 1 : 0);
            if(per < 1) per = 1;
            chainPools.push_back(new ThreadPool(per));
        }
        chainDecision = 1;
        for(int i = 0; i < numModels; i++)
            lastEnd[i] = (indices[i] == 0) ? 0 : (indices[i] == numModels - 1 ? 1 : -1);
        betas[0] = 1.0;
        for(int k = 1; k < numModels; k++)
            betas[k] = 1.0 / (1.0 + deltaT * k);
    }
    RandomVariable::setActiveInstance(&swapRng);
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    coldModelIdx = -1;
    for(int i = 0; i < numModels; i++)
        if(indices[i] == 0){ coldModelIdx = i; break; }

    unsigned long target = gen + nGens;
    if(gen == 0 && tuning == false)
        sample(0);

    bool canParallel = (numModels > 1 && threadPool.size() > 1);
    while(gen < target){
        unsigned long blockLen = std::min(resampleEvery, target - gen);
        unsigned long gen0 = gen;

        coldModelIdx = -1;
        for(int i = 0; i < numModels; i++)
            if(indices[i] == 0){ coldModelIdx = i; break; }

        auto runBlock = [this, blockLen, gen0](int i){
            RandomVariable::setActiveInstance(models[i]->getRng());
            ThreadPool::setCurrent(chainPools[i]);
            RandomVariable* r = models[i]->getRng();
            double heat = calcHeating(indices[i]);
            bool cold = (i == coldModelIdx);
            for(unsigned long b = 1; b <= blockLen; b++){
                double lnProp = models[i]->update();
                if(lnProp == -INFINITY){
                    models[i]->updateForRejection();
                }else{
                    double nL = models[i]->lnLikelihood();
                    double nP = models[i]->lnPriorProbability();
                    double lnAcc = heat * ((nL - currLnL[i]) + (nP - currLnP[i])) + lnProp;
                    if(log(r->uniformRv()) < lnAcc){
                        currLnL[i] = nL;
                        currLnP[i] = nP;
                        models[i]->updateForAcceptance();
                    }else{
                        models[i]->updateForRejection();
                    }
                }
                if(cold && tuning == false && (gen0 + b) % (unsigned long)thinning == 0)
                    writeColdSample(gen0 + b);
            }
            ThreadPool::setCurrent(nullptr);
        };

        if(canParallel){
            std::vector<std::function<void()>> tasks;
            tasks.reserve(numModels);
            for(int i = 0; i < numModels; i++)
                tasks.push_back([&runBlock, i](){ runBlock(i); });
            threadPool.parallelTasks(tasks);
        }else{
            for(int i = 0; i < numModels; i++)
                runBlock(i);
        }
        gen += blockLen;

        if(verbose && tuning == false && gen % (unsigned long)thinning == 0){
            const size_t numAccepted = std::count(recentAcceptRej.begin(), recentAcceptRej.end(), true);
            const double acceptanceRate = recentAcceptRej.empty() ? 0.0 : static_cast<double>(numAccepted) / recentAcceptRej.size();
            std::ostringstream os;
            os << std::fixed << std::setprecision(2) << "chain " << runLabel << "  " << gen << " -- posterior "
               << (currLnL[coldModelIdx] + currLnP[coldModelIdx]) << " likelihood " << currLnL[coldModelIdx]
               << " prior " << currLnP[coldModelIdx] << " | swap accept " << acceptanceRate << "\n";
            ChainRunner::logLine(os.str());
        }

        std::vector<int> slotOf(numModels);
        for(int s = 0; s < numModels; s++)
            slotOf[indices[s]] = s;
        int parity = (int)(numSwapSweeps & 1ul);
        for(int k = parity; k + 1 < numModels; k += 2){
            int a = slotOf[k], b = slotOf[k + 1];
            double lnAcc = (calcHeating(k) - calcHeating(k + 1)) *
                           ((currLnL[b] + currLnP[b]) - (currLnL[a] + currLnP[a]));
            bool acc = std::log(rng.uniformRv()) < lnAcc;
            if(acc){
                indices[a] = k + 1;
                indices[b] = k;
            }
            emaAcc[k] = 0.98 * emaAcc[k] + 0.02 * (acc ? 1.0 : 0.0);
            recentAcceptRej.push_back(acc);
            if(recentAcceptRej.size() > 10000)
                recentAcceptRej.pop_front();
            adaptLadder(acc);
        }
        numSwapSweeps++;
        int topRung = numModels - 1;
        for(int s = 0; s < numModels; s++){
            if(indices[s] == topRung && lastEnd[s] == 0)
                lastEnd[s] = 1;
            else if(indices[s] == 0){
                if(lastEnd[s] == 1)
                    roundTrips++;
                lastEnd[s] = 0;
            }
        }

        if(tuning == false)
            writeCheckpoint();
    }
}

void MetropolisCoupledMcmc::sample(unsigned long n) {
    if(n == 0){
        params.addFilepath(paramOut, true);
        std::vector<std::string> cn = {"n", "posterior", "likelihood", "prior"};
        std::vector<std::string> headStr = models[coldModelIdx]->getParameterNames();
        cn.insert( cn.end(), headStr.begin(), headStr.end() );
        params.addColumnNamesTSV(cn);

        std::vector<std::string> latNames = models[coldModelIdx]->getLatentNames();
        if(latNames.empty() == false && UserSettings::userSettings().getWriteLatentLog()){
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
        fixedMask.assign(3, false);
        for(bool fb : models[coldModelIdx]->getParameterFixedMask())
            fixedMask.push_back(fb);
    }

    writeColdSample(n);
    writeCheckpoint();
}

void MetropolisCoupledMcmc::writeColdSample(unsigned long n) {
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
}

void MetropolisCoupledMcmc::writeCheckpoint(void) {
    std::string path = paramOut + ".ckp";
    std::string tmp = path + ".tmp";
    std::ofstream os(tmp);
    os << std::setprecision(17);
    os << gen << ' ' << thinning << ' ' << deltaT << ' ' << numModels << ' ' << swapAdaptCount << ' ' << swapAdaptAcc << ' ' << swapAdaptAtt << '\n';
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
    is >> gen >> storedSf >> deltaT >> nm >> swapAdaptCount >> swapAdaptAcc >> swapAdaptAtt;
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

void MetropolisCoupledMcmc::adaptLadder(bool swapAccepted) {
    swapAdaptAtt++;
    if(swapAccepted)
        swapAdaptAcc++;
    if(swapAdaptAtt >= 100){
        double ar = (double)swapAdaptAcc / (double)swapAdaptAtt;
        swapAdaptCount++;
        double gain = 1.0 / std::sqrt((double)swapAdaptCount);
        deltaT *= std::exp(gain * (ar - 0.40));
        if(deltaT < 1e-4) deltaT = 1e-4;
        if(deltaT > 10.0)  deltaT = 10.0;
        rebuildLadder();
        swapAdaptAcc = 0;
        swapAdaptAtt = 0;
    }
}

void MetropolisCoupledMcmc::rebuildLadder(void){
    if(numModels <= 2)
        return;
    const int G = numModels - 1;
    const double betaFloor = 1.0 / (1.0 + deltaT * G);
    std::vector<double> oldB = betas;
    std::vector<double> cum(numModels, 0.0);
    for(int k = 0; k < G; k++)
        cum[k + 1] = cum[k] + (1.0 - emaAcc[k]);
    const double Lam = cum[G];
    betas[0] = 1.0;
    betas[G] = betaFloor;
    if(Lam <= 1e-9){
        for(int k = 1; k < G; k++)
            betas[k] = 1.0 - (1.0 - betaFloor) * (double)k / (double)G;
        return;
    }
    int m = 0;
    for(int k = 1; k < G; k++){
        const double want = Lam * (double)k / (double)G;
        while(m < G - 1 && cum[m + 1] < want)
            m++;
        const double seg = cum[m + 1] - cum[m];
        const double frac = (seg > 1e-12) ? (want - cum[m]) / seg : 0.0;
        double b = oldB[m] + (oldB[m + 1] - oldB[m]) * frac;
        const double hi = betas[k - 1] - 1e-9;
        const double lo = betaFloor + 1e-9 * (double)(G - k);
        betas[k] = (b > hi) ? hi : (b < lo ? lo : b);
    }
}
