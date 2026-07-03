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
#include <future>
#include <iomanip>
#include <iostream>
#include <sstream>

MetropolisCoupledMcmc::MetropolisCoupledMcmc(unsigned long ng, int pf, int sf, std::vector<PhylogeneticModel*> m, unsigned int masterSeed) : numCycles(ng), printFrequency(pf),
    sampleFrequency(sf),
    models(m),
    numModels(models.size()),
    coldModelIdx(-1),
    deltaT(0.2),
    threadPool(1){
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
    writeTrees = (models[0]->getTree()->getNumBackbone() > 0);
    gen = 0;

    ThreadPool::shared().setChainCap(std::max(1, settings.getNumCores() / threadPool.size()));
}

double MetropolisCoupledMcmc::calcHeating(int idx){
    if(idx == 0)
        return 1.0;
    return (1 / (1 + deltaT * idx));
}

void MetropolisCoupledMcmc::run(void) {
    init();
    advance(numCycles);
    finalize();
}

void MetropolisCoupledMcmc::init(void) {
    if(numCycles >= std::numeric_limits<unsigned long>::max())
        Msg::error("Numer of cycles requested is greater than largest possible value");
    RandomVariable::setActiveInstance(&swapRng);
    int idx = 0;
    for(PhylogeneticModel* m : models){
        currLnL.push_back(m->lnLikelihood());
        currLnP.push_back(m->lnPriorProbability());
        indices.push_back(idx);
        idx++;
    }
    gen = 0;
}

void MetropolisCoupledMcmc::finalize(void) {
    params.closeTSV();
    trees.closeTSV();
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
        
        if (n % printFrequency == 0){
            std::cout << std::fixed << std::setprecision(2);
            const size_t numAccepted = std::count(recentAcceptRej.begin(), recentAcceptRej.end(), true);
            const double acceptanceRate = static_cast<double>(numAccepted) / recentAcceptRej.size();
            std::cout << n << " -- "
                << currLnL[coldModelIdx] << " -> " << newLnL[coldModelIdx]
                << " | Chain swap acceptance ratio: " << acceptanceRate << "\n";
            models[coldModelIdx]->print();
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
        
        if (n == 1 || n == numCycles || n % sampleFrequency == 0 )
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

        if(writeTrees)
            trees.addFilepath(treeOut, true); // no CN for tree file

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

    if(writeTrees)
        trees.appendDataTSV(models[coldModelIdx]->getTree()->getBackboneNewickString());

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
    os << gen << ' ' << sampleFrequency << ' ' << deltaT << ' ' << numModels << '\n';
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
    std::ifstream is(path);
    if(is.is_open() == false)
        Msg::error("could not open checkpoint file for -resume: " + path);
    int nm, storedSf;
    is >> gen >> storedSf >> deltaT >> nm;
    if(storedSf != sampleFrequency){
        Msg::warning("-s " + std::to_string(sampleFrequency) + " differs from the pre-resume thinning " + std::to_string(storedSf) + "; forcing " + std::to_string(storedSf) + " to keep the log spacing uniform.");
        sampleFrequency = storedSf;
    }
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
    if(is.fail())
        Msg::error("checkpoint file is truncated or corrupt: " + path);
    return true;
}

void MetropolisCoupledMcmc::resumeOutputs(void) {
    RandomVariable::setActiveInstance(&swapRng);

    std::vector<std::string> headStr = models[coldModelIdx]->getParameterNames();
    traceNms.clear();
    traceNms.push_back("posterior");
    traceNms.push_back("likelihood");
    traceNms.push_back("prior");
    for(const std::string& s : headStr)
        traceNms.push_back(s);
    traceCols.assign(traceNms.size(), std::vector<double>());

    std::ifstream pin(paramOut);
    std::string header;
    std::vector<std::string> keep;
    if(pin.is_open()){
        std::getline(pin, header);
        std::string line;
        while(std::getline(pin, line)){
            if(line.empty())
                continue;
            std::stringstream ss(line);
            double nval;
            ss >> nval;
            if((unsigned long)nval > gen)
                break;
            keep.push_back(line);
        }
        pin.close();
    }

    for(const std::string& line : keep){
        std::stringstream ss(line);
        double nval;
        ss >> nval;
        double v;
        for(size_t j = 0; j < traceCols.size() && (ss >> v); j++)
            traceCols[j].push_back(v);
    }

    std::ofstream pout(paramOut, std::ios::out | std::ios::trunc);
    pout << header << '\n';
    for(const std::string& line : keep)
        pout << line << '\n';
    pout.close();

    params.addFilepath(paramOut, false);
    if(writeTrees){
        std::ifstream tin(treeOut);
        std::vector<std::string> tkeep;
        std::string tline;
        while(std::getline(tin, tline) && tkeep.size() < keep.size())
            tkeep.push_back(tline);
        tin.close();
        std::ofstream tout(treeOut, std::ios::out | std::ios::trunc);
        for(const std::string& l : tkeep)
            tout << l << '\n';
        tout.close();
        trees.addFilepath(treeOut, false);
    }
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
