#include "Mcmc.hpp"
#include "Msg.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"
#include "WriteTSV.hpp"

#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>


Mcmc::Mcmc(int ng, int pf, int sf, PhylogeneticModel* m) : model(m), numCycles(ng), printFrequency(pf), sampleFrequency(sf), gen(0), curLnL(0.0), curLnP(0.0) {
    UserSettings& settings = UserSettings::userSettings();
    treeOut = settings.getTreeOutput();
    paramOut = settings.getParamOutput();
    writeTrees = (model->getTree()->getNumBackbone() > 0);
}

void Mcmc::run(void) {
    init();
    advance((unsigned long)numCycles);
    finalize();
}

void Mcmc::init(void) {
    RandomVariable::setActiveInstance(model->getRng());
    curLnL = model->lnLikelihood();
    curLnP = model->lnPriorProbability();
    gen = 0;
}

void Mcmc::advance(unsigned long nGens) {
    RandomVariable::setActiveInstance(model->getRng());
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    unsigned long target = gen + nGens;
    while (gen < target) {
        gen++;
        unsigned long n = gen;

        double lnProposalRatio = model->update();
        double newLnL = model->lnLikelihood();
        double lnLikelihoodRatio = newLnL - curLnL;
        double newLnP = model->lnPriorProbability();
        double lnPriorRatio = newLnP - curLnP;
        double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;

        bool acceptMove = (std::log(rng.uniformRv()) < lnR);

        if (n % printFrequency == 0)
            {
            std::cout << std::fixed << std::setprecision(2);
            std::cout << n << " -- " << curLnL << " -> " << newLnL << "\n";
            model->print();
            }

        if (acceptMove == true)
            {
            curLnL = newLnL;
            curLnP = newLnP;
            model->updateForAcceptance();
            }
        else
            {
            model->updateForRejection();
            }

        if (n == 1 || n == (unsigned long)numCycles || n % sampleFrequency == 0)
            sample(n, curLnL, curLnP);
    }
}

void Mcmc::finalize(void) {
    params.closeTSV();
    trees.closeTSV();
}

void Mcmc::sample(unsigned long n, double lnL, double lnP) {
    bool cpuTime = UserSettings::userSettings().getCpuTime();
    if(n == 1){
        params.addFilepath(paramOut, true);
        std::vector<std::string> cn = {"n", "posterior", "likelihood", "prior"};
        std::vector<std::string> headStr = model->getParameterNames();
        cn.insert( cn.end(), headStr.begin(), headStr.end() );
        if(cpuTime)
            cn.push_back("cpu_s");
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

    std::vector<double> dat = {(double)n, lnL + lnP, lnL, lnP};
    std::vector<double> parmStr = model->getParameterString();
    dat.insert( dat.end(), parmStr.begin(), parmStr.end() );
    if(cpuTime)
        dat.push_back((double)std::clock() / CLOCKS_PER_SEC);
    params.appendDataTSV(dat);

    if(writeTrees)
        trees.appendDataTSV(model->getTree()->getBackboneNewickString());

    std::vector<double> tv = {lnL + lnP, lnL, lnP};
    tv.insert(tv.end(), parmStr.begin(), parmStr.end());
    for(size_t j = 0; j < traceCols.size() && j < tv.size(); j++)
        traceCols[j].push_back(tv[j]);

    writeCheckpoint();
}

void Mcmc::writeCheckpoint(void) {
    std::string path = paramOut + ".ckp";
    std::string tmp = path + ".tmp";
    std::ofstream os(tmp);
    os << std::setprecision(17);
    os << gen << ' ' << sampleFrequency << '\n';
    os << curLnL << ' ' << curLnP << '\n';
    model->getRng()->writeState(os);
    model->writeState(os);
    os.flush();
    os.close();
    std::rename(tmp.c_str(), path.c_str());
}

bool Mcmc::loadCheckpoint(void) {
    std::string path = paramOut + ".ckp";
    std::ifstream is(path);
    if(is.is_open() == false)
        Msg::error("could not open checkpoint file for -resume: " + path);
    int storedSf;
    is >> gen >> storedSf;
    if(storedSf != sampleFrequency){
        Msg::warning("-s " + std::to_string(sampleFrequency) + " differs from the pre-resume thinning " + std::to_string(storedSf) + "; forcing " + std::to_string(storedSf));
        sampleFrequency = storedSf;
    }
    is >> curLnL >> curLnP;
    model->getRng()->readState(is);
    model->readState(is);
    if(is.fail())
        Msg::error("checkpoint file is truncated or corrupt: " + path);
    return true;
}

void Mcmc::resumeOutputs(void) {
    RandomVariable::setActiveInstance(model->getRng());

    std::vector<std::string> headStr = model->getParameterNames();
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
