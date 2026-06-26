#include "Mcmc.hpp"
#include "Msg.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"
#include "WriteTSV.hpp"

#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>


Mcmc::Mcmc(int ng, int pf, int sf, PhylogeneticModel* m) : model(m), numCycles(ng), printFrequency(pf), sampleFrequency(sf), gen(0), curLnL(0.0), curLnP(0.0) {
    UserSettings& settings = UserSettings::userSettings();
    treeOut = settings.getTreeOutput();
    paramOut = settings.getParamOutput();
}

void Mcmc::run(void) {
    init();
    advance((unsigned long)numCycles);
    finalize();
}

void Mcmc::init(void) {
    if(numCycles >= std::numeric_limits<unsigned long>::max())
        Msg::error("numCycles requested in greater than largest possible value");
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

    trees.appendDataTSV(model->getTree()->getNewickString());

    std::vector<double> tv = {lnL + lnP, lnL, lnP};
    tv.insert(tv.end(), parmStr.begin(), parmStr.end());
    for(size_t j = 0; j < traceCols.size() && j < tv.size(); j++)
        traceCols[j].push_back(tv[j]);
}
