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


Mcmc::Mcmc(int ng, int thin, PhylogeneticModel* m) : model(m), numCycles(ng), thinning(thin), curLnL(0.0), curLnP(0.0) {
    UserSettings& settings = UserSettings::userSettings();
    gen = 0;
    treeOut = settings.getTreeOutput();
    paramOut = settings.getParamOutput();
    writeTrees = (treeOut.empty() == false && model->getTree()->getNumBackbone() > 0);
}

void Mcmc::printMoveDiagnostics(int rep) {
    std::cout << "  run " << rep << "  ";
    model->print();
}

Tree* Mcmc::getTree(void) { return model->getTree(); }
bool Mcmc::treeHasFossils(void) { return model->treeIncludesFossils(); }

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

        if (verbose && n % thinning == 0) {
            std::ostringstream os;
            os << std::fixed << std::setprecision(2) << "chain " << runLabel << "  " << n << " -- posterior " << (curLnL + curLnP) << " likelihood " << curLnL << "\n";
            ChainRunner::logLine(os.str());
        }

        if (n == 1 || n == (unsigned long)numCycles || n % thinning == 0)
            sample(n, curLnL, curLnP);
    }
}

void Mcmc::finalize(void) {
    params.closeTSV();
    if(writeTrees)
        trees.appendDataTSV("END;");
    trees.closeTSV();
    if(writeLatent)
        latent.closeTSV();
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

        std::vector<std::string> latNames = model->getLatentNames();
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

    std::vector<double> dat = {(double)n, lnL + lnP, lnL, lnP};
    std::vector<double> parmStr = model->getParameterString();
    dat.insert( dat.end(), parmStr.begin(), parmStr.end() );
    if(cpuTime)
        dat.push_back((double)std::clock() / CLOCKS_PER_SEC);
    params.appendDataTSV(dat);

    if(writeLatent){
        std::vector<double> ls = model->getLatentString();
        std::vector<double> ldat = {(double)n};
        ldat.insert(ldat.end(), ls.begin(), ls.end());
        latent.appendDataTSV(ldat);
        for(size_t j = 0; j < latentCols.size() && j < ls.size(); j++)
            latentCols[j].push_back(ls[j]);
    }

    if(writeTrees)
        trees.appendDataTSV("\tTREE STATE_" + std::to_string(n) + " = [&R] " + model->getTree()->getBackboneNewickString(model->treeIncludesFossils()));

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
    os << gen << ' ' << thinning << '\n';
    os << curLnL << ' ' << curLnP << '\n';
    model->getRng()->writeState(os);
    model->writeState(os);
    os.flush();
    os.close();
    std::rename(tmp.c_str(), path.c_str());
}

bool Mcmc::loadCheckpoint(void) {
    std::string path = paramOut + ".ckp";
    std::ifstream is = openCheckpoint(path);
    int storedSf;
    is >> gen >> storedSf;
    reconcileThinning(storedSf, thinning);
    is >> curLnL >> curLnP;
    model->getRng()->readState(is);
    model->readState(is);
    requireCheckpointIntact(is, path);
    return true;
}

RandomVariable* Mcmc::resumeRng(void) {
    return model->getRng();
}

std::vector<std::string> Mcmc::resumeParameterNames(void) {
    return model->getParameterNames();
}
