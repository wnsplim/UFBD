#include "ChainRunner.hpp"
#include "ConvergenceRunner.hpp"
#include "FBDInput.hpp"
#include "FBDTreeModel.hpp"
#include "Mcmc.hpp"
#include "MetropolisCoupledMcmc.hpp"
#include "Msg.hpp"
#include "ParameterBranchRates.hpp"
#include "RandomVariable.hpp"
#include "RelaxedClockTreeModel.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#include <fstream>
#include <iostream>
#include <random>
#include <execinfo.h>
#include <csignal>
#include <cstdio>
#include <cstdlib>

static void crashBacktrace(int sig){
    void* frames[64];
    int n = backtrace(frames, 64);
    fprintf(stderr, "\n=== fatal signal %d — backtrace (%d frames) ===\n", sig, n);
    backtrace_symbols_fd(frames, n, 2);
    _exit(128 + sig);
}

int main(int argc, const char* argv[]) {

    signal(SIGSEGV, crashBacktrace);
    signal(SIGABRT, crashBacktrace);

    UserSettings& settings = UserSettings::userSettings();
    settings.initializeSettings(argc, argv);
    settings.print();

    unsigned int masterSeed = settings.getSeed();

    FBDInput input(settings.getTreeFile(), settings.getCladesFile(), settings.getFossilFile());

    Tree* pt = input.getTree();

    bool seq = settings.getSequenceFile().empty() == false;
    bool hessian = settings.getHessianFile().empty() == false;
    if(seq && hessian)
        Msg::error("-seq and -hessian are mutually exclusive.");
    if(settings.clockOrCtmcConfigured() && seq == false && hessian == false)
        Msg::warning("ㅜo -seq or -hessian provided, running pure FBD model.");
    ClockModel cm = ClockModel::UCLN;
    std::string cn = settings.getClockModelName();
    if(cn == "gbm")  cm = ClockModel::GBM;
    // WN + GBMC clock: halt — detached, not selectable

    int numCoupledChains = settings.getNumCoupledChains();
    int numRuns = settings.getNumRuns();
    bool autoStop = settings.getAutoChainLength();
    int thin = settings.getThinning();

    auto makeModel = [&](unsigned int sd) -> PhylogeneticModel* {
        if(seq)
            return new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), settings.getSequenceFile(), settings.getPartitionFile(), settings.getModelNStates(), settings.getNumCats(), cm, settings.getRgeneGamma(), settings.getSigma2Gamma(), sd);
        if(hessian)
            return new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), settings.getHessianFile(), settings.getTreeFile(), settings.getModelNStates(), cm, settings.getRgeneGamma(), settings.getSigma2Gamma(), sd);
        return new FBDTreeModel(pt, input.getClades(), input.getFossils(), sd);
    };

    bool resume = settings.getResume();

    if(numRuns == 1 && autoStop == false){
        ChainRunner* mcmc;
        if(numCoupledChains > 1){
            std::vector<PhylogeneticModel*> models(numCoupledChains);
            for(int c = 0; c < numCoupledChains; c++)
                models[c] = makeModel(masterSeed + (unsigned int)c);
            mcmc = new MetropolisCoupledMcmc(settings.getChainLength(), thin, models, masterSeed);
        }else
            mcmc = new Mcmc((int)settings.getChainLength(), thin, makeModel(masterSeed));
        mcmc->setVerbose(true);
        mcmc->setLabel(0);
        if(resume){
            mcmc->loadCheckpoint();
            mcmc->resumeOutputs();
            unsigned long g = mcmc->currentGen();
            if(g < settings.getChainLength())
                mcmc->advance(settings.getChainLength() - g);
            mcmc->finalize();
        }else{
            mcmc->init();
            mcmc->advance(settings.getChainLength());
            mcmc->finalize();
        }
        std::string sbase = settings.getParamOutput();
        size_t sdp = sbase.rfind(".log"); if(sdp != std::string::npos) sbase = sbase.substr(0, sdp);
        if(sbase.empty()) sbase = "out";
        if(mcmc->getTree() != nullptr)
            writeSummaryTree(mcmc->getTree(), mcmc->traceNames(), mcmc->traceColumns(), mcmc->latentNames(), mcmc->latentColumns(), settings.getBurninFraction(), sbase + ".tree", mcmc->treeHasFossils());
        delete mcmc;
    }else{
        unsigned long ncyc = autoStop ? settings.getMaxGen() : settings.getChainLength();
        std::vector<ChainRunner*> reps;
        for(int r = 0; r < numRuns; r++){
            unsigned int base = masterSeed + (unsigned int)(r * (numCoupledChains + 1));
            if(numCoupledChains > 1){
                std::vector<PhylogeneticModel*> models(numCoupledChains);
                for(int c = 0; c < numCoupledChains; c++)
                    models[c] = makeModel(base + (unsigned int)c);
                reps.push_back(new MetropolisCoupledMcmc(ncyc, thin, models, base));
            }else{
                int ng = (ncyc > 2000000000UL) ? 2000000000 : (int)ncyc;
                reps.push_back(new Mcmc(ng, thin, makeModel(base)));
            }
        }
        ConvergenceRunner cr(reps, settings.getParamOutput(), settings.getTreeOutput());
        cr.setEmitSummaryTree(true);
        cr.setVerbose(true);
        cr.run();
        for(ChainRunner* c : reps)
            delete c;
    }

    std::string base = settings.getParamOutput();
    size_t dp = base.rfind(".log"); if(dp != std::string::npos) base = base.substr(0, dp);
    auto wrote = [](const std::string& p){ if(p.empty()) return false; std::ifstream f(p); return f.good(); };
    std::cout << "-----------------------------------------------------------------------\n";
    if(wrote(base + "_bulk.log"))             std::cout << "MCMC log               -> " << base << "_bulk.log\n";
    else if(wrote(settings.getParamOutput())) std::cout << "MCMC log               -> " << settings.getParamOutput() << "\n";
    if(wrote(base + "_bulk.trees"))           std::cout << "tree log               -> " << base << "_bulk.trees\n";
    else if(wrote(settings.getTreeOutput()))  std::cout << "tree log               -> " << settings.getTreeOutput() << "\n";
    if(wrote(base + ".tree"))                 std::cout << "posterior summary tree -> " << base << ".tree\n";

    return 0;
}
