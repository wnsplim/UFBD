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

    unsigned int masterSeed = settings.getSeedSet() ? settings.getSeed() : std::random_device{}();

    FBDInput input(settings.getTreeFile(), settings.getCladesFile(), settings.getFossilFile());

    Tree* pt = input.getTree();
    pt->print();

    bool seq = settings.getSequenceFile().empty() == false;
    bool hessian = settings.getHessianFile().empty() == false;
    if(seq == false && hessian == false)
        Msg::warning("No sequence of Hessian file supplied: running the FBD model without a molecular clock.");
    ClockModel cm = ClockModel::UCLN;
    std::string cn = settings.getClockModelName();
    if(cn == "wn")        cm = ClockModel::WN;
    else if(cn == "gbm")  cm = ClockModel::GBM;
    // GBMC clock: halt — detached, not selectable

    int numCoupledChains = settings.getNumCoupledChains();
    int numRuns = settings.getNumRuns();
    bool autoStop = settings.getAutoChainLength();
    int pf = settings.getPrintFrequency();
    int sf = settings.getSampleFrequency();

    auto makeModel = [&](unsigned int sd) -> PhylogeneticModel* {
        if(seq)
            return new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), settings.getSequenceFile(), settings.getPartitionFile(), settings.getModelNStates(), settings.getNumCats(), cm, settings.getRgeneGamma(), settings.getSigma2Gamma(), sd);
        if(hessian)
            return new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), settings.getHessianFile(), settings.getTreeFile(), settings.getModelNStates(), cm, settings.getRgeneGamma(), settings.getSigma2Gamma(), sd);
        return new FBDTreeModel(pt, input.getClades(), input.getFossils(), sd);
    };

    bool resume = settings.getResume();

    if(numRuns == 1 && autoStop == false){
        if(numCoupledChains > 1){
            std::cout << "Running Metropolis-coupled MCMC with " << numCoupledChains << " coupled chains across " << settings.getNumThreads() << " threads\n";
            std::cout << "-----------------------------------------------------------------------" << std::endl;
            std::vector<PhylogeneticModel*> models(numCoupledChains);
            for(int c = 0; c < numCoupledChains; c++)
                models[c] = makeModel(masterSeed + (unsigned int)c);
            MetropolisCoupledMcmc mcmc(settings.getChainLength(), pf, sf, models, masterSeed);
            if(resume){
                mcmc.loadCheckpoint();
                mcmc.resumeOutputs();
                unsigned long g = mcmc.currentGen();
                if(g < settings.getChainLength())
                    mcmc.advance(settings.getChainLength() - g);
                mcmc.finalize();
            }else
                mcmc.run();
        }else{
            std::cout << "Running standard MCMC\n";
            std::cout << "-----------------------------------------------------------------------" << std::endl;
            Mcmc mcmc((int)settings.getChainLength(), pf, sf, makeModel(masterSeed));
            if(resume){
                mcmc.loadCheckpoint();
                mcmc.resumeOutputs();
                unsigned long g = mcmc.currentGen();
                if(g < settings.getChainLength())
                    mcmc.advance(settings.getChainLength() - g);
                mcmc.finalize();
            }else
                mcmc.run();
        }
    }else{
        unsigned long ncyc = autoStop ? settings.getMaxGen() : settings.getChainLength();
        std::cout << "Running " << numRuns << " independent runs";
        if(numCoupledChains > 1) std::cout << " of " << numCoupledChains << "-chain MC^3";
        if(autoStop) std::cout << ", auto-stopping on R-hat/ESS";
        std::cout << "\n-----------------------------------------------------------------------" << std::endl;
        std::vector<ChainRunner*> reps;
        for(int r = 0; r < numRuns; r++){
            unsigned int base = masterSeed + (unsigned int)(r * (numCoupledChains + 1));
            if(numCoupledChains > 1){
                std::vector<PhylogeneticModel*> models(numCoupledChains);
                for(int c = 0; c < numCoupledChains; c++)
                    models[c] = makeModel(base + (unsigned int)c);
                reps.push_back(new MetropolisCoupledMcmc(ncyc, pf, sf, models, base));
            }else{
                int ng = (ncyc > 2000000000UL) ? 2000000000 : (int)ncyc;
                reps.push_back(new Mcmc(ng, pf, sf, makeModel(base)));
            }
        }
        ConvergenceRunner cr(reps, settings.getParamOutput(), settings.getTreeOutput());
        cr.run();
        for(ChainRunner* c : reps)
            delete c;
    }

    return 0;
}
