#include "FBDInput.hpp"
#include "FBDTreeModel.hpp"
#include "Mcmc.hpp"
#include "MetropolisCoupledMcmc.hpp"
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
    bool dating = seq || (settings.getHessianFile().empty() == false);
    ClockModel cm = ClockModel::UCLN;
    std::string cn = settings.getClockModelName();
    if(cn == "wn")        cm = ClockModel::WN;
    else if(cn == "gbm")  cm = ClockModel::GBM;
    else if(cn == "gbmc") cm = ClockModel::GBMC;

    int numChains = settings.getNumChains();
    if(numChains > 1){
        std::cout << "Running Metropolis-coupled MCMC with " << numChains << " chains parallelized across " << settings.getNumThreads() << " threads \n";
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::vector<PhylogeneticModel*> models;
        models.resize(numChains);
        for(int i = 0; i < numChains; i++)
            models[i] = seq
                ? (PhylogeneticModel*)new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), settings.getSequenceFile(), settings.getPartitionFile(), settings.getModelNStates(), settings.getNumCats(), cm, settings.getRgeneGamma(), settings.getSigma2Gamma(), masterSeed + i)
                : dating
                    ? (PhylogeneticModel*)new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), settings.getHessianFile(), settings.getTreeFile(), settings.getModelNStates(), cm, settings.getRgeneGamma(), settings.getSigma2Gamma(), masterSeed + i)
                    : (PhylogeneticModel*)new FBDTreeModel(pt, input.getClades(), input.getFossils(), masterSeed + i);
        MetropolisCoupledMcmc mcmc(settings.getChainLength(), settings.getPrintFrequency(), settings.getSampleFrequency(), models, masterSeed);
        mcmc.run();
    }else if (numChains == 1){
        std::cout << "Running standard MCMC \n";
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        PhylogeneticModel* model = seq
            ? (PhylogeneticModel*)new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), settings.getSequenceFile(), settings.getPartitionFile(), settings.getModelNStates(), settings.getNumCats(), cm, settings.getRgeneGamma(), settings.getSigma2Gamma(), masterSeed)
            : dating
                ? (PhylogeneticModel*)new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), settings.getHessianFile(), settings.getTreeFile(), settings.getModelNStates(), cm, settings.getRgeneGamma(), settings.getSigma2Gamma(), masterSeed)
                : (PhylogeneticModel*)new FBDTreeModel(pt, input.getClades(), input.getFossils(), masterSeed);
        Mcmc mcmc(settings.getChainLength(), settings.getPrintFrequency(), settings.getSampleFrequency(), model);
        mcmc.run();
    }

    return 0;
}
