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
#include "PathModel.hpp"
#include "BridgeSampler.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <streambuf>
#include <execinfo.h>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#ifdef __linux__
#include <sys/prctl.h>
#endif

static void crashBacktrace(int sig){
    void* frames[64];
    int n = backtrace(frames, 64);
    fprintf(stderr, "\n=== fatal signal %d — backtrace (%d frames) ===\n", sig, n);
    backtrace_symbols_fd(frames, n, 2);
    _exit(128 + sig);
}

class TeeBuf : public std::streambuf {
    public:
        TeeBuf(std::streambuf* a, std::streambuf* b) : bufA(a), bufB(b) {}
        void setB(std::streambuf* b) { bufB = b; }
    protected:
        int overflow(int c) override {
            if(c == traits_type::eof()) return c;
            int r1 = bufA->sputc((char)c);
            int r2 = bufB->sputc((char)c);
            return (r1 == traits_type::eof() || r2 == traits_type::eof()) ? traits_type::eof() : c;
        }
        int sync() override { return (bufA->pubsync() == 0 && bufB->pubsync() == 0) ? 0 : -1; }
    private:
        std::streambuf* bufA;
        std::streambuf* bufB;
};

static std::streambuf* gOrigCout = nullptr;
static std::streambuf* gOrigCerr = nullptr;
static void restoreConsoleStreams(){
    if(gOrigCout) std::cout.rdbuf(gOrigCout);
    if(gOrigCerr) std::cerr.rdbuf(gOrigCerr);
}

static std::string clockString(std::chrono::system_clock::time_point tp){
    std::time_t t = std::chrono::system_clock::to_time_t(tp);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    return std::string(buf);
}

static std::string elapsedString(double secs){
    long s = (long)(secs + 0.5);
    long h = s / 3600; s %= 3600;
    long m = s / 60;   s %= 60;
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%ldh %02ldm %02lds", h, m, s);
    return std::string(buf);
}

int main(int argc, const char* argv[]) {

#ifdef __linux__
    if(const char* parent = std::getenv("UFBD_BF_PARENT")){
        prctl(PR_SET_PDEATHSIG, SIGKILL);
        if(getppid() != (pid_t)std::atol(parent))
            return 1;
    }
#endif
    signal(SIGSEGV, crashBacktrace);
    signal(SIGABRT, crashBacktrace);

    std::ostringstream preBuf;
    static std::ofstream consoleFile;
    static TeeBuf teeOut(std::cout.rdbuf(), preBuf.rdbuf());
    static TeeBuf teeErr(std::cerr.rdbuf(), preBuf.rdbuf());
    gOrigCout = std::cout.rdbuf(&teeOut);
    gOrigCerr = std::cerr.rdbuf(&teeErr);
    std::atexit(restoreConsoleStreams);

    auto runStart = std::chrono::system_clock::now();
    bool wantsHelp = false;
    for(int i = 1; i < argc; i++){
        std::string a = argv[i];
        if(a == "-h" || a == "-help"){ wantsHelp = true; break; }
    }
    if(wantsHelp == false)
        std::cout << "Run started:  " << clockString(runStart) << "\n";

    UserSettings& settings = UserSettings::userSettings();
    Msg::setDeferWarnings(true);
    settings.initializeSettings(argc, argv);

    std::string logBase = settings.getParamOutput();
    std::string consolePath = "console.txt";
    std::string bfBase;
    bool ladder = settings.getModelBSettings() != nullptr && settings.getBfLadder().size() != 1;
    if(logBase.empty() == false){
        size_t p = logBase.rfind(".log");
        std::string stem = (p != std::string::npos ? logBase.substr(0, p) : logBase);
        bfBase = stem;
        consolePath = stem + (ladder ? "_bf" : "") + ".console.txt";
    }
    for(const std::string& pth : {settings.getParamOutput(), settings.getTreeOutput(), consolePath}){
        if(pth.empty()) continue;
        std::filesystem::path d = std::filesystem::path(pth).parent_path();
        if(d.empty() == false) std::filesystem::create_directories(d);
    }
    consoleFile.open(consolePath);
    if(consoleFile.is_open()){
        consoleFile << preBuf.str();
        teeOut.setB(consoleFile.rdbuf());
        teeErr.setB(consoleFile.rdbuf());
    }else{
        restoreConsoleStreams();
    }

    FBDInput input(settings.getTreeFile(), settings.getCladesFile(), settings.getFossilFile(), settings.getUnsampledExtantsFile());
    Tree* pt = input.getTree();

    if(settings.getPathPilot())
        settings.setChainLength(BridgeSampler::pilotGenerations());
    settings.print();
    if(Msg::hasDeferredWarnings()){
        Msg::flushWarnings();
        std::cout << "-----------------------------------------------------------------------\n";
    }
    Msg::setDeferWarnings(false);

    bool seq = settings.getSequenceFile().empty() == false;
    bool hessian = settings.getHessianFile().empty() == false;
    if(seq && hessian)
        Msg::error("-sequence and -hessian flags are mutually exclusive.");
    if(settings.clockOrCtmcConfigured() && seq == false && hessian == false)
        Msg::warning("no -sequence or -hessian given, so the clock and substitution settings are ignored.");
    std::vector<ClockModel> cm;
    for(const std::string& cn : settings.getClockModelNames())
        cm.push_back(cn == "gbm" ? ClockModel::GBM : ClockModel::UCLN);
    // WN + GBMC clock: halt — detached, not selectable

    int numCoupledChains = settings.getNumCoupledChains();
    int numRuns = settings.getNumRuns();
    bool autoStop = settings.getAutoChainLength();
    int thin = settings.getThinning();

    auto makeModel = [&](unsigned int sd, int run) -> PhylogeneticModel* {
        PhylogeneticModel* model;
        if(seq)
            model = new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), settings.getSequenceFile(), settings.getPartitionFile(), settings.getModelNStates(), settings.getNumCats(), cm, settings.getRgeneGamma(), settings.getSigma2Gamma(), sd);
        else if(hessian)
            model = new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), settings.getHessianFile(), settings.getTreeFile(), settings.getModelNStates(), cm, settings.getRgeneGamma(), settings.getSigma2Gamma(), sd);
        else
            model = new FBDTreeModel(pt, input.getClades(), input.getFossils(), sd);
        static bool printedMap = false;
        if(printedMap == false){
            printedMap = true;
            std::string rm;
            if(FBDTreeModel* fm = dynamic_cast<FBDTreeModel*>(model))              rm = fm->getRateMap();
            else if(RelaxedClockTreeModel* rc = dynamic_cast<RelaxedClockTreeModel*>(model)) rm = rc->getRateMap();
            if(rm.empty() == false) std::cout << "Rate intervals:\n" << rm;
        }
        if(settings.getModelBSettings() != nullptr){
            UserSettings* bs = settings.getModelBSettings();
            unsigned int sdB = sd + 104729u;
            std::vector<ClockModel> cmB;
            for(const std::string& cn : bs->getClockModelNames())
                cmB.push_back(cn == "gbm" ? ClockModel::GBM : ClockModel::UCLN);
            PhylogeneticModel* mb;
            if(seq)
                mb = new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), bs->getSequenceFile(), bs->getPartitionFile(), bs->getModelNStates(), bs->getNumCats(), cmB, bs->getRgeneGamma(), bs->getSigma2Gamma(), sdB, bs);
            else if(hessian)
                mb = new RelaxedClockTreeModel(pt, input.getClades(), input.getFossils(), bs->getHessianFile(), bs->getTreeFile(), bs->getModelNStates(), cmB, bs->getRgeneGamma(), bs->getSigma2Gamma(), sdB, bs);
            else
                mb = new FBDTreeModel(pt, input.getClades(), input.getFossils(), sdB, bs);
            PathModel* path = new PathModel(model, mb, input.getFossils(), sd + 7919u);
            path->setStart(run);
            model = path;
        }
        return model;
    };

    auto runOnce = [&](){
    bool resume = settings.getResume();
    unsigned int masterSeed = settings.getSeed() + (unsigned int)std::llround(settings.getPathT() * 1e6) * 1000003u + (settings.getPathPilot() ? 999983u : 0u);

    if(numRuns == 1 && autoStop == false){
        ChainRunner* mcmc;
        if(numCoupledChains > 1){
            std::vector<PhylogeneticModel*> models(numCoupledChains);
            for(int c = 0; c < numCoupledChains; c++)
                models[c] = makeModel(masterSeed + (unsigned int)c, 0);
            mcmc = new MetropolisCoupledMcmc(settings.getChainLength(), thin, models, masterSeed);
        }else
            mcmc = new Mcmc((int)settings.getChainLength(), thin, makeModel(masterSeed, 0));
        mcmc->setVerbose(true);
        mcmc->setLabel(0);
        if(resume && std::ifstream(settings.getParamOutput() + ".ckp").good() == false && std::ifstream(ChainRunner::tunePath(settings.getParamOutput()) + ".ckp").good() == false)
            Msg::error("-resume: nothing to resume from.");
        if(resume && std::ifstream(settings.getParamOutput() + ".ckp").good()){
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
        if(settings.getArLog())
            mcmc->printMoveDiagnostics(0);
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
                    models[c] = makeModel(base + (unsigned int)c, r);
                reps.push_back(new MetropolisCoupledMcmc(ncyc, thin, models, base));
            }else{
                int ng = (ncyc > 2000000000UL) ? 2000000000 : (int)ncyc;
                reps.push_back(new Mcmc(ng, thin, makeModel(base, r)));
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
    if(wrote(settings.getParamOutput())) std::cout << "MCMC log               -> " << settings.getParamOutput() << "\n";
    if(wrote(settings.getTreeOutput()))  std::cout << "tree log               -> " << settings.getTreeOutput() << "\n";
    if(wrote(base + ".tree"))            std::cout << "posterior summary tree -> " << base << ".tree\n";
    };

    if(ladder)
        BridgeSampler(bfBase).run();
    else
        runOnce();

    auto runEnd = std::chrono::system_clock::now();
    std::cout << "Run finished: " << clockString(runEnd)
              << "   elapsed " << elapsedString(std::chrono::duration<double>(runEnd - runStart).count()) << "\n";

    return 0;
}
