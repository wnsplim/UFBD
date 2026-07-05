#include "ConvergenceRunner.hpp"
#include "ChainRunner.hpp"
#include "Convergence.hpp"
#include "Msg.hpp"
#include "ThreadPool.hpp"
#include "UserSettings.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>

ConvergenceRunner::ConvergenceRunner(std::vector<ChainRunner*> reps, const std::string& po, const std::string& to) : replicates(reps), paramBase(po), treeBase(to) {
    int M = (int)replicates.size();
    for(int r = 0; r < M; r++){
        std::string suffix = (M > 1) ? (".run" + std::to_string(r)) : "";
        std::string pf = paramBase + suffix;
        std::string tf = treeBase + suffix;
        repParamFiles.push_back(pf);
        repTreeFiles.push_back(tf);
        replicates[r]->setOutputPaths(pf, tf);
    }
}

bool ConvergenceRunner::run(void){
    UserSettings& s = UserSettings::userSettings();
    bool autoStop = s.getAutoChainLength();
    unsigned long maxGen = autoStop ? s.getMaxGen() : s.getChainLength();
    long blockGens = (long)s.getEssThreshold() * (long)s.getThinning();
    if(blockGens < 1) blockGens = 1;

    unsigned long gen = 0;
    if(s.getResume()){
        for(ChainRunner* c : replicates){
            c->loadCheckpoint();
            c->resumeOutputs();
            gen = std::max(gen, c->currentGen());
        }
        for(ChainRunner* c : replicates){
            unsigned long cg = c->currentGen();
            if(cg < gen)
                c->advance(gen - cg);
        }
    }else{
        for(ChainRunner* c : replicates)
            c->init();
    }

    int M = (int)replicates.size();
    int nCores = UserSettings::userSettings().getNumCores();
    bool parallelRuns = (M > 1 && UserSettings::userSettings().getNumCoupledChains() == 1);
    int nWorkers = parallelRuns ? std::min(nCores, M) : 0;
    std::vector<ThreadPool*> runPools;
    for(int t = 0; t < nWorkers; t++){
        int per = nCores / nWorkers + (t < nCores % nWorkers ? 1 : 0);
        if(per < 1) per = 1;
        runPools.push_back(new ThreadPool(per));
    }

    bool stoppedEarly = false;
    while(gen < maxGen){
        unsigned long step = (unsigned long)blockGens;
        if(gen + step > maxGen) step = maxGen - gen;
        if(parallelRuns){
            std::vector<std::thread> workers;
            for(int t = 0; t < nWorkers; t++){
                workers.emplace_back([this, t, nWorkers, step, &runPools](){
                    ThreadPool::setCurrent(runPools[t]);
                    for(int r = t; r < (int)replicates.size(); r += nWorkers)
                        replicates[r]->advance(step);
                    ThreadPool::setCurrent(nullptr);
                });
            }
            for(std::thread& w : workers)
                w.join();
        }else{
            for(ChainRunner* c : replicates)
                c->advance(step);
        }
        gen += step;
        bool met = report(gen, false);
        if(autoStop && met){ stoppedEarly = true; break; }
    }

    for(ChainRunner* c : replicates)
        c->finalize();

    for(ThreadPool* p : runPools)
        delete p;

    std::cout << "-----------------------------------------------------------------------\n";
    bool met = report(gen, true);
    writeMerged();

    std::string crit = (replicates.size() >= 2) ? "R-hat/ESS" : "ESS";
    if(autoStop){
        if(stoppedEarly)
            std::cout << "Auto-stop at generation " << gen << ": " << crit << " stopping thresholds reached for all logged quantities.\n";
        else
            Msg::warning("reached the -max-gen limit (" + std::to_string(maxGen) + ") with quantities still below the " + crit + " thresholds.");
    }else if(met == false){
        Msg::warning("at the requested -n (" + std::to_string(maxGen) + ") some quantities remain below the " + crit + " thresholds.");
    }
    return met;
}

bool ConvergenceRunner::report(unsigned long gen, bool finalPass){
    UserSettings& s = UserSettings::userSettings();
    double burn = s.getBurninFraction();
    double rhatMax = s.getRhatThreshold();
    double essMin = s.getEssThreshold();

    const std::vector<std::string>& names = replicates[0]->traceNames();
    if(names.empty())
        return false;
    int nP = (int)names.size();

    size_t minPost = (size_t)-1;
    for(ChainRunner* c : replicates){
        const std::vector<std::vector<double>>& cols = c->traceColumns();
        size_t n = cols.empty() ? 0 : cols[0].size();
        size_t post = n - (size_t)(n * burn);
        if(post < minPost) minPost = post;
    }
    if(minPost < (size_t)essMin){
        std::cout << "gen " << gen << "  (accumulating: " << minPost << "/" << (int)essMin << " post-burnin samples)\n";
        return false;
    }

    bool multiChain = replicates.size() >= 2;
    double worstRhat = 1.0;
    std::string worstRhatName;
    double minEss = -1.0;
    std::string minEssName;
    int nBelow = 0;
    std::vector<std::string> bad;

    for(int p = 0; p < nP; p++){
        std::vector<std::vector<double>> chains;
        for(ChainRunner* c : replicates){
            const std::vector<double>& col = c->traceColumns()[p];
            int b = (int)(col.size() * burn);
            chains.push_back(std::vector<double>(col.begin() + b, col.end()));
        }
        Convergence::Diagnostic d = Convergence::diagnose(chains);
        if(multiChain && std::isnan(d.rhat) == false && d.rhat > worstRhat){ worstRhat = d.rhat; worstRhatName = names[p]; }
        double worstChainEss = -1.0;
        for(const std::vector<double>& ch : chains){
            Convergence::Diagnostic dc = Convergence::diagnose(std::vector<std::vector<double>>(1, ch));
            double ec = std::fmin(dc.bulkEss, dc.tailEss);
            if(std::isnan(ec))
                continue;
            if(worstChainEss < 0.0 || ec < worstChainEss)
                worstChainEss = ec;
        }
        if(worstChainEss >= 0.0 && (minEss < 0.0 || worstChainEss < minEss)){ minEss = worstChainEss; minEssName = names[p]; }
        bool rhatOk = (multiChain == false) || std::isnan(d.rhat) || d.rhat <= rhatMax;
        bool essOk = (worstChainEss < 0.0) || worstChainEss >= essMin;
        if(rhatOk == false || essOk == false){
            nBelow++;
            if((int)bad.size() < 8) bad.push_back(names[p]);
        }
    }

    std::cout << std::fixed;
    std::cout << "gen " << gen;
    if(multiChain)
        std::cout << "  max R-hat " << std::setprecision(4) << worstRhat << " (" << worstRhatName << ")";
    std::cout << "  min per-chain ESS " << std::setprecision(0) << minEss << " (" << minEssName << ")\n";

    if(finalPass && nBelow > 0){
        std::string crit = multiChain ? ("R-hat > " + std::to_string(rhatMax) + " or per-chain ESS < " + std::to_string((int)essMin))
                                      : ("per-chain ESS < " + std::to_string((int)essMin));
        std::string msg = std::to_string(nBelow) + " of " + std::to_string(nP) + " quantities have " + crit + " : ";
        for(size_t i = 0; i < bad.size(); i++){ msg += bad[i]; if(i + 1 < bad.size()) msg += ", "; }
        if((int)bad.size() < nBelow) msg += ", ...";
        Msg::warning(msg);
    }

    return nBelow == 0;
}

void ConvergenceRunner::writeMerged(void){
    if(replicates.size() < 2)
        return;
    double burn = UserSettings::userSettings().getBurninFraction();

    std::ofstream mp(paramBase);
    bool wroteHeader = false;
    for(const std::string& fn : repParamFiles){
        std::ifstream in(fn);
        if(in.is_open() == false) continue;
        std::vector<std::string> lines;
        std::string line;
        while(std::getline(in, line)) lines.push_back(line);
        in.close();
        if(lines.empty()) continue;
        if(wroteHeader == false){ mp << lines[0] << "\n"; wroteHeader = true; }
        int nData = (int)lines.size() - 1;
        int b = (int)(nData * burn);
        for(int i = 1 + b; i < (int)lines.size(); i++)
            mp << lines[i] << "\n";
    }
    mp.close();

    bool anyTree = false;
    for(const std::string& fn : repTreeFiles){
        std::ifstream in(fn);
        if(in.is_open()){ anyTree = true; in.close(); break; }
    }
    if(anyTree == false)
        return;

    std::ofstream mt(treeBase);
    for(const std::string& fn : repTreeFiles){
        std::ifstream in(fn);
        if(in.is_open() == false) continue;
        std::vector<std::string> lines;
        std::string line;
        while(std::getline(in, line)) lines.push_back(line);
        in.close();
        int b = (int)(lines.size() * burn);
        for(int i = b; i < (int)lines.size(); i++)
            mt << lines[i] << "\n";
    }
    mt.close();
}
