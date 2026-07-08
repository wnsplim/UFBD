#include "ConvergenceRunner.hpp"
#include "ChainRunner.hpp"
#include "Convergence.hpp"
#include "Tree.hpp"
#include "Msg.hpp"
#include "ThreadPool.hpp"
#include "UserSettings.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>

static std::string bulkPath(const std::string& p){
    size_t d = p.rfind('.');
    return (d == std::string::npos) ? p + "_bulk" : p.substr(0, d) + "_bulk" + p.substr(d);
}

ConvergenceRunner::ConvergenceRunner(std::vector<ChainRunner*> reps, const std::string& po, const std::string& to) : replicates(reps), paramBase(po), treeBase(to) {
    int M = (int)replicates.size();
    for(int r = 0; r < M; r++){
        std::string pf = paramBase, tf = treeBase;
        if(M > 1){
            std::string tag = "_chain" + std::to_string(r);
            size_t dp = paramBase.rfind('.');
            pf = (dp == std::string::npos) ? paramBase + tag : paramBase.substr(0, dp) + tag + paramBase.substr(dp);
            size_t dt = treeBase.rfind('.');
            tf = (dt == std::string::npos) ? treeBase + tag : treeBase.substr(0, dt) + tag + treeBase.substr(dt);
        }
        repParamFiles.push_back(pf);
        repTreeFiles.push_back(tf);
        replicates[r]->setOutputPaths(pf, tf);
    }
}

bool ConvergenceRunner::run(void){
    UserSettings& s = UserSettings::userSettings();
    if(verbose)
        for(int r = 0; r < (int)replicates.size(); r++){
            replicates[r]->setVerbose(true);
            replicates[r]->setLabel(r);
        }
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
        if(autoStop && report(gen, false)){ stoppedEarly = true; break; }
    }

    for(ChainRunner* c : replicates)
        c->finalize();

    for(ThreadPool* p : runPools)
        delete p;

    std::cout << "-----------------------------------------------------------------------\n";
    report(gen, true);
    writeMerged();

    if(emitSummary && replicates.empty() == false && replicates[0]->getTree() != nullptr){
        const std::vector<std::string>& names = replicates[0]->traceNames();
        std::vector<std::vector<double>> pooled(names.size());
        double burn = s.getBurninFraction();
        for(ChainRunner* c : replicates){
            const std::vector<std::vector<double>>& cols = c->traceColumns();
            for(size_t j = 0; j < cols.size() && j < pooled.size(); j++){
                size_t b = (size_t)(cols[j].size() * burn);
                pooled[j].insert(pooled[j].end(), cols[j].begin() + b, cols[j].end());
            }
        }
        std::string base = paramBase;
        size_t dp = base.rfind(".log");
        if(dp != std::string::npos) base = base.substr(0, dp);
        if(base.empty()) base = "out";
        writeSummaryTree(replicates[0]->getTree(), names, pooled, 0.0, base + ".tree");
    }

    if(autoStop){
        std::string crit = (replicates.size() >= 2) ? "R-hat/ESS" : "ESS";
        if(stoppedEarly)
            std::cout << "Auto-stop at generation " << gen << ": " << crit << " stopping thresholds reached for all logged quantities.\n";
        else
            Msg::warning("reached the -max-gen limit (" + std::to_string(maxGen) + ") with quantities still below the " + crit + " thresholds.");
    }
    return stoppedEarly;
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
    if(finalPass == false && minPost < (size_t)essMin)
        return false;

    bool multiChain = replicates.size() >= 2;
    double worstRhat = 1.0;
    std::string worstRhatName;
    double minEss = -1.0;
    std::string minEssName;
    double minBulkEss = -1.0;
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
        bool isDensity = (names[p] == "posterior" || names[p] == "likelihood" || names[p] == "prior");
        if(multiChain && isDensity == false && std::isnan(d.rhat) == false && d.rhat > worstRhat){ worstRhat = d.rhat; worstRhatName = names[p]; }
        if(std::isnan(d.bulkEss) == false && (minBulkEss < 0.0 || d.bulkEss < minBulkEss)) minBulkEss = d.bulkEss;
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
        bool rhatOk = (multiChain == false) || isDensity || std::isnan(d.rhat) || d.rhat <= rhatMax;
        bool essOk = (worstChainEss < 0.0) || worstChainEss >= essMin;
        if(rhatOk == false || essOk == false){
            nBelow++;
            if((int)bad.size() < 8) bad.push_back(names[p]);
        }
    }

    lastMaxRhat = worstRhat; lastMinChainEss = minEss; lastMinBulkEss = minBulkEss;

    if((long)gen != lastReportedGen){
        std::cout << std::fixed;
        std::cout << "gen " << gen;
        if(multiChain)
            std::cout << "  max R-hat " << std::setprecision(4) << worstRhat << " (" << worstRhatName << ")";
        std::cout << "  min per-chain ESS " << std::setprecision(0) << minEss << " (" << minEssName << ")\n";
        lastReportedGen = (long)gen;
    }

    if(finalPass && nBelow > 0 && s.getAutoChainLength()){
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
    int thin = UserSettings::userSettings().getThinning();
    long g = 0;

    std::ofstream mp(bulkPath(paramBase));
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
        for(int i = 1 + b; i < (int)lines.size(); i++){
            g += thin;
            size_t tab = lines[i].find('\t');
            if(tab == std::string::npos) mp << lines[i] << "\n";
            else mp << g << lines[i].substr(tab) << "\n";
        }
    }
    mp.close();

    bool anyTree = false;
    for(const std::string& fn : repTreeFiles){
        std::ifstream in(fn);
        if(in.is_open()){ anyTree = true; in.close(); break; }
    }
    if(anyTree == false)
        return;

    std::ofstream mt(bulkPath(treeBase));
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
