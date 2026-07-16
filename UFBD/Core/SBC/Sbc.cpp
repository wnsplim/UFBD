#include "Sbc.hpp"
#include "FBDInput.hpp"
#include "FBDTreeModel.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"
#include "ChainRunner.hpp"
#include "Mcmc.hpp"
#include "ConvergenceRunner.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace {

double drawPrior(const Probability::PriorSpec& s, RandomVariable* rng){
    using namespace Probability;
    switch(s.family){
        case PriorFamily::FIXED:            return s.p1;
        case PriorFamily::UNIFORM:          return Uniform::rv(rng, s.p1, s.p2);
        case PriorFamily::LOGNORMAL:        return std::exp(s.p1 + s.p2 * Normal::rv(rng));
        case PriorFamily::EXPONENTIAL:      return Exponential::rv(rng, s.p1);
        case PriorFamily::GAMMA:            return Gamma::rv(rng, s.p1, s.p2);
        case PriorFamily::TRUNCATED_NORMAL: return TruncatedNormal::rv(rng, s.p1, s.p2, 0.0, std::numeric_limits<double>::infinity());
        case PriorFamily::IMPROPER:         break;
    }
    return 0.0;
}

double drawChunk(const std::vector<int>& groups, const std::map<int,Probability::PriorSpec>& groupPrior,
                 const Probability::PriorSpec& base, int bin, std::map<int,double>& chunkVal, RandomVariable* rng){
    int gid = groups.empty() ? bin : groups[bin];
    std::map<int,double>::iterator it = chunkVal.find(gid);
    if(it != chunkVal.end())
        return it->second;
    std::map<int,Probability::PriorSpec>::const_iterator pit = groupPrior.find(gid);
    const Probability::PriorSpec& sp = (pit != groupPrior.end() && pit->second.set) ? pit->second : base;
    double v = drawPrior(sp, rng);
    chunkVal[gid] = v;
    return v;
}

double truthForName(const std::string& name, const SimParams& p, bool origin, const std::vector<std::string>& psiTypeNames){
    auto idxOf = [](const std::string& rest) -> int { return rest.empty() ? 0 : std::stoi(rest); };
    if(name.rfind("lambda", 0) == 0)  return p.lambda[idxOf(name.substr(6))];
    if(name.rfind("mu", 0) == 0)      return p.mu[idxOf(name.substr(2))];
    if(name.rfind("psi_", 0) == 0){
        std::string rest = name.substr(4);
        for(size_t ti = 0; ti < psiTypeNames.size(); ti++){
            const std::string& tn = psiTypeNames[ti];
            if(rest == tn)                       return p.psi[ti][0];
            if(rest.rfind(tn + "_", 0) == 0)     return p.psi[ti][std::stoi(rest.substr(tn.size() + 1))];
        }
        return std::numeric_limits<double>::quiet_NaN();
    }
    if(name.rfind("psi", 0) == 0)     return p.psi[0][idxOf(name.substr(3))];
    if(name == "originAge" && origin) return p.startAge;
    return std::numeric_limits<double>::quiet_NaN();
}

void initAges(Tree* tree, std::vector<Fossil>& fossils, bool origin, const Probability::PriorSpec& xp){
    double maxFoss = 0.0;
    for(Fossil& f : fossils)
        if(f.getMaxAge() > maxFoss) maxFoss = f.getMaxAge();
    int numInternal = tree->getNumNodes() - tree->getNumBackbone();
    double unit = (maxFoss > 0.0) ? maxFoss / (numInternal + 1) : 1.0;
    std::map<Node*,double> minAges;
    if(maxFoss > 0.0)
        minAges[tree->getRoot()] = maxFoss * 1.05;
    tree->assignStartingAges(minAges, unit);

    double mean = Probability::priorMean(xp.family, xp.p1, xp.p2, xp.p3);
    double target = origin ? 0.9 * mean : mean;
    if(origin == false && target < maxFoss)
        target = maxFoss * 1.05;
    double cur = tree->getRoot()->getTime();
    if(cur > 0.0)
        tree->scaleInternalAges(target / cur);
}

double ksPvalue(double d, long n){
    double t = std::sqrt((double)n) * d;
    if(t < 1e-12) return 1.0;
    double sum = 0.0;
    for(int k = 1; k <= 200; k++)
        sum += (k % 2 == 1 ? 1.0 : -1.0) * std::exp(-2.0 * k * k * t * t);
    double p = 2.0 * sum;
    if(p < 0.0) p = 0.0;
    if(p > 1.0) p = 1.0;
    return p;
}

double quantile(const std::vector<double>& s, double q){
    double idx = q * (s.size() - 1);
    size_t lo = (size_t)idx;
    double frac = idx - lo;
    if(lo + 1 < s.size())
        return s[lo] * (1.0 - frac) + s[lo + 1] * frac;
    return s[lo];
}

}

SimParams Sbc::drawParams(void){
    SimParams p;
    p.intervalStart = cfg.intervalStart;
    p.rho = cfg.rho;
    p.bb = cfg.bb;
    p.originConditioning = cfg.originConditioning;
    p.condEvent = cfg.condEvent;
    p.startAge = drawPrior(cfg.startAgePrior, rng);
    int nT = cfg.numPsiTypes;
    p.psi.assign(nT, std::vector<double>());
    p.psiIdx.assign(nT, std::vector<int>());
    int lPrev = -1, mPrev = -1;
    std::vector<int> pPrev(nT, -1);
    std::map<int,double> lamChunk, muChunk;
    std::vector<std::map<int,double>> psiChunk(nT);
    static const std::vector<int> noGroups;
    static const std::map<int,Probability::PriorSpec> noGroupPrior;
    for(double s : cfg.intervalStart){
        int li = 0, mi = 0;
        for(int k = 0; k < (int)cfg.lambdaTimes.size(); k++) if(cfg.lambdaTimes[k] <= s) li = k;
        for(int k = 0; k < (int)cfg.muTimes.size(); k++)     if(cfg.muTimes[k] <= s)     mi = k;
        if(li > lPrev){ p.lambda.push_back(drawChunk(cfg.lambdaGroups, cfg.lambdaGroupPrior, cfg.lambdaPrior, li, lamChunk, rng)); lPrev = li; }
        if(mi > mPrev){ p.mu.push_back(drawChunk(cfg.muGroups, cfg.muGroupPrior, cfg.muPrior, mi, muChunk, rng)); mPrev = mi; }
        p.lambdaIdx.push_back(li);
        p.muIdx.push_back(mi);
        for(int t = 0; t < nT; t++){
            int pi = 0;
            for(int k = 0; k < (int)cfg.psiTimes[t].size(); k++) if(cfg.psiTimes[t][k] <= s) pi = k;
            if(pi > pPrev[t]){
                const std::vector<int>& g = (t < (int)cfg.psiGroups.size()) ? cfg.psiGroups[t] : noGroups;
                const std::map<int,Probability::PriorSpec>& gp = (t < (int)cfg.psiGroupPriors.size()) ? cfg.psiGroupPriors[t] : noGroupPrior;
                p.psi[t].push_back(drawChunk(g, gp, cfg.psiPriors[t], pi, psiChunk[t], rng));
                pPrev[t] = pi;
            }
            p.psiIdx[t].push_back(pi);
        }
    }
    return p;
}

void Sbc::run(void){
    if(cfg.emitFiles)
        runEmit();
    else if(cfg.simulateOnly)
        runSimulateOnly();
    else
        runInference();
}

void Sbc::runEmit(void){
    ForwardSimulator sim(rng);
    const char* asg = cfg.originConditioning ? "TOTAL" : "CROWN";
    for(int rep = 0; rep < cfg.numReps; rep++){
        SimParams truth = drawParams();
        SimResult r = sim.simulate(truth);
        std::string base = cfg.dumpPrefix + "_rep" + std::to_string(rep);

        std::ofstream tf(base + ".tree");
        tf << r.backboneNewick << "\n";
        tf.close();

        std::ofstream cf(base + ".clades");
        cf << "whole\t";
        for(int i = 1; i <= r.numBackbone; i++)
            cf << (i > 1 ? "," : "") << "T" << i;
        cf << "\n";
        cf.close();

        std::ofstream ff(base + ".fossils");
        for(size_t i = 0; i < r.fossilAges.size(); i++){
            ff << "F" << (i + 1) << '\t' << r.fossilAges[i] << '\t' << r.fossilAges[i] << "\twhole\t" << asg;
            if(cfg.numPsiTypes > 1) ff << '\t' << cfg.psiTypeNames[r.fossilTypes[i]];
            ff << '\n';
        }
        for(int i = 0; i < r.numUE; i++)
            ff << "U" << (i + 1) << "\t0\t0\twhole\t" << asg << '\n';
        ff.close();

        std::ofstream xf(base + ".truth");
        xf << "lambda0\t" << truth.lambda[0] << "\nmu0\t" << truth.mu[0];
        for(int t = 0; t < cfg.numPsiTypes; t++)
            for(size_t b = 0; b < truth.psi[t].size(); b++){
                std::string lab = (cfg.numPsiTypes > 1)
                    ? ("psi_" + cfg.psiTypeNames[t] + (truth.psi[t].size() > 1 ? "_" + std::to_string(b) : ""))
                    : ("psi" + std::to_string(b));
                xf << "\n" << lab << "\t" << truth.psi[t][b];
            }
        xf << "\nx\t" << truth.startAge << '\n';
        xf.close();

        printf("emit rep %d: %d backbone, %zu fossil, %d UE -> %s.{tree,clades,fossils,truth}\n",
               rep, r.numBackbone, r.fossilAges.size(), r.numUE, base.c_str());
    }
}

void Sbc::runSimulateOnly(void){
    ForwardSimulator sim(rng);
    double sExt = 0, sFoss = 0, sBB = 0, sUE = 0, sLam = 0, sMu = 0, sPsi = 0, sX = 0;
    long minE = std::numeric_limits<long>::max(), maxE = 0;
    long minF = std::numeric_limits<long>::max(), maxF = 0;
    long minB = std::numeric_limits<long>::max(), maxB = 0;
    long nLt2 = 0, nEq1 = 0;
    for(int i = 0; i < cfg.numReps; i++){
        SimParams p = drawParams();
        SimResult r = sim.simulate(p);
        sExt += r.numExtantSampled; sFoss += r.numFossils; sBB += r.numBackbone; sUE += r.numUE;
        sLam += p.lambda[0]; sMu += p.mu[0]; sPsi += p.psi[0][0]; sX += p.startAge;
        if(r.numBackbone < minB) minB = r.numBackbone;
        if(r.numBackbone > maxB) maxB = r.numBackbone;
        if(r.numExtantSampled < minE) minE = r.numExtantSampled;
        if(r.numExtantSampled > maxE) maxE = r.numExtantSampled;
        if(r.numFossils < minF) minF = r.numFossils;
        if(r.numFossils > maxF) maxF = r.numFossils;
        if(r.numBackbone < 2) nLt2++;
        if(r.numBackbone == 1) nEq1++;
    }
    double n = cfg.numReps;
    printf("SBC simulate-only: %d replicates, %s conditioning, rho=%.3g, bb=%.3g, %zu interval(s)\n",
           cfg.numReps, cfg.originConditioning ? "origin" : "crown", cfg.rho, cfg.bb, cfg.intervalStart.size());
    printf("  drawn means : lambda0 %.5g  mu0 %.5g  psi0 %.5g  x %.5g\n", sLam/n, sMu/n, sPsi/n, sX/n);
    printf("  #extant(tot): mean %.3f  [%ld..%ld]\n", sExt/n, minE, maxE);
    printf("  #backbone   : mean %.3f  [%ld..%ld]\n", sBB/n, minB, maxB);
    printf("  #UE         : mean %.3f\n", sUE/n);
    printf("  #fossils    : mean %.3f  [%ld..%ld]\n", sFoss/n, minF, maxF);
}

void Sbc::runInference(void){
    ForwardSimulator sim(rng);
    UserSettings& us = UserSettings::userSettings();
    int nRuns = us.getNumRuns();
    int thin = us.getThinning();
    bool autoStop = us.getAutoChainLength();
    unsigned long ncyc = autoStop ? us.getMaxGen() : us.getChainLength();
    int ng = (ncyc > 2000000000UL) ? 2000000000 : (int)ncyc;
    double burnFrac = us.getBurninFraction();
    std::map<std::string, std::vector<double>> ranks;
    std::map<std::string, long> cov50, cov90;
    int nUnconverged = 0;
    int progStep = (cfg.numReps + 9) / 10;
    if(progStep < 1) progStep = 1;
    std::ofstream rankOut;
    std::vector<std::string> outCols;
    bool liveHeader = false;

    for(int rep = 0; rep < cfg.numReps; rep++){
        SimParams truth = drawParams();
        SimResult r = sim.simulate(truth);

        std::string repBase = cfg.dumpPrefix + "_rep" + std::to_string(rep);

        Tree* tree = new Tree(r.backboneNewick);
        std::vector<std::string> taxa;
        for(int i = 1; i <= r.numBackbone; i++)
            taxa.push_back("T" + std::to_string(i));
        Node* crown;
        if(r.numBackbone > 0){
            tree->validateBackbone();
            crown = tree->getMRCA(taxa);
        }else{
            crown = tree->getCrown();
        }
        std::vector<Clade> clades;
        clades.push_back(Clade("whole", taxa, crown, crown->getAncestor()));
        Assignment asg = cfg.originConditioning ? Assignment::TOTAL : Assignment::CROWN;
        std::vector<Fossil> fossils;
        for(size_t i = 0; i < r.fossilAges.size(); i++){
            std::string ty = (cfg.numPsiTypes > 1) ? cfg.psiTypeNames[r.fossilTypes[i]] : "";
            fossils.push_back(Fossil("F" + std::to_string(i + 1), r.fossilAges[i], r.fossilAges[i], "whole", asg, ty));
        }
        for(int i = 0; i < r.numUE; i++)
            fossils.push_back(Fossil("U" + std::to_string(i + 1), 0.0, 0.0, "whole", asg));
        if(r.numBackbone > 0)
            initAges(tree, fossils, cfg.originConditioning, cfg.startAgePrior);

        std::vector<ChainRunner*> chains;
        std::vector<FBDTreeModel*> models;
        for(int m = 0; m < nRuns; m++){
            unsigned int modelSeed = (unsigned int)(rng->uniformRv() * 4294967295.0);
            FBDTreeModel* model = new FBDTreeModel(tree, clades, fossils, modelSeed);
            models.push_back(model);
            chains.push_back(new Mcmc(ng, thin, model));
        }
        ConvergenceRunner cr(chains, repBase + ".log", repBase + ".trees");
        if(cr.run() == false)
            nUnconverged++;
        double repMaxRhat = cr.getMaxRhat(), repMinChainEss = cr.getMinChainEss(), repBulkEss = cr.getMinBulkEss();
        if((rep + 1) % progStep == 0 || rep + 1 == cfg.numReps){
            printf("rep %d/%d\n", rep + 1, cfg.numReps);
            fflush(stdout);
        }

        const std::vector<std::string>& names = chains[0]->traceNames();
        std::map<std::string, double> thisRep;
        for(size_t c = 0; c < names.size(); c++){
            double t = truthForName(names[c], truth, cfg.originConditioning, cfg.psiTypeNames);
            if(std::isnan(t))
                continue;
            std::vector<double> s;
            for(ChainRunner* ch : chains){
                const std::vector<double>& col = ch->traceColumns()[c];
                size_t bIdx = (size_t)(burnFrac * col.size());
                s.insert(s.end(), col.begin() + bIdx, col.end());
            }
            if(s.empty())
                continue;
            std::sort(s.begin(), s.end());
            long below = 0, equal = 0;
            for(double x : s){
                if(x < t) below++;
                else if(x == t) equal++;
            }
            double rank = ((double)below + rng->uniformRv() * (double)equal) / (double)s.size();
            ranks[names[c]].push_back(rank);
            thisRep[names[c]] = rank;
            if(t >= quantile(s, 0.25) && t <= quantile(s, 0.75)) cov50[names[c]]++;
            if(t >= quantile(s, 0.05) && t <= quantile(s, 0.95)) cov90[names[c]]++;
        }
        if(outCols.empty())
            for(size_t c = 0; c < names.size(); c++)
                if(thisRep.count(names[c])) outCols.push_back(names[c]);
        if(cfg.dumpPrefix.empty() == false){
            if(liveHeader == false){
                rankOut.open(cfg.dumpPrefix + "_ranks.tsv");
                for(size_t c = 0; c < outCols.size(); c++)
                    rankOut << (c ? "\t" : "") << outCols[c];
                rankOut << "\tminChainESS\tbulkESS\tmaxRhat\n";
                liveHeader = true;
            }
            for(size_t c = 0; c < outCols.size(); c++){
                std::map<std::string, double>::iterator f = thisRep.find(outCols[c]);
                rankOut << (c ? "\t" : "") << (f != thisRep.end() ? f->second : std::numeric_limits<double>::quiet_NaN());
            }
            rankOut << "\t" << repMinChainEss << "\t" << repBulkEss << "\t" << repMaxRhat << "\n";
            rankOut.flush();
        }
        for(ChainRunner* ch : chains) delete ch;
        for(FBDTreeModel* mo : models) delete mo;
        delete tree;
    }

    printf("SBC inference: %d reps, %d chains/rep, conditioning=%s, rho=%.3g, bb=%.3g\n",
           cfg.numReps, nRuns, cfg.originConditioning ? "origin" : "crown", cfg.rho, cfg.bb);
    if(nUnconverged > 0)
        printf("  WARNING: %d of %d reps did not reach the R-hat/ESS thresholds.\n",
               nUnconverged, cfg.numReps);
    printf("  %-12s %8s %8s %7s %7s\n", "param", "KS_D", "KS_p", "cov50", "cov90");
    for(const std::string& nm : outCols){
        std::vector<double> v = ranks[nm];
        long R = (long)v.size();
        if(R == 0) continue;
        std::sort(v.begin(), v.end());
        double D = 0.0;
        for(long i = 0; i < R; i++){
            double lo = (double)i / R, hi = (double)(i + 1) / R;
            D = std::max(D, std::max(v[i] - lo, hi - v[i]));
        }
        double ksp = ksPvalue(D, R);
        double c50 = (double)cov50[nm] / R;
        double c90 = (double)cov90[nm] / R;
        printf("  %-12s %8.4f %8.4f %7.4f %7.4f\n", nm.c_str(), D, ksp, c50, c90);
    }

    if(cfg.dumpPrefix.empty() == false && liveHeader)
        printf("  ranks saved in %s_ranks.tsv\n", cfg.dumpPrefix.c_str());
}
