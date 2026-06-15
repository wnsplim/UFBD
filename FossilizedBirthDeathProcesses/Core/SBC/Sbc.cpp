#include "Sbc.hpp"
#include "FBDInput.hpp"
#include "FBDTreeModel.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

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

double truthForName(const std::string& name, const SimParams& p, bool origin){
    auto idxOf = [](const std::string& rest) -> int { return rest.empty() ? 0 : std::stoi(rest); };
    if(name.rfind("lambda", 0) == 0)  return p.lambda[idxOf(name.substr(6))];
    if(name.rfind("psi", 0) == 0)     return p.psi[idxOf(name.substr(3))];
    if(name.rfind("mu", 0) == 0)      return p.mu[idxOf(name.substr(2))];
    if(name == "originAge" && origin) return p.startAge;
    return std::numeric_limits<double>::quiet_NaN();
}

void initAges(Tree* tree, std::vector<Fossil>& fossils, bool origin, const Probability::PriorSpec& xp){
    double maxFoss = 0.0;
    for(Fossil& f : fossils)
        if(f.getMaxAge() > maxFoss) maxFoss = f.getMaxAge();
    int numInternal = tree->getNumNodes() - tree->getNumTaxa();
    double unit = (maxFoss > 0.0) ? maxFoss / (numInternal + 1) : 1.0;
    std::map<Node*,double> minAges;
    if(maxFoss > 0.0)
        minAges[tree->getCrown()] = maxFoss * 1.05;
    tree->assignStartingAges(minAges, unit);

    double mean = Probability::priorMean(xp.family, xp.p1, xp.p2);
    double target = origin ? 0.9 * mean : mean;
    if(origin == false && target < maxFoss)
        target = maxFoss * 1.05;
    double cur = tree->getCrown()->getTime();
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
    p.originConditioning = cfg.originConditioning;
    p.startAge = drawPrior(cfg.startAgePrior, rng);
    for(size_t i = 0; i < cfg.intervalStart.size(); i++){
        p.lambda.push_back(drawPrior(cfg.lambdaPrior, rng));
        p.mu.push_back(drawPrior(cfg.muPrior, rng));
        p.psi.push_back(drawPrior(cfg.psiPrior, rng));
    }
    return p;
}

void Sbc::run(void){
    if(cfg.simulateOnly)
        runSimulateOnly();
    else
        runInference();
}

void Sbc::runSimulateOnly(void){
    ForwardSimulator sim(rng);
    double sExt = 0, sFoss = 0, sLam = 0, sMu = 0, sPsi = 0, sX = 0;
    long minE = std::numeric_limits<long>::max(), maxE = 0;
    long minF = std::numeric_limits<long>::max(), maxF = 0;
    long nLt2 = 0, nEq1 = 0;
    for(int i = 0; i < cfg.numReps; i++){
        SimParams p = drawParams();
        SimResult r = sim.simulate(p);
        sExt += r.numExtantSampled; sFoss += r.numFossils;
        sLam += p.lambda[0]; sMu += p.mu[0]; sPsi += p.psi[0]; sX += p.startAge;
        if(r.numExtantSampled < minE) minE = r.numExtantSampled;
        if(r.numExtantSampled > maxE) maxE = r.numExtantSampled;
        if(r.numFossils < minF) minF = r.numFossils;
        if(r.numFossils > maxF) maxF = r.numFossils;
        if(r.numExtantSampled < 2) nLt2++;
        if(r.numExtantSampled == 1) nEq1++;
    }
    double n = cfg.numReps;
    printf("SBC simulate-only: %d replicates, %s conditioning, rho=%.3g, %zu interval(s)\n",
           cfg.numReps, cfg.originConditioning ? "origin" : "crown", cfg.rho, cfg.intervalStart.size());
    printf("  drawn means : lambda0 %.5g  mu0 %.5g  psi0 %.5g  x %.5g\n", sLam/n, sMu/n, sPsi/n, sX/n);
    printf("  #extant     : mean %.3f  [%ld..%ld]\n", sExt/n, minE, maxE);
    printf("  #fossils    : mean %.3f  [%ld..%ld]\n", sFoss/n, minF, maxF);
    printf("  P(nExt<2)=%.4f  P(nExt=1)=%.4f\n", nLt2/n, nEq1/n);
}

void Sbc::runInference(void){
    ForwardSimulator sim(rng);
    long burnin = (long)(cfg.burninFraction * cfg.mcmcGen);
    std::map<std::string, std::vector<double>> ranks;
    std::map<std::string, long> cov50, cov90;
    std::vector<int> repNExt, repNFoss;

    for(int rep = 0; rep < cfg.numReps; rep++){
        SimParams truth = drawParams();
        SimResult r = sim.simulate(truth);
        repNExt.push_back(r.numExtantSampled);
        repNFoss.push_back(r.numFossils);

        Tree* tree = new Tree(r.backboneNewick);
        tree->validateBackbone();
        std::vector<std::string> taxa;
        for(int i = 1; i <= r.numExtantSampled; i++)
            taxa.push_back("T" + std::to_string(i));
        Node* crown = tree->getMRCA(taxa);
        std::vector<Clade> clades;
        clades.push_back(Clade("whole", taxa, crown, crown->getAncestor()));
        Assignment asg = cfg.originConditioning ? Assignment::TOTAL : Assignment::CROWN;
        std::vector<Fossil> fossils;
        for(size_t i = 0; i < r.fossilAges.size(); i++)
            fossils.push_back(Fossil("F" + std::to_string(i + 1), r.fossilAges[i], r.fossilAges[i], "whole", asg));
        initAges(tree, fossils, cfg.originConditioning, cfg.startAgePrior);

        unsigned int modelSeed = (unsigned int)(rng->uniformRv() * 4294967295.0);
        FBDTreeModel model(tree, clades, fossils, modelSeed);

        RandomVariable::setActiveInstance(model.getRng());
        RandomVariable& mrng = RandomVariable::randomVariableInstance();

        std::vector<std::string> names = model.getParameterNames();
        std::vector<std::vector<double>> samp(names.size());

        double curLnL = model.lnLikelihood();
        double curLnP = model.lnPriorProbability();
        for(long n = 1; n <= cfg.mcmcGen; n++){
            double lpr = model.update();
            double newLnL = model.lnLikelihood();
            double newLnP = model.lnPriorProbability();
            double lnR = (newLnL - curLnL) + (newLnP - curLnP) + lpr;
            if(std::log(mrng.uniformRv()) < lnR){
                curLnL = newLnL; curLnP = newLnP;
                model.updateForAcceptance();
            }else{
                model.updateForRejection();
            }
            if(n > burnin && (n % cfg.mcmcThin == 0)){
                std::vector<double> v = model.getParameterString();
                for(size_t c = 0; c < names.size(); c++)
                    samp[c].push_back(v[c]);
            }
        }

        if(samp[0].empty())
            Msg::error("SBC: no post-burnin samples collected; mcmcGen too small relative to burnin/thin");
        for(size_t c = 0; c < names.size(); c++){
            double t = truthForName(names[c], truth, cfg.originConditioning);
            if(std::isnan(t))
                continue;
            std::vector<double>& s = samp[c];
            std::sort(s.begin(), s.end());
            long below = 0;
            for(double x : s)
                if(x < t) below++;
            ranks[names[c]].push_back((double)below / (double)s.size());
            if(t >= quantile(s, 0.25) && t <= quantile(s, 0.75)) cov50[names[c]]++;
            if(t >= quantile(s, 0.05) && t <= quantile(s, 0.95)) cov90[names[c]]++;
        }
    }

    printf("SBC inference: %d reps, %s, rho=%.3g, %ld gen/rep (burnin %.2f, thin %d)\n",
           cfg.numReps, cfg.originConditioning ? "origin" : "crown", cfg.rho,
           cfg.mcmcGen, cfg.burninFraction, cfg.mcmcThin);
    printf("  %-12s %8s %8s %10s %7s %7s %7s\n", "param", "KS_D", "KS_p", "chi2", "chi2_p", "cov50", "cov90");
    for(std::map<std::string, std::vector<double>>::iterator it = ranks.begin(); it != ranks.end(); ++it){
        std::vector<double> v = it->second;
        long R = (long)v.size();
        std::sort(v.begin(), v.end());
        double D = 0.0;
        for(long i = 0; i < R; i++){
            double lo = (double)i / R, hi = (double)(i + 1) / R;
            D = std::max(D, std::max(v[i] - lo, hi - v[i]));
        }
        double ksp = ksPvalue(D, R);
        std::vector<long> bin(cfg.rankBins, 0);
        for(double x : v){
            int b = (int)(x * cfg.rankBins);
            if(b >= cfg.rankBins) b = cfg.rankBins - 1;
            if(b < 0) b = 0;
            bin[b]++;
        }
        double expected = (double)R / cfg.rankBins;
        double chi2 = 0.0;
        for(int b = 0; b < cfg.rankBins; b++){
            double d = bin[b] - expected;
            chi2 += d * d / expected;
        }
        double chip = 1.0 - Probability::ChiSquare::cdf(cfg.rankBins - 1, chi2);
        double c50 = (double)cov50[it->first] / R;
        double c90 = (double)cov90[it->first] / R;
        printf("  %-12s %8.4f %8.4f %10.3f %7.4f %7.4f %7.4f\n", it->first.c_str(), D, ksp, chi2, chip, c50, c90);
    }

    if(cfg.dumpPrefix.empty() == false){
        std::ofstream out(cfg.dumpPrefix + "_ranks.tsv");
        std::vector<std::string> cols;
        for(std::map<std::string, std::vector<double>>::iterator it = ranks.begin(); it != ranks.end(); ++it)
            cols.push_back(it->first);
        for(size_t c = 0; c < cols.size(); c++)
            out << (c ? "\t" : "") << cols[c];
        out << "\n";
        for(int rep = 0; rep < cfg.numReps; rep++){
            for(size_t c = 0; c < cols.size(); c++)
                out << (c ? "\t" : "") << ranks[cols[c]][rep];
            out << "\n";
        }
        printf("  wrote ranks to %s_ranks.tsv\n", cfg.dumpPrefix.c_str());
    }
}
