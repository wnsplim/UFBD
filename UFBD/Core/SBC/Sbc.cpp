#include "Sbc.hpp"
#include "FBDInput.hpp"
#include "FBDTreeModel.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterOUField.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"
#include "ChainRunner.hpp"
#include "Mcmc.hpp"
#include "ConvergenceRunner.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

double drawPrior(const Probability::PriorSpec& s, RandomVariable* rng){
    return Probability::priorRv(s.family, s.p1, s.p2, rng);
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

std::vector<double> drawOU(const std::vector<double>& times, const OUHyperSpec& ou, std::vector<double>& hyper, RandomVariable* rng){
    ParameterOUField f(1.0, nullptr, (int)times.size(), times, 0.0, nullptr, ou);
    f.drawFromPrior(rng);
    hyper = { f.getTheta(), f.getSdEq(), f.getNu() };
    std::vector<double> r;
    for(int i = 0; i < f.getNumBins(); i++)
        r.push_back(f.getRate(i));
    return r;
}

bool isRankedName(const std::string& nm){
    return nm != "originAge" && (nm[0] != 'x' || std::isdigit((unsigned char)nm[1]) == 0);
}

void tipNames(Node* n, std::vector<std::string>& out){
    if(n->getIsTip()){
        out.push_back(n->getName());
        return;
    }
    for(Node* d : n->getDescendants())
        tipNames(d, out);
}

const char* assignmentName(Assignment a){
    return (a == Assignment::CROWN) ? "CROWN" : (a == Assignment::STEM) ? "STEM" : "TOTAL";
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

double lnBinomTail(long n, double p, long k, bool upper){
    if(upper)
        return (k <= 0) ? 0.0 : Probability::Helper::lnIncompleteBeta((double)k, (double)(n - k + 1), p);
    return (k >= n) ? 0.0 : Probability::Helper::lnIncompleteBeta((double)(n - k), (double)(k + 1), 1.0 - p);
}

double ecdfLnGamma(std::vector<double> u, const std::vector<double>& z){
    long n = (long)u.size();
    std::sort(u.begin(), u.end());
    double g = 0.0;
    for(double zj : z){
        long k = (long)(std::upper_bound(u.begin(), u.end(), zj) - u.begin());
        double lo = lnBinomTail(n, zj, k, false);
        double hi = lnBinomTail(n, zj, k, true);
        double t = std::log(2.0) + std::min(lo, hi);
        if(t < g) g = t;
    }
    return g;
}

std::vector<double> nullLnGammas(long n, const std::vector<double>& z, int nNull, RandomVariable* rng){
    std::vector<double> g(nNull);
    std::vector<double> u((size_t)n);
    for(int b = 0; b < nNull; b++){
        for(long i = 0; i < n; i++) u[i] = rng->uniformRv();
        g[b] = ecdfLnGamma(u, z);
    }
    std::sort(g.begin(), g.end());
    return g;
}

double quantile(const std::vector<double>& s, double q){
    double idx = q * (s.size() - 1);
    size_t lo = (size_t)idx;
    double frac = idx - lo;
    if(lo + 1 < s.size())
        return s[lo] * (1.0 - frac) + s[lo + 1] * frac;
    return s[lo];
}

double lnQuantile(const std::vector<double>& s, double q){
    double idx = q * (s.size() - 1);
    size_t lo = (size_t)idx;
    double frac = idx - lo;
    if(lo + 1 < s.size())
        return std::log(std::exp(s[lo]) * (1.0 - frac) + std::exp(s[lo + 1]) * frac);
    return s[lo];
}

std::vector<double> ecdfGrid(long n){
    int nGrid = (int)std::min(n, (long)200);
    std::vector<double> z;
    for(int j = 1; j < nGrid; j++)
        z.push_back((double)j / (double)nGrid);
    return z;
}

void printRankSummary(const std::vector<std::string>& names, std::map<std::string, std::vector<double>>& ranks,
                      const SbcConfig& cfg, RandomVariable* rng,
                      std::map<std::string, long>* cov50, std::map<std::string, long>* cov90){
    struct Row { std::string nm; long n; double lg, lg5, lgBH, lgBY, lgp, p; };
    std::vector<Row> rows;
    for(const std::string& nm : names)
        if((long)ranks[nm].size() >= 2)
            rows.push_back({ nm, (long)ranks[nm].size(), 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 });
    long m = (long)rows.size();
    double harmonic = 0.0;
    for(long j = 1; j <= m; j++) harmonic += 1.0 / (double)j;
    double levelBH = cfg.fdrQ / (double)m, levelBY = cfg.fdrQ / ((double)m * harmonic);
    double smallest = (cfg.level > 0.0) ? std::min(levelBY, cfg.level) : levelBY;
    int nNull = (cfg.numNull > 0) ? cfg.numNull : std::max(4000, (int)std::ceil(20.0 / smallest));
    std::map<long, std::vector<double>> nullByN;
    for(Row& r : rows){
        std::vector<double> z = ecdfGrid(r.n);
        if(nullByN.count(r.n) == 0)
            nullByN[r.n] = nullLnGammas(r.n, z, nNull, rng);
        const std::vector<double>& nul = nullByN[r.n];
        r.lg = ecdfLnGamma(ranks[r.nm], z);
        long k = (long)(std::upper_bound(nul.begin(), nul.end(), r.lg) - nul.begin());
        r.lg5 = lnQuantile(nul, 0.05);
        r.lgBH = lnQuantile(nul, levelBH);
        r.lgBY = lnQuantile(nul, levelBY);
        if(cfg.level > 0.0)
            r.lgp = lnQuantile(nul, cfg.level);
        r.p = (double)(k + 1) / (double)(nNull + 1);
    }
    std::vector<size_t> order(m);
    for(long i = 0; i < m; i++) order[i] = (size_t)i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b){ return rows[a].p < rows[b].p; });
    std::vector<double> pBonf(m), pBH(m), pBY(m);
    double runBH = 1.0, runBY = 1.0;
    for(long i = m - 1; i >= 0; i--){
        double p = rows[order[i]].p;
        pBonf[order[i]] = std::min(1.0, p * (double)m);
        runBH = std::min(runBH, p * (double)m / (double)(i + 1));
        runBY = std::min(runBY, p * (double)m * harmonic / (double)(i + 1));
        pBH[order[i]] = runBH;
        pBY[order[i]] = runBY;
    }
    printf("  q = %g, m = %ld: BH p = %.4g, BY p = %.4g", cfg.fdrQ, m, levelBH, levelBY);
    if(cfg.level > 0.0) printf(", fixed p = %.4g", cfg.level);
    printf("; %d null rank sets\n", nNull);
    printf("  %-14s %6s %14s %14s %14s", "quantity", "n", "log(g/g_5)", "log(g/g_BH)", "log(g/g_BY)");
    if(cfg.level > 0.0) printf(" %14s", "log(g/g_p)");
    printf(" %9s %9s %9s %9s", "p", "p_Bonf", "p_BH", "p_BY");
    if(cov50 != nullptr) printf(" %7s %7s", "cov50", "cov90");
    printf("\n");
    for(long i = 0; i < m; i++){
        const Row& r = rows[i];
        printf("  %-14s %6ld %14.3f %14.3f %14.3f", r.nm.c_str(), r.n, r.lg - r.lg5, r.lg - r.lgBH, r.lg - r.lgBY);
        if(cfg.level > 0.0) printf(" %14.3f", r.lg - r.lgp);
        printf(" %9.4g %9.4g %9.4g %9.4g", r.p, pBonf[i], pBH[i], pBY[i]);
        if(cov50 != nullptr && cov50->count(r.nm) > 0)
            printf(" %7.4f %7.4f", (double)(*cov50)[r.nm] / r.n, (double)(*cov90)[r.nm] / r.n);
        else if(cov50 != nullptr)
            printf(" %7s %7s", "-", "-");
        printf("\n");
    }
}

void addNodeQuantities(const std::vector<std::vector<std::pair<std::string, double>>>& perRep,
                       std::map<std::string, std::vector<double>>& ranks, std::vector<std::string>& names){
    RandomVariable pick(1u);
    RandomVariable* rng = &pick;
    for(const char* nm : { "x0", "x1", "random_node" })
        names.push_back(nm);
    for(const std::vector<std::pair<std::string, double>>& rep : perRep){
        std::vector<double> internal;
        for(const std::pair<std::string, double>& nr : rep){
            if(nr.first == "x0") ranks["x0"].push_back(nr.second);
            else internal.push_back(nr.second);
            if(nr.first == "x1") ranks["x1"].push_back(nr.second);
        }
        if(internal.empty() == false)
            ranks["random_node"].push_back(internal[(size_t)(rng->uniformRv() * (double)internal.size())]);
    }
}

}

SimParams Sbc::drawParams(void){
    SimParams p;
    p.intervalStart = cfg.intervalStart;
    p.rho = cfg.rho;
    p.ue = cfg.ue;
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
    std::vector<double> lamOU, muOU;
    std::vector<std::vector<double>> psiOU(nT);
    p.psiOU.assign(nT, std::vector<double>());
    if(cfg.lambdaMode == RateMode::OU) lamOU = drawOU(cfg.lambdaTimes, cfg.lambdaOU, p.lambdaOU, rng);
    if(cfg.muMode == RateMode::OU)     muOU = drawOU(cfg.muTimes, cfg.muOU, p.muOU, rng);
    for(int t = 0; t < nT; t++)
        if(cfg.psiModes[t] == RateMode::OU) psiOU[t] = drawOU(cfg.psiTimes[t], cfg.psiOU[t], p.psiOU[t], rng);
    for(double s : cfg.intervalStart){
        int li = 0, mi = 0;
        for(int k = 0; k < (int)cfg.lambdaTimes.size(); k++) if(cfg.lambdaTimes[k] <= s) li = k;
        for(int k = 0; k < (int)cfg.muTimes.size(); k++)     if(cfg.muTimes[k] <= s)     mi = k;
        if(li > lPrev){ p.lambda.push_back(cfg.lambdaMode == RateMode::OU ? lamOU[li] : drawChunk(cfg.lambdaGroups, cfg.lambdaGroupPrior, cfg.lambdaPrior, li, lamChunk, rng)); lPrev = li; }
        if(mi > mPrev){ p.mu.push_back(cfg.muMode == RateMode::OU ? muOU[mi] : drawChunk(cfg.muGroups, cfg.muGroupPrior, cfg.muPrior, mi, muChunk, rng)); mPrev = mi; }
        p.lambdaIdx.push_back(li);
        p.muIdx.push_back(mi);
        for(int t = 0; t < nT; t++){
            int pi = 0;
            for(int k = 0; k < (int)cfg.psiTimes[t].size(); k++) if(cfg.psiTimes[t][k] <= s) pi = k;
            if(pi > pPrev[t]){
                const std::vector<int>& g = (t < (int)cfg.psiGroups.size()) ? cfg.psiGroups[t] : noGroups;
                const std::map<int,Probability::PriorSpec>& gp = (t < (int)cfg.psiGroupPriors.size()) ? cfg.psiGroupPriors[t] : noGroupPrior;
                p.psi[t].push_back(cfg.psiModes[t] == RateMode::OU ? psiOU[t][pi] : drawChunk(g, gp, cfg.psiPriors[t], pi, psiChunk[t], rng));
                pPrev[t] = pi;
            }
            p.psiIdx[t].push_back(pi);
        }
    }
    return p;
}

void Sbc::buildInputs(const SimResult& r, Tree* tree, std::vector<Clade>& clades, std::vector<Fossil>& fossils, std::vector<int>& zoneNode){
    std::vector<std::string> taxa;
    for(int i = 1; i <= r.numBackbone; i++)
        taxa.push_back("T" + std::to_string(i));
    Node* crown = tree->getCrown();
    if(r.numBackbone > 0){
        tree->validateBackbone();
        crown = tree->getMRCA(taxa);
    }
    clades.push_back(Clade("whole", taxa, crown, crown->getAncestor()));
    std::map<Node*, std::string> cladeOf;
    cladeOf[crown] = "whole";
    Assignment asg = cfg.originConditioning ? Assignment::TOTAL : Assignment::CROWN;
    for(size_t i = 0; i < r.fossilAges.size(); i++){
        std::string ty = (cfg.numPsiTypes > 1) ? cfg.psiTypeNames[r.fossilTypes[i]] : "";
        double lo = r.fossilAges[i], hi = r.fossilAges[i];
        if(cfg.ageRange > 0.0){
            hi = r.fossilAges[i] + rng->uniformRv() * cfg.ageRange;
            lo = std::max(0.0, hi - cfg.ageRange);
        }
        Node* x = crown;
        Assignment a = asg;
        if(cfg.zone >= 0.0){
            std::vector<std::string> et;
            for(int l : r.fossilEdgeTaxa[i])
                et.push_back("T" + std::to_string(l));
            Node* u = tree->getMRCA(et);
            std::vector<Node*> path;
            for(Node* p = u; ; p = p->getAncestor()){
                path.push_back(p);
                if(p == crown) break;
            }
            for(int j = (int)path.size() - 1; j >= 0; j--){
                x = path[j];
                if(x->getIsTip() || rng->uniformRv() < cfg.zone){ a = Assignment::TOTAL; break; }
                if(x == u){ a = Assignment::STEM; break; }
            }
            if(cladeOf.count(x) == 0){
                std::vector<std::string> sub;
                tipNames(x, sub);
                cladeOf[x] = "n" + std::to_string(x->getOffset());
                clades.push_back(Clade(cladeOf[x], sub, x, x->getAncestor()));
            }
            zoneNode.push_back(u->getOffset());
        }
        fossils.push_back(Fossil("F" + std::to_string(i + 1), lo, hi, cladeOf[x], a, ty));
    }
    for(int i = 0; i < r.numUE; i++)
        fossils.push_back(Fossil("U" + std::to_string(i + 1), 0.0, 0.0, "whole", asg));
}

void Sbc::run(void){
    if(cfg.summarize)
        runSummarize();
    else if(cfg.emitFiles)
        runEmit();
    else if(cfg.simulateOnly)
        runSimulateOnly();
    else
        runInference();
}

void Sbc::runSummarize(void){
    for(size_t f = 0; f < cfg.ranksPaths.size(); f++){
        std::map<std::string, std::vector<double>> ranks;
        std::vector<std::string> names;
        const std::string& path = cfg.ranksPaths[f];
        std::string dir = path.substr(0, path.find_last_of('/'));
        size_t cut = dir.find_last_of('/');
        cut = (cut == std::string::npos || cut == 0) ? std::string::npos : dir.find_last_of('/', cut - 1);
        std::string label = (cut == std::string::npos) ? dir : dir.substr(cut + 1);
        std::ifstream rf(path);
        if(rf.is_open() == false)
            Msg::error("cannot open " + path);
        std::string line, cell;
        std::getline(rf, line);
        std::vector<std::string> header;
        std::stringstream hs(line);
        while(std::getline(hs, cell, '\t'))
            header.push_back(cell);
        std::vector<bool> use(header.size());
        for(size_t c = 0; c < header.size(); c++){
            const std::string& h = header[c];
            use[c] = h != "bulkESS" && h != "tailESS" && h != "maxRhat" && h != "simAttempts" && h.rfind("t_", 0) != 0;
            if(use[c])
                names.push_back(h);
        }
        while(std::getline(rf, line)){
            std::stringstream ls(line);
            for(size_t c = 0; c < header.size() && std::getline(ls, cell, '\t'); c++){
                double v = std::stod(cell);
                if(use[c] && std::isnan(v) == false)
                    ranks[header[c]].push_back(v);
            }
        }
        if(f < cfg.nodeRanksPaths.size()){
            std::ifstream nf(cfg.nodeRanksPaths[f]);
            if(nf.is_open() == false)
                Msg::error("cannot open " + cfg.nodeRanksPaths[f]);
            std::getline(nf, line);
            bool chunked = line.rfind("chunk", 0) == 0;
            std::vector<std::vector<std::pair<std::string, double>>> perRep;
            std::string lastKey;
            while(std::getline(nf, line)){
                std::stringstream ls(line);
                std::string chunk, rep, node, rank;
                if(chunked)
                    std::getline(ls, chunk, '\t');
                std::getline(ls, rep, '\t');
                std::getline(ls, node, '\t');
                std::getline(ls, rank, '\t');
                std::string key = chunk + ":" + rep;
                if(perRep.empty() || key != lastKey)
                    perRep.push_back({});
                lastKey = key;
                perRep.back().push_back({ node, std::stod(rank) });
            }
            addNodeQuantities(perRep, ranks, names);
        }
        printf("%s\n", label.c_str());
        RandomVariable fileRng(cfg.seed);
        printRankSummary(names, ranks, cfg, &fileRng, nullptr, nullptr);
    }
}

void Sbc::runEmit(void){
    ForwardSimulator sim(rng);
    for(int rep = 0; rep < cfg.numReps; rep++){
        SimParams truth = drawParams();
        SimResult r = sim.simulate(truth);
        std::string base = cfg.dumpPrefix + "_rep" + std::to_string(rep);

        std::ofstream tf(base + ".tree");
        tf << r.backboneNewick << "\n";
        tf.close();

        Tree tree(r.backboneNewick);
        std::vector<Clade> clades;
        std::vector<Fossil> fossils;
        std::vector<int> zoneNode;
        buildInputs(r, &tree, clades, fossils, zoneNode);

        std::ofstream cf(base + ".clades");
        for(Clade& c : clades){
            cf << c.getName() << '\t';
            for(size_t i = 0; i < c.getTaxa().size(); i++)
                cf << (i > 0 ? "," : "") << c.getTaxa()[i];
            cf << "\n";
        }
        cf.close();

        std::ofstream ff(base + ".fossils");
        for(Fossil& f : fossils){
            ff << f.getTaxon() << '\t' << f.getMinAge() << '\t' << f.getMaxAge() << '\t' << f.getClade() << '\t' << assignmentName(f.getAssignment());
            if(f.getType().empty() == false) ff << '\t' << f.getType();
            ff << '\n';
        }
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
    printf("SBC simulate-only: %d replicates, %s conditioning, rho=%.3g, ue=%.3g, %zu interval(s)\n",
           cfg.numReps, cfg.originConditioning ? "origin" : "crown", cfg.rho, cfg.ue, cfg.intervalStart.size());
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
    std::vector<std::vector<std::pair<std::string, double>>> nodeRanksPerRep;
    int nUnconverged = 0;
    int progStep = (cfg.numReps + 9) / 10;
    if(progStep < 1) progStep = 1;
    std::ofstream rankOut, nodeOut;
    std::vector<std::string> outCols;
    bool liveHeader = false;

    for(int rep = 0; rep < cfg.numReps; rep++){
        SimParams truth = drawParams();
        SimResult r = sim.simulate(truth);

        std::string repBase = cfg.dumpPrefix + "_rep" + std::to_string(rep);

        Tree* tree = new Tree(r.backboneNewick);
        std::vector<Clade> clades;
        std::vector<Fossil> fossils;
        std::vector<int> zoneNode;
        buildInputs(r, tree, clades, fossils, zoneNode);
        if(r.numBackbone > 0)
            initAges(tree, fossils, cfg.originConditioning, cfg.startAgePrior);

        std::vector<double> yTrue = r.fossilAges;
        yTrue.resize(fossils.size(), 0.0);
        std::vector<double> zTrue;
        std::vector<char> saTrue;
        for(size_t i = 0; i < r.fossilAttachAge.size(); i++){
            zTrue.push_back(r.fossilAttachAge[i]);
            saTrue.push_back(r.fossilIsSA[i]);
        }
        for(size_t i = 0; i < r.ueAttachAge.size(); i++){
            zTrue.push_back(r.ueAttachAge[i]);
            saTrue.push_back(0);
        }
        double truthLnL = 0.0;
        double truthLnPrior = std::numeric_limits<double>::quiet_NaN();
        rng->uniformRv(); // keeps each seed's replicate datasets identical to runs that built a separate truth model
        std::vector<ChainRunner*> chains;
        std::vector<FBDTreeModel*> models;
        for(int m = 0; m < nRuns; m++){
            unsigned int modelSeed = (unsigned int)(rng->uniformRv() * 4294967295.0);
            FBDTreeModel* model = new FBDTreeModel(tree, clades, fossils, modelSeed);
            truthLnL = model->evaluateAtGeneratingState(truth.lambda, truth.mu, truth.psi,
                                                        truth.lambdaOU, truth.muOU, truth.psiOU, truth.startAge,
                                                        r.backboneNodeAges, yTrue, zTrue, saTrue, zoneNode, truthLnPrior);
            models.push_back(model);
            chains.push_back(new Mcmc(ng, thin, model));
            chains.back()->setWriteCheckpoints(false);
        }
        std::map<std::string, double> truthVal;
        std::vector<std::string> truthNames = models[0]->getParameterNames();
        std::vector<double> truthRow = models[0]->getParameterString();
        for(size_t i = 0; i < truthNames.size(); i++)
            truthVal[truthNames[i]] = truthRow[i];
        ConvergenceRunner cr(chains, repBase + ".log", repBase + ".trees");
        if(cr.run() == false)
            nUnconverged++;
        double repMaxRhat = cr.getMaxRhat(), repBulkEss = cr.getMinBulkEss(), repTailEss = cr.getMinTailEss();
        if((rep + 1) % progStep == 0 || rep + 1 == cfg.numReps){
            printf("rep %d/%d\n", rep + 1, cfg.numReps);
            fflush(stdout);
        }

        const std::vector<std::string>& names = chains[0]->traceNames();
        std::vector<std::vector<double>> post(names.size());
        for(size_t c = 0; c < names.size(); c++)
            for(ChainRunner* ch : chains){
                const std::vector<double>& col = ch->traceColumns()[c];
                size_t bIdx = std::max((size_t)1, (size_t)(burnFrac * col.size()));
                post[c].insert(post[c].end(), col.begin() + bIdx, col.end());
            }

        std::map<std::string, double> thisRep;
        auto rankOf = [&](const std::vector<double>& s, double t){
            long below = 0, equal = 0;
            for(double x : s){
                if(x < t) below++;
                else if(x == t) equal++;
            }
            return ((double)below + rng->uniformRv() * (double)equal) / (double)s.size();
        };
        auto addQuantity = [&](const std::string& nm, std::vector<double> s, double t){
            if(std::isnan(t))
                return;
            std::sort(s.begin(), s.end());
            double rank = rankOf(s, t);
            ranks[nm].push_back(rank);
            thisRep[nm] = rank;
            if(t >= quantile(s, 0.25) && t <= quantile(s, 0.75)) cov50[nm]++;
            if(t >= quantile(s, 0.05) && t <= quantile(s, 0.95)) cov90[nm]++;
            if(rep == 0) outCols.push_back(nm);
        };

        int iLik = -1, iPri = -1;
        std::map<std::string, int> colOf;
        for(size_t c = 0; c < names.size(); c++){
            if(names[c] == "likelihood") iLik = (int)c;
            if(names[c] == "prior")      iPri = (int)c;
            colOf[names[c]] = (int)c;
        }
        for(size_t c = 0; c < names.size(); c++){
            std::map<std::string, double>::iterator t = truthVal.find(names[c]);
            if(t != truthVal.end() && isRankedName(names[c]))
                addQuantity(names[c], post[c], t->second);
        }
        addQuantity("logLik", post[iLik], truthLnL);
        addQuantity("logPrior", post[iPri], truthLnPrior);
        std::vector<std::pair<std::string, double>> nodeRank;
        if(truthVal.count("originAge") > 0)
            nodeRank.push_back({ "x0", rankOf(post[colOf["originAge"]], truthVal["originAge"]) });
        for(int i = 1; truthVal.count("x" + std::to_string(i)) > 0; i++){
            std::string xi = "x" + std::to_string(i);
            nodeRank.push_back({ xi, rankOf(post[colOf[xi]], truthVal[xi]) });
        }
        nodeRanksPerRep.push_back(nodeRank);
        if(cfg.dumpPrefix.empty() == false){
            if(liveHeader == false){
                rankOut.open(cfg.dumpPrefix + "_ranks.tsv");
                rankOut.precision(17);
                for(size_t c = 0; c < outCols.size(); c++)
                    rankOut << (c ? "\t" : "") << outCols[c];
                rankOut << "\tbulkESS\ttailESS\tmaxRhat\tt_lambda\tt_mu\tt_psi\tt_logLik\tt_logPrior\tsimAttempts\n";
                nodeOut.open(cfg.dumpPrefix + "_noderanks.tsv");
                nodeOut.precision(17);
                nodeOut << "rep\tnode\trank\n";
                liveHeader = true;
            }
            for(const std::pair<std::string, double>& nr : nodeRank)
                nodeOut << rep << "\t" << nr.first << "\t" << nr.second << "\n";
            nodeOut.flush();
            for(size_t c = 0; c < outCols.size(); c++){
                std::map<std::string, double>::iterator f = thisRep.find(outCols[c]);
                rankOut << (c ? "\t" : "") << (f != thisRep.end() ? f->second : std::numeric_limits<double>::quiet_NaN());
            }
            rankOut << "\t" << repBulkEss << "\t" << repTailEss << "\t" << repMaxRhat
                    << "\t" << truth.lambda[0] << "\t" << truth.mu[0] << "\t" << truth.psi[0][0]
                    << "\t" << truthLnL << "\t" << truthLnPrior << "\t" << r.attempts << "\n";
            rankOut.flush();
        }
        for(ChainRunner* ch : chains) delete ch;
        for(FBDTreeModel* mo : models) delete mo;
        delete tree;
    }

    printf("SBC inference: %d reps, %d chains/rep, conditioning=%s, rho=%.3g, ue=%.3g\n",
           cfg.numReps, nRuns, cfg.originConditioning ? "origin" : "crown", cfg.rho, cfg.ue);
    if(nUnconverged > 0)
        printf("  WARNING: %d of %d reps did not reach the R-hat/ESS thresholds.\n",
               nUnconverged, cfg.numReps);
    std::vector<std::string> names = outCols;
    addNodeQuantities(nodeRanksPerRep, ranks, names);
    printRankSummary(names, ranks, cfg, rng, &cov50, &cov90);

    if(cfg.dumpPrefix.empty() == false && liveHeader)
        printf("  ranks saved in %s_ranks.tsv\n", cfg.dumpPrefix.c_str());
}
