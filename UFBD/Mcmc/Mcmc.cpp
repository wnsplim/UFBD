#include "Mcmc.hpp"
#include "Msg.hpp"
#include "PhylogeneticModel.hpp"
#include "SequenceLikelihood.hpp"
#include "RandomVariable.hpp"
#include "Node.hpp"
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
    if(std::isfinite(curLnL) == false || std::isfinite(curLnP) == false){
        std::ostringstream os;
        os << std::setprecision(6) << "the initial state has lnLikelihood " << curLnL
           << " and lnPrior " << curLnP;
        Msg::error(os.str());
    }
    gen = 0;

    UserSettings& settings = UserSettings::userSettings();
    if(settings.clockPresent() && settings.getSigma2Param() == Sigma2Param::PNCP && settings.getPncpTuningGens() > 0){
        model->setChainLabel(runLabel);
        tuning = true;
        advance(settings.getPncpTuningGens());
        model->freezePncpTuning();
        tuning = false;
        gen = 0;
    }
}

void Mcmc::advance(unsigned long nGens) {
    RandomVariable::setActiveInstance(model->getRng());
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    unsigned long target = gen + nGens;
    if(gen == 0 && tuning == false)
        sample(0, curLnL, curLnP);
    while (gen < target) {
        gen++;
        unsigned long n = gen;

        static const bool chkState = (getenv("FBD_CHK_STATE") != nullptr);
        if(chkState && gen > 1){
            static long nBad = 0;
            double rl = model->lnLikelihood(), rp = model->lnPriorProbability();
            if(std::fabs(rl - curLnL) > 1e-6 || std::fabs(rp - curLnP) > 1e-6){
                nBad++;
                if(nBad <= 3)
                    std::cout << "[state] VIOLATION gen " << gen << "  lnL " << curLnL << " vs " << rl
                              << "  lnP " << curLnP << " vs " << rp << "\n" << std::flush;
            }
            if(gen % 20000 == 0)
                std::cout << "[state] gen " << gen << "  violations " << nBad << "\n" << std::flush;
        }

        static const bool chkCache = (getenv("FBD_CHK_CACHE") != nullptr);
        if(chkCache && gen > 1){
            static long nBad = 0;
            static double worst = 0.0;
            double cached = model->lnLikelihood();
            model->invalidateLikelihoodCache();
            double fresh = model->lnLikelihood();
            double e = std::fabs(cached - fresh);
            if(std::isfinite(cached) == false && std::isfinite(fresh) == false)
                e = 0.0;
            if(e > worst) worst = e;
            if(e > 1e-6){
                nBad++;
                if(nBad <= 3)
                    std::cout << "[cache] CORRUPT gen " << gen << "  cached " << cached
                              << "  fresh " << fresh << "  diff " << e << "\n" << std::flush;
            }
            if(gen % 20000 == 0)
                std::cout << "[cache] gen " << gen << "  corrupt " << nBad << " (worst " << worst
                          << ")\n" << std::flush;
        }

        static const bool chkPrior = (getenv("FBD_CHK_PRIOR") != nullptr);
        if(chkPrior && gen > 1){
            static long nBad = 0;
            double cached = model->lnPriorProbability();
            model->invalidatePriorCache();
            double fresh = model->lnPriorProbability();
            double e = std::fabs(cached - fresh);
            if(std::isfinite(cached) == false && std::isfinite(fresh) == false) e = 0.0;
            if(e > 1e-6 || (std::isfinite(cached) != std::isfinite(fresh))){
                nBad++;
                if(nBad <= 3)
                    std::cout << "[prior] STALE gen " << gen << "  cached " << cached
                              << "  fresh " << fresh << "\n" << std::flush;
            }
            if(gen % 20000 == 0)
                std::cout << "[prior] gen " << gen << "  stale " << nBad << "\n" << std::flush;
        }

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

        static const bool trc = (getenv("FBD_TRACE") != nullptr);
        if(trc)
            fprintf(stderr, "[t] %lu mv=%d h=%.17g curL=%.17g curP=%.17g a=%d\n", gen,
                    model->getLastMoveType(), lnProposalRatio, curLnL, curLnP, (int)acceptMove);

        if (verbose && tuning == false && n % thinning == 0) {
            std::ostringstream os;
            os << std::fixed << std::setprecision(2) << "chain " << runLabel << "  " << n << " -- posterior " << (curLnL + curLnP) << " likelihood " << curLnL << " prior " << curLnP << "\n";
            ChainRunner::logLine(os.str());
        }

        if (tuning == false && n % thinning == 0)
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
    if(n == 0){
        params.addFilepath(paramOut, true);
        std::vector<std::string> cn = {"n", "posterior", "likelihood", "prior"};
        std::vector<std::string> headStr = model->getParameterNames();
        cn.insert( cn.end(), headStr.begin(), headStr.end() );
        if(cpuTime)
            cn.push_back("cpu_s");
        params.addColumnNamesTSV(cn);

        std::vector<std::string> latNames = model->getLatentNames();
        if(latNames.empty() == false && UserSettings::userSettings().getWriteLatentLog()){
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
            std::vector<std::string> zleg = model->getZoneLegend();
            if(zleg.empty() == false){
                std::string zp = paramOut;
                size_t zdp = zp.rfind(".log");
                if(zdp != std::string::npos) zp = zp.substr(0, zdp);
                zp += "_zones.tsv";
                std::ofstream zf(zp);
                for(const std::string& line : zleg) zf << line << "\n";
            }
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
        fixedMask.assign(3, false);
        for(bool fb : model->getParameterFixedMask())
            fixedMask.push_back(fb);
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
    static const bool chkCkp = (getenv("FBD_CHK_CKP") != nullptr);
    if(chkCkp)
        checkCheckpointRoundTrip(path);
    return true;
}

void Mcmc::checkCheckpointRoundTrip(const std::string& path) {
    std::ostringstream os;
    os << std::setprecision(17);
    os << gen << ' ' << thinning << '\n';
    os << curLnL << ' ' << curLnP << '\n';
    model->getRng()->writeState(os);
    model->writeState(os);

    std::vector<std::string> onDisk, roundTrip;
    std::ifstream f(path);
    for(std::string t; f >> t; )
        onDisk.push_back(t);
    std::istringstream rt(os.str());
    for(std::string t; rt >> t; )
        roundTrip.push_back(t);

    size_t n = (onDisk.size() < roundTrip.size()) ? onDisk.size() : roundTrip.size();
    size_t i = 0;
    while(i < n && onDisk[i] == roundTrip[i])
        i++;
    if(i == n && onDisk.size() == roundTrip.size()){
        std::cout << "[CHK_CKP] writeState(readState(x)) == x  (" << onDisk.size() << " tokens)\n";
        return;
    }
    std::cout << "[CHK_CKP] STATE LOST ON RESUME: first mismatch at token " << i
              << "  (on-disk " << onDisk.size() << " tokens, round-trip " << roundTrip.size() << ")\n";
    size_t lo = (i > 6) ? i - 6 : 0;
    size_t hi = (i + 7 < n) ? i + 7 : n;
    std::cout << "[CHK_CKP]   on-disk   :";
    for(size_t j = lo; j < hi; j++)
        std::cout << (j == i ? " >>" : " ") << onDisk[j];
    std::cout << "\n[CHK_CKP]   round-trip:";
    for(size_t j = lo; j < hi; j++)
        std::cout << (j == i ? " >>" : " ") << roundTrip[j];
    std::cout << '\n';
}

RandomVariable* Mcmc::resumeRng(void) {
    return model->getRng();
}

std::vector<std::string> Mcmc::resumeParameterNames(void) {
    return model->getParameterNames();
}

std::vector<std::string> Mcmc::resumeLatentNames(void) {
    return model->getLatentNames();
}
