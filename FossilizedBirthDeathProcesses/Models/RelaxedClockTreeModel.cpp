#include "RelaxedClockTreeModel.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "ApproxBranchLengthLikelihood.hpp"
#include "FBDTreeModel.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

RelaxedClockTreeModel::RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& hessianFile, const std::string& mlTreeFile, int nStates, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param, unsigned int seed){
    rng.setSeed(seed);
    if(UserSettings::userSettings().getModel() == Model::BD){
        double cur = t->getCrown()->getTime();
        if(cur > 0.0)
            t->scaleInternalAges(1.0 / cur);
    }
    fbd = new FBDTreeModel(t, clades, fossils, seed);
    if(const char* ap = std::getenv("INIT_AGE_POW")){
        double pw = std::atof(ap);
        Tree* tt = fbd->getTree();
        double ca = tt->getCrown()->getTime();
        for(Node* n : tt->getDownPassSequence())
            if(n != tt->getCrown() && n->getIsTip() == false)
                n->setTime(ca * std::pow(n->getTime() / ca, pw));
    }
    std::vector<std::string> rogue;
    for(Fossil& f : fossils)
        rogue.push_back(f.getTaxon());
    lik = new ApproxBranchLengthLikelihood(hessianFile, mlTreeFile, rogue, 0, nStates);
    clock = new ParameterBranchRates(1.0, this, fbd->getTree(), lik->getNumPartitions(), clockModel, rgeneParam, sigma2Param);
    parameters.push_back(fbd->getParameterTree());
    lastMoveType = 2;
}

double RelaxedClockTreeModel::lnLikelihood(void){
    return lik->computeLnL(fbd->getTree(), clock->getAbsoluteRates());
}

double RelaxedClockTreeModel::lnPriorProbability(void){
    double bd = fbd->lnLikelihood();
    double fp = fbd->lnPriorProbability();
    double cl = clock->lnProbability();
    if(std::getenv("DUMP_COMP")){
        static long cnt = 0;
        if((cnt++ % 100000) == 0)
            std::fprintf(stderr, "COMP n=%ld approxBL=%.2f bdTime=%.2f fbdPrior=%.2f gbm=%.2f s2prior=%.2f muprior=%.2f s2=%.3f\n",
                         cnt, lik->computeLnL(fbd->getTree(), clock->getAbsoluteRates()), bd, fp,
                         clock->dbgGbm(), clock->dbgSigma2Prior(), clock->dbgMuPrior(), clock->getLocusSigma2(0));
    }
    return bd + fp + cl;
}

double RelaxedClockTreeModel::update(void){
    RandomVariable& r = RandomVariable::randomVariableInstance();
    double u = r.uniformRv();
    if(u < 0.25){ lastMoveType = 0; return clock->update(); }
    if(u < 0.85){ lastMoveType = 1; return clock->constantDistanceMove(); }
    if(u < 0.95){
        lastMoveType = 3;
        double cc = std::exp(0.02 * (r.uniformRv() - 0.5));
        int kAge = fbd->getTree()->scaleInternalAges(cc);
        clock->scaleAll(1.0 / cc);
        int nRate = clock->getNumLoci() * (1 + clock->getNumBranchNodes());
        return ((double)kAge - (double)nRate) * std::log(cc);
    }
    lastMoveType = 2; return fbd->update();
}

void RelaxedClockTreeModel::updateForAcceptance(void){
    if(lastMoveType == 0)
        clock->updateForAcceptance();
    else if(lastMoveType == 1){
        clock->updateForAcceptance();
        fbd->getParameterTree()->updateForAcceptance();
    }else if(lastMoveType == 3){
        clock->commitAll();
        fbd->getParameterTree()->updateForAcceptance();
    }else
        fbd->updateForAcceptance();
}

void RelaxedClockTreeModel::updateForRejection(void){
    if(lastMoveType == 0)
        clock->updateForRejection();
    else if(lastMoveType == 1){
        clock->updateForRejection();
        fbd->getParameterTree()->updateForRejection();
    }else if(lastMoveType == 3){
        clock->restoreAll();
        fbd->getParameterTree()->updateForRejection();
    }else
        fbd->updateForRejection();
}

std::vector<std::string> RelaxedClockTreeModel::getParameterNames(void){
    std::vector<std::string> n;
    n.push_back("crownAge");
    for(const std::string& s : fbd->getParameterNames())
        n.push_back(s);
    for(int p = 0; p < clock->getNumLoci(); p++){
        std::string suf = (clock->getNumLoci() > 1) ? std::to_string(p) : "";
        n.push_back("clockMean" + suf);
        n.push_back("clockSigma2" + suf);
    }
    return n;
}

std::vector<double> RelaxedClockTreeModel::getParameterString(void){
    std::vector<double> v;
    v.push_back(fbd->getTree()->getCrown()->getTime());
    for(double x : fbd->getParameterString())
        v.push_back(x);
    for(int p = 0; p < clock->getNumLoci(); p++){
        v.push_back(clock->getLocusRate(p));
        v.push_back(clock->getLocusSigma2(p));
    }
    return v;
}

void RelaxedClockTreeModel::print(void){
    fbd->print();
    clock->print();
}
