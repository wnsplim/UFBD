#include "RelaxedClockTreeModel.hpp"

#include <cmath>

#include "ApproxBranchLengthLikelihood.hpp"
#include "FBDTreeModel.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

RelaxedClockTreeModel::RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& hessianFile, const std::string& mlTreeFile, int nStates, ClockModel clockModel, const double* rgenePara, const double* sigma2Para, unsigned int seed){
    rng.setSeed(seed);
    UserSettings& us = UserSettings::userSettings();
    fixAge = us.getFixAge();
    fixedRoot = (fixAge > 0.0);
    isOrigin = (us.getConditioning() == Conditioning::ORIGIN);
    fbd = nullptr;
    parameterTree = nullptr;
    std::vector<std::string> rogue;

    if(fixedRoot){
        double cur = t->getRoot()->getTime();
        double target = isOrigin ? 0.5 * fixAge : fixAge;
        if(cur > 0.0)
            t->scaleInternalAges(target / cur);
        parameterTree = new ParameterTree(1.0, this);
        parameterTree->setTree(t);
        if(isOrigin == false)
            parameterTree->setFixRoot(true);
        lik = new ApproxBranchLengthLikelihood(hessianFile, mlTreeFile, rogue, 0, nStates);
        clock = new ParameterBranchRates(1.0, this, parameterTree->getTree(), lik->getNumPartitions(), clockModel, rgenePara, sigma2Para);
        parameters.push_back(parameterTree);
    }else{
        fbd = new FBDTreeModel(t, clades, fossils, seed);
        for(Fossil& f : fossils)
            rogue.push_back(f.getTaxon());
        lik = new ApproxBranchLengthLikelihood(hessianFile, mlTreeFile, rogue, 0, nStates);
        clock = new ParameterBranchRates(1.0, this, fbd->getTree(), lik->getNumPartitions(), clockModel, rgenePara, sigma2Para);
        parameters.push_back(fbd->getParameterTree());
    }
    clockWeight = 0.5;
    lastMoveType = 2;
}

double RelaxedClockTreeModel::bdPrior(void){
    double lam = 1.0, mu = 0.999, rho = 0.5;
    double d = lam - mu;
    Tree* tr = parameterTree->getTree();
    if(isOrigin && tr->getRoot()->getTime() >= fixAge)
        return -INFINITY;
    double lnP = 0.0;
    for(Node* n : tr->getDownPassSequence()){
        if(n->getIsTip())
            continue;
        double e = std::exp(-d * n->getTime());
        double denom = rho * lam + (lam * (1.0 - rho) - mu) * e;
        double p1 = rho * d * d * e / (denom * denom);
        lnP += std::log(lam * p1);
    }
    return lnP;
}

double RelaxedClockTreeModel::lnLikelihood(void){
    Tree* tr = fixedRoot ? parameterTree->getTree() : fbd->getTree();
    return lik->computeLnL(tr, clock->getAbsoluteRates());
}

double RelaxedClockTreeModel::lnPriorProbability(void){
    if(fixedRoot)
        return bdPrior() + clock->lnProbability();
    return fbd->lnLikelihood() + fbd->lnPriorProbability() + clock->lnProbability();
}

double RelaxedClockTreeModel::update(void){
    RandomVariable& r = RandomVariable::randomVariableInstance();
    double u = r.uniformRv();
    if(fixedRoot){
        if(u < 0.4){ lastMoveType = 0; return clock->update(); }
        if(u < 0.7){ lastMoveType = 1; return clock->constantDistanceMove(); }
        lastMoveType = 2; return parameterTree->update();
    }
    if(u < clockWeight){ lastMoveType = 0; return clock->update(); }
    lastMoveType = 2; return fbd->update();
}

void RelaxedClockTreeModel::updateForAcceptance(void){
    if(lastMoveType == 0)
        clock->updateForAcceptance();
    else if(lastMoveType == 1){
        clock->updateForAcceptance();
        parameterTree->updateForAcceptance();
    }else if(fixedRoot)
        parameterTree->updateForAcceptance();
    else
        fbd->updateForAcceptance();
}

void RelaxedClockTreeModel::updateForRejection(void){
    if(lastMoveType == 0)
        clock->updateForRejection();
    else if(lastMoveType == 1){
        clock->updateForRejection();
        parameterTree->updateForRejection();
    }else if(fixedRoot)
        parameterTree->updateForRejection();
    else
        fbd->updateForRejection();
}

std::vector<std::string> RelaxedClockTreeModel::getParameterNames(void){
    std::vector<std::string> n;
    if(fixedRoot)
        n.push_back("rootAge");
    else
        n = fbd->getParameterNames();
    for(int p = 0; p < clock->getNumLoci(); p++){
        std::string suf = (clock->getNumLoci() > 1) ? std::to_string(p) : "";
        n.push_back("clockMean" + suf);
        n.push_back("clockSigma2" + suf);
    }
    return n;
}

std::vector<double> RelaxedClockTreeModel::getParameterString(void){
    std::vector<double> v;
    if(fixedRoot)
        v.push_back(parameterTree->getTree()->getRoot()->getTime());
    else
        v = fbd->getParameterString();
    for(int p = 0; p < clock->getNumLoci(); p++){
        v.push_back(clock->getLocusRate(p));
        v.push_back(clock->getLocusSigma2(p));
    }
    return v;
}

void RelaxedClockTreeModel::print(void){
    if(fixedRoot == false)
        fbd->print();
    clock->print();
}
