#include "RelaxedClockTreeModel.hpp"

#include "ApproxBranchLengthLikelihood.hpp"
#include "FBDTreeModel.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

RelaxedClockTreeModel::RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& hessianFile, const std::string& mlTreeFile, int nStates, ClockModel clockModel, const double* rgenePara, const double* sigma2Para, unsigned int seed){
    rng.setSeed(seed);
    fbd = new FBDTreeModel(t, clades, fossils, seed);
    std::vector<std::string> rogue;
    for(Fossil& f : fossils)
        rogue.push_back(f.getTaxon());
    lik = new ApproxBranchLengthLikelihood(hessianFile, mlTreeFile, rogue, 0, nStates);
    clock = new ParameterBranchRates(1.0, this, fbd->getTree(), lik->getNumPartitions(), clockModel, rgenePara, sigma2Para);
    parameters.push_back(fbd->getParameterTree());
    clockWeight = 0.5;
    lastWasClock = false;
}

double RelaxedClockTreeModel::lnLikelihood(void){
    return lik->computeLnL(fbd->getTree(), clock->getAbsoluteRates());
}

double RelaxedClockTreeModel::lnPriorProbability(void){
    return fbd->lnLikelihood() + fbd->lnPriorProbability() + clock->lnProbability();
}

double RelaxedClockTreeModel::update(void){
    RandomVariable& r = RandomVariable::randomVariableInstance();
    if(r.uniformRv() < clockWeight){
        lastWasClock = true;
        return clock->update();
    }
    lastWasClock = false;
    return fbd->update();
}

void RelaxedClockTreeModel::updateForAcceptance(void){
    if(lastWasClock)
        clock->updateForAcceptance();
    else
        fbd->updateForAcceptance();
}

void RelaxedClockTreeModel::updateForRejection(void){
    if(lastWasClock)
        clock->updateForRejection();
    else
        fbd->updateForRejection();
}

std::vector<std::string> RelaxedClockTreeModel::getParameterNames(void){
    std::vector<std::string> n = fbd->getParameterNames();
    for(int p = 0; p < clock->getNumLoci(); p++){
        std::string suf = (clock->getNumLoci() > 1) ? std::to_string(p) : "";
        n.push_back("clockMean" + suf);
        n.push_back("clockSigma2" + suf);
    }
    return n;
}

std::vector<double> RelaxedClockTreeModel::getParameterString(void){
    std::vector<double> v = fbd->getParameterString();
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
