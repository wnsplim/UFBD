#include "RelaxedClockTreeModel.hpp"

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
    return fbd->lnLikelihood() + fbd->lnPriorProbability() + clock->lnProbability();
}

double RelaxedClockTreeModel::update(void){
    RandomVariable& r = RandomVariable::randomVariableInstance();
    double u = r.uniformRv();
    if(u < 0.4){ lastMoveType = 0; return clock->update(); }
    if(u < 0.55){ lastMoveType = 1; return clock->constantDistanceMove(); }
    lastMoveType = 2; return fbd->update();
}

void RelaxedClockTreeModel::updateForAcceptance(void){
    if(lastMoveType == 0)
        clock->updateForAcceptance();
    else if(lastMoveType == 1){
        clock->updateForAcceptance();
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
