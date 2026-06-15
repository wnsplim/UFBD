#include "RelaxedClockTreeModel.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "ApproxBranchLengthLikelihood.hpp"
#include "FBDTreeModel.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "ParameterUnresolvedFossils.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

RelaxedClockTreeModel::RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& hessianFile, const std::string& mlTreeFile, int nStates, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param, unsigned int seed){
    rng.setSeed(seed);
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);
    fbd = new FBDTreeModel(t, clades, fossils, seed);
    std::vector<std::string> rogue;
    for(Fossil& f : fossils)
        rogue.push_back(f.getTaxon());
    lik = new ApproxBranchLengthLikelihood(hessianFile, mlTreeFile, rogue, 0, nStates);
    if(clockModel == ClockModel::CIR)
        clock = new ParameterBranchRatesCIR(1.0, this, fbd->getTree(), lik->getNumPartitions(), rgeneParam, sigma2Param);
    else
        clock = new ParameterBranchRates(1.0, this, fbd->getTree(), lik->getNumPartitions(), clockModel, rgeneParam, sigma2Param);
    clock->setUnresolvedFossils(fbd->getUnresolvedFossils());
    parameters.push_back(fbd->getParameterTree());
    lastMoveType = 2;
    RandomVariable::setActiveInstance(prevRng);
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
    if(u < 0.25){ lastMoveType = 0; return clock->update(); }
    if(u < 0.60){ lastMoveType = 1; return clock->constantDistanceMove(); }
    if(u < 0.75){ lastMoveType = 5; return fbd->doClockNodeAge(); }
    if(u < 0.85){
        lastMoveType = 4;
        return clock->rateAgeSubtreeMove();
    }
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
    }else if(lastMoveType == 4){
        clock->updateForAcceptance();
        fbd->getParameterTree()->updateForAcceptance();
        if(fbd->getUnresolvedFossils() != nullptr)
            fbd->getUnresolvedFossils()->updateForAcceptance();
    }else if(lastMoveType == 5){
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
    }else if(lastMoveType == 4){
        clock->updateForRejection();
        fbd->getParameterTree()->updateForRejection();
        if(fbd->getUnresolvedFossils() != nullptr)
            fbd->getUnresolvedFossils()->updateForRejection();
    }else if(lastMoveType == 5){
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
