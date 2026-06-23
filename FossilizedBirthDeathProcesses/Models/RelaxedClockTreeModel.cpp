#include "RelaxedClockTreeModel.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "AlignmentReader.hpp"
#include "ApproxBranchLengthLikelihood.hpp"
#include "FBDTreeModel.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterDouble.hpp"
#include "ParameterSimplex.hpp"
#include "ParameterTree.hpp"
#include "ParameterUnresolvedFossils.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "SequenceLikelihood.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

void RelaxedClockTreeModel::buildClock(ClockModel clockModel, const double* rgeneParam, const double* sigma2Param){
    int nLoci = (lik != nullptr) ? lik->getNumPartitions() : seqLik->getNumPartitions();
    if(clockModel == ClockModel::CIR)
        clock = new ParameterBranchRatesCIR(1.0, this, fbd->getTree(), nLoci, rgeneParam, sigma2Param);
    else
        clock = new ParameterBranchRates(1.0, this, fbd->getTree(), nLoci, clockModel, rgeneParam, sigma2Param);
    clock->setUnresolvedFossils(fbd->getUnresolvedFossils());
}

RelaxedClockTreeModel::RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& hessianFile, const std::string& mlTreeFile, int nStates, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param, unsigned int seed){
    rng.setSeed(seed);
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);
    fbd = new FBDTreeModel(t, clades, fossils, seed);
    seqLik = nullptr;
    std::vector<std::string> rogue;
    for(Fossil& f : fossils)
        rogue.push_back(f.getTaxon());
    lik = new ApproxBranchLengthLikelihood(hessianFile, mlTreeFile, rogue, 0, nStates);
    buildClock(clockModel, rgeneParam, sigma2Param);
    parameters.push_back(fbd->getParameterTree());
    lastSubstParm = nullptr;
    lastMoveType = 2;
    RandomVariable::setActiveInstance(prevRng);
}

RelaxedClockTreeModel::RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& sequenceFile, const std::string& partitionFile, int nStates, int numCats, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param, unsigned int seed){
    rng.setSeed(seed);
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);
    fbd = new FBDTreeModel(t, clades, fossils, seed);
    lik = nullptr;
    seqLik = new SequenceLikelihood(nStates, numCats);

    AlignmentReader aln(sequenceFile, partitionFile, nStates);
    for(int p = 0; p < aln.getNumPartitions(); p++)
        seqLik->addPartition(aln.getTaxa(), aln.getPatterns(p), aln.getWeights(p));

    buildClock(clockModel, rgeneParam, sigma2Param);

    int nLoci = seqLik->getNumPartitions();
    int nExch = nStates * (nStates - 1) / 2;
    for(int p = 0; p < nLoci; p++){
        std::string suf = (nLoci > 1) ? std::to_string(p) : "";
        exch.push_back(new ParameterSimplex(1.0, this, "exch" + suf, nExch, 1.0, 500.0));
        freq.push_back(new ParameterSimplex(1.0, this, "freq" + suf, nStates, 1.0, 300.0));
        ParameterDouble* a = new ParameterDouble(1.0, this, "alpha" + suf, 0.0, 1.0e8);
        a->setPrior(Probability::PriorFamily::UNIFORM, 0.0, 1.0e8);
        a->setValue(0.5);
        ParameterDouble* pv = new ParameterDouble(1.0, this, "pinv" + suf, 0.0, 1.0);
        pv->setPrior(Probability::PriorFamily::UNIFORM, 0.0, 1.0);
        pv->setValue(0.1);
        alpha.push_back(a);
        pinv.push_back(pv);
    }

    parameters.push_back(fbd->getParameterTree());
    lastSubstParm = nullptr;
    lastMoveType = 2;
    RandomVariable::setActiveInstance(prevRng);
}

double RelaxedClockTreeModel::lnLikelihood(void){
    if(seqLik == nullptr)
        return lik->computeLnL(fbd->getTree(), clock->getAbsoluteRates());
    int nLoci = seqLik->getNumPartitions();
    std::vector<std::vector<double>> exchV(nLoci), freqV(nLoci);
    std::vector<double> alphaV(nLoci), pinvV(nLoci);
    for(int p = 0; p < nLoci; p++){
        exchV[p] = exch[p]->getValue();
        freqV[p] = freq[p]->getValue();
        alphaV[p] = alpha[p]->getValue();
        pinvV[p] = pinv[p]->getValue();
    }
    return seqLik->computeLnL(fbd->getTree(), clock->getAbsoluteRates(), exchV, freqV, alphaV, pinvV);
}

double RelaxedClockTreeModel::substitutionUpdate(void){
    RandomVariable& r = RandomVariable::randomVariableInstance();
    int p = (int)(r.uniformRv() * seqLik->getNumPartitions());
    double v = r.uniformRv();
    if(v < 0.40)       lastSubstParm = exch[p];
    else if(v < 0.65)  lastSubstParm = freq[p];
    else if(v < 0.82)  lastSubstParm = alpha[p];
    else               lastSubstParm = pinv[p];
    return lastSubstParm->update();
}

double RelaxedClockTreeModel::lnPriorProbability(void){
    double lnp = fbd->lnLikelihood() + fbd->lnPriorProbability() + clock->lnProbability();
    if(seqLik != nullptr)
        for(int p = 0; p < seqLik->getNumPartitions(); p++)
            lnp += exch[p]->lnProbability() + freq[p]->lnProbability() + alpha[p]->lnProbability() + pinv[p]->lnProbability();
    return lnp;
}

double RelaxedClockTreeModel::update(void){
    RandomVariable& r = RandomVariable::randomVariableInstance();
    if(seqLik != nullptr && r.uniformRv() < 0.25){ lastMoveType = 6; return substitutionUpdate(); }
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
    if(lastMoveType == 6)
        lastSubstParm->updateForAcceptance();
    else if(lastMoveType == 0)
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
    if(lastMoveType == 6)
        lastSubstParm->updateForRejection();
    else if(lastMoveType == 0)
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
    if(seqLik != nullptr)
        for(int p = 0; p < seqLik->getNumPartitions(); p++){
            std::string suf = (seqLik->getNumPartitions() > 1) ? std::to_string(p) : "";
            for(int i = 0; i < (int)exch[p]->getValue().size(); i++)
                n.push_back("exch" + suf + "_" + std::to_string(i));
            for(int i = 0; i < (int)freq[p]->getValue().size(); i++)
                n.push_back("freq" + suf + "_" + std::to_string(i));
            n.push_back("alpha" + suf);
            n.push_back("pinv" + suf);
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
    if(seqLik != nullptr)
        for(int p = 0; p < seqLik->getNumPartitions(); p++){
            for(double x : exch[p]->getValue())
                v.push_back(x);
            for(double x : freq[p]->getValue())
                v.push_back(x);
            v.push_back(alpha[p]->getValue());
            v.push_back(pinv[p]->getValue());
        }
    return v;
}

void RelaxedClockTreeModel::print(void){
    fbd->print();
    clock->print();
}
