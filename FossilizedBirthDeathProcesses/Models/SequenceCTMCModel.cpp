#include "SequenceCTMCModel.hpp"

#include <string>

#include "AlignmentReader.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "ParameterSimplex.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "SequenceLikelihood.hpp"

SequenceCTMCModel::SequenceCTMCModel(PhylogeneticModel* owner, const std::string& sequenceFile, const std::string& partitionFile, int nStates, int numCats){
    this->nStates = nStates;
    seqLik = new SequenceLikelihood(nStates, numCats);

    AlignmentReader aln(sequenceFile, partitionFile, nStates);
    for(int p = 0; p < aln.getNumPartitions(); p++)
        seqLik->addPartition(aln.getTaxa(), aln.getPatterns(p), aln.getWeights(p));

    int nLoci = seqLik->getNumPartitions();
    int nExch = nStates * (nStates - 1) / 2;
    for(int p = 0; p < nLoci; p++){
        std::string suf = (nLoci > 1) ? std::to_string(p) : "";
        exch.push_back(new ParameterSimplex(1.0, owner, "exch" + suf, nExch, 1.0, 500.0));
        freq.push_back(new ParameterSimplex(1.0, owner, "freq" + suf, nStates, 1.0, 300.0));
        ParameterDouble* a = new ParameterDouble(1.0, owner, "alpha" + suf, 0.0, 1.0e8);
        a->setPrior(Probability::PriorFamily::UNIFORM, 0.0, 1.0e8);
        a->setValue(0.5);
        ParameterDouble* pv = new ParameterDouble(1.0, owner, "pinv" + suf, 0.0, 1.0);
        pv->setPrior(Probability::PriorFamily::UNIFORM, 0.0, 1.0);
        pv->setValue(0.1);
        alpha.push_back(a);
        pinv.push_back(pv);
    }
    lastSubstParm = nullptr;
}

int SequenceCTMCModel::getNumPartitions(void) const {
    return seqLik->getNumPartitions();
}

double SequenceCTMCModel::computeLnL(Tree* tree, const std::vector<std::vector<double> >& branchRates){
    int nLoci = seqLik->getNumPartitions();
    std::vector<std::vector<double> > exchV(nLoci), freqV(nLoci);
    std::vector<double> alphaV(nLoci), pinvV(nLoci);
    for(int p = 0; p < nLoci; p++){
        exchV[p] = exch[p]->getValue();
        freqV[p] = freq[p]->getValue();
        alphaV[p] = alpha[p]->getValue();
        pinvV[p] = pinv[p]->getValue();
    }
    return seqLik->computeLnL(tree, branchRates, exchV, freqV, alphaV, pinvV);
}

double SequenceCTMCModel::lnPrior(void){
    double lnp = 0.0;
    for(int p = 0; p < seqLik->getNumPartitions(); p++)
        lnp += exch[p]->lnProbability() + freq[p]->lnProbability() + alpha[p]->lnProbability() + pinv[p]->lnProbability();
    return lnp;
}

double SequenceCTMCModel::update(void){
    RandomVariable& r = RandomVariable::randomVariableInstance();
    int p = (int)(r.uniformRv() * seqLik->getNumPartitions());
    double v = r.uniformRv();
    if(v < 0.40)       lastSubstParm = exch[p];
    else if(v < 0.65)  lastSubstParm = freq[p];
    else if(v < 0.82)  lastSubstParm = alpha[p];
    else               lastSubstParm = pinv[p];
    return lastSubstParm->update();
}

void SequenceCTMCModel::updateForAcceptance(void){
    lastSubstParm->updateForAcceptance();
}

void SequenceCTMCModel::updateForRejection(void){
    lastSubstParm->updateForRejection();
}

void SequenceCTMCModel::appendParameterNames(std::vector<std::string>& names){
    for(int p = 0; p < seqLik->getNumPartitions(); p++){
        std::string suf = (seqLik->getNumPartitions() > 1) ? std::to_string(p) : "";
        for(int i = 0; i < (int)exch[p]->getValue().size(); i++)
            names.push_back("exch" + suf + "_" + std::to_string(i));
        for(int i = 0; i < (int)freq[p]->getValue().size(); i++)
            names.push_back("freq" + suf + "_" + std::to_string(i));
        names.push_back("alpha" + suf);
        names.push_back("pinv" + suf);
    }
}

void SequenceCTMCModel::appendParameterValues(std::vector<double>& values){
    for(int p = 0; p < seqLik->getNumPartitions(); p++){
        for(double x : exch[p]->getValue())
            values.push_back(x);
        for(double x : freq[p]->getValue())
            values.push_back(x);
        values.push_back(alpha[p]->getValue());
        values.push_back(pinv[p]->getValue());
    }
}
