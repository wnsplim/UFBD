#include <string>

#include "AlignmentReader.hpp"
#include "EmpiricalRateMatrix.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "ParameterSimplex.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "SequenceLikelihood.hpp"
#include "SequenceCTMCModel.hpp"
#include "UserSettings.hpp"

SequenceCTMCModel::SequenceCTMCModel(PhylogeneticModel* owner, const std::string& sequenceFile, const std::string& partitionFile, int nStates, int numCats){
    this->owner = owner;
    this->nStates = nStates;
    this->numCats = numCats;
    seqLik = new SequenceLikelihood(nStates, numCats);

    AlignmentReader aln(sequenceFile, partitionFile, nStates);
    for(int p = 0; p < aln.getNumPartitions(); p++){
        seqLik->addPartition(aln.getTaxa(), aln.getPatterns(p), aln.getWeights(p));
        const std::vector<std::vector<int> >& pat = aln.getPatterns(p);
        const std::vector<int>& wt = aln.getWeights(p);
        std::vector<double> cnt(nStates, 0.0);
        for(size_t t = 0; t < pat.size(); t++)
            for(size_t h = 0; h < wt.size(); h++)
                for(int s = 0; s < nStates; s++)
                    if(pat[t][h] == (1 << s)) cnt[s] += wt[h];
        double tot = 0.0; for(double c : cnt) tot += c;
        if(tot > 0.0) for(double& c : cnt) c /= tot;
        else          for(double& c : cnt) c = 1.0 / nStates;
        observedFreq.push_back(cnt);
    }

    lastSubstParm = nullptr;
}

void SequenceCTMCModel::buildParameters(void){
    UserSettings& us = UserSettings::userSettings();
    std::string fm = us.getFreqMode();
    empirical     = (us.getSubstModel() != "gtr");
    useInvariant  = us.getUseInvariant();
    useGammaHet   = (numCats > 1);
    bool freqObserved = (fm == "empirical");
    freqEstimated = (freqObserved == false) && (empirical == false || fm == "estimate");

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
        pv->setValue(useInvariant ? 0.1 : 0.0);
        alpha.push_back(a);
        pinv.push_back(pv);
    }

    std::vector<double> datFreq;
    if(empirical){
        EmpiricalRateMatrix em(us.getSubstModel(), nStates);
        std::vector<double> e = em.getExchangeabilities();
        double se = 0.0; for(double x : e) se += x; for(double& x : e) x /= se;
        datFreq = em.getFrequencies();
        double sf = 0.0; for(double x : datFreq) sf += x; for(double& x : datFreq) x /= sf;
        for(int p = 0; p < nLoci; p++)
            exch[p]->setValue(e);
    }
    for(int p = 0; p < nLoci; p++){
        if(freqObserved)
            freq[p]->setValue(observedFreq[p]);
        else if(empirical)
            freq[p]->setValue(datFreq);
    }
}

int SequenceCTMCModel::getNumPartitions(void) const {
    return seqLik->getNumPartitions();
}

double SequenceCTMCModel::computeLnL(Tree* tree, const std::vector<std::vector<double> >& branchRates, const std::vector<std::vector<BranchMGF> >& branchMGF){
    int nLoci = seqLik->getNumPartitions();
    std::vector<std::vector<double> > exchV(nLoci), freqV(nLoci);
    std::vector<double> alphaV(nLoci), pinvV(nLoci);
    for(int p = 0; p < nLoci; p++){
        exchV[p] = exch[p]->getValue();
        freqV[p] = freq[p]->getValue();
        alphaV[p] = alpha[p]->getValue();
        pinvV[p] = pinv[p]->getValue();
    }
    return seqLik->computeLnL(tree, branchRates, exchV, freqV, alphaV, pinvV, branchMGF);
}

double SequenceCTMCModel::lnPrior(void){
    double lnp = 0.0;
    for(int p = 0; p < seqLik->getNumPartitions(); p++){
        if(empirical == false)  lnp += exch[p]->lnProbability();
        if(freqEstimated)       lnp += freq[p]->lnProbability();
        if(useGammaHet)         lnp += alpha[p]->lnProbability();
        if(useInvariant)        lnp += pinv[p]->lnProbability();
    }
    return lnp;
}

bool SequenceCTMCModel::hasMovableParams(void) const {
    return (empirical == false) || freqEstimated || useGammaHet || useInvariant;
}

double SequenceCTMCModel::update(void){
    RandomVariable& r = RandomVariable::randomVariableInstance();
    int p = (int)(r.uniformRv() * (int)exch.size());
    Parameter* cand[4]; double w[4]; int nc = 0;
    if(empirical == false)  { cand[nc] = exch[p]; w[nc++] = 0.40; }
    if(freqEstimated)       { cand[nc] = freq[p]; w[nc++] = 0.25; }
    if(useGammaHet)         { cand[nc] = alpha[p]; w[nc++] = 0.17; }
    if(useInvariant)        { cand[nc] = pinv[p]; w[nc++] = 0.18; }
    double tot = 0.0; for(int i = 0; i < nc; i++) tot += w[i];
    double v = r.uniformRv() * tot;
    double acc = 0.0; int idx = nc - 1;
    for(int i = 0; i < nc; i++){ acc += w[i]; if(v < acc){ idx = i; break; } }
    lastSubstParm = cand[idx];
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

void SequenceCTMCModel::writeState(std::ostream& os){
    for(int p = 0; p < seqLik->getNumPartitions(); p++){
        exch[p]->writeState(os);
        freq[p]->writeState(os);
        alpha[p]->writeState(os);
        pinv[p]->writeState(os);
    }
}

void SequenceCTMCModel::readState(std::istream& is){
    for(int p = 0; p < seqLik->getNumPartitions(); p++){
        exch[p]->readState(is);
        freq[p]->readState(is);
        alpha[p]->readState(is);
        pinv[p]->readState(is);
    }
}
