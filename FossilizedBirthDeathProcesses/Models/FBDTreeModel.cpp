#include "Node.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "ParameterShrinkageField.hpp"
#include "FBDTreeModel.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"
#include "Probability.hpp"
#include "Msg.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>

FBDTreeModel::FBDTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, unsigned int seed) :
    PhylogeneticModel(){

    lastWasJointScale = false;
    lastWasUpDown = false;
    lastWasRateVec = false;
    lastRateVec = nullptr;
    rateVecStep = 0.2;
    shrinkStep = 0.2;
    rvAccW = rvAttW = seAccW = seAttW = 0;
    upDownStep = 0.1;
    upDownTotal = 0;
    cacheInit = false;
    lambdaField = nullptr;
    muField = nullptr;
    psiField = nullptr;
    rng.setSeed(seed);
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);

    parameterTree = new ParameterTree(1.0, this);
    isFBD = (UserSettings::userSettings().getModel() == Model::FBD);

    int numBackbone = t->getNumBackbone();
    int numUE = 0;
    for(Fossil& f : fossils)
        if(f.getMaxAge() == 0.0) numUE++;
    Conditioning condPoint = UserSettings::userSettings().getConditioning();
    ConditioningEvent condEvent = UserSettings::userSettings().getConditioningEvent();
    if(condPoint == Conditioning::CROWN && numBackbone < 2)
        Msg::error("crown conditioning requires at least 2 backbone tips, but the backbone tree has " + std::to_string(numBackbone) + ".");
    if(condEvent == ConditioningEvent::SURVIVAL && (numBackbone + numUE) < 1)
        Msg::error("survival conditioning requires at least 1 extant taxon.");
    if(condEvent == ConditioningEvent::EXTINCT && (numBackbone + numUE) > 0)
        Msg::error("extinct conditioning requires 0 extant taxa, but the data has " + std::to_string(numBackbone + numUE) + ".");

    originAge = nullptr;
    if(UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN){
        double x0init = t->getCrown()->getTime();
        for(Fossil& f : fossils)
            if(f.getMaxAge() > x0init)
                x0init = f.getMaxAge();
        x0init *= 1.05;
        originAge = new ParameterDouble(1.0, this, "originAge", 0.0, std::numeric_limits<double>::max());
        UserSettings& us = UserSettings::userSettings();
        if(us.getConditionAgePriorSet()){
            originAge->setPrior(us.getConditionAgePrior(), us.getConditionAgePriorP1(), us.getConditionAgePriorP2());
            double pm = Probability::priorMean(us.getConditionAgePrior(), us.getConditionAgePriorP1(), us.getConditionAgePriorP2());
            if(pm > x0init)
                x0init = pm;
        }else{
            originAge->setPrior(Probability::PriorFamily::IMPROPER, 0.0, 0.0);
        }
        originAge->setValue(x0init);
    }

    if(isFBD){
        Tree* working = new Tree(*t);
        if(originAge != nullptr)
            working->addOriginNode(originAge->getValue());
        resolveFossils(working, clades, fossils);
        parameterTree->setTree(working);
        delete working;
    }else{
        parameterTree->setTree(t);
    }
    parameters.push_back(parameterTree);
    if(originAge != nullptr)
        parameters.push_back(originAge);
    
    intervalStart.push_back(0.0);
    for(double t : UserSettings::userSettings().getSkylineTimes())
        intervalStart.push_back(t);
    Probability::PriorSpec lp = UserSettings::userSettings().getLambdaPrior();
    Probability::PriorSpec mp = UserSettings::userSettings().getMuPrior();
    Probability::PriorSpec pp = UserSettings::userSettings().getPsiPrior();
    UserSettings& rateUs = UserSettings::userSettings();
    int nB = (int)intervalStart.size();
    bool lamSmooth = (rateUs.getLambdaMode() == RateMode::SMOOTH);
    bool muSmooth  = (rateUs.getMuMode() == RateMode::SMOOTH);
    bool psiSmooth = (rateUs.getPsiMode() == RateMode::SMOOTH);
    if(lamSmooth && nB < 2){ Msg::warning("Smoothing (HSMRF) set for speciation rate (lambda) but only single rate interval."); lamSmooth = false; }
    if(muSmooth && nB < 2){ Msg::warning("Smoothting (HSMRF) set for extinction rate (mu) but only single rate interval."); muSmooth = false; }
    if(psiSmooth && nB < 2){ Msg::warning("Smoothing (HSMRF) set for sampling rate (psi) but only single rate interval."); psiSmooth = false; }
    double nShifts = rateUs.getHsmrfShifts();
    double shiftSize = rateUs.getHsmrfShiftSize();
    if(lamSmooth || muSmooth || psiSmooth){
        if(lamSmooth){
            lambdaField = new ParameterShrinkageField(1.0, this, nB, lp, nShifts, shiftSize);
            parameters.push_back(lambdaField);
        }else{
            for(int i = 0; i < nB; i++){
                std::string suf = (nB > 1) ? std::to_string(i) : "";
                ParameterDouble* l = new ParameterDouble(1.0, this, "lambda" + suf, 0.0, std::numeric_limits<double>::max());
                if(lp.set) l->setPrior(lp.family, lp.p1, lp.p2);
                lambda.push_back(l);
                parameters.push_back(l);
            }
        }
        if(muSmooth){
            muField = new ParameterShrinkageField(1.0, this, nB, mp, nShifts, shiftSize);
            parameters.push_back(muField);
        }else{
            for(int i = 0; i < nB; i++){
                std::string suf = (nB > 1) ? std::to_string(i) : "";
                ParameterDouble* m = new ParameterDouble(1.0, this, "mu" + suf, 0.0, std::numeric_limits<double>::max());
                if(mp.set) m->setPrior(mp.family, mp.p1, mp.p2);
                mu.push_back(m);
                parameters.push_back(m);
            }
        }
        if(psiSmooth){
            psiField = new ParameterShrinkageField(1.0, this, nB, pp, nShifts, shiftSize);
            parameters.push_back(psiField);
        }else{
            for(int i = 0; i < nB; i++){
                std::string suf = (nB > 1) ? std::to_string(i) : "";
                ParameterDouble* p = new ParameterDouble(1.0, this, "psi" + suf, 0.0, std::numeric_limits<double>::max());
                if(pp.set) p->setPrior(pp.family, pp.p1, pp.p2);
                psi.push_back(p);
                parameters.push_back(p);
            }
        }
    }else{
        for(size_t i = 0; i < intervalStart.size(); i++){
            std::string suf = (intervalStart.size() > 1) ? std::to_string(i) : "";
            ParameterDouble* l = new ParameterDouble(1.0, this, "lambda" + suf, 0.0, std::numeric_limits<double>::max());
            ParameterDouble* m = new ParameterDouble(1.0, this, "mu" + suf, 0.0, std::numeric_limits<double>::max());
            if(lp.set) l->setPrior(lp.family, lp.p1, lp.p2);
            if(mp.set) m->setPrior(mp.family, mp.p1, mp.p2);
            lambda.push_back(l);
            mu.push_back(m);
            parameters.push_back(l);
            parameters.push_back(m);
            ParameterDouble* p = new ParameterDouble(1.0, this, "psi" + suf, 0.0, std::numeric_limits<double>::max());
            if(pp.set) p->setPrior(pp.family, pp.p1, pp.p2);
            psi.push_back(p);
            parameters.push_back(p);
        }
    }
    for(size_t i = 1; i < lambda.size(); i++) lambda[i]->setValue(lambda[0]->getValue());
    for(size_t i = 1; i < mu.size(); i++)     mu[i]->setValue(mu[0]->getValue());
    for(size_t i = 1; i < psi.size(); i++)    psi[i]->setValue(psi[0]->getValue());
    rho = UserSettings::userSettings().getRho();

    unresolvedFossils = nullptr;
    if(isFBD){
        Tree* wt = parameterTree->getTree();
        for(Fossil& f : fossils){
            Clade* clade = nullptr;
            for(Clade& c : clades)
                if(c.getName() == f.getClade()){
                    clade = &c;
                    break;
                }
            fossilCrown.push_back(wt->getMRCA(clade->getTaxa()));
            bool originCond = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);
            fossilOrigin.push_back(originCond ? wt->getCrown() : wt->getNodeByOffset(clade->getOrigin()->getOffset()));
        }
    }else{
        unresolvedFossils = new ParameterUnresolvedFossils(1.0, this, parameterTree->getTree(), clades, fossils, originAge);
        parameters.push_back(unresolvedFossils);
    }

    if(isFBD){
        parameterTree->setProposalProbability(78.0);
        for(ParameterDouble* l : lambda) l->setProposalProbability(15.0);
        for(ParameterDouble* m : mu)     m->setProposalProbability(15.0);
        for(ParameterDouble* p : psi)    p->setProposalProbability(15.0);
    }
    double fieldWeight = (isFBD ? 15.0 : 1.0) * (double)nB;
    if(lambdaField) lambdaField->setProposalProbability(fieldWeight);
    if(muField)     muField->setProposalProbability(fieldWeight);
    if(psiField)    psiField->setProposalProbability(fieldWeight);

    double sum = 0.0;
    for(Parameter* p : parameters)
        sum += p->getProposalProbability();
    for(Parameter* p : parameters)
        p->setProposalProbability(p->getProposalProbability() / sum);

    RandomVariable::setActiveInstance(prevRng);
}

double FBDTreeModel::lambdaAt(int i){
    if(lambdaField != nullptr)
        return lambdaField->getRate(i);
    ParameterDouble* d = lambda[i];
    return d->getValue();
}

double FBDTreeModel::muAt(int i){
    if(muField != nullptr)
        return muField->getRate(i);
    ParameterDouble* d = mu[i];
    return d->getValue();
}

double FBDTreeModel::psiAt(int i){
    if(psiField != nullptr)
        return psiField->getRate(i);
    ParameterDouble* d = psi[i];
    return d->getValue();
}

static bool nodeInSubtree(Node* node, Node* subtreeCrown){
    for(Node* p = node; ; p = p->getAncestor()){
        if(p == subtreeCrown)
            return true;
        if(p->getAncestor() == p) // reached the crown self-loop
            return false;
    }
}

static bool nodeOnStalk(Node* n, Node* crown, Node* origin){
    for(Node* p = crown; p != origin; p = p->getAncestor()){
        if(p == n)
            return true;
        if(p->getAncestor() == p)
            break;
    }
    return false;
}

std::vector<std::string> FBDTreeModel::getParameterNames(void){
    std::vector<std::string> names;
    for(Parameter* p : parameters)
        if( p != parameterTree){ //by convention, exclude parameter tree from these getters
            ParameterShrinkageField* sf = dynamic_cast<ParameterShrinkageField*>(p);
            if(sf != nullptr){
                std::string base = (p == lambdaField) ? "lambda" : (p == muField ? "mu" : "psi");
                for(int i = 0; i < sf->getNumBins(); i++)
                    names.push_back(base + std::to_string(i));
                continue;
            }
            ParameterUnresolvedFossils* uf = dynamic_cast<ParameterUnresolvedFossils*>(p);
            names.push_back(uf != nullptr ? "nSA" : p->getName());
        }
    if(isFBD)
        names.push_back("nSA");
    return names;
}

int FBDTreeModel::countResolvedSA(void){
    Tree* tree = parameterTree->getTree();
    int s = 0;
    for(Node* n : tree->getDownPassSequence())
        if(n->getIsTip() && n->getIsFossil() && n->getAncestor()->getTime() == n->getTime())
            s++;
    return s;
}

std::vector<double> FBDTreeModel::getParameterString(void){
    std::vector<double> vals;
    for(Parameter* p : parameters)
        if( p != parameterTree){ //by convention, exclude parameter tree from these getters
            ParameterShrinkageField* sf = dynamic_cast<ParameterShrinkageField*>(p);
            if(sf != nullptr){
                for(int i = 0; i < sf->getNumBins(); i++)
                    vals.push_back(sf->getRate(i));
                continue;
            }
            ParameterDouble* pd = dynamic_cast<ParameterDouble*>(p);
            if(pd != nullptr){
                vals.push_back(pd->getValue());
                continue;
            }
            ParameterUnresolvedFossils* uf = dynamic_cast<ParameterUnresolvedFossils*>(p);
            if(uf != nullptr)
                vals.push_back((double)uf->getNumSampledAncestors());
        }
    if(isFBD)
        vals.push_back((double)countResolvedSA());

    return vals;
}

double FBDTreeModel::lnLikelihood(void){
    double v = calculateFBDProbability();
    if(std::isfinite(v) == false)
        return -INFINITY;
    return v;
}

double FBDTreeModel::lnPriorProbability(void){
    double lnP = 0.0;
    
    for(Parameter* p : parameters)
        if( p != parameterTree)
            lnP += p->lnProbability();

    UserSettings& us = UserSettings::userSettings();
    if(us.getConditioning() == Conditioning::CROWN){
        double crownAge = parameterTree->getTree()->getCrown()->getTime();
        lnP += Probability::priorLnPdf(us.getConditionAgePrior(), us.getConditionAgePriorP1(), us.getConditionAgePriorP2(), crownAge, 0.0, std::numeric_limits<double>::max());
    }

    return lnP;
}

void FBDTreeModel::print(void){
    for(Parameter* p : parameters){
        if(p != parameterTree && p->getParmPrintConsole() == true)
            std::cout << p->getName() << " (A/R): " << p->getAcceptanceRatio() << "\t";
    }
    std::cout << "tree (A/R): " << parameterTree->getAcceptanceRatio() << "\tscaleLambda: " << parameterTree->getScaleLambda() << "\n";
}

std::vector<ParameterDouble*>* FBDTreeModel::pickIidRateVector(void){
    std::vector<std::vector<ParameterDouble*>*> cands;
    if(lambdaField == nullptr && lambda.size() >= 2) cands.push_back(&lambda);
    if(muField == nullptr && mu.size() >= 2)         cands.push_back(&mu);
    if(psiField == nullptr && psi.size() >= 2)       cands.push_back(&psi);
    if(cands.empty())
        return nullptr;
    return cands[(int)(rng.uniformRv() * cands.size())];
}

double FBDTreeModel::doRateVectorScale(void){
    lastWasRateVec = true;
    lastRateVecScale = true;
    lastRateVec = pickIidRateVector();
    double mB = 0.95;
    double d = mB + Probability::Normal::rv(&rng) * std::sqrt(1.0 - mB * mB);
    if(rng.uniformRv() < 0.5) d = -d;
    double c = std::exp(rateVecStep * d);
    for(ParameterDouble* p : *lastRateVec)
        p->scaleProposed(c);
    rvAttW++;
    if(rvAttW >= 200){
        rateVecStep *= std::exp((double)rvAccW / rvAttW - 0.3);
        if(rateVecStep < 1e-3) rateVecStep = 1e-3;
        if(rateVecStep > 10.0)  rateVecStep = 10.0;
        rvAccW = 0;
        rvAttW = 0;
    }
    return (double)lastRateVec->size() * std::log(c);
}

double FBDTreeModel::doRateShrinkExpand(void){
    lastWasRateVec = true;
    lastRateVecScale = false;
    lastRateVec = pickIidRateVector();
    int n = (int)lastRateVec->size();
    double logMean = 0.0;
    for(ParameterDouble* p : *lastRateVec)
        logMean += std::log(p->getValue());
    logMean /= (double)n;
    double mB = 0.95;
    double d = mB + Probability::Normal::rv(&rng) * std::sqrt(1.0 - mB * mB);
    if(rng.uniformRv() < 0.5) d = -d;
    double a = std::exp(shrinkStep * d);
    for(ParameterDouble* p : *lastRateVec){
        double cur = p->getValue();
        double target = std::exp(logMean + a * (std::log(cur) - logMean));
        p->scaleProposed(target / cur);
    }
    seAttW++;
    if(seAttW >= 200){
        shrinkStep *= std::exp((double)seAccW / seAttW - 0.3);
        if(shrinkStep < 1e-3) shrinkStep = 1e-3;
        if(shrinkStep > 10.0)  shrinkStep = 10.0;
        seAccW = 0;
        seAttW = 0;
    }
    return (double)(n - 1) * std::log(a);
}

double FBDTreeModel::update(void){
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);

    lastWasJointScale = false;
    lastWasUpDown = false;
    lastWasRateVec = false;
    bool haveIid = (lambdaField == nullptr && lambda.size() >= 2)
                || (muField == nullptr && mu.size() >= 2)
                || (psiField == nullptr && psi.size() >= 2);
    if(haveIid && rng.uniformRv() < 0.20){
        double r = (rng.uniformRv() < 0.5) ? doRateVectorScale() : doRateShrinkExpand();
        RandomVariable::setActiveInstance(prevRng);
        return r;
    }

    double u = rng.uniformRv();
    double sum = 0.0;
    for(Parameter* p  : parameters){
        sum += p->getProposalProbability();
        if(u < sum){
            updatedParameter = p;
            break;
        }
    }
    if(updatedParameter == parameterTree && isFBD == false && unresolvedFossils != nullptr){
        Tree* t = parameterTree->getTree();
        int numSlideable = 0;
        for(Node* n : t->getDownPassSequence())
            if(n != t->getCrown() && n->getIsTip() == false)
                numSlideable++;
        double fixedWeight = 3.0;
        double slideAndCrown = numSlideable + fixedWeight;
        double uMove = rng.uniformRv() * (slideAndCrown + 3.0 * fixedWeight);
        if(uMove >= slideAndCrown + 2.0 * fixedWeight){
            double r = doUpDownScale();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove >= slideAndCrown){
            lastWasJointScale = true;
            double r = (uMove < slideAndCrown + fixedWeight) ? doJointScale() : doSubtreeScale();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        std::map<Node*,double> floors;
        computeAgeFloors(floors);
        t->setAgeFloors(floors);
    }
    else if(updatedParameter == parameterTree && isFBD){
        parameterTree->getTree()->setLastUpdateWasScale(false);
        double wNE = 15.0;
        double wWB = 10.0;
        double wWE = 10.0;
        double wTreeScale = 3.0;
        double wSA = 15.0;
        double wUpDown = 10.0;
        double wAge = 40.0;
        double uMove = rng.uniformRv() * (wNE + wWB + wWE + wTreeScale + wSA + wUpDown + wAge);
        if(uMove < wNE){
            double r = doNarrowExchange();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove < wNE + wWB){
            double r = doWilsonBalding();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove < wNE + wWB + wWE){
            double r = doWideExchange();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove < wNE + wWB + wWE + wTreeScale){
            double r = doTreeScale();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove < wNE + wWB + wWE + wTreeScale + wSA){
            double r = doSARJMCMC();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove < wNE + wWB + wWE + wTreeScale + wSA + wUpDown){
            double r = doUpDownScale();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
    }
    double ratio = updatedParameter->update();

    RandomVariable::setActiveInstance(prevRng);
    return ratio;
}

void FBDTreeModel::updateForAcceptance(void){
    if(lastWasRateVec){
        for(ParameterDouble* p : *lastRateVec)
            p->commitProposed();
        if(lastRateVecScale) rvAccW++; else seAccW++;
        return;
    }
    if(lastWasUpDown){
        for(ParameterDouble* l : lambda) l->commitProposed();
        for(ParameterDouble* m : mu) m->commitProposed();
        for(ParameterDouble* p : psi) p->commitProposed();
        if(originAge != nullptr) originAge->commitProposed();
        parameterTree->updateForAcceptance();
        if(unresolvedFossils != nullptr) unresolvedFossils->updateForAcceptance();
        upDownTotal++;
        upDownRecent.push_back(true);
        if(upDownRecent.size() > 1000) upDownRecent.pop_front();
    }else if(lastWasJointScale){
        parameterTree->updateForAcceptance();
        unresolvedFossils->updateForAcceptance();
    }else{
        updatedParameter->updateForAcceptance();
    }
    if(isFBD && originAge != nullptr)
        parameterTree->getTree()->getCrown()->setTime(originAge->getValue());
}

void FBDTreeModel::updateForRejection(void){
    if(lastWasRateVec){
        for(ParameterDouble* p : *lastRateVec)
            p->restoreProposed();
        return;
    }
    if(lastWasUpDown){
        for(ParameterDouble* l : lambda) l->restoreProposed();
        for(ParameterDouble* m : mu) m->restoreProposed();
        for(ParameterDouble* p : psi) p->restoreProposed();
        if(originAge != nullptr) originAge->restoreProposed();
        parameterTree->updateForRejection();
        if(unresolvedFossils != nullptr) unresolvedFossils->updateForRejection();
        upDownTotal++;
        upDownRecent.push_back(false);
        if(upDownRecent.size() > 1000) upDownRecent.pop_front();
    }else if(lastWasJointScale){
        parameterTree->updateForRejection();
        unresolvedFossils->updateForRejection();
    }else{
        updatedParameter->updateForRejection();
    }
    if(isFBD && originAge != nullptr)
        parameterTree->getTree()->getCrown()->setTime(originAge->getValue());
}

double FBDTreeModel::calculateFBDProbability(void){
    Tree* tree = parameterTree->getTree();

    if(isFBD && originAge != nullptr)
        tree->getCrown()->setTime(originAge->getValue());

    int numInternalNodes = tree->getNumNodes() - tree->getNumBackbone();
    double crownAge = tree->getCrown()->getTime();
    std::vector<Node*> dpseq = tree->getDownPassSequence();
    
    lambdaVal = lambdaAt(0);
    muVal = muAt(0);
    rhoVal = rho;
    psiVal = psiAt(0);
    double log4LambdaRho = std::log(4*lambdaVal*rhoVal);
    
    prepareIntervals();

    if(isFBD)
        return calculateResolvedFBD();

    bool useOrigin = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);
    double fbdProb = 0.0;

    //term 1: conditioning
    if(useOrigin){
        double x0 = originAge->getValue();
        if(x0 < crownAge)
            return -INFINITY;
        for(int i = 0; i < unresolvedFossils->getNumFossils(); i++){
            if(unresolvedFossils->getFossilAge(i) > x0)
                return -INFINITY;
            if(i != unresolvedFossils->getSpineIdx() && unresolvedFossils->getAttachAge(i) > x0)
                return -INFINITY;
        }
        double lx0 = lambdaAt(findIndex(x0));
        fbdProb -= std::log(lx0);
        fbdProb -= calculateLnConditioning(x0);
        fbdProb += std::log(4 * lx0 * rhoVal);
        fbdProb += lnD(x0) - std::log(4.0);
    }else{
        double lr = lambdaAt(findIndex(crownAge));
        fbdProb -= 2 * (std::log(lr) + calculateLnSurvival(crownAge));
        fbdProb += std::log(4 * lr * rhoVal);
        fbdProb += lnD(crownAge) - std::log(4.0);
    }

    //term 2: main body
    for(Node* n : dpseq){
        if(n->getIsTip())
            continue;
        bool hasChild = false;
        for(Node* c : n->getNeighbors())
            if(c != n->getAncestor()){ hasChild = true; break; }
        if(hasChild)
            fbdProb += std::log(lambdaAt(findIndex(n->getTime())) * rhoVal) + lnD(n->getTime());
    }

    //term 3: fossil attachment
    int numFossils = unresolvedFossils->getNumFossils();
    if(originAge != nullptr)
        unresolvedFossils->syncSpine(originAge->getValue());
    updateGammaCache();
    for(int i = 0; i < numFossils; i++){
        if(unresolvedFossils->isSA(i)){
            fbdProb += std::log(psiAt(findIndex(unresolvedFossils->getFossilAge(i)))) + cachedGammaLn[i];
            continue;
        }
        if(unresolvedFossils->isUE(i)){
            fbdProb += uePqLn(unresolvedFossils->getAttachAge(i)) + cachedGammaLn[i];
            continue;
        }
        if(i == unresolvedFossils->getSpineIdx()){
            double ys = unresolvedFossils->getFossilAge(i);
            fbdProb += std::log(psiAt(findIndex(ys))) + std::log(calculateP0(ys)) - lnD(ys);
            continue;
        }
        fbdProb += fossilPqLn(unresolvedFossils->getFossilAge(i), unresolvedFossils->getAttachAge(i)) + cachedGammaLn[i];
    }
    return fbdProb;
}

double FBDTreeModel::calculateResolvedFBD(void){
    Tree* tree = parameterTree->getTree();
    Node* crown = tree->getCrown();
    double crownAge = crown->getTime();
    bool useOrigin = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);

    double lnP;
    if(useOrigin){
        double x0 = originAge->getValue();
        for(Node* c : crown->getNeighbors())
            if(c != crown->getAncestor() && c->getTime() >= x0)
                return -INFINITY;
        lnP = -calculateLnConditioning(x0);
    }else{
        lnP = -2.0 * calculateLnSurvival(crownAge);
    }

    for(Node* n : tree->getDownPassSequence()){
        if(n != crown)
            lnP += lnD(n->getAncestor()->getTime()) - lnD(n->getTime());
        if(n->getIsTip()){
            if(n->getIsFossil() == false)
                lnP += std::log(rhoVal);
            else if(n->getAncestor()->getTime() == n->getTime())
                lnP += std::log(psiAt(findIndex(n->getTime())));
            else
                lnP += std::log(psiAt(findIndex(n->getTime()))) + std::log(calculateP0(n->getTime()));
        }
        else if(n != crown){
            bool fakeSplit = false;
            for(Node* c : n->getNeighbors())
                if(c != n->getAncestor() && c->getIsTip() && c->getIsFossil() && c->getTime() == n->getTime()){ fakeSplit = true; break; }
            if(fakeSplit == false)
                lnP += std::log(2.0 * lambdaAt(findIndex(n->getTime())));
        }
    }
    return lnP;
}

double FBDTreeModel::lnD(double t){
    return (t <= 0.0) ? 0.0 : lnDPrev[findIndex(t)] + std::log(4.0) - calculateLnQtAt(findIndex(t), t);
}

double FBDTreeModel::fossilPqLn(double y, double z){
    return std::log(psiAt(findIndex(y))) + std::log(2*lambdaAt(findIndex(z))) + std::log(calculateP0(y)) + lnD(z) - lnD(y);
}

double FBDTreeModel::uePqLn(double z){
    return std::log(rhoVal) + std::log(2*lambdaAt(findIndex(z))) + lnD(z);
}

double FBDTreeModel::calculateLnSurvival(double t){
    if(intervalStart.size() == 1){
        double a = lambdaVal - muVal;
        double B = lambdaVal * (1.0 - rhoVal) - muVal;
        double lnAbsNum = std::log(std::abs(rhoVal * a));
        double lnAbsDenom;
        if(-a * t > 0.0)
            lnAbsDenom = (-a * t) + std::log(std::abs(lambdaVal * rhoVal * std::exp(a * t) + B));
        else
            lnAbsDenom = std::log(std::abs(lambdaVal * rhoVal + B * std::exp(-a * t)));
        return lnAbsNum - lnAbsDenom;
    }
    return std::log(1.0 - calculateP0Hat(t));
}

double FBDTreeModel::calculateLnConditioning(double t){
    switch(UserSettings::userSettings().getConditioningEvent()){
        case ConditioningEvent::SURVIVAL:
            return calculateLnSurvival(t);
        case ConditioningEvent::ANYSAMPLE:
            return std::log(1.0 - calculateP0(t));
        case ConditioningEvent::EXTINCT: {
            double pHat = 1.0 - std::exp(calculateLnSurvival(t));
            double d = pHat - calculateP0(t);
            return (d > 0.0) ? std::log(d) : -INFINITY;
        }
    }
    return -INFINITY;
}

void FBDTreeModel::prepareIntervals(void){
    size_t n = intervalStart.size();
    c1Vec.assign(n, 0.0);
    c2Vec.assign(n, 0.0);
    ePrev.assign(n, 0.0);
    lnDPrev.assign(n, 0.0);
    for(size_t i = 0; i < n; i++){
        double li = lambdaAt(i);
        double mi = muAt(i);
        double pi = psiAt(i);
        c1Vec[i] = std::abs(std::sqrt(std::pow(li - mi - pi, 2) + 4*li*pi));
        if(i == 0){
            ePrev[0] = 1.0;
            lnDPrev[0] = 0.0;
            c2Vec[0] = (-li + mi + 2*li * rhoVal + pi) / c1Vec[0];
        }else{
            double s = intervalStart[i];
            ePrev[i] = calculateP0At((int)i - 1, s);
            lnDPrev[i] = lnDPrev[i-1] + std::log(4.0) - calculateLnQtAt((int)i - 1, s);
            c2Vec[i] = ((1.0 - 2.0 * ePrev[i]) * li + mi + pi) / c1Vec[i];
        }
    }
    if(n > 1){
        c1HatVec.assign(n, 0.0);
        c2HatVec.assign(n, 0.0);
        ePrevHat.assign(n, 0.0);
        for(size_t i = 0; i < n; i++){
            double li = lambdaAt(i);
            double mi = muAt(i);
            c1HatVec[i] = std::abs(li - mi);
            if(i == 0){
                ePrevHat[0] = 1.0;
                c2HatVec[0] = (-li + mi + 2*li * rhoVal) / c1HatVec[0];
            }else{
                ePrevHat[i] = calculateP0HatAt((int)i - 1, intervalStart[i]);
                c2HatVec[i] = ((1.0 - 2.0 * ePrevHat[i]) * li + mi) / c1HatVec[i];
            }
        }
    }
}

int FBDTreeModel::findIndex(double t){
    int i = 0;
    for(size_t j = 1; j < intervalStart.size(); j++)
        if(t >= intervalStart[j])
            i = (int)j;
    return i;
}

double FBDTreeModel::calculateLnQtAt(int i, double t){
    double tau = t - intervalStart[i];
    double a = c1Vec[i] * tau;
    double c2 = c2Vec[i];
    double bracket = std::pow(1.0 + c2, 2) + std::exp(-a) * 2.0 * (1.0 - c2 * c2) + std::exp(-2.0 * a) * std::pow(1.0 - c2, 2);
    return a + std::log(bracket);
}

double FBDTreeModel::calculateP0At(int i, double t){
    double tau = t - intervalStart[i];
    double li = lambdaAt(i);
    double tmp = -li + muAt(i) + psiAt(i);
    tmp += c1Vec[i] * (std::exp(-c1Vec[i] * tau) * (1 - c2Vec[i]) - (1+c2Vec[i]) ) / ( std::exp(-c1Vec[i] * tau) * (1 - c2Vec[i]) + (1+c2Vec[i])  );
    tmp /= 2*li;
    return 1 + tmp;
}

double FBDTreeModel::calculateP0(double t){
    return calculateP0At(findIndex(t), t);
}

double FBDTreeModel::calculateP0HatAt(int i, double t){
    double tau = t - intervalStart[i];
    double li = lambdaAt(i);
    double tmp = -li + muAt(i);
    tmp += c1HatVec[i] * (std::exp(-c1HatVec[i] * tau) * (1 - c2HatVec[i]) - (1+c2HatVec[i]) ) / ( std::exp(-c1HatVec[i] * tau) * (1 - c2HatVec[i]) + (1+c2HatVec[i])  );
    tmp /= 2*li;
    return 1 + tmp;
}

double FBDTreeModel::calculateP0Hat(double t){
    return calculateP0HatAt(findIndex(t), t);
}

void FBDTreeModel::buildEulerIndex(void){
    Tree* tree = parameterTree->getTree();
    int n = tree->getNumNodes();
    subPre.assign(n, 0);
    subSize.assign(n, 0);
    nodesByPre.assign(n, nullptr);
    int counter = 0;
    std::vector<std::pair<Node*,bool> > stk;
    stk.push_back(std::make_pair(tree->getCrown(), false));
    while(stk.empty() == false){
        std::pair<Node*,bool> top = stk.back();
        stk.pop_back();
        Node* nd = top.first;
        if(top.second == false){
            subPre[nd->getOffset()] = counter;
            nodesByPre[counter] = nd;
            counter++;
            stk.push_back(std::make_pair(nd, true));
            for(Node* c : nd->getDescendants())
                stk.push_back(std::make_pair(c, false));
        }else{
            int s = 1;
            for(Node* c : nd->getDescendants())
                s += subSize[c->getOffset()];
            subSize[nd->getOffset()] = s;
        }
    }
    eulerBuilt = true;
}

bool FBDTreeModel::inSub(Node* node, Node* subtreeCrown){
    int pc = subPre[subtreeCrown->getOffset()];
    int pn = subPre[node->getOffset()];
    return pc <= pn && pn < pc + subSize[subtreeCrown->getOffset()];
}

double FBDTreeModel::computeGamma(double z, int i){
    Tree* tree = parameterTree->getTree();
    if(eulerBuilt == false)
        buildEulerIndex();
    Node* crown = unresolvedFossils->getCrownNode(i);
    bool stem = unresolvedFossils->getIsStem(i);
    bool total = (unresolvedFossils->getIsCrown(i) == false && stem == false);
    double count = 0.0;
    if(stem == false && crown == tree->getCrown()){
        int younger = (int)(std::lower_bound(sortedYounger.begin(), sortedYounger.end(), z) - sortedYounger.begin());
        int older   = (int)(std::upper_bound(sortedOlder.begin(), sortedOlder.end(), z) - sortedOlder.begin());
        count += (double)(younger - older);
    }else{
        int lo = subPre[crown->getOffset()];
        int hi = lo + subSize[crown->getOffset()];
        for(int k = lo; k < hi; k++){
            Node* n = nodesByPre[k];
            if(n == tree->getCrown())
                continue;
            Node* anc = n->getAncestor();
            if(n->getTime() < z && z < anc->getTime()){
                bool inZone;
                if(stem){
                    inZone = (n == crown);
                }else{
                    inZone = inSub(anc, crown);
                    if(inZone == false && total && n == crown)
                        inZone = true;
                }
                if(inZone)
                    count++;
            }
        }
    }
    if((total || stem) && crown == tree->getCrown() && originAge != nullptr){
        double x0 = originAge->getValue();
        if(tree->getNumBackbone() == 0){
            if(z >= x0)
                count++;
        }else{
            if(tree->getCrown()->getTime() < z && z < x0)
                count++;
        }
    }
    bool SymmetryCorrection = (UserSettings::userSettings().getModel() == Model::UFBD);
    bool focalIsTip = (unresolvedFossils->isSA(i) == false);
    int numFossils = unresolvedFossils->getNumFossils();
    for(int j = 0; j < numFossils; j++){
        if(j == i)
            continue;
        if(unresolvedFossils->isSA(j))
            continue;
        double yj = unresolvedFossils->getFossilAge(j);
        double zj = unresolvedFossils->getAttachAge(j);
        if(yj >= z || z >= zj)
            continue;
        bool reciprocal;
        if(stem){
            if(unresolvedFossils->getCrownNode(j) != crown)
                continue;
            bool jStem = unresolvedFossils->getIsStem(j);
            bool jTotalOnStem = (unresolvedFossils->getIsCrown(j) == false) && (jStem == false) && (zj >= crown->getTime());
            if(jStem == false && jTotalOnStem == false)
                continue;
            reciprocal = true;
        }else{
            if(inSub(unresolvedFossils->getCrownNode(j), crown) == false)
                continue;
            if(total == false && zj > crown->getTime())
                continue;
            Node* crownJ = unresolvedFossils->getCrownNode(j);
            bool jStem = unresolvedFossils->getIsStem(j);
            bool jTotal = (unresolvedFossils->getIsCrown(j) == false) && (jStem == false);
            bool iPendReachesRj = jTotal ? true
                                : jStem  ? (z >= crownJ->getTime())
                                :          (z <= crownJ->getTime());
            reciprocal = inSub(crown, crownJ) && iPendReachesRj;
        }
        double w = 1.0;
        if(SymmetryCorrection && focalIsTip && j != unresolvedFossils->getSpineIdx() && reciprocal)
            w = 0.5;
        count += w;
    }
    return count;
}

void FBDTreeModel::updateGammaCache(void){
    Tree* tree = parameterTree->getTree();
    int nf = unresolvedFossils->getNumFossils();
    std::vector<Node*>& dpseq = tree->getDownPassSequence();

    sortedYounger.clear();
    sortedOlder.clear();
    for(Node* nd : dpseq){
        if(nd == tree->getCrown())
            continue;
        sortedYounger.push_back(nd->getTime());
        sortedOlder.push_back(nd->getAncestor()->getTime());
    }
    std::sort(sortedYounger.begin(), sortedYounger.end());
    std::sort(sortedOlder.begin(), sortedOlder.end());

    if(cacheInit == false){
        cachedGammaLn.assign(nf, 0.0);
        gammaStale.assign(nf, 1);
        prevY.assign(nf, -1.0);
        prevZ.assign(nf, -1.0);
        prevSA.assign(nf, -1);
        prevNodeAge.assign(tree->getNumNodes(), -1.0);
        prevX0 = -1.0;
        cacheInit = true;
    }

    if(originAge != nullptr){
        double x0 = originAge->getValue();
        if(x0 != prevX0){
            Node* crown = tree->getCrown();
            for(int i = 0; i < nf; i++)
                if(unresolvedFossils->getCrownNode(i) == crown && unresolvedFossils->getIsCrown(i) == false)
                    gammaStale[i] = 1;
            prevX0 = x0;
        }
    }

    for(int i = 0; i < nf; i++){
        double yi = unresolvedFossils->getFossilAge(i);
        double zi = unresolvedFossils->getAttachAge(i);
        int sai = unresolvedFossils->isSA(i) ? 1 : 0;
        if(yi == prevY[i] && zi == prevZ[i] && sai == prevSA[i])
            continue;
        gammaStale[i] = 1;
        bool wasTerm = (prevSA[i] == 0);
        bool isTerm = (sai == 0);
        if((wasTerm || isTerm) && prevY[i] >= 0.0){
            double lo = std::min(prevY[i], yi);
            double hi = -INFINITY;
            if(wasTerm) hi = std::max(hi, prevZ[i]);
            if(isTerm)  hi = std::max(hi, zi);
            for(int j = 0; j < nf; j++){
                if(gammaStale[j]) continue;
                double tj = unresolvedFossils->isSA(j) ? unresolvedFossils->getFossilAge(j) : unresolvedFossils->getAttachAge(j);
                if(tj > lo && tj < hi) gammaStale[j] = 1;
            }
        }
    }

    std::vector<std::pair<double,double> > changedIntervals;
    for(Node* n : dpseq){
        if(n == tree->getCrown()) continue;
        Node* anc = n->getAncestor();
        double pc = prevNodeAge[n->getOffset()];
        double pp = prevNodeAge[anc->getOffset()];
        if(pc < 0.0 || pp < 0.0) continue;
        if(n->getTime() != pc || anc->getTime() != pp)
            changedIntervals.push_back(std::make_pair(std::min(n->getTime(), pc), std::max(anc->getTime(), pp)));
    }
    if(changedIntervals.empty() == false){
        for(int i = 0; i < nf; i++){
            if(gammaStale[i]) continue;
            double ti = unresolvedFossils->isSA(i) ? unresolvedFossils->getFossilAge(i) : unresolvedFossils->getAttachAge(i);
            for(std::pair<double,double>& iv : changedIntervals)
                if(ti > iv.first && ti < iv.second){ gammaStale[i] = 1; break; }
        }
    }

    for(int i = 0; i < nf; i++){
        if(gammaStale[i] == 0) continue;
        double g = computeGamma(unresolvedFossils->getAttachAge(i), i);
        cachedGammaLn[i] = (g > 0.0) ? std::log(g) : -INFINITY;
        gammaStale[i] = 0;
    }

    for(int i = 0; i < nf; i++){
        prevY[i] = unresolvedFossils->getFossilAge(i);
        prevZ[i] = unresolvedFossils->getAttachAge(i);
        prevSA[i] = unresolvedFossils->isSA(i) ? 1 : 0;
    }
    for(Node* n : dpseq)
        prevNodeAge[n->getOffset()] = n->getTime();
}

void FBDTreeModel::computeAgeFloors(std::map<Node*,double>& floors){
    int numFossils = unresolvedFossils->getNumFossils();
    for(int i = 0; i < numFossils; i++){
        if(unresolvedFossils->isSA(i))
            continue;
        Node* node = unresolvedFossils->getMaxAttachNode(i);
        double bound = unresolvedFossils->getAttachAge(i);
        if(unresolvedFossils->getFossilAge(i) > bound)
            bound = unresolvedFossils->getFossilAge(i);
        std::map<Node*,double>::iterator it = floors.find(node);
        if(it == floors.end() || bound > it->second)
            floors[node] = bound;
    }
}

void FBDTreeModel::resolveFossils(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils){
    RandomVariable& rv = RandomVariable::randomVariableInstance();
    for(Fossil& f : fossils){
        Clade* clade = nullptr;
        for(Clade& c : clades)
            if(c.getName() == f.getClade()){
                clade = &c;
                break;
            }
        if(clade == nullptr)
            Msg::error("fossil '" + f.getTaxon() + "' assigned to undefined clade '" + f.getClade() + "'");
        Node* crown = t->getMRCA(clade->getTaxa());
        bool originCond = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);
        Node* origin = originCond ? t->getCrown() : t->getNodeByOffset(clade->getOrigin()->getOffset());
        bool isCrown = (f.getAssignment() == Assignment::CROWN);
        if(crown->getAncestor() == crown)
            isCrown = true;
        double y = 0.5 * (f.getMinAge() + f.getMaxAge());
        std::vector<Node*> hosts;
        std::vector<double> los;
        std::vector<double> his;
        enumerateFossilHosts(t, crown, origin, isCrown, y, hosts, los, his);
        int k = (int)(rv.uniformRv() * hosts.size());
        double z = los[k] + rv.uniformRv() * (his[k] - los[k]);
        t->insertFossilTip(hosts[k], f.getTaxon(), y, z);

        fossilName.push_back(f.getTaxon());
        fossilIsCrown.push_back(isCrown);
        fossilY.push_back(y);
    }
}

void FBDTreeModel::enumerateFossilHosts(Tree* t, Node* crown, Node* origin, bool isCrown, double y, std::vector<Node*>& hosts, std::vector<double>& los, std::vector<double>& his){
    double ceiling = isCrown ? crown->getTime() : origin->getTime();
    for(Node* n : t->getDownPassSequence()){
        if(n == t->getCrown())
            continue;
        Node* anc = n->getAncestor();
        bool inZone = nodeInSubtree(anc, crown);
        if(inZone == false && isCrown == false && nodeOnStalk(n, crown, origin))
            inZone = true;
        if(inZone == false)
            continue;
        double lo = std::max(y, n->getTime());
        double hi = std::min(ceiling, anc->getTime());
        if(lo < hi){
            hosts.push_back(n);
            los.push_back(lo);
            his.push_back(hi);
        }
    }
}

int FBDTreeModel::fossilIndexByName(const std::string& nm){
    for(size_t i = 0; i < fossilName.size(); i++)
        if(fossilName[i] == nm)
            return (int)i;
    return -1;
}

double FBDTreeModel::doWilsonBalding(void){
    Tree* tree = parameterTree->getTree();
    std::vector<Node*> prunable;
    enumeratePrunableRoots(tree, prunable);
    if(prunable.empty())
        return -INFINITY;
    double pFwd = (double)prunable.size();
    Node* r = prunable[(int)(rng.uniformRv() * pFwd)];
    Node* split = r->getAncestor();
    Node* hostParent = split->getAncestor();
    Node* hostChild = nullptr;
    for(Node* nb : split->getNeighbors())
        if(nb != hostParent && nb != r)
            hostChild = nb;

    std::vector<Node*> crowns;
    std::vector<char> isCrowns;
    std::vector<Node*> origins;
    double ceilingS = std::numeric_limits<double>::max();
    std::vector<Node*> stack;
    stack.push_back(r);
    while(stack.empty() == false){
        Node* x = stack.back();
        stack.pop_back();
        if(x->getIsTip()){
            int idx = fossilIndexByName(x->getName());
            crowns.push_back(fossilCrown[idx]);
            isCrowns.push_back(fossilIsCrown[idx] ? 1 : 0);
            origins.push_back(fossilOrigin[idx]);
            double c = fossilIsCrown[idx] ? fossilCrown[idx]->getTime() : fossilOrigin[idx]->getTime();
            if(c < ceilingS)
                ceilingS = c;
        }else{
            for(Node* nb : x->getNeighbors())
                if(nb != x->getAncestor())
                    stack.push_back(nb);
        }
    }

    double rAge = r->getTime();
    double oldRange = tree->isSATip(hostChild) ? 1.0 : std::min(ceilingS, hostParent->getTime()) - std::max(rAge, hostChild->getTime());

    hostParent->removeNeighbor(split);
    split->removeNeighbor(hostParent);
    hostChild->removeNeighbor(split);
    split->removeNeighbor(hostChild);
    hostParent->addNeighbor(hostChild);
    hostChild->addNeighbor(hostParent);
    hostChild->setAncestor(hostParent);
    tree->initializeDownPassSequence();

    std::vector<Node*> hosts;
    std::vector<double> los;
    std::vector<double> his;
    enumerateSubtreeHosts(tree, crowns, isCrowns, origins, rAge, ceilingS, hosts, los, his);
    int k = (int)(rng.uniformRv() * (double)hosts.size());
    Node* newChild = hosts[k];
    Node* newParent = newChild->getAncestor();
    double newRange = (los[k] == his[k]) ? 1.0 : his[k] - los[k];
    double z = (los[k] == his[k]) ? los[k] : los[k] + rng.uniformRv() * (his[k] - los[k]);

    newParent->removeNeighbor(newChild);
    newChild->removeNeighbor(newParent);
    newParent->addNeighbor(split);
    split->addNeighbor(newParent);
    split->setAncestor(newParent);
    newChild->addNeighbor(split);
    split->addNeighbor(newChild);
    newChild->setAncestor(split);
    split->setTime(z);
    tree->initializeDownPassSequence();
    tree->reindexNodes();

    std::vector<Node*> prunable2;
    enumeratePrunableRoots(tree, prunable2);
    double pRev = (double)prunable2.size();

    return std::log(newRange / oldRange) + std::log(pFwd / pRev);
}

double FBDTreeModel::doNarrowExchange(void){
    Tree* tree = parameterTree->getTree();
    const std::vector<Node*>& dp = tree->getDownPassSequence();
    Node* i = dp[(int)(rng.uniformRv() * (double)dp.size())];
    Node* parent = i->getAncestor();
    if(i == tree->getCrown() || parent == tree->getCrown())
        return -INFINITY;
    Node* grandparent = parent->getAncestor();
    Node* uncle = nullptr;
    for(Node* nb : grandparent->getNeighbors())
        if(nb != parent && nb != grandparent->getAncestor())
            uncle = nb;
    if(uncle == nullptr)
        return -INFINITY;
    if(uncle->getTime() >= parent->getTime())
        return -INFINITY;
    if(tree->isSATip(i) || tree->isSATip(uncle))
        return -INFINITY;
    if(subtreeAllFossil(i) == false || subtreeAllFossil(uncle) == false)
        return -INFINITY;
    if(subtreeFossilsValidAt(tree, i, grandparent) == false)
        return -INFINITY;
    if(subtreeFossilsValidAt(tree, uncle, parent) == false)
        return -INFINITY;

    grandparent->removeNeighbor(uncle);
    uncle->removeNeighbor(grandparent);
    parent->removeNeighbor(i);
    i->removeNeighbor(parent);
    grandparent->addNeighbor(i);
    i->addNeighbor(grandparent);
    i->setAncestor(grandparent);
    parent->addNeighbor(uncle);
    uncle->addNeighbor(parent);
    uncle->setAncestor(parent);
    tree->initializeDownPassSequence();
    tree->reindexNodes();
    return 0.0;
}

double FBDTreeModel::doWideExchange(void){
    Tree* tree = parameterTree->getTree();
    const std::vector<Node*>& dp = tree->getDownPassSequence();
    double n = (double)dp.size();
    Node* i = dp[(int)(rng.uniformRv() * n)];
    Node* j = dp[(int)(rng.uniformRv() * n)];
    if(i == j)
        return -INFINITY;
    Node* pi = i->getAncestor();
    Node* pj = j->getAncestor();
    if(i == tree->getCrown() || j == tree->getCrown() || pi == tree->getCrown() || pj == tree->getCrown())
        return -INFINITY;
    if(pi == pj)
        return -INFINITY;
    if(nodeInSubtree(i, j) || nodeInSubtree(j, i))
        return -INFINITY;
    if(i->getTime() >= pj->getTime() || j->getTime() >= pi->getTime())
        return -INFINITY;
    if(tree->isSATip(i) || tree->isSATip(j))
        return -INFINITY;
    if(subtreeAllFossil(i) == false || subtreeAllFossil(j) == false)
        return -INFINITY;
    if(subtreeFossilsValidAt(tree, i, pj) == false || subtreeFossilsValidAt(tree, j, pi) == false)
        return -INFINITY;

    pi->removeNeighbor(i);
    i->removeNeighbor(pi);
    pj->removeNeighbor(j);
    j->removeNeighbor(pj);
    pi->addNeighbor(j);
    j->addNeighbor(pi);
    j->setAncestor(pi);
    pj->addNeighbor(i);
    i->addNeighbor(pj);
    i->setAncestor(pj);
    tree->initializeDownPassSequence();
    tree->reindexNodes();
    return 0.0;
}

double FBDTreeModel::doTreeScale(void){
    Tree* tree = parameterTree->getTree();
    double m = std::exp(parameterTree->getScaleLambda() * (rng.uniformRv() - 0.5));
    int numScaled = tree->scaleInternalAges(m);
    for(Node* n : tree->getDownPassSequence())
        if(n != tree->getCrown() && tree->isSATip(n) == false && n->getTime() >= n->getAncestor()->getTime())
            return -INFINITY;
    return numScaled * std::log(m);
}

double FBDTreeModel::doSARJMCMC(void){
    Tree* tree = parameterTree->getTree();
    std::vector<Node*> fossils;
    for(Node* n : tree->getDownPassSequence())
        if(n->getIsTip() && n->getIsFossil())
            fossils.push_back(n);
    if(fossils.empty())
        return -INFINITY;
    Node* f = fossils[(int)(rng.uniformRv() * fossils.size())];
    Node* sp = f->getAncestor();
    bool useOrigin = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);
    bool spIsCrown = (sp == tree->getCrown());
    if(spIsCrown && useOrigin == false)
        return -INFINITY;
    Node* gp = sp->getAncestor();
    Node* sib = nullptr;
    for(Node* c : sp->getNeighbors())
        if(c != gp && c != f)
            sib = c;
    double y = f->getTime();
    double maxAge = spIsCrown ? originAge->getValue() : gp->getTime();
    double range = maxAge - y;
    if(range <= 0.0)
        return -INFINITY;
    if(sp->getTime() == y){
        sp->setTime(y + rng.uniformRv() * range);
        return std::log(range);
    }
    if(y < sib->getTime())
        return -INFINITY;
    if(sib->getIsTip() && sib->getIsFossil() && sib->getTime() == y)
        return -INFINITY;
    sp->setTime(y);
    return -std::log(range);
}

double FBDTreeModel::doUpDownScale(void){
    lastWasUpDown = true;
    double ar = 0.0;
    for(bool b : upDownRecent)
        if(b)
            ar++;
    if(upDownRecent.empty() == false)
        ar /= upDownRecent.size();
    if(upDownTotal > 0 && upDownTotal % 100 == 0)
        upDownStep *= std::exp((1.0 / std::sqrt((double)(upDownTotal / 100))) * (ar - 0.3));

    double mBact = 0.95;
    double sBact = std::sqrt(1.0 - mBact * mBact);
    double delta = mBact + Probability::Normal::rv(&rng) * sBact;
    if(Probability::Uniform::rv(&rng, 0.0, 1.0) < 0.5)
        delta = -delta;
    double c = std::exp(upDownStep * delta);
    double invc = 1.0 / c;
    double lnc = std::log(c);

    int nUp = 0;
    for(ParameterDouble* l : lambda){ l->scaleProposed(c); nUp++; }
    for(ParameterDouble* m : mu){ m->scaleProposed(c); nUp++; }
    for(ParameterDouble* p : psi){ p->scaleProposed(c); nUp++; }
    double h = nUp * lnc;

    Tree* t = parameterTree->getTree();
    t->setLastUpdateWasScale(true);
    if(unresolvedFossils != nullptr){
        double zJac = unresolvedFossils->scaleAllAttachAges(invc);
        if(zJac == -INFINITY)
            return -INFINITY;
        h += zJac;
    }
    int numScaledDown = t->scaleInternalAges(invc);
    h += numScaledDown * std::log(invc);
    if(originAge != nullptr){
        originAge->scaleProposed(invc);
        h += std::log(invc);
    }
    for(Node* n : t->getDownPassSequence())
        if(n != t->getCrown() && n->getAncestor()->getTime() < n->getTime())
            return -INFINITY;
    return h;
}

double FBDTreeModel::doClockNodeAge(void){
    setupNodeAgeFloors();
    return parameterTree->getTree()->updateNodeAge();
}

void FBDTreeModel::setupNodeAgeFloors(void){
    std::map<Node*,double> floors;
    computeAgeFloors(floors);
    parameterTree->getTree()->setAgeFloors(floors);
    parameterTree->getTree()->setLastUpdateWasScale(false);
}

double FBDTreeModel::doJointScale(void){
    double m = std::exp(parameterTree->getScaleLambda() * (rng.uniformRv() - 0.5));
    parameterTree->getTree()->setLastUpdateWasScale(true);
    double zJac = unresolvedFossils->scaleAllAttachAges(m);
    if(zJac == -INFINITY)
        return -INFINITY;
    int numScaled = parameterTree->getTree()->scaleInternalAges(m);
    return numScaled * std::log(m) + zJac;
}

double FBDTreeModel::doSubtreeScale(void){
    Tree* tree = parameterTree->getTree();
    tree->setLastUpdateWasScale(false);

    std::vector<Node*> candidates;
    for(Node* n : tree->getDownPassSequence())
        if(n != tree->getCrown() && n->getIsTip() == false)
            candidates.push_back(n);
    if(candidates.empty()){
        unresolvedFossils->scaleAttachAges(std::vector<int>(), 1.0);
        return 0.0;
    }

    Node* node = candidates[(int)(rng.uniformRv() * candidates.size())];
    double oldAge = node->getTime();
    double parentAge = node->getAncestor()->getTime();
    double oldestTip = 0.0;
    for(Node* d : tree->getAllDescendants(node))
        if(d->getIsTip() && d->getTime() > oldestTip)
            oldestTip = d->getTime();
    double sf = (oldestTip + rng.uniformRv() * (parentAge - oldestTip)) / oldAge;

    std::vector<int> insideZ;
    int nf = unresolvedFossils->getNumFossils();
    for(int i = 0; i < nf; i++){
        if(unresolvedFossils->isSA(i))
            continue;
        if(nodeInSubtree(unresolvedFossils->getMaxAttachNode(i), node))
            insideZ.push_back(i);
    }

    double zJac = unresolvedFossils->scaleAttachAges(insideZ, sf);
    if(zJac == -INFINITY)
        return -INFINITY;
    int numScaled = tree->scaleSubtreeAges(node, sf);
    return (numScaled - 1) * std::log(sf) + zJac;
}

void FBDTreeModel::enumeratePrunableRoots(Tree* t, std::vector<Node*>& roots){
    Node* treeCrown = t->getCrown();
    std::set<Node*> allFossil;
    for(Node* n : t->getDownPassSequence()){
        bool af;
        if(n->getIsTip()){
            af = n->getIsFossil();
        }else{
            af = true;
            for(Node* nb : n->getNeighbors()){
                if(nb == n->getAncestor())
                    continue;
                if(allFossil.find(nb) == allFossil.end()){
                    af = false;
                    break;
                }
            }
        }
        if(af){
            allFossil.insert(n);
            if(n->getAncestor() != treeCrown && t->isSATip(n) == false)
                roots.push_back(n);
        }
    }
}

void FBDTreeModel::enumerateSubtreeHosts(Tree* t, std::vector<Node*>& crowns, std::vector<char>& isCrowns, std::vector<Node*>& origins, double rAge, double ceilingS, std::vector<Node*>& hosts, std::vector<double>& los, std::vector<double>& his){
    for(Node* n : t->getDownPassSequence()){
        if(n == t->getCrown())
            continue;
        Node* anc = n->getAncestor();
        bool allInZone = true;
        for(size_t f = 0; f < crowns.size(); f++){
            bool inZone = nodeInSubtree(anc, crowns[f]);
            if(inZone == false && isCrowns[f] == 0 && nodeOnStalk(n, crowns[f], origins[f]))
                inZone = true;
            if(inZone == false){
                allInZone = false;
                break;
            }
        }
        if(allInZone == false)
            continue;
        double lo = std::max(rAge, n->getTime());
        double hi = std::min(ceilingS, anc->getTime());
        if(lo < hi){
            hosts.push_back(n);
            los.push_back(lo);
            his.push_back(hi);
        }
        if(n->getIsTip() && n->getIsFossil() && t->isSATip(n) == false && rAge < n->getTime() && n->getTime() < ceilingS){
            hosts.push_back(n);
            los.push_back(n->getTime());
            his.push_back(n->getTime());
        }
    }
}

bool FBDTreeModel::subtreeFossilsValidAt(Tree* t, Node* s, Node* g){
    std::vector<Node*> stack;
    stack.push_back(s);
    while(stack.empty() == false){
        Node* x = stack.back();
        stack.pop_back();
        if(x->getIsTip()){
            int idx = fossilIndexByName(x->getName());
            Node* crown = fossilCrown[idx];
            bool inZone = nodeInSubtree(g, crown);
            if(inZone == false && fossilIsCrown[idx] == false && nodeOnStalk(g, crown, fossilOrigin[idx]))
                inZone = true;
            if(inZone == false)
                return false;
        }else{
            for(Node* nb : x->getNeighbors())
                if(nb != x->getAncestor())
                    stack.push_back(nb);
        }
    }
    return true;
}

bool FBDTreeModel::subtreeAllFossil(Node* n){
    if(n->getIsTip())
        return n->getIsFossil();
    for(Node* nb : n->getNeighbors())
        if(nb != n->getAncestor() && subtreeAllFossil(nb) == false)
            return false;
    return true;
}
