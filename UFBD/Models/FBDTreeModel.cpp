#include "Node.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "ParameterShrinkageField.hpp"
#include "FBDTreeModel.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "Serialize.hpp"
#include "ThreadPool.hpp"
#include "UserSettings.hpp"
#include "Probability.hpp"
#include "Msg.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>

FBDTreeModel::FBDTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, unsigned int seed) :
    PhylogeneticModel(){

    lastWasJointScale = false;
    lastWasUpDown = false;
    lastWasRateVec = false;
    lastRateVec = nullptr;
    rateVecStep = 0.2;
    shrinkStep = 0.2;
    rvAccW = rvAttW = seAccW = seAttW = 0;
    lastWasFbdRate = false;
    lastWasJointRate = false;
    jointRateParam = nullptr;
    jointRateStep = 0.5;
    jrAttW = 0;
    jrAccW = 0;
    turnoverStep = 0.1;
    frAccW = frAttW = 0;
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
    isResolved = (UserSettings::userSettings().getModel() == Model::RFBD);

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

    if(isResolved){
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
    
    UserSettings& rateUs = UserSettings::userSettings();
    std::vector<double> lambdaTimes, muTimes, psiTimes;
    lambdaTimes.push_back(0.0);
    for(double t : rateUs.getLambdaSkylineTimes()) lambdaTimes.push_back(t);
    muTimes.push_back(0.0);
    for(double t : rateUs.getMuSkylineTimes())     muTimes.push_back(t);
    psiTimes.push_back(0.0);
    for(double t : rateUs.getPsiSkylineTimes())    psiTimes.push_back(t);
    intervalStart.push_back(0.0);
    for(double t : rateUs.getSkylineTimes())
        intervalStart.push_back(t);
    for(double s : intervalStart){
        int li = 0, mi = 0, pi = 0;
        for(int k = 0; k < (int)lambdaTimes.size(); k++) if(lambdaTimes[k] <= s) li = k;
        for(int k = 0; k < (int)muTimes.size(); k++)     if(muTimes[k] <= s)     mi = k;
        for(int k = 0; k < (int)psiTimes.size(); k++)    if(psiTimes[k] <= s)    pi = k;
        lambdaIdx.push_back(li);
        muIdx.push_back(mi);
        psiIdx.push_back(pi);
    }
    Probability::PriorSpec lp = rateUs.getLambdaPrior();
    Probability::PriorSpec mp = rateUs.getMuPrior();
    Probability::PriorSpec pp = rateUs.getPsiPrior();
    Probability::PriorSpec defRate{true, Probability::PriorFamily::EXPONENTIAL, 10.0, 1.0};
    if(!lp.set) lp = defRate;
    if(!mp.set) mp = defRate;
    if(!pp.set) pp = defRate;
    int nLambda = (int)lambdaTimes.size();
    int nMu = (int)muTimes.size();
    int nPsi = (int)psiTimes.size();
    bool lamSmooth = (rateUs.getLambdaMode() == RateMode::SMOOTH);
    bool muSmooth  = (rateUs.getMuMode() == RateMode::SMOOTH);
    bool psiSmooth = (rateUs.getPsiMode() == RateMode::SMOOTH);
    if(lamSmooth && nLambda < 2){ Msg::warning("Smoothing (HSMRF) set for speciation rate (lambda) but only single rate interval."); lamSmooth = false; }
    if(muSmooth && nMu < 2){ Msg::warning("Smoothing (HSMRF) set for extinction rate (mu) but only single rate interval."); muSmooth = false; }
    if(psiSmooth && nPsi < 2){ Msg::warning("Smoothing (HSMRF) set for sampling rate (psi) but only single rate interval."); psiSmooth = false; }
    double nShifts = rateUs.getHsmrfShifts();
    double shiftSize = rateUs.getHsmrfShiftSize();
    double lam0 = Probability::priorMean(lp.family, lp.p1, lp.p2);
    if(!(lam0 > 0.0 && std::isfinite(lam0))) lam0 = 0.1;
    double mu0 = Probability::priorMean(mp.family, mp.p1, mp.p2);
    if(!(mu0 > 0.0 && std::isfinite(mu0))) mu0 = 0.1;
    if(mu0 == lam0) mu0 = 0.5 * lam0;
    double psi0 = Probability::priorMean(pp.family, pp.p1, pp.p2);
    if(!(psi0 > 0.0 && std::isfinite(psi0))) psi0 = 0.1;
    if(psi0 == lam0) psi0 = 0.5 * lam0;
    if(lamSmooth){
        lambdaField = new ParameterShrinkageField(1.0, this, nLambda, lp, nShifts, shiftSize, lam0);
        parameters.push_back(lambdaField);
    }else{
        for(int i = 0; i < nLambda; i++){
            std::string suf = (nLambda > 1) ? std::to_string(i) : "";
            ParameterDouble* l = new ParameterDouble(1.0, this, "lambda" + suf, 0.0, std::numeric_limits<double>::max());
            l->setPrior(lp.family, lp.p1, lp.p2);
            l->setValue(lam0);
            lambda.push_back(l);
            parameters.push_back(l);
        }
    }
    if(muSmooth){
        muField = new ParameterShrinkageField(1.0, this, nMu, mp, nShifts, shiftSize, mu0);
        parameters.push_back(muField);
    }else{
        for(int i = 0; i < nMu; i++){
            std::string suf = (nMu > 1) ? std::to_string(i) : "";
            ParameterDouble* m = new ParameterDouble(1.0, this, "mu" + suf, 0.0, std::numeric_limits<double>::max());
            m->setPrior(mp.family, mp.p1, mp.p2);
            m->setValue(mu0);
            mu.push_back(m);
            parameters.push_back(m);
        }
    }
    if(psiSmooth){
        psiField = new ParameterShrinkageField(1.0, this, nPsi, pp, nShifts, shiftSize, psi0);
        parameters.push_back(psiField);
    }else{
        for(int i = 0; i < nPsi; i++){
            std::string suf = (nPsi > 1) ? std::to_string(i) : "";
            ParameterDouble* p = new ParameterDouble(1.0, this, "psi" + suf, 0.0, std::numeric_limits<double>::max());
            p->setPrior(pp.family, pp.p1, pp.p2);
            p->setValue(psi0);
            psi.push_back(p);
            parameters.push_back(p);
        }
    }
    for(size_t i = 1; i < lambda.size(); i++) lambda[i]->setValue(lambda[0]->getValue());
    for(size_t i = 1; i < mu.size(); i++)     mu[i]->setValue(mu[0]->getValue());
    for(size_t i = 1; i < psi.size(); i++)    psi[i]->setValue(psi[0]->getValue());
    rho = UserSettings::userSettings().getRho();

    unresolvedFossils = nullptr;
    if(isResolved){
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
            fossilOrigin.push_back(originCond ? wt->getRoot() : wt->getNodeByOffset(clade->getOrigin()->getOffset()));
        }
    }else{
        unresolvedFossils = new ParameterUnresolvedFossils(1.0, this, parameterTree->getTree(), clades, fossils, originAge);
        parameters.push_back(unresolvedFossils);
    }

    if(isResolved){
        parameterTree->setProposalProbability(78.0);
        for(ParameterDouble* l : lambda) l->setProposalProbability(15.0);
        for(ParameterDouble* m : mu)     m->setProposalProbability(15.0);
        for(ParameterDouble* p : psi)    p->setProposalProbability(15.0);
    }
    double fieldBase = (isResolved ? 15.0 : 1.0);
    if(lambdaField) lambdaField->setProposalProbability(fieldBase * (double)nLambda);
    if(muField)     muField->setProposalProbability(fieldBase * (double)nMu);
    if(psiField)    psiField->setProposalProbability(fieldBase * (double)nPsi);

    double sum = 0.0;
    for(Parameter* p : parameters)
        sum += p->getProposalProbability();
    for(Parameter* p : parameters)
        p->setProposalProbability(p->getProposalProbability() / sum);

    RandomVariable::setActiveInstance(prevRng);
}

double FBDTreeModel::lambdaAt(int i){
    int j = lambdaIdx[i];
    if(lambdaField != nullptr)
        return lambdaField->getRate(j);
    return lambda[j]->getValue();
}

double FBDTreeModel::muAt(int i){
    int j = muIdx[i];
    if(muField != nullptr)
        return muField->getRate(j);
    return mu[j]->getValue();
}

double FBDTreeModel::psiAt(int i){
    int j = psiIdx[i];
    if(psiField != nullptr)
        return psiField->getRate(j);
    return psi[j]->getValue();
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
    std::vector<Node*> bbNodes = parameterTree->getTree()->getBackboneAgeNodes();
    for(size_t i = 0; i < bbNodes.size(); i++)
        names.push_back("x" + std::to_string(i + 1));
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
    if(isResolved)
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
    for(Node* n : parameterTree->getTree()->getBackboneAgeNodes())
        vals.push_back(n->getTime());
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
    if(isResolved)
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
        if(p != parameterTree)
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

double FBDTreeModel::doTurnoverMove(void){
    lastWasFbdRate = true;
    double lnH = 0.0;
    for(size_t u = 0; u < intervalStart.size(); u++){
        int i = lambdaIdx[u];
        int j = muIdx[u];
        double lam = lambda[i]->getValue();
        double muv = mu[j]->getValue();
        double d = lam - muv;
        if(d <= 0.0)
            return -std::numeric_limits<double>::infinity();
        double t = muv / lam;
        double tNew = t + turnoverStep * (rng.uniformRv() - 0.5);
        while(tNew <= 0.0 || tNew >= 1.0){
            if(tNew <= 0.0) tNew = -tNew;
            if(tNew >= 1.0) tNew = 2.0 - tNew;
        }
        double lamNew = d / (1.0 - tNew);
        double muNew = tNew * lamNew;
        lambda[i]->scaleProposed(lamNew / lam);
        mu[j]->scaleProposed(muNew / muv);
        lnH += 2.0 * std::log(lamNew / lam);
    }
    frAttW++;
    if(frAttW >= 200){
        turnoverStep *= std::exp((double)frAccW / frAttW - 0.3);
        if(turnoverStep < 1e-3) turnoverStep = 1e-3;
        if(turnoverStep > 0.9)  turnoverStep = 0.9;
        frAccW = 0;
        frAttW = 0;
    }
    return lnH;
}

// straddle count: intervals (lo<hi) with lo<zq<hi = #{lo<zq} - #{hi<=zq}
namespace {
int countStraddling(const std::vector<double>& los, const std::vector<double>& his, double zq){
    int below = (int)(std::lower_bound(los.begin(), los.end(), zq) - los.begin());
    int above = (int)(std::upper_bound(his.begin(), his.end(), zq) - his.begin());
    return below - above;
}
// pendants as (z_attach,y) sorted by z_attach; y passed separately (sorted)
int countStraddling(const std::vector<double>& ySorted, const std::vector<std::pair<double,double>>& zy, double zq){
    int below = (int)(std::lower_bound(ySorted.begin(), ySorted.end(), zq) - ySorted.begin());
    int above = (int)(std::upper_bound(zy.begin(), zy.end(), zq,
                     [](double v, const std::pair<double,double>& p){ return v < p.first; }) - zy.begin());
    return below - above;
}
// straddlers with z_attach rootward of T (leave a crown/sub zone)
double countAboveCrown(const std::vector<std::pair<double,double>>& zy, double T, double zq){
    double s = 0.0;
    for(std::vector<std::pair<double,double>>::const_iterator it = std::upper_bound(zy.begin(), zy.end(), std::pair<double,double>(T, INFINITY));
             it != zy.end(); ++it)
        if(it->second < zq && zq < it->first)
            s += 1.0;
    return s;
}
}

// backbone lineages crossing z in fossil i's zone R_i
double FBDTreeModel::cladeBackboneLineages(int i, double z){
    Tree* tree = parameterTree->getTree();
    if(eulerBuilt == false)
        buildEulerIndex();
    Node* crown = unresolvedFossils->getCrownNode(i);
    bool stem = unresolvedFossils->getIsStem(i);
    bool total = (unresolvedFossils->getIsCrown(i) == false && stem == false);
    double count = 0.0;
    if(stem == false && crown == tree->getCrown()){
        count += (double)countStraddling(sortedYounger, sortedOlder, z);   // all edges in-zone
    }else{
        int lo = subPre[crown->getOffset()];
        int hi = lo + subSize[crown->getOffset()];
        for(int k = lo; k < hi; k++){
            Node* n = nodesByPre[k];
            if(n == tree->getRoot())
                continue;
            Node* anc = n->getAncestor();
            if(n->getTime() < z && z < anc->getTime()){
                bool inZone;
                if(stem){
                    inZone = (n == crown);              // stem edge only
                }else{
                    inZone = inSub(anc, crown);         // crown group
                    if(inZone == false && total && n == crown)
                        inZone = true;                  // total owns clade stem
                }
                if(inZone)
                    count++;
            }
        }
    }
    if((total || stem) && crown == tree->getCrown() && originAge != nullptr){   // trunk crown->origin
        double x0 = originAge->getValue();
        if(tree->getNumBackbone() == 0){
            if(z >= x0)
                count++;                                // empty backbone: origin spine
        }else{
            if(tree->getCrown()->getTime() < z && z < x0)
                count++;
        }
    }
    return count;
}

double FBDTreeModel::slabAntideriv(double z, int k){
    double w = (1.0 + c2Vec[k]) * std::exp(c1Vec[k]*(z - intervalStart[k])) + (1.0 - c2Vec[k]);
    return -8.0 * lambdaAt(k) * std::exp(lnDPrev[k]) / (c1Vec[k]*(1.0 + c2Vec[k])) * (1.0/w);
}

double FBDTreeModel::jointRateFossilProposal(int i, bool doDraw, bool& saOut, double& zOut){
    double y = unresolvedFossils->getFossilAge(i);
    double lo = unresolvedFossils->getMinAttachAge(i);
    double hi = unresolvedFossils->getMaxAttachAge(i);
    std::vector<double> bpts;
    bpts.push_back(lo); bpts.push_back(hi);
    for(double a : sortedYounger) if(a > lo && a < hi) bpts.push_back(a);
    for(double a : sortedOlder)   if(a > lo && a < hi) bpts.push_back(a);
    for(double s : intervalStart) if(s > lo && s < hi) bpts.push_back(s);
    std::sort(bpts.begin(), bpts.end());
    bpts.erase(std::unique(bpts.begin(), bpts.end()), bpts.end());

    double pfac = unresolvedFossils->isUE(i) ? 1.0 : calculateP0(y) / std::exp(lnD(y));
    std::vector<double> mass, pa, pb;
    std::vector<int> pk;
    double sumMass = 0.0;
    for(size_t j = 0; j + 1 < bpts.size(); j++){
        double a = bpts[j], b = bpts[j+1];
        if(b - a <= 0.0) continue;
        int k = findIndex(0.5*(a+b));
        double g = cladeBackboneLineages(i, 0.5*(a+b));
        double m = g * (slabAntideriv(b,k) - slabAntideriv(a,k));
        if(m < 0.0) m = 0.0;
        mass.push_back(m); pa.push_back(a); pb.push_back(b); pk.push_back(k);
        sumMass += m;
    }
    double atom = (unresolvedFossils->isUE(i) || lo > y) ? 0.0 : cladeBackboneLineages(i, y);
    double slab = pfac * sumMass;
    double tot = atom + slab;

    if(doDraw == false){
        if(unresolvedFossils->isSA(i))
            return std::log(atom / tot);
        double z = unresolvedFossils->getAttachAge(i);
        double dens = pfac * 2.0 * lambdaAt(findIndex(z)) * std::exp(lnD(z)) * cladeBackboneLineages(i, z);
        return std::log(dens / tot);
    }
    if(rng.uniformRv() < atom / tot){
        saOut = true; zOut = y;
        return std::log(atom / tot);
    }
    double u = rng.uniformRv() * sumMass;
    double acc = 0.0;
    int sel = (int)mass.size() - 1;
    for(size_t j = 0; j < mass.size(); j++){ acc += mass[j]; if(u <= acc){ sel = (int)j; break; } }
    int k = pk[sel];
    double Fa = slabAntideriv(pa[sel], k), Fb = slabAntideriv(pb[sel], k);
    double Ft = Fa + rng.uniformRv() * (Fb - Fa);
    double C = 8.0 * lambdaAt(k) * std::exp(lnDPrev[k]) / (c1Vec[k]*(1.0 + c2Vec[k]));
    double w = -C / Ft;
    double z = intervalStart[k] + std::log((w - (1.0 - c2Vec[k]))/(1.0 + c2Vec[k])) / c1Vec[k];
    if(z <= pa[sel]) z = pa[sel] + 1e-9;
    if(z >= pb[sel]) z = pb[sel] - 1e-9;
    saOut = false; zOut = z;
    double dens = pfac * 2.0 * lambdaAt(findIndex(z)) * std::exp(lnD(z)) * cladeBackboneLineages(i, z);
    return std::log(dens / tot);
}

void FBDTreeModel::adaptJointRate(void){
    if(jrAttW < 200)
        return;
    double ar = (double)jrAccW / (double)jrAttW;
    jointRateStep *= std::exp(ar - 0.3);
    if(jointRateStep < 1e-3) jointRateStep = 1e-3;
    if(jointRateStep > 10.0) jointRateStep = 10.0;
    jrAccW = 0;
    jrAttW = 0;
}

double FBDTreeModel::jointRateFossilMove(void){
    std::vector<ParameterDouble*> allr;
    for(ParameterDouble* p : lambda) allr.push_back(p);
    for(ParameterDouble* p : mu)     allr.push_back(p);
    for(ParameterDouble* p : psi)    allr.push_back(p);
    if(allr.empty())
        return -std::numeric_limits<double>::infinity();
    prepareIntervals();
    updateGammaCache();
    int nf = unresolvedFossils->getNumFossils();
    bool s; double z;
    int spine = unresolvedFossils->getSpineIdx();
    double logqOld = 0.0;
    for(int i = 0; i < nf; i++)
        if(i != spine)
            logqOld += jointRateFossilProposal(i, false, s, z);
    ParameterDouble* p = allr[(int)(rng.uniformRv() * allr.size())];
    double c = std::exp(jointRateStep * (rng.uniformRv() - 0.5));
    p->scaleProposed(c);
    prepareIntervals();
    unresolvedFossils->beginBulkMove();
    double logqNew = 0.0;
    for(int i = 0; i < nf; i++){
        if(i == spine) continue;
        logqNew += jointRateFossilProposal(i, true, s, z);
        unresolvedFossils->setProposedAttach(i, s ? unresolvedFossils->getFossilAge(i) : z);
    }
    lastWasJointRate = true;
    jointRateParam = p;
    return std::log(c) + (logqOld - logqNew);
}

double FBDTreeModel::update(void){
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);

    lastWasJointScale = false;
    lastWasUpDown = false;
    lastWasRateVec = false;
    lastWasJointRate = false;
    bool haveIid = (lambdaField == nullptr && lambda.size() >= 2)
                || (muField == nullptr && mu.size() >= 2)
                || (psiField == nullptr && psi.size() >= 2);
    static int jointRateGate = [](){ const char* e = getenv("FBD_JOINT_RATE"); return e ? atoi(e) : 0; }();
    if(jointRateGate && haveIid && unresolvedFossils != nullptr && rng.uniformRv() < 0.05 + 0.25 * std::min(jointRateStep / 0.3, 1.0)){
        double r = jointRateFossilMove();
        RandomVariable::setActiveInstance(prevRng);
        return r;
    }
    if(haveIid && rng.uniformRv() < 0.20){
        double r = (rng.uniformRv() < 0.5) ? doRateVectorScale() : doRateShrinkExpand();
        RandomVariable::setActiveInstance(prevRng);
        return r;
    }
    lastWasFbdRate = false;
    if(lambdaField == nullptr && muField == nullptr && lambda.size() >= 1 && mu.size() >= 1 && rng.uniformRv() < 0.10){
        double r = doTurnoverMove();
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
    if(updatedParameter == parameterTree && isResolved == false && unresolvedFossils != nullptr){
        Tree* t = parameterTree->getTree();
        int numSlideable = 0;
        for(Node* n : t->getDownPassSequence())
            if(n != t->getRoot() && n->getIsTip() == false)
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
    else if(updatedParameter == parameterTree && isResolved){
        parameterTree->getTree()->setLastUpdateWasScale(false);
        double wNE = 10.0;
        double wWB = 10.0;
        double wWE = 10.0;
        double wTreeScale = 3.0;
        double wSA = 10.0;
        double wUpDown = 10.0;
        double wAge = 20.0;
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
    if(lastWasJointRate){
        jointRateParam->commitProposed();
        unresolvedFossils->updateForAcceptance();
        jrAccW++; jrAttW++; adaptJointRate();
        return;
    }
    if(lastWasFbdRate){
        for(ParameterDouble* l : lambda) l->commitProposed();
        for(ParameterDouble* m : mu) m->commitProposed();
        frAccW++;
        return;
    }
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
    if(isResolved && originAge != nullptr)
        parameterTree->getTree()->getRoot()->setTime(originAge->getValue());
}

void FBDTreeModel::updateForRejection(void){
    if(lastWasJointRate){
        jointRateParam->restoreProposed();
        unresolvedFossils->updateForRejection();
        jrAttW++; adaptJointRate();
        return;
    }
    if(lastWasFbdRate){
        for(ParameterDouble* l : lambda) l->restoreProposed();
        for(ParameterDouble* m : mu) m->restoreProposed();
        return;
    }
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
    if(isResolved && originAge != nullptr)
        parameterTree->getTree()->getRoot()->setTime(originAge->getValue());
}

void FBDTreeModel::writeState(std::ostream& os){
    for(Parameter* p : parameters)
        p->writeState(os);
    os << rateVecStep << ' ' << shrinkStep << ' ' << turnoverStep << ' ' << upDownStep << ' ' << upDownTotal << '\n';
    os << rvAccW << ' ' << rvAttW << ' ' << seAccW << ' ' << seAttW << ' ' << frAccW << ' ' << frAttW << '\n';
    Serialize::writeBoolDeque(os, upDownRecent);
}

void FBDTreeModel::readState(std::istream& is){
    for(Parameter* p : parameters)
        p->readState(is);
    is >> rateVecStep >> shrinkStep >> turnoverStep >> upDownStep >> upDownTotal;
    is >> rvAccW >> rvAttW >> seAccW >> seAttW >> frAccW >> frAttW;
    Serialize::readBoolDeque(is, upDownRecent);
    cacheInit = false;
}

double FBDTreeModel::calculateFBDProbability(void){
    Tree* tree = parameterTree->getTree();

    if(isResolved && originAge != nullptr)
        tree->getRoot()->setTime(originAge->getValue());

    int numInternalNodes = tree->getNumNodes() - tree->getNumBackbone();
    double crownAge = tree->getCrown()->getTime();
    std::vector<Node*> dpseq = tree->getDownPassSequence();
    
    rhoVal = rho;

    prepareIntervals();

    if(isResolved)
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
    int nDp = (int)dpseq.size();
    std::vector<double> termNode(nDp, 0.0);
    ThreadPool::current().parallelFor(OP_FBD, nDp, [&](int lo, int hi){
        for(int idx = lo; idx < hi; idx++){
            Node* n = dpseq[idx];
            if(n->getIsTip())
                continue;
            bool hasChild = false;
            for(Node* c : n->getNeighbors())
                if(c != n->getAncestor()){ hasChild = true; break; }
            if(hasChild)
                termNode[idx] = std::log(lambdaAt(findIndex(n->getTime())) * rhoVal) + lnD(n->getTime());
        }
    });
    for(int idx = 0; idx < nDp; idx++)
        fbdProb += termNode[idx];

    //term 3: fossil attachment
    int numFossils = unresolvedFossils->getNumFossils();
    if(originAge != nullptr)
        unresolvedFossils->syncSpine(originAge->getValue());
    updateGammaCache();
    int spineIdx = unresolvedFossils->getSpineIdx();
    std::vector<double> termFoss(numFossils, 0.0);
    ThreadPool::current().parallelFor(OP_FBD, numFossils, [&](int lo, int hi){
        for(int i = lo; i < hi; i++){
            if(unresolvedFossils->isSA(i)){
                termFoss[i] = std::log(psiAt(findIndex(unresolvedFossils->getFossilAge(i)))) + cachedGammaLn[i];
                continue;
            }
            if(i == spineIdx && unresolvedFossils->isUE(i)){
                termFoss[i] = 0.0;
                continue;
            }
            if(unresolvedFossils->isUE(i)){
                termFoss[i] = uePqLn(unresolvedFossils->getAttachAge(i)) + cachedGammaLn[i];
                continue;
            }
            if(i == spineIdx){
                double ys = unresolvedFossils->getFossilAge(i);
                termFoss[i] = std::log(psiAt(findIndex(ys))) + std::log(calculateP0(ys)) - lnD(ys);
                continue;
            }
            termFoss[i] = fossilPqLn(unresolvedFossils->getFossilAge(i), unresolvedFossils->getAttachAge(i)) + cachedGammaLn[i];
        }
    });
    for(int i = 0; i < numFossils; i++)
        fbdProb += termFoss[i];
    return fbdProb;
}

double FBDTreeModel::calculateResolvedFBD(void){
    Tree* tree = parameterTree->getTree();
    Node* root = tree->getRoot();
    double rootAge = root->getTime();
    bool useOrigin = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);

    for(Node* n : tree->getDownPassSequence())
        if(n != root && tree->isSATip(n) == false && n->getTime() >= n->getAncestor()->getTime())
            return -INFINITY;

    double lnP;
    if(useOrigin){
        double x0 = originAge->getValue();
        for(Node* c : root->getNeighbors())
            if(c != root->getAncestor() && c->getTime() >= x0)
                return -INFINITY;
        lnP = -calculateLnConditioning(x0);
    }else{
        lnP = -2.0 * calculateLnSurvival(rootAge);
    }

    for(Node* n : tree->getDownPassSequence()){
        if(n != root)
            lnP += lnD(n->getAncestor()->getTime()) - lnD(n->getTime());
        if(n->getIsTip()){
            if(n->getIsFossil() == false)
                lnP += std::log(rhoVal);
            else if(n->getAncestor()->getTime() == n->getTime())
                lnP += std::log(psiAt(findIndex(n->getTime())));
            else
                lnP += std::log(psiAt(findIndex(n->getTime()))) + std::log(calculateP0(n->getTime()));
        }
        else if(n != root){
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
    int k = findIndex(t);
    double lam = lambdaAt(k), mu = muAt(k), c1 = c1HatVec[k], c2 = c2HatVec[k];
    double E = std::exp(-c1 * (t - intervalStart[k]));
    double N  = E * (1.0 - c2) * ((lam - mu) - c1) + (1.0 + c2) * ((lam - mu) + c1);
    double Dp = E * (1.0 - c2) + (1.0 + c2);
    return std::log(N) - std::log(2.0 * lam) - std::log(Dp);
}

double FBDTreeModel::calculateLnAnySample(double t){
    int k = findIndex(t);
    double lam = lambdaAt(k), mu = muAt(k), psi = psiAt(k), c1 = c1Vec[k], c2 = c2Vec[k];
    double beta = lam - mu - psi, bpc, bmc;
    if(beta >= 0.0){ bpc = beta + c1; bmc = -4.0*lam*psi/bpc; }
    else           { bmc = beta - c1; bpc = -4.0*lam*psi/bmc; }
    double E = std::exp(-c1 * (t - intervalStart[k]));
    double N  = E * (1.0 - c2) * bmc + (1.0 + c2) * bpc;
    double Dp = E * (1.0 - c2) + (1.0 + c2);
    return std::log(N) - std::log(2.0 * lam) - std::log(Dp);
}

double FBDTreeModel::calculateLnConditioning(double t){
    switch(UserSettings::userSettings().getConditioningEvent()){
        case ConditioningEvent::SURVIVAL:
            return calculateLnSurvival(t);
        case ConditioningEvent::ANYSAMPLE:
            return calculateLnAnySample(t);
        case ConditioningEvent::EXTINCT: {
            double d = std::exp(calculateLnAnySample(t)) - std::exp(calculateLnSurvival(t));
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
    stk.push_back(std::make_pair(tree->getRoot(), false));
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

// gamma_i(z) = backbone lineages crossing z in R_i + sum_j w_ij [pendant j crosses z in R_i]
double FBDTreeModel::computeGamma(double z, int i){
    if(eulerBuilt == false)
        buildEulerIndex();
    Node* crown = unresolvedFossils->getCrownNode(i);
    bool stem = unresolvedFossils->getIsStem(i);
    bool total = (unresolvedFossils->getIsCrown(i) == false && stem == false);
    double count = cladeBackboneLineages(i, z);

    bool halfFix = (UserSettings::userSettings().getModel() == Model::UFBD);   // w=1/2 on reciprocal pairs
    bool focalIsTip = (unresolvedFossils->isSA(i) == false);
    int numFossils = unresolvedFossils->getNumFossils();

    if(wholeTreeTotalFast == 1){                       // one clade, all TOTAL or all CROWN, no stem
        if(total){
            double crossing = (double)countStraddling(sortedFossilY, sortedFossilZ, z);   // focal self-cancels
            if(halfFix && focalIsTip){
                double spineCrosses = 0.0;             // spine keeps full weight (owns x0)
                int sp = unresolvedFossils->getSpineIdx();
                if(i != sp && sp >= 0 && unresolvedFossils->isSA(sp) == false){
                    double ys = unresolvedFossils->getFossilAge(sp), zs = unresolvedFossils->getAttachAge(sp);
                    if(ys < z && z < zs)
                        spineCrosses = 1.0;
                }
                count += 0.5 * crossing + 0.5 * spineCrosses;
            }else{
                count += crossing;
            }
        }else{
            double crossing = (double)countStraddling(sortedFossilY, sortedZY, z);
            crossing -= countAboveCrown(sortedZY, crown->getTime(), z);        // crown-cap
            if(halfFix && focalIsTip)
                count += 0.5 * crossing;
            else
                count += crossing;
        }
        return count;
    }
    if(multiCladeFast == 1){                            // mixed CROWN/TOTAL over clades, no stem
        const CladeGammaIndex& g = cladeGamma.find(crown)->second;
        double T = crown->getTime();
        double nested = (double)countStraddling(g.subY, g.subZY, z);
        if(total == false)
            nested -= countAboveCrown(g.subZY, T, z);   // crown-cap
        if(halfFix && focalIsTip){
            double reciprocal;                          // pendants owned as TOTAL/CROWN of this clade
            if(total == false)
                reciprocal = (double)countStraddling(g.totY, g.totZY, z) + (double)countStraddling(g.crY, g.crZY, z)
                             - countAboveCrown(g.totZY, T, z) - countAboveCrown(g.crZY, T, z);
            else{
                reciprocal = (double)countStraddling(g.totY, g.totZY, z);
                if(z <= T)
                    reciprocal += (double)countStraddling(g.crY, g.crZY, z);
            }
            count += nested - 0.5 * reciprocal;
        }else{
            count += nested;
        }
        return count;
    }
    for(int j = 0; j < numFossils; j++){               // general reference: zone host law per fossil
        if(j == i)
            continue;
        if(unresolvedFossils->isSA(j))
            continue;
        double yj = unresolvedFossils->getFossilAge(j);
        double zj = unresolvedFossils->getAttachAge(j);
        if(yj >= z || z >= zj)
            continue;                                   // j must straddle z
        bool reciprocal;                                // both orientations admissible -> 1/2 eligible
        if(stem){
            if(unresolvedFossils->getCrownNode(j) != crown)
                continue;
            bool jStem = unresolvedFossils->getIsStem(j);
            bool jTotalOnStem = (unresolvedFossils->getIsCrown(j) == false) && (jStem == false) && (zj >= crown->getTime());
            if(jStem == false && jTotalOnStem == false)
                continue;                               // stem-i reaches only same-clade stem
            reciprocal = true;
        }else{
            if(inSub(unresolvedFossils->getCrownNode(j), crown) == false)
                continue;
            if(total == false && zj > crown->getTime())
                continue;                               // crown-cap
            Node* crownJ = unresolvedFossils->getCrownNode(j);
            bool jStem = unresolvedFossils->getIsStem(j);
            bool jTotal = (unresolvedFossils->getIsCrown(j) == false) && (jStem == false);
            bool iPendReachesRj = jTotal ? true
                                : jStem  ? (z >= crownJ->getTime())
                                :          (z <= crownJ->getTime());
            reciprocal = inSub(crown, crownJ) && iPendReachesRj;
        }
        double w = 1.0;
        if(halfFix && focalIsTip && j != unresolvedFossils->getSpineIdx() && reciprocal)
            w = 0.5;
        count += w;
    }
    return count;
}

void FBDTreeModel::updateGammaCache(void){
    Tree* tree = parameterTree->getTree();
    int nf = unresolvedFossils->getNumFossils();
    std::vector<Node*>& dpseq = tree->getDownPassSequence();

    Node* root = tree->getRoot();
    int budget = 0;
    for(int x = (int)sortedYounger.size(); x > 1; x >>= 1) budget++;
    int edits = (cacheInit && sortedYounger.empty() == false) ? 0 : -1;
    if(edits == 0){
        for(Node* nd : dpseq)
            if(nd->getTime() != prevNodeAge[nd->getOffset()]){
                edits += (nd != root ? 1 : 0) + (int)nd->getDescendants().size();
                if(edits > budget){ edits = -1; break; }
            }
    }
    if(edits >= 0){
        for(Node* nd : dpseq){
            double oldA = prevNodeAge[nd->getOffset()];
            double newA = nd->getTime();
            if(oldA == newA)
                continue;
            if(nd != root){
                sortedYounger.erase(std::lower_bound(sortedYounger.begin(), sortedYounger.end(), oldA));
                sortedYounger.insert(std::lower_bound(sortedYounger.begin(), sortedYounger.end(), newA), newA);
            }
            int nc = (int)nd->getDescendants().size();
            for(int c = 0; c < nc; c++){
                sortedOlder.erase(std::lower_bound(sortedOlder.begin(), sortedOlder.end(), oldA));
                sortedOlder.insert(std::lower_bound(sortedOlder.begin(), sortedOlder.end(), newA), newA);
            }
        }
    }else{
        sortedYounger.clear();
        sortedOlder.clear();
        for(Node* nd : dpseq){
            if(nd == root)
                continue;
            sortedYounger.push_back(nd->getTime());
            sortedOlder.push_back(nd->getAncestor()->getTime());
        }
        std::sort(sortedYounger.begin(), sortedYounger.end());
        std::sort(sortedOlder.begin(), sortedOlder.end());
    }

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
            double xlo = std::min(prevX0, x0), xhi = std::max(prevX0, x0);
            for(int i = 0; i < nf; i++){
                if(unresolvedFossils->getCrownNode(i) != root || unresolvedFossils->getIsCrown(i))
                    continue;
                double zi2 = unresolvedFossils->getAttachAge(i);
                if(zi2 >= xlo && zi2 < xhi)
                    gammaStale[i] = 1;
            }
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
            bool pureZ = (yi == prevY[i] && wasTerm && isTerm);
            double lo, hi;
            if(pureZ){
                lo = std::min(prevZ[i], zi);
                hi = std::max(prevZ[i], zi);
            }else{
                lo = std::min(prevY[i], yi);
                hi = -INFINITY;
                if(wasTerm) hi = std::max(hi, prevZ[i]);
                if(isTerm)  hi = std::max(hi, zi);
            }
            for(int j = 0; j < nf; j++){
                if(gammaStale[j]) continue;
                double tj = unresolvedFossils->isSA(j) ? unresolvedFossils->getFossilAge(j) : unresolvedFossils->getAttachAge(j);
                bool aff = pureZ ? (tj >= lo && tj < hi) : (tj > lo && tj < hi);
                if(pureZ && aff == false && (unresolvedFossils->getIsCrown(j) || unresolvedFossils->getIsStem(j))){
                    double ct = unresolvedFossils->getCrownNode(j)->getTime();
                    aff = (ct >= lo && ct < hi);
                }
                if(aff)
                    gammaStale[j] = 1;
            }
        }
    }

    std::vector<std::pair<double,double> > changedIntervals;
    for(Node* n : dpseq){
        if(n == tree->getRoot()) continue;
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

    if(wholeTreeTotalFast < 0){
        wholeTreeTotalFast = 1;
        Node* commonClade = (nf > 0) ? unresolvedFossils->getCrownNode(0) : tree->getCrown();
        bool anyCrown = false, anyTotal = false;
        for(int i = 0; i < nf; i++){
            if(unresolvedFossils->getCrownNode(i) != commonClade || unresolvedFossils->getIsStem(i)){
                wholeTreeTotalFast = 0;
                break;
            }
            if(unresolvedFossils->getIsCrown(i)) anyCrown = true; else anyTotal = true;
        }
        if(wholeTreeTotalFast == 1 && anyCrown && anyTotal)
            wholeTreeTotalFast = 0;
        fastIsCrown = (wholeTreeTotalFast == 1 && anyCrown) ? 1 : 0;
    }
    if(wholeTreeTotalFast == 1){
        sortedFossilY.clear();
        for(int i = 0; i < nf; i++)
            if(unresolvedFossils->isSA(i) == false)
                sortedFossilY.push_back(unresolvedFossils->getFossilAge(i));
        std::sort(sortedFossilY.begin(), sortedFossilY.end());
        if(fastIsCrown){
            sortedZY.clear();
            for(int i = 0; i < nf; i++)
                if(unresolvedFossils->isSA(i) == false)
                    sortedZY.push_back(std::make_pair(unresolvedFossils->getAttachAge(i), unresolvedFossils->getFossilAge(i)));
            std::sort(sortedZY.begin(), sortedZY.end());
        }else{
            sortedFossilZ.clear();
            for(int i = 0; i < nf; i++)
                if(unresolvedFossils->isSA(i) == false)
                    sortedFossilZ.push_back(unresolvedFossils->getAttachAge(i));
            std::sort(sortedFossilZ.begin(), sortedFossilZ.end());
        }
    }

    if(multiCladeFast < 0){
        multiCladeFast = 0;
        if(wholeTreeTotalFast == 0 && tree->getNumBackbone() > 0){
            bool anyStem = false;
            for(int i = 0; i < nf; i++)
                if(unresolvedFossils->getIsStem(i)){ anyStem = true; break; }
            if(anyStem == false){
                multiCladeFast = 1;
                activeClades.clear();
                for(int i = 0; i < nf; i++){
                    Node* c = unresolvedFossils->getCrownNode(i);
                    bool found = false;
                    for(Node* a : activeClades) if(a == c){ found = true; break; }
                    if(found == false) activeClades.push_back(c);
                }
            }
        }
    }
    if(multiCladeFast == 1){
        if(eulerBuilt == false)
            buildEulerIndex();
        cladeGamma.clear();
        for(Node* c : activeClades)
            cladeGamma[c];
        for(int i = 0; i < nf; i++){
            if(unresolvedFossils->isSA(i))
                continue;
            Node* cj = unresolvedFossils->getCrownNode(i);
            double yj = unresolvedFossils->getFossilAge(i);
            double zj = unresolvedFossils->getAttachAge(i);
            CladeGammaIndex& own = cladeGamma[cj];
            if(unresolvedFossils->getIsCrown(i)){
                own.crY.push_back(yj);
                own.crZY.push_back(std::make_pair(zj, yj));
            }else{
                own.totY.push_back(yj);
                own.totZY.push_back(std::make_pair(zj, yj));
            }
            for(Node* c : activeClades)
                if(inSub(cj, c)){
                    CladeGammaIndex& g = cladeGamma[c];
                    g.subY.push_back(yj);
                    g.subZY.push_back(std::make_pair(zj, yj));
                }
        }
        for(Node* c : activeClades){
            CladeGammaIndex& g = cladeGamma[c];
            std::sort(g.subY.begin(), g.subY.end());
            std::sort(g.subZY.begin(), g.subZY.end());
            std::sort(g.totY.begin(), g.totY.end());
            std::sort(g.totZY.begin(), g.totZY.end());
            std::sort(g.crY.begin(), g.crY.end());
            std::sort(g.crZY.begin(), g.crZY.end());
        }
    }

    std::vector<int> staleIdx;
    for(int i = 0; i < nf; i++)
        if(gammaStale[i])
            staleIdx.push_back(i);
    if(eulerBuilt == false)
        buildEulerIndex();
    ThreadPool::current().parallelFor(OP_GAMMA, (int)staleIdx.size(), [&](int a, int b){
        for(int k = a; k < b; k++){
            int i = staleIdx[k];
            double g = computeGamma(unresolvedFossils->getAttachAge(i), i);
            cachedGammaLn[i] = (g > 0.0) ? std::log(g) : -INFINITY;
            gammaStale[i] = 0;
        }
    });

    for(int i = 0; i < nf; i++){
        prevY[i] = unresolvedFossils->getFossilAge(i);
        prevZ[i] = unresolvedFossils->getAttachAge(i);
        prevSA[i] = unresolvedFossils->isSA(i) ? 1 : 0;
    }
    for(Node* n : dpseq)
        prevNodeAge[n->getOffset()] = n->getTime();
}

void FBDTreeModel::computeAgeFloors(std::map<Node*,double>& floors){
    if(isResolved)
        return;
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
        Node* crown = t->getMRCA(clade->getTaxa());
        bool originCond = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);
        Node* origin = originCond ? t->getRoot() : t->getNodeByOffset(clade->getOrigin()->getOffset());
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
        if(n == t->getRoot())
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
    if(i == tree->getRoot() || parent == tree->getRoot())
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
    if(i == tree->getRoot() || j == tree->getRoot() || pi == tree->getRoot() || pj == tree->getRoot())
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
        if(n != tree->getRoot() && tree->isSATip(n) == false && n->getTime() >= n->getAncestor()->getTime())
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
    bool spIsRoot = (sp == tree->getRoot());
    if(spIsRoot && useOrigin == false)
        return -INFINITY;
    Node* gp = sp->getAncestor();
    Node* sib = nullptr;
    for(Node* c : sp->getNeighbors())
        if(c != gp && c != f)
            sib = c;
    double y = f->getTime();
    double maxAge = spIsRoot ? originAge->getValue() : gp->getTime();
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
        if(n != t->getRoot() && n->getAncestor()->getTime() < n->getTime())
            return -INFINITY;
    return h;
}

double FBDTreeModel::getOriginAgeValue(void){
    return originAge->getValue();
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
        if(n != tree->getRoot() && n->getIsTip() == false)
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
    Node* treeRoot = t->getRoot();
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
            if(n->getAncestor() != treeRoot && t->isSATip(n) == false)
                roots.push_back(n);
        }
    }
}

void FBDTreeModel::enumerateSubtreeHosts(Tree* t, std::vector<Node*>& crowns, std::vector<char>& isCrowns, std::vector<Node*>& origins, double rAge, double ceilingS, std::vector<Node*>& hosts, std::vector<double>& los, std::vector<double>& his){
    for(Node* n : t->getDownPassSequence()){
        if(n == t->getRoot())
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
