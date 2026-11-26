#include <map>
#include <functional>
#include <iomanip>
#include "Node.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "ParameterOUField.hpp"
#include "FBDTreeModel.hpp"
#include "PathModel.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "Serialize.hpp"
#include "ThreadPool.hpp"
#include "UserSettings.hpp"
#include "Probability.hpp"
#include "Msg.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>

static double offsetFromLambda(double rate0, double lam0, const Probability::PriorSpec& sp){
    if(rate0 != lam0)
        return rate0;
    const double hi = std::numeric_limits<double>::max();
    const double cand[4] = {0.5 * rate0, 1.5 * rate0, rate0 * (1.0 - 1e-3), rate0 * (1.0 + 1e-3)};
    for(int i = 0; i < 4; i++)
        if(cand[i] > 0.0 && cand[i] != lam0 &&
           std::isfinite(Probability::priorLnPdf(sp.family, sp.p1, sp.p2, cand[i], 0.0, hi, sp.p3)))
            return cand[i];
    return rate0;
}

static int numFree(const std::vector<ParameterDouble*>& v);

FBDTreeModel::FBDTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, unsigned int seed, UserSettings* rateSettings) :
    PhylogeneticModel(){

    lastMoveKind = MK_PARAM;
    lastRateVec = nullptr;
    lastRateField = nullptr;
    rateVecStep = 0.2;
    shrinkStep = 0.2;
    rvAccW = rvAttW = rvAtt = seAccW = seAttW = seAtt = 0;
    lastTreeMove = TM_NONE;
    for(int i = 0; i < TM_COUNT; i++) tmAcc[i] = tmAtt[i] = 0;
    azAcc = 0;
    azAtt = 0;
    gammaStaleTot = 0;
    gammaStaleCalls = 0;
    shiftStep = 0.005;
    saBatch = 0;
    saBatchF = 1.0;
    rsAccW = rsAttW = 0;
    rsAcc = rsTot = rsAdapt = 0;
    upDownStep = 0.1;
    upDownTotal = 0;
    cacheInit = false;
    zoneInit = false;
    rng.setSeed(seed);
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);

    parameterTree = new ParameterTree(1.0, this);
    isResolved = (UserSettings::userSettings().getModel() == Model::RFBD);

    auto resolveType = [&](Fossil& f, const std::vector<std::string>& names) -> int {
        if(names.size() <= 1 || f.getMaxAge() <= 0.0) return 0;
        if(f.getType().empty())
            Msg::error("fossil '" + f.getTaxon() + "' has no preservation type, but -psi_types declares multiple types.");
        std::vector<std::string>::const_iterator it = std::find(names.begin(), names.end(), f.getType());
        if(it == names.end())
            Msg::error("fossil '" + f.getTaxon() + "' has preservation type '" + f.getType() + "' not listed in -psi_types.");
        return (int)(it - names.begin());
    };

    std::vector<Fossil> unrFossils, backboneDefs;
    if(isResolved == false){
        for(Fossil& f : fossils){
            Node* tip = t->getTaxonNode(f.getTaxon());
            if(tip == nullptr || tip->getIsTip() == false){ unrFossils.push_back(f); continue; }
            tip->setIsFossil(true);
            tip->setTime(0.5 * (f.getMinAge() + f.getMaxAge()));
            tip->setFossilAgeRange(f.getMinAge(), f.getMaxAge());
            backboneFossils.push_back({tip, f.getMinAge(), f.getMaxAge()});
            backboneDefs.push_back(f);
        }
        if(backboneFossils.empty() == false){
            t->liftInternalAgesAboveChildren();
            static bool warnedNotYoungest = false;
            if(warnedNotYoungest == false){
                warnedNotYoungest = true;
                std::map<std::string, Node*> cladeCrown;
                for(Clade& c : clades)
                    cladeCrown[c.getName()] = c.getCrown();
                auto isUnder = [](Node* a, Node* anc) -> bool {
                    for(Node* p = a; ; p = p->getAncestor()){
                        if(p == anc) return true;
                        if(p->getAncestor() == p) return false;
                    }
                };
                for(size_t k = 0; k < backboneFossils.size(); k++){
                    Node* fTip = backboneFossils[k].tip;
                    std::string fTaxon = fTip->getName();
                    double minCompat = std::numeric_limits<double>::infinity();
                    for(Fossil& g : fossils){
                        if(g.getInBackbone() || g.getAssignment() == Assignment::STEM) continue;
                        if(isUnder(fTip, cladeCrown[g.getClade()]) == false) continue;
                        if(g.getMinAge() < minCompat) minCompat = g.getMinAge();
                    }
                    if(backboneFossils[k].yMax >= minCompat)
                        Msg::warning("backbone fossil '" + fTaxon + "' is not the absolute youngest among clade-compatible fossils; it is blocked from being a sampled ancestor.");
                }
            }
        }
    }else{
        unrFossils = fossils;
    }

    int numUE = 0;
    for(Fossil& f : fossils)
        if(f.getMaxAge() == 0.0) numUE++;
    numExtantTips = t->getNumExtant();
    int numExtant = numExtantTips;
    Conditioning condPoint = UserSettings::userSettings().getConditioning();
    ConditioningEvent condEvent = UserSettings::userSettings().getConditioningEvent();
    if(condPoint == Conditioning::CROWN && t->getNumBackbone() < 2)
        Msg::error("crown conditioning requires at least 2 backbone tips, but the backbone tree has " + std::to_string(t->getNumBackbone()) + ".");
    if(condEvent == ConditioningEvent::SURVIVAL && (numExtant + numUE) < 1)
        Msg::error("survival conditioning requires at least 1 extant taxon.");
    if(condEvent == ConditioningEvent::EXTINCT && (numExtant + numUE) > 0)
        Msg::error("extinct conditioning requires 0 extant taxa, but the data has " + std::to_string(numExtant + numUE) + ".");

    originAge = nullptr;
    if(UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN){
        double floor = t->getCrown()->getTime();
        double oldestFossil = 0.0;
        for(Fossil& f : fossils){
            if(f.getMaxAge() > floor)
                floor = f.getMaxAge();
            if(f.getMaxAge() > oldestFossil)
                oldestFossil = f.getMaxAge();
        }
        double x0init = floor * 1.05;
        originAge = new ParameterDouble(1.0, this, "originAge", 0.0, std::numeric_limits<double>::max());
        UserSettings& us = UserSettings::userSettings();
        if(us.getConditionAgePriorSet()){
            double off = us.getConditionAgePriorP3();
            originAge->setPrior(us.getConditionAgePrior(), us.getConditionAgePriorP1(), us.getConditionAgePriorP2(), off);
            double pm = Probability::priorMean(us.getConditionAgePrior(), us.getConditionAgePriorP1(), us.getConditionAgePriorP2(), off);
            if(pm > x0init)
                x0init = pm;
            bool bounded = (us.getConditionAgePrior() == Probability::PriorFamily::UNIFORM ||
                            us.getConditionAgePrior() == Probability::PriorFamily::EMPIRICAL);
            double hi = us.getConditionAgePriorP2();
            if(bounded && x0init >= hi)
                x0init = 0.5 * (floor + hi);
            if(x0init <= off)
                x0init = off + 0.05 * (floor > 0.0 ? floor : 1.0);
        }else{
            originAge->setPrior(Probability::PriorFamily::IMPROPER, 0.0, 0.0);
        }
        if(originAge->getProposalProbability() > 0.0)
            originAge->setValue(x0init);
        else if(originAge->getValue() < oldestFossil)
            Msg::error("fixed origin age " + std::to_string(originAge->getValue()) + " is below the oldest fossil age " + std::to_string(oldestFossil) + ".");
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
    if(UserSettings::userSettings().getConditioning() == Conditioning::CROWN &&
       UserSettings::userSettings().getConditionAgePrior() == Probability::PriorFamily::FIXED){
        double fixedCrown = UserSettings::userSettings().getConditionAgePriorP1();
        double oldestFossil = 0.0;
        for(Fossil& f : fossils)
            if(f.getMaxAge() > oldestFossil)
                oldestFossil = f.getMaxAge();
        if(fixedCrown < oldestFossil)
            Msg::error("fixed crown age " + std::to_string(fixedCrown) + " is below the oldest fossil age " + std::to_string(oldestFossil) + ".");
        Tree* crownTree = parameterTree->getTree();
        double curCrown = crownTree->getCrown()->getTime();
        if(curCrown > 0.0 && curCrown != fixedCrown)
            crownTree->scaleInternalAges(fixedCrown / curCrown);
        crownTree->setCrownFixed(true);
    }
    for(BackboneFossil& bf : backboneFossils)
        bf.tip = parameterTree->getTree()->getNodeByOffset(bf.tip->getOffset());
    parameters.push_back(parameterTree);
    if(originAge != nullptr)
        parameters.push_back(originAge);

    rho = UserSettings::userSettings().getRho();
    buildRates(rateSettings != nullptr ? *rateSettings : UserSettings::userSettings(), a);
    for(Fossil& f : backboneDefs)
        a.backboneType.push_back(resolveType(f, a.psiTypeNames));
    for(Fossil& f : unrFossils){
        a.fossilType.push_back(resolveType(f, a.psiTypeNames));
        a.fossilTypeByName[f.getTaxon()] = a.fossilType.back();
        unrFossilName.push_back(f.getTaxon());
    }

    unresolvedFossils = nullptr;
    if(isResolved){
        Tree* wt = parameterTree->getTree();
        for(Fossil& f : unrFossils){
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
        unresolvedFossils = new ParameterUnresolvedFossils(1.0, this, parameterTree->getTree(), clades, unrFossils, originAge);
        parameters.push_back(unresolvedFossils);
        for(int i = 0; i < unresolvedFossils->getNumFossils(); i++)
            if(unresolvedFossils->isUE(i))
                ueFossils.push_back(i);
    }

    if(isResolved){
        parameterTree->setProposalProbability(78.0);
        for(ParameterDouble* l : a.lambda) l->setProposalProbability(15.0);
        for(ParameterDouble* m : a.mu)     m->setProposalProbability(15.0);
        for(auto& pv : a.psi) for(ParameterDouble* p : pv) p->setProposalProbability(15.0);
    }
    double fieldBase = (isResolved ? 15.0 : 1.0);
    if(a.lambdaField) a.lambdaField->setProposalProbability(fieldBase * (double)a.lambdaField->getNumBins());
    if(a.muField)     a.muField->setProposalProbability(fieldBase * (double)a.muField->getNumBins());
    for(int tp = 0; tp < a.numPsiTypes; tp++) if(a.psiField[tp]) a.psiField[tp]->setProposalProbability(fieldBase * (double)a.psiField[tp]->getNumBins());

    double sum = 0.0;
    for(Parameter* p : parameters)
        sum += p->getProposalProbability();
    for(Parameter* p : parameters)
        p->setProposalProbability(p->getProposalProbability() / sum);

    haveIid = (a.lambdaField != nullptr || numFree(a.lambda) >= 2) || (a.muField != nullptr || numFree(a.mu) >= 2);
    for(int tp = 0; tp < a.numPsiTypes; tp++)
        if(a.psiField[tp] != nullptr || numFree(a.psi[tp]) >= 2) haveIid = true;
    haveLMField = (a.lambdaField != nullptr || a.muField != nullptr);

    UserSettings& us = UserSettings::userSettings();
    psiActive.assign(a.numPsiTypes, 1);
    if(rateSettings == nullptr && us.getPsiRjMain() >= 0){
        psiRjMain = us.getPsiRjMain();
        fossilTypeOrig = a.fossilType;
        backboneTypeOrig = a.backboneType;
        psiRjLab = PathModel::lnTypeLabels(a.psiTypeNames, fossils);
        psiRjCheck = us.getResampleEvery();
        psiRjWeight = -psiRjLab;
        size_t nm = psiShiftCoords().size();
        psiRjShift.assign(nm, 0.0);
        psiRjSum[0].assign(nm, 0.0);
        psiRjSum[1].assign(nm, 0.0);
        psiRjShiftable.assign(nm, 1);
        if(a.psiField[psiRjMain] == nullptr)
            for(size_t b = 0; b < nm; b++)
                psiRjShiftable[b] = a.psi[psiRjMain][b]->getProposalProbability() > 0.0;
    }

    RandomVariable::setActiveInstance(prevRng);
}

void FBDTreeModel::applyPsiModel(int m){
    psiModel = m;
    for(size_t i = 0; i < fossilTypeOrig.size(); i++){
        a.fossilType[i] = (m == 1) ? psiRjMain : fossilTypeOrig[i];
        a.fossilTypeByName[unrFossilName[i]] = a.fossilType[i];
    }
    for(size_t k = 0; k < backboneTypeOrig.size(); k++)
        a.backboneType[k] = (m == 1) ? psiRjMain : backboneTypeOrig[k];
    for(int t = 0; t < a.numPsiTypes; t++)
        psiActive[t] = (m == 0 || t == psiRjMain);
    a.fbdMemoInit = false;
}

std::vector<double> FBDTreeModel::psiMainValues(void){
    std::vector<double> v;
    if(ParameterOUField* f = a.psiField[psiRjMain]){
        for(int b = 0; b < f->getNumBins(); b++)
            v.push_back(f->getRate(b));
        v.push_back(f->getTheta());
    }else
        for(ParameterDouble* p : a.psi[psiRjMain])
            v.push_back(p->getValue());
    return v;
}

void FBDTreeModel::setPsiMainValues(const std::vector<double>& v){
    if(ParameterOUField* f = a.psiField[psiRjMain])
        f->setState(std::vector<double>(v.begin(), v.end() - 1), v.back(), f->getSdEq(), f->getNu());
    else
        for(size_t b = 0; b < v.size(); b++)
            a.psi[psiRjMain][b]->setValue(v[b]);
}

std::vector<double> FBDTreeModel::psiShiftCoords(void){
    std::vector<double> v;
    if(ParameterOUField* f = a.psiField[psiRjMain]){
        if(f->isThetaFixed() == false)
            v.push_back(f->getTheta());
    }else
        for(ParameterDouble* p : a.psi[psiRjMain])
            v.push_back(std::log(p->getValue()));
    return v;
}

double FBDTreeModel::psiRjJacobian(void){
    if(ParameterOUField* f = a.psiField[psiRjMain])
        return psiRjShift.empty() ? 0.0 : f->getNumBins() * psiRjShift[0];
    double j = 0.0;
    for(double s : psiRjShift)
        j += s;
    return j;
}

double FBDTreeModel::lnPsiMainPrior(void){
    if(psiRjShift.empty())
        return 0.0;
    if(a.psiField[psiRjMain] != nullptr)
        return a.psiField[psiRjMain]->lnProbability();
    double lp = 0.0;
    for(ParameterDouble* p : a.psi[psiRjMain])
        lp += p->lnProbability();
    return lp;
}

void FBDTreeModel::scalePsiMain(double sgn){
    psiRjSaved = psiMainValues();
    std::vector<double> v = psiRjSaved;
    if(a.psiField[psiRjMain] != nullptr){
        double s = psiRjShift.empty() ? 0.0 : psiRjShift[0];
        for(size_t b = 0; b + 1 < v.size(); b++)
            v[b] *= std::exp(sgn * s);
        v.back() += sgn * s;
    }else
        for(size_t b = 0; b < v.size(); b++)
            v[b] *= std::exp(sgn * psiRjShift[b]);
    setPsiMainValues(v);
}

double FBDTreeModel::lnPsiRjTerm(void){
    return (psiModel == 0) ? 0.0 : psiRjWeight;
}

double FBDTreeModel::psiRjAdapt(unsigned long k){
    if(psiRjTuning() == false)
        return 0.0;
    psiRjVisits[psiModel]++;
    psiRjN[psiModel]++;
    std::vector<double> v = psiShiftCoords();
    for(size_t b = 0; b < v.size(); b++)
        if(psiRjShiftable[b])
            psiRjSum[psiModel][b] += v[b];
    double d = (psiModel == 0) ? psiRjStep : -psiRjStep;
    psiRjWeight += d;
    if(k % psiRjCheck == 0 && std::min(psiRjVisits[0], psiRjVisits[1]) >= 0.4 * (psiRjVisits[0] + psiRjVisits[1])){
        for(size_t b = 0; b < psiRjShift.size(); b++)
            psiRjShift[b] = psiRjShiftable[b] ? psiRjSum[1][b] / psiRjN[1] - psiRjSum[0][b] / psiRjN[0] : 0.0;
        psiRjShiftSum = psiRjJacobian();
        if(psiRjStage == 0)
            for(int m = 0; m < 2; m++){
                std::fill(psiRjSum[m].begin(), psiRjSum[m].end(), 0.0);
                psiRjN[m] = 0.0;
            }
        psiRjStage++;
        psiRjStep = (psiRjStep * (psiRjVisits[0] + psiRjVisits[1]) <= std::log(2.0)) ? 0.0 : 0.5 * psiRjStep;
        psiRjVisits[0] = psiRjVisits[1] = 0.0;
    }
    return (psiModel == 1) ? d : 0.0;
}

std::vector<double> FBDTreeModel::getPsiRjState(void){
    std::vector<double> s = {psiRjWeight, psiRjStep, psiRjVisits[0], psiRjVisits[1], (double)psiRjStage, psiRjN[0], psiRjN[1]};
    s.insert(s.end(), psiRjShift.begin(), psiRjShift.end());
    s.insert(s.end(), psiRjSum[0].begin(), psiRjSum[0].end());
    s.insert(s.end(), psiRjSum[1].begin(), psiRjSum[1].end());
    return s;
}

double FBDTreeModel::setPsiRjState(const std::vector<double>& s){
    double d = (psiModel == 1) ? s[0] - psiRjWeight : 0.0;
    psiRjWeight = s[0];
    psiRjStep = s[1];
    psiRjVisits[0] = s[2];
    psiRjVisits[1] = s[3];
    psiRjStage = (int)s[4];
    psiRjN[0] = s[5];
    psiRjN[1] = s[6];
    size_t nm = psiRjShift.size();
    for(size_t b = 0; b < nm; b++){
        psiRjShift[b] = s[7 + b];
        psiRjSum[0][b] = s[7 + nm + b];
        psiRjSum[1][b] = s[7 + 2 * nm + b];
    }
    psiRjShiftSum = psiRjJacobian();
    return d;
}

void FBDTreeModel::buildRates(UserSettings& us, SkylineRates& s){
    std::vector<double> lambdaTimes, muTimes;
    lambdaTimes.push_back(0.0);
    for(double t : us.getLambdaSkylineTimes()) lambdaTimes.push_back(t);
    muTimes.push_back(0.0);
    for(double t : us.getMuSkylineTimes())     muTimes.push_back(t);
    s.numPsiTypes = us.getNumPsiTypes();
    s.psiTypeNames = us.getPsiTypeNames();
    std::vector<std::vector<double>> psiTimes(s.numPsiTypes);
    for(int tp = 0; tp < s.numPsiTypes; tp++){
        psiTimes[tp].push_back(0.0);
        for(double t : us.getPsiSkylineTimes(tp)) psiTimes[tp].push_back(t);
    }
    s.intervalStart.push_back(0.0);
    for(double t : us.getSkylineTimes())
        s.intervalStart.push_back(t);
    s.psiIdx.assign(s.numPsiTypes, std::vector<int>());
    for(double st : s.intervalStart){
        int li = 0, mi = 0;
        for(int k = 0; k < (int)lambdaTimes.size(); k++) if(lambdaTimes[k] <= st) li = k;
        for(int k = 0; k < (int)muTimes.size(); k++)     if(muTimes[k] <= st)     mi = k;
        s.lambdaIdx.push_back(li);
        s.muIdx.push_back(mi);
        for(int tp = 0; tp < s.numPsiTypes; tp++){
            int pi = 0;
            for(int k = 0; k < (int)psiTimes[tp].size(); k++) if(psiTimes[tp][k] <= st) pi = k;
            s.psiIdx[tp].push_back(pi);
        }
    }
    Probability::PriorSpec lp = us.getLambdaPrior();
    Probability::PriorSpec mp = us.getMuPrior();
    if(us.getLambdaMode() == RateMode::OU && (lp.set || us.getLambdaGroupPrior().empty() == false)) Msg::error("-lambda_prior is only used under mode=indep; under mode=ou set the OU level with -lambda_ou_theta.");
    if(us.getMuMode() == RateMode::OU && (mp.set || us.getMuGroupPrior().empty() == false)) Msg::error("-mu_prior is only used under mode=indep; under mode=ou set the OU level with -mu_ou_theta.");
    Probability::PriorSpec defRate{true, Probability::PriorFamily::EXPONENTIAL, 5.0, 1.0};
    if(!lp.set) lp = defRate;
    if(!mp.set) mp = defRate;
    int nLambda = (int)lambdaTimes.size();
    int nMu = (int)muTimes.size();
    double lam0 = Probability::priorMean(lp.family, lp.p1, lp.p2, lp.p3);
    if(!(lam0 > 0.0 && std::isfinite(lam0))) lam0 = 0.1;
    double mu0 = Probability::priorMean(mp.family, mp.p1, mp.p2, mp.p3);
    if(!(mu0 > 0.0 && std::isfinite(mu0))) mu0 = 0.1;
    {
        std::vector<int> b2c = buildSkylineRates("lambda", "", nLambda, lambdaTimes, us.getLambdaMode(), us.getLambdaOU(), lp, lam0, std::numeric_limits<double>::quiet_NaN(), us.getLambdaGroups(), us.getLambdaGroupNames(), us.getLambdaGroupPrior(), s.lambda, s.lambdaField, s.lambdaName);
        for(int& u : s.lambdaIdx) u = b2c[u];
        appendRateMap(lambdaTimes, b2c, s.lambdaName);
    }
    {
        std::vector<int> b2c = buildSkylineRates("mu", "", nMu, muTimes, us.getMuMode(), us.getMuOU(), mp, mu0, lam0, us.getMuGroups(), us.getMuGroupNames(), us.getMuGroupPrior(), s.mu, s.muField, s.muName);
        for(int& u : s.muIdx) u = b2c[u];
        appendRateMap(muTimes, b2c, s.muName);
    }
    s.psi.assign(s.numPsiTypes, std::vector<ParameterDouble*>());
    s.psiField.assign(s.numPsiTypes, nullptr);
    s.psiName.assign(s.numPsiTypes, std::vector<std::string>());
    for(int tp = 0; tp < s.numPsiTypes; tp++){
        Probability::PriorSpec pp = us.getPsiPrior(tp);
        if(us.getPsiMode(tp) == RateMode::OU && (pp.set || us.getPsiGroupPrior(tp).empty() == false)) Msg::error("-psi_prior is only used under mode=indep; under mode=ou set the OU level with -psi_ou_theta.");
        if(!pp.set) pp = defRate;
        double psi0 = Probability::priorMean(pp.family, pp.p1, pp.p2, pp.p3);
        if(!(psi0 > 0.0 && std::isfinite(psi0))) psi0 = 0.1;
        int nPsi = (int)psiTimes[tp].size();
        std::string prefix = (s.numPsiTypes > 1) ? ("psi_" + s.psiTypeNames[tp]) : "psi";
        std::string sep = (s.numPsiTypes > 1) ? "_" : "";
        std::vector<int> b2c = buildSkylineRates(prefix, sep, nPsi, psiTimes[tp], us.getPsiMode(tp), us.getPsiOU(tp), pp, psi0, lam0, us.getPsiGroups(tp), us.getPsiGroupNames(tp), us.getPsiGroupPrior(tp), s.psi[tp], s.psiField[tp], s.psiName[tp]);
        for(int& u : s.psiIdx[tp]) u = b2c[u];
        appendRateMap(psiTimes[tp], b2c, s.psiName[tp]);
    }
}

std::vector<int> FBDTreeModel::buildSkylineRates(const std::string& prefix, const std::string& sep, int nBins, const std::vector<double>& times, RateMode mode, const OUHyperSpec& ou, const Probability::PriorSpec& basePrior, double rate0, double avoid, const std::vector<int>& groupIds, const std::vector<std::string>& groupLabels, const std::map<int,Probability::PriorSpec>& groupPrior, std::vector<ParameterDouble*>& outVec, ParameterOUField*& outField, std::vector<std::string>& outNames){
    rate0 = offsetFromLambda(rate0, avoid, basePrior);
    std::vector<int> binToChunk(nBins, 0);
    std::map<int,int> gidToChunk;
    std::vector<int> chunkGid, chunkMinBin;
    for(int i = 0; i < nBins; i++){
        int gid = groupIds.empty() ? i : groupIds[i];
        std::map<int,int>::iterator it = gidToChunk.find(gid);
        if(it == gidToChunk.end()){ int c = (int)chunkGid.size(); gidToChunk[gid] = c; chunkGid.push_back(gid); chunkMinBin.push_back(i); binToChunk[i] = c; }
        else binToChunk[i] = it->second;
    }
    int nChunks = (int)chunkGid.size();
    outNames.clear();
    for(int c = 0; c < nChunks; c++){
        if(nChunks == 1){ outNames.push_back(prefix); continue; }
        if(chunkGid[c] >= 0 && chunkGid[c] < (int)groupLabels.size())
            outNames.push_back(prefix + "_" + groupLabels[chunkGid[c]]);
        else
            outNames.push_back(prefix + sep + std::to_string(chunkMinBin[c]));
    }

    bool useOu = (mode == RateMode::OU);
    if(useOu && nChunks < 2){
        Msg::warning("prior_mode=ou set for " + prefix + " but there is only one rate interval; using an independent per-bin prior instead.");
        useOu = false;
    }
    if(useOu){
        std::vector<int> chunkLastBin(nChunks, -1), chunkCount(nChunks, 0);
        for(int i = 0; i < nBins; i++){ chunkLastBin[binToChunk[i]] = i; chunkCount[binToChunk[i]]++; }
        for(int c = 0; c < nChunks; c++)
            if(chunkLastBin[c] - chunkMinBin[c] + 1 != chunkCount[c])
                Msg::error("prior_mode=ou requires contiguous bins; " + prefix + " has a union bin.");
        std::vector<double> loEdges(nChunks);
        for(int c = 0; c < nChunks; c++)
            loEdges[c] = times[chunkMinBin[c]];
        outField = new ParameterOUField(1.0, this, nChunks, loEdges, rate0, originAge, ou);
        parameters.push_back(outField);
        return binToChunk;
    }

    for(int c = 0; c < nChunks; c++){
        ParameterDouble* p = new ParameterDouble(1.0, this, outNames[c], 0.0, std::numeric_limits<double>::max());
        std::map<int,Probability::PriorSpec>::const_iterator pit = groupPrior.find(chunkGid[c]);
        bool ov = (pit != groupPrior.end() && pit->second.set);
        const Probability::PriorSpec& sp = ov ? pit->second : basePrior;
        p->setPrior(sp.family, sp.p1, sp.p2);
        p->setValue(offsetFromLambda(ov ? Probability::priorMean(sp.family, sp.p1, sp.p2, sp.p3) : rate0, avoid, sp));
        outVec.push_back(p);
        parameters.push_back(p);
    }
    return binToChunk;
}

void FBDTreeModel::appendRateMap(const std::vector<double>& times, const std::vector<int>& binToChunk, const std::vector<std::string>& names){
    int nBins = (int)binToChunk.size();
    double off = UserSettings::userSettings().getAgeOffset();
    for(int c = 0; c < (int)names.size(); c++){
        std::string iv;
        bool first = true;
        for(int i = 0; i < nBins; i++){
            if(binToChunk[i] != c) continue;
            char buf[64];
            if(i + 1 < nBins) std::snprintf(buf, sizeof(buf), "[%g, %g)", times[i] + off, times[i + 1] + off);
            else              std::snprintf(buf, sizeof(buf), "[%g, inf)", times[i] + off);
            iv += (first ? "" : " + ");
            iv += buf;
            first = false;
        }
        rateMapRows.push_back(std::make_pair(names[c], iv));
    }
}

std::string FBDTreeModel::getRateMap(void){
    size_t w = 0;
    for(const std::pair<std::string,std::string>& r : rateMapRows)
        if(r.first.size() > w) w = r.first.size();
    std::string s;
    for(const std::pair<std::string,std::string>& r : rateMapRows){
        std::string nm = r.first;
        nm.resize(w, ' ');
        s += "    " + nm + "  " + r.second + "\n";
    }
    return s;
}

double FBDTreeModel::lambdaAt(int i){
    int j = a.lambdaIdx[i];
    return (a.lambdaField != nullptr) ? a.lambdaField->getRate(j) : a.lambda[j]->getValue();
}

double FBDTreeModel::muAt(int i){
    int j = a.muIdx[i];
    if(a.muField != nullptr)
        return a.muField->getRate(j);
    return a.mu[j]->getValue();
}

double FBDTreeModel::psiOfTypeAt(int type, int i){
    int j = a.psiIdx[type][i];
    return (a.psiField[type] != nullptr) ? a.psiField[type]->getRate(j) : a.psi[type][j]->getValue();
}

double FBDTreeModel::psiTotalAt(int i){
    double s = 0.0;
    for(int t = 0; t < a.numPsiTypes; t++)
        if(psiActive[t])
            s += psiOfTypeAt(t, i);
    return s;
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

std::vector<std::string> FBDTreeModel::getRateNames(const SkylineRates& s){
    std::vector<std::string> names;
    if(s.lambdaField != nullptr){
        for(const std::string& n : s.lambdaName) names.push_back(n);
        names.push_back("ou_lambda_theta"); names.push_back("ou_lambda_sdEq"); names.push_back("ou_lambda_nu");
    }else
        for(ParameterDouble* p : s.lambda) names.push_back(p->getName());
    if(s.muField != nullptr){
        for(const std::string& n : s.muName) names.push_back(n);
        names.push_back("ou_mu_theta"); names.push_back("ou_mu_sdEq"); names.push_back("ou_mu_nu");
    }else
        for(ParameterDouble* p : s.mu) names.push_back(p->getName());
    for(int tp = 0; tp < s.numPsiTypes; tp++){
        if(s.psiField[tp] != nullptr){
            for(const std::string& n : s.psiName[tp]) names.push_back(n);
            std::string pf = "ou_psi" + (s.numPsiTypes > 1 ? "_" + s.psiTypeNames[tp] : std::string(""));
            names.push_back(pf + "_theta"); names.push_back(pf + "_sdEq"); names.push_back(pf + "_nu");
        }else
            for(ParameterDouble* p : s.psi[tp]) names.push_back(p->getName());
    }
    return names;
}

std::vector<double> FBDTreeModel::getRateString(const SkylineRates& s){
    std::vector<double> vals;
    if(s.lambdaField != nullptr){
        for(int i = 0; i < s.lambdaField->getNumBins(); i++) vals.push_back(s.lambdaField->getRate(i));
        vals.push_back(s.lambdaField->getTheta()); vals.push_back(s.lambdaField->getSdEq()); vals.push_back(s.lambdaField->getNu());
    }else
        for(ParameterDouble* p : s.lambda) vals.push_back(p->getValue());
    if(s.muField != nullptr){
        for(int i = 0; i < s.muField->getNumBins(); i++) vals.push_back(s.muField->getRate(i));
        vals.push_back(s.muField->getTheta()); vals.push_back(s.muField->getSdEq()); vals.push_back(s.muField->getNu());
    }else
        for(ParameterDouble* p : s.mu) vals.push_back(p->getValue());
    for(int tp = 0; tp < s.numPsiTypes; tp++){
        if(s.psiField[tp] != nullptr){
            for(int i = 0; i < s.psiField[tp]->getNumBins(); i++) vals.push_back(s.psiField[tp]->getRate(i));
            vals.push_back(s.psiField[tp]->getTheta()); vals.push_back(s.psiField[tp]->getSdEq()); vals.push_back(s.psiField[tp]->getNu());
        }else
            for(ParameterDouble* p : s.psi[tp]) vals.push_back(p->getValue());
    }
    return vals;
}

std::vector<std::pair<double,double>> FBDTreeModel::getRateSupports(const SkylineRates& s){
    std::vector<std::pair<double,double>> sup;
    auto add = [&](ParameterOUField* f, const std::vector<ParameterDouble*>& v){
        if(f != nullptr){
            sup.insert(sup.end(), f->getNumBins() + 1, {-INFINITY, INFINITY});
            sup.push_back(f->sdSupport());
            sup.push_back(f->nuSupport());
        }else
            for(ParameterDouble* p : v) sup.push_back(p->support());
    };
    add(s.lambdaField, s.lambda);
    add(s.muField, s.mu);
    for(int tp = 0; tp < s.numPsiTypes; tp++)
        add(s.psiField[tp], s.psi[tp]);
    return sup;
}

std::vector<size_t> FBDTreeModel::rateBlockSizes(const SkylineRates& s){
    std::vector<size_t> b = {s.lambdaField != nullptr ? s.lambdaName.size() + 3 : s.lambda.size(), s.muField != nullptr ? s.muName.size() + 3 : s.mu.size()};
    for(int tp = 0; tp < s.numPsiTypes; tp++)
        b.push_back(s.psiField[tp] != nullptr ? s.psiName[tp].size() + 3 : s.psi[tp].size());
    return b;
}

std::vector<std::string> FBDTreeModel::getParameterNames(void){
    std::vector<std::string> names = getRateNames(a);
    if(psiRjMain >= 0)
        names.insert(names.end(), {"model", "delta", "log_prior_odds"});
    if(unresolvedFossils != nullptr || isResolved)
        names.push_back("nSA");
    if(originAge != nullptr)
        names.push_back("originAge");
    size_t numAgeNodes = getAgeLogNodes().size();
    for(size_t i = 0; i < numAgeNodes; i++)
        names.push_back("x" + std::to_string(i + 1));
    return names;
}

void FBDTreeModel::appendFieldMask(std::vector<bool>& m, ParameterOUField* f, size_t nRates){
    m.insert(m.end(), nRates, false);
    m.push_back(f->isThetaFixed());
    m.push_back(f->isSdFixed());
    m.push_back(f->isNuFixed());
}

void FBDTreeModel::appendRateMask(std::vector<bool>& m, const SkylineRates& s){
    if(s.lambdaField != nullptr) appendFieldMask(m, s.lambdaField, s.lambdaName.size());
    else for(ParameterDouble* p : s.lambda) m.push_back(p->getProposalProbability() == 0.0);
    if(s.muField != nullptr) appendFieldMask(m, s.muField, s.muName.size());
    else for(ParameterDouble* p : s.mu) m.push_back(p->getProposalProbability() == 0.0);
    for(int tp = 0; tp < s.numPsiTypes; tp++){
        if(s.psiField[tp] != nullptr) appendFieldMask(m, s.psiField[tp], s.psiName[tp].size());
        else for(ParameterDouble* p : s.psi[tp]) m.push_back(p->getProposalProbability() == 0.0);
    }
}

std::vector<bool> FBDTreeModel::getParameterFixedMask(void){
    std::vector<bool> m;
    appendRateMask(m, a);
    if(psiRjMain >= 0)
        m.insert(m.end(), {false, true, true});
    if(unresolvedFossils != nullptr){
        int nNonUE = 0;
        for(int i = 0; i < unresolvedFossils->getNumFossils(); i++)
            if(unresolvedFossils->isUE(i) == false) nNonUE++;
        m.push_back(nNonUE <= (unresolvedFossils->getSpineIdx() >= 0 ? 1 : 0));
    }else if(isResolved)
        m.push_back(false);
    if(originAge != nullptr) m.push_back(originAge->getProposalProbability() == 0.0);
    std::vector<Node*> logNodes = getAgeLogNodes();
    Node* crownNode = parameterTree->getTree()->getCrown();
    bool cf = parameterTree->getTree()->getCrownFixed();
    for(size_t i = 0; i < logNodes.size(); i++) m.push_back(cf && logNodes[i] == crownNode);
    return m;
}

int FBDTreeModel::countResolvedSA(void){
    Tree* tree = parameterTree->getTree();
    int s = 0;
    for(Node* n : tree->getDownPassSequence())
        if(n->getIsTip() && n->getIsFossil() && n->getIsSA())
            s++;
    return s;
}

std::vector<double> FBDTreeModel::getParameterString(void){
    std::vector<double> vals = getRateString(a);
    if(psiRjMain >= 0){
        int cur = psiModel;
        auto at0 = [&](){ applyPsiModel(0); return lnLikelihood() + lnPsiMainPrior(); };
        auto at1 = [&](){ applyPsiModel(1); return lnLikelihood() + lnPsiRjTerm() - psiRjWeight + lnPsiMainPrior() + psiRjShiftSum; };
        double j0, j1;
        if(cur == 0){
            j0 = at0();
            scalePsiMain(1.0);
            j1 = at1();
        }else{
            j1 = at1();
            scalePsiMain(-1.0);
            j0 = at0();
        }
        setPsiMainValues(psiRjSaved);
        applyPsiModel(cur);
        vals.insert(vals.end(), {(double)cur, j1 - j0, psiRjWeight});
    }
    if(unresolvedFossils != nullptr)
        vals.push_back((double)unresolvedFossils->getNumSampledAncestors());
    else if(isResolved)
        vals.push_back((double)countResolvedSA());
    double off = UserSettings::userSettings().getAgeOffset();
    if(originAge != nullptr)
        vals.push_back(originAge->getValue() + off);
    for(Node* n : getAgeLogNodes())
        vals.push_back(n->getTime() + off);
    return vals;
}

std::vector<std::string> FBDTreeModel::getLatentNames(void){
    if(unresolvedFossils != nullptr && zoneIndexBuilt == false)
        buildZoneIndex();
    std::vector<std::string> names;
    auto unrName = [&](int i){ return (i < (int)unrFossilName.size()) ? unrFossilName[i] : std::to_string(i); };
    int nUnr = (unresolvedFossils != nullptr) ? unresolvedFossils->getNumFossils() : 0;
    for(int i = 0; i < nUnr; i++)
        if(unresolvedFossils->isUE(i) == false && unresolvedFossils->getYMin(i) < unresolvedFossils->getYMax(i))
            names.push_back("y_" + unrName(i));
    for(const BackboneFossil& bf : backboneFossils)
        if(bf.yMin < bf.yMax)
            names.push_back("y_" + bf.tip->getName());
    for(int i = 0; i < nUnr; i++)
        names.push_back("z_" + unrName(i));
    for(int i = 0; i < nUnr; i++)
        if(unresolvedFossils->isUE(i) == false)
            names.push_back("sa_" + unrName(i));
    for(int i = 0; i < nUnr; i++)
        if(unresolvedFossils->getAttachmentZoneDomain(i).size() > 1)
            names.push_back("az_" + unrName(i));
    return names;
}

std::vector<double> FBDTreeModel::getLatentString(void){
    std::vector<double> vals;
    double off = UserSettings::userSettings().getAgeOffset();
    int nUnr = (unresolvedFossils != nullptr) ? unresolvedFossils->getNumFossils() : 0;
    for(int i = 0; i < nUnr; i++)
        if(unresolvedFossils->isUE(i) == false && unresolvedFossils->getYMin(i) < unresolvedFossils->getYMax(i))
            vals.push_back(unresolvedFossils->getFossilAge(i) + off);
    for(const BackboneFossil& bf : backboneFossils)
        if(bf.yMin < bf.yMax)
            vals.push_back(bf.tip->getTime() + off);
    for(int i = 0; i < nUnr; i++)
        vals.push_back(unresolvedFossils->getAttachAge(i) + off);
    for(int i = 0; i < nUnr; i++)
        if(unresolvedFossils->isUE(i) == false)
            vals.push_back(unresolvedFossils->isSA(i) ? 1.0 : 0.0);
    for(int i = 0; i < nUnr; i++)
        if(unresolvedFossils->getAttachmentZoneDomain(i).size() > 1)
            vals.push_back((double)unresolvedFossils->getAttachmentZone(i));
    return vals;
}

std::vector<std::string> FBDTreeModel::getZoneLegend(void){
    std::vector<std::string> lines;
    if(unresolvedFossils == nullptr || zoneCladeName.empty())
        return lines;
    int nUnr = unresolvedFossils->getNumFossils();
    bool anyVar = false;
    for(int i = 0; i < nUnr; i++)
        if(unresolvedFossils->getAttachmentZoneDomain(i).size() > 1){ anyVar = true; break; }
    if(anyVar == false)
        return lines;
    lines.push_back("zone_index\tclade\ttype");
    for(int z = 0; z < (int)zoneCladeName.size(); z++){
        const char* ty = (zoneType[z] == 'C') ? "crown" : (zoneType[z] == 'S' ? "stem" : "total");
        lines.push_back(std::to_string(z) + "\t" + zoneCladeName[z] + "\t" + ty);
    }
    return lines;
}

double FBDTreeModel::lnLikelihood(void){
    static const bool noLik = (std::getenv("FBD_NO_LIK") != nullptr);
    if(noLik)
        return lnSupport();
    double v = calculateFBDProbability();
    if(std::isfinite(v) == false)
        return -INFINITY;
    return (psiModel == 1) ? v + psiRjLab : v;
}

double FBDTreeModel::lnSupport(void){
    if(isResolved)
        return 0.0;
    int nf = unresolvedFossils->getNumFossils();
    if(UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN){
        double x0 = originAge->getValue();
        if(x0 < parameterTree->getTree()->getCrown()->getTime())
            return -INFINITY;
        for(int i = 0; i < nf; i++)
            if(unresolvedFossils->getFossilAge(i) > x0 || unresolvedFossils->getAttachAge(i) > x0)
                return -INFINITY;
    }
    updateGammaCache();
    int spineIdx = unresolvedFossils->getSpineIdx();
    for(int i = 0; i < nf; i++)
        if(i != spineIdx && std::isfinite(cachedGammaLn[i]) == false)
            return -INFINITY;
    return 0.0;
}

double FBDTreeModel::evaluateAtGeneratingState(const std::vector<double>& lam, const std::vector<double>& muV,
                                               const std::vector<std::vector<double>>& psiV,
                                               const std::vector<double>& lamOU, const std::vector<double>& muOU,
                                               const std::vector<std::vector<double>>& psiOU, double x0,
                                               const std::vector<double>& nodeAges, const std::vector<double>& y,
                                               const std::vector<double>& z, const std::vector<char>& sa,
                                               const std::vector<int>& zoneNode, double& lnPriorOut){
    for(size_t i = 0; i < a.lambda.size(); i++)
        a.lambda[i]->setValue(lam[i]);
    for(size_t i = 0; i < a.mu.size(); i++)
        a.mu[i]->setValue(muV[i]);
    for(size_t t = 0; t < a.psi.size(); t++)
        for(size_t i = 0; i < a.psi[t].size(); i++)
            a.psi[t][i]->setValue(psiV[t][i]);
    if(a.lambdaField != nullptr)
        a.lambdaField->setState(lam, lamOU[0], lamOU[1], lamOU[2]);
    if(a.muField != nullptr)
        a.muField->setState(muV, muOU[0], muOU[1], muOU[2]);
    for(size_t t = 0; t < a.psiField.size(); t++)
        if(a.psiField[t] != nullptr)
            a.psiField[t]->setState(psiV[t], psiOU[t][0], psiOU[t][1], psiOU[t][2]);

    Tree* tr = parameterTree->getTree();
    if(originAge != nullptr){
        originAge->setValue(x0);
        Node* rt = tr->getRoot();
        if(rt->getDescendants().size() == 1)
            rt->setTime(x0);
    }
    std::vector<Node*> stack, order;
    stack.push_back(tr->getRoot());
    while(stack.empty() == false){
        Node* n = stack.back();
        stack.pop_back();
        std::vector<Node*>& ch = n->getDescendants();
        if(ch.size() >= 2)
            order.push_back(n);
        for(int k = (int)ch.size() - 1; k >= 0; k--)
            stack.push_back(ch[k]);
    }
    for(size_t i = 0; i < order.size(); i++)
        order[i]->setTime(nodeAges[i]);
    parameterTree->commitTree();

    if(unresolvedFossils != nullptr){
        unresolvedFossils->setFossilAges(y);
        int nf = unresolvedFossils->getNumFossils();
        for(int i = 0; i < nf; i++){
            if(sa[i])
                unresolvedFossils->flipSA(i, true);
            else
                unresolvedFossils->setAttachAge(i, z[i]);
            unresolvedFossils->acceptFossil(i);
        }
    }
    invalidateGammaCache();
    zoneInit = false;
    zoneIndexBuilt = false;
    eulerBuilt = false;
    if(zoneNode.empty() == false){
        buildZoneIndex();
        int rootOff = tr->getRoot()->getOffset();
        for(size_t i = 0; i < zoneNode.size(); i++)
            unresolvedFossils->initAttachmentZone((int)i, zoneNode[i] == rootOff ? trunkMinZone : minZoneOfNode[zoneNode[i]]);
        zoneInit = true;
    }

    lnPriorOut = lnPriorProbability();
    return lnLikelihood();
}

double FBDTreeModel::lnPriorProbability(void){
    double lnP = 0.0;
    for(Parameter* p : parameters)
        if( p != parameterTree)
            lnP += p->lnProbability();

    UserSettings& us = UserSettings::userSettings();
    if(us.getConditioning() == Conditioning::CROWN){
        double crownAge = parameterTree->getTree()->getCrown()->getTime();
        lnP += Probability::priorLnPdf(us.getConditionAgePrior(), us.getConditionAgePriorP1(), us.getConditionAgePriorP2(), crownAge, 0.0, std::numeric_limits<double>::max(), us.getConditionAgePriorP3());
    }
    if(psiRjMain >= 0)
        lnP += lnPsiRjTerm();

    return lnP;
}

void FBDTreeModel::print(void){
    if(UserSettings::userSettings().getArLog() == false)
        return;
    std::ostream& os = std::cout;
    std::ios_base::fmtflags f = os.flags();
    std::streamsize pr = os.precision();
    os << std::fixed << std::setprecision(3);
    for(Parameter* p : parameters){
        if(p != parameterTree){
            std::cout << p->getName() << " (A/R): " << p->getAcceptanceRatio();
            double st = p->getStep();
            if(st >= 0.0) std::cout << " step: " << st;
            std::cout << "\t";
        }
    }
    if(a.lambdaField) a.lambdaField->printPerBinAccept(std::cout, "lambda");
    if(a.muField) a.muField->printPerBinAccept(std::cout, "mu");
    for(size_t t = 0; t < a.psiField.size(); t++)
        if(a.psiField[t]) a.psiField[t]->printPerBinAccept(std::cout, "psi");
    if(unresolvedFossils != nullptr)
        unresolvedFossils->print();
    static const char* tmName[TM_COUNT] = {"NE","WB","WE","treeScale","SARJ","upDown","jointScale","subtree","nodeAge","crown"};
    std::cout << "tree[";
    for(int i = 0; i < TM_COUNT; i++)
        if(tmAtt[i] > 0)
            std::cout << tmName[i] << ":" << (double)tmAcc[i] / (double)tmAtt[i] << "[" << tmAcc[i] << "/" << tmAtt[i] << "] ";
    std::cout << "] treeScaleStep: " << parameterTree->getScaleLambda()
              << " subtreeStep: " << subtreeStep << " jointScaleStep: " << jointScaleStep;
    if(rsTot > 0)
        std::cout << " rateShift[A/R:" << (double)rsAcc / (double)rsTot
                  << " step:" << std::setprecision(5) << shiftStep << std::setprecision(3)
                  << " batch:" << saBatch << "]";
    if(rvAtt > 0)
        std::cout << " rateVecScale[A/R:" << (double)rvAcc / (double)rvAtt << " " << rvAcc << "/" << rvAtt << " step:" << rateVecStep << "]";
    if(seAtt > 0)
        std::cout << " rateShrinkExpand[A/R:" << (double)seAcc / (double)seAtt << " " << seAcc << "/" << seAtt << " step:" << shrinkStep << "]";
    if(azAtt > 0)
        std::cout << " azGibbs[A/R:" << (double)azAcc / (double)azAtt << " " << azAcc << "/" << azAtt << "]";
    if(rjAtt > 0)
        std::cout << " psiRj[A/R:" << (double)rjAcc / (double)rjAtt << " " << rjAcc << "/" << rjAtt << "]";
    if(gammaStaleCalls > 0)
        std::cout << " gammaStale[mean:" << (double)gammaStaleTot / (double)gammaStaleCalls
                  << " of " << unresolvedFossils->getNumFossils() << "]";
    std::cout << "\n";
    os.flags(f);
    os.precision(pr);
}

static int numFree(const std::vector<ParameterDouble*>& v){
    int k = 0;
    for(ParameterDouble* p : v)
        if(p->getProposalProbability() > 0.0)
            k++;
    return k;
}

void FBDTreeModel::pickRateTarget(std::vector<ParameterDouble*>*& outVec, ParameterOUField*& outFld){
    outVec = nullptr;
    outFld = nullptr;
    std::vector<std::vector<ParameterDouble*>*> vcands;
    std::vector<ParameterOUField*> fcands;
    if(a.lambdaField == nullptr && numFree(a.lambda) >= 2) vcands.push_back(&a.lambda);
    if(a.muField == nullptr && numFree(a.mu) >= 2) vcands.push_back(&a.mu);
    for(int tp = 0; tp < a.numPsiTypes; tp++){
        if(a.psiField[tp] == nullptr && numFree(a.psi[tp]) >= 2) vcands.push_back(&a.psi[tp]);
    }
    if(a.lambdaField != nullptr) fcands.push_back(a.lambdaField);
    if(a.muField != nullptr) fcands.push_back(a.muField);
    for(int tp = 0; tp < a.numPsiTypes; tp++)
        if(a.psiField[tp] != nullptr) fcands.push_back(a.psiField[tp]);
    int n = (int)(vcands.size() + fcands.size());
    if(n == 0) return;
    int idx = (int)(rng.uniformRv() * n);
    if(idx < (int)vcands.size()) outVec = vcands[idx];
    else outFld = fcands[idx - (int)vcands.size()];
}

double FBDTreeModel::doRateVectorScale(void){
    lastMoveKind = MK_RATEVEC;
    lastRateVecScale = true;
    pickRateTarget(lastRateVec, lastRateField);
    double mB = 0.95;
    double d = mB + Probability::Normal::rv(&rng) * std::sqrt(1.0 - mB * mB);
    if(rng.uniformRv() < 0.5) d = -d;
    double c = std::exp(rateVecStep * d);
    double h;
    if(lastRateField != nullptr){
        h = lastRateField->scaleAllProposed(c, true);
    }else{
        int nFree = 0;
        for(ParameterDouble* p : *lastRateVec)
            if(p->getProposalProbability() > 0.0){
                p->scaleProposed(c);
                nFree++;
            }
        h = (double)nFree * std::log(c);
    }
    rvAttW++;
    rvAtt++;
    if(rvAttW >= 200){
        double gain = 1.0 / std::sqrt((double)(rvAtt / 200));
        rateVecStep *= std::exp(gain * ((double)rvAccW / rvAttW - 0.3));
        if(rateVecStep < 1e-3) rateVecStep = 1e-3;
        if(rateVecStep > 10.0)  rateVecStep = 10.0;
        rvAccW = 0;
        rvAttW = 0;
    }
    return h;
}

double FBDTreeModel::doRateShrinkExpand(void){
    lastMoveKind = MK_RATEVEC;
    lastRateVecScale = false;
    pickRateTarget(lastRateVec, lastRateField);
    int n = 0;
    double logMean = 0.0;
    if(lastRateField == nullptr){
        n = numFree(*lastRateVec);
        for(ParameterDouble* p : *lastRateVec)
            if(p->getProposalProbability() > 0.0)
                logMean += std::log(p->getValue());
        logMean /= (double)n;
    }
    double mB = 0.95;
    double d = mB + Probability::Normal::rv(&rng) * std::sqrt(1.0 - mB * mB);
    if(rng.uniformRv() < 0.5) d = -d;
    double a = std::exp(shrinkStep * d);
    double h;
    if(lastRateField != nullptr){
        h = lastRateField->shrinkExpandProposed(a);
    }else{
        for(ParameterDouble* p : *lastRateVec){
            if(p->getProposalProbability() == 0.0)
                continue;
            double cur = p->getValue();
            double target = std::exp(logMean + a * (std::log(cur) - logMean));
            p->scaleProposed(target / cur);
        }
        h = (double)(n - 1) * std::log(a);
    }
    seAttW++;
    seAtt++;
    if(seAttW >= 200){
        double gain = 1.0 / std::sqrt((double)(seAtt / 200));
        shrinkStep *= std::exp(gain * ((double)seAccW / seAttW - 0.3));
        if(shrinkStep < 1e-3) shrinkStep = 1e-3;
        if(shrinkStep > 10.0)  shrinkStep = 10.0;
        seAccW = 0;
        seAttW = 0;
    }
    return h;
}

double FBDTreeModel::doRateShift(void){
    lastMoveKind = MK_RATESHIFT;
    if(unresolvedFossils != nullptr)
        unresolvedFossils->beginBatchMove();

    double m = 0.95;
    double b = m + Probability::Normal::rv(&rng) * std::sqrt(1.0 - m * m);
    if(rng.uniformRv() < 0.5)
        b = -b;
    double d = shiftStep * b;

    rsAttW++;
    rsTot++;
    if(rsAttW >= 200){
        double ar = (double)rsAccW / (double)rsAttW;
        rsAdapt++;
        double gain = 1.0 / std::sqrt((double)rsAdapt);
        shiftStep *= std::exp(gain * (ar - 0.3));
        if(shiftStep < 1e-5) shiftStep = 1e-5;
        if(shiftStep > 0.5)  shiftStep = 0.5;
        saBatchF *= std::exp(gain * (ar - 0.3));
        if(saBatchF < 1.0) saBatchF = 1.0;
        saBatch = (int)(saBatchF + 0.5);
        rsAccW = 0;
        rsAttW = 0;
    }

    double lnH = 0.0;
    if(a.lambdaField != nullptr){
        double j = a.lambdaField->shiftRates(d);
        if(std::isinf(j)) return -std::numeric_limits<double>::infinity();
        lnH += j;
    }else{
        for(ParameterDouble* l : a.lambda)
            if(l->getValue() + d <= 0.0)
                return -std::numeric_limits<double>::infinity();
        for(ParameterDouble* l : a.lambda)
            l->shiftProposed(d);
    }
    if(a.muField != nullptr){
        double j = a.muField->shiftRates(d);
        if(std::isinf(j)) return -std::numeric_limits<double>::infinity();
        lnH += j;
    }else{
        for(ParameterDouble* mv : a.mu)
            if(mv->getValue() + d <= 0.0)
                return -std::numeric_limits<double>::infinity();
        for(ParameterDouble* mv : a.mu)
            mv->shiftProposed(d);
    }

    if(unresolvedFossils == nullptr)
        return lnH;

    bool toSA = (d < 0.0);
    std::vector<int> pool;
    int nOther = 0;
    int nf = unresolvedFossils->getNumFossils();
    for(int i = 0; i < nf; i++){
        if(unresolvedFossils->saEligible(i) == false)
            continue;
        if(unresolvedFossils->isSA(i) == toSA) nOther++;
        else                                   pool.push_back(i);
    }
    int B = saBatch;
    if(B > (int)pool.size())
        return -std::numeric_limits<double>::infinity();

    for(int k = 0; k < B; k++){
        int r = k + (int)(rng.uniformRv() * (double)(pool.size() - k));
        std::swap(pool[k], pool[r]);
        lnH += unresolvedFossils->flipSA(pool[k], toSA);
    }
    lnH += lnChoose((int)pool.size(), B) - lnChoose(nOther + B, B);
    rebuildStalkIndex();
    return lnH;
}

double FBDTreeModel::lnChoose(int n, int k){
    return std::lgamma((double)n + 1.0) - std::lgamma((double)k + 1.0) - std::lgamma((double)(n - k) + 1.0);
}

// crossing count: intervals (lo<hi) with lo<zq<hi = #{lo<zq} - #{hi<=zq}
namespace {
void eraseStalkZ(std::vector<double>& sz, double z){
    sz.erase(std::lower_bound(sz.begin(), sz.end(), z));
}
void insertStalkZ(std::vector<double>& sz, double z){
    sz.insert(std::lower_bound(sz.begin(), sz.end(), z), z);
}
void eraseStalk(std::vector<double>& sy, std::vector<double>& sz, double y, double z){
    sy.erase(std::lower_bound(sy.begin(), sy.end(), y));
    eraseStalkZ(sz, z);
}
void insertStalk(std::vector<double>& sy, std::vector<double>& sz, double y, double z){
    sy.insert(std::lower_bound(sy.begin(), sy.end(), y), y);
    insertStalkZ(sz, z);
}
void eraseZIdx(std::vector<std::pair<double,int> >& v, double z, int i){
    v.erase(std::lower_bound(v.begin(), v.end(), std::make_pair(z, i)));
}
void insertZIdx(std::vector<std::pair<double,int> >& v, double z, int i){
    v.insert(std::lower_bound(v.begin(), v.end(), std::make_pair(z, i)), std::make_pair(z, i));
}
void stalkChangeBands(bool saO, double yO, double zO, bool saN, double yN, double zN,
                      double& b1lo, double& b1hi, double& b2lo, double& b2hi){
    b1lo = b1hi = b2lo = b2hi = 0.0;
    if(saO && saN)
        return;
    if(saO){ b1lo = yN; b1hi = zN; return; }
    if(saN){ b1lo = yO; b1hi = zO; return; }
    b1lo = std::min(yO, yN); b1hi = std::max(yO, yN);
    b2lo = std::min(zO, zN); b2hi = std::max(zO, zN);
}
void collectInBand(const std::vector<std::pair<double,int> >& v, double lo, double hi, std::vector<int>& out){
    if(hi <= lo)
        return;
    std::vector<std::pair<double,int> >::const_iterator b =
        std::lower_bound(v.begin(), v.end(), std::make_pair(lo, std::numeric_limits<int>::min()));
    std::vector<std::pair<double,int> >::const_iterator e =
        std::upper_bound(v.begin(), v.end(), std::make_pair(hi, std::numeric_limits<int>::max()));
    for(; b != e; ++b)
        out.push_back(b->second);
}
int countStraddling(const std::vector<double>& los, const std::vector<double>& his, double zq){
    int below = (int)(std::lower_bound(los.begin(), los.end(), zq) - los.begin());
    int above = (int)(std::upper_bound(his.begin(), his.end(), zq) - his.begin());
    return below - above;
}
}

double FBDTreeModel::zoneBackboneEdges(int mz, double z){
    double count = (double)countStraddling(zoneEdges[mz].yng, zoneEdges[mz].old, z);
    if(mz == trunkMinZone){
        Tree* tree = parameterTree->getTree();
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

double FBDTreeModel::update(void){
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);

    if(isResolved == false && unresolvedFossils != nullptr)
        updateGammaCache();

    lastMoveKind = MK_PARAM;
    lastTreeMove = TM_NONE;
    if(psiRjMain >= 0 && rng.uniformRv() < 0.05){
        lastMoveKind = MK_PSIRJ;
        prevPsiModel = psiModel;
        applyPsiModel(1 - psiModel);
        double sgn = (prevPsiModel == 0) ? 1.0 : -1.0;
        scalePsiMain(sgn);
        rjAtt++;
        RandomVariable::setActiveInstance(prevRng);
        return sgn * psiRjShiftSum;
    }
    if(haveIid && rng.uniformRv() < 0.20){
        double r = (rng.uniformRv() < 0.5) ? doRateVectorScale() : doRateShrinkExpand();
        RandomVariable::setActiveInstance(prevRng);
        return r;
    }
    if(haveLMField && rng.uniformRv() < 0.10){
        double r = doRateShift();
        RandomVariable::setActiveInstance(prevRng);
        return r;
    }

    if(azGibbsIdx.empty() == false && rng.uniformRv() < 0.15){
        lastMoveKind = MK_AZGIBBS;
        azAtt++;
        double r = doAttachmentZoneGibbs();
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
        double blockUnit = 0.1 * numSlideable;
        if(blockUnit < fixedWeight) blockUnit = fixedWeight;
        double slideAndCrown = numSlideable + fixedWeight;
        double subtreeW = (numSlideable > 0) ? blockUnit : 0.0;
        double uMove = rng.uniformRv() * (slideAndCrown + 2.0 * blockUnit + subtreeW);
        if(uMove >= slideAndCrown + blockUnit + subtreeW){
            lastTreeMove = TM_UPDOWN;
            double r = doUpDownScale();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove >= slideAndCrown){
            lastMoveKind = MK_JOINTSCALE;
            double r;
            if(uMove < slideAndCrown + blockUnit){ lastTreeMove = TM_JOINTSCALE; r = doJointScale(); }
            else { lastTreeMove = TM_SUBTREE; r = doSubtreeScale(); }
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        std::map<Node*,double> floors;
        computeAgeFloors(floors);
        parameterTree->setAgeFloors(floors);
    }
    else if(updatedParameter == parameterTree && isResolved){ // Not deployed
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
            lastTreeMove = TM_NE;
            double r = doNarrowExchange();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove < wNE + wWB){
            lastTreeMove = TM_WB;
            double r = doWilsonBalding();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove < wNE + wWB + wWE){
            lastTreeMove = TM_WE;
            double r = doWideExchange();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove < wNE + wWB + wWE + wTreeScale){
            lastTreeMove = TM_TREESCALE;
            double r = doTreeScale();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove < wNE + wWB + wWE + wTreeScale + wSA){
            lastTreeMove = TM_SARJ;
            double r = doSARJMCMC();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
        if(uMove < wNE + wWB + wWE + wTreeScale + wSA + wUpDown){
            lastTreeMove = TM_UPDOWN;
            double r = doUpDownScale();
            RandomVariable::setActiveInstance(prevRng);
            return r;
        }
    }
    double ratio = updatedParameter->update();
    if(updatedParameter == parameterTree)
        lastTreeMove = parameterTree->getTree()->getLastUpdateWasScale() ? TM_CROWN : TM_NODEAGE;
    if((lastTreeMove == TM_NODEAGE || lastTreeMove == TM_CROWN) && isResolved == false && unresolvedFossils != nullptr){
        Tree* t = parameterTree->getTree();
        if(t->getLastAgeNode() != nullptr)
            ratio += relabelAcrossNode(t->getLastAgeNode(), t->getLastAgeNodeOldAge());
    }

    RandomVariable::setActiveInstance(prevRng);
    return ratio;
}

void FBDTreeModel::updateForAcceptance(void){
    if(lastMoveKind == MK_PSIRJ){
        rjAcc++;
        return;
    }
    if(lastTreeMove != TM_NONE){ tmAtt[lastTreeMove]++; tmAcc[lastTreeMove]++; }
    if(lastTreeMove == TM_SUBTREE) subtreeAccW++;
    else if(lastTreeMove == TM_JOINTSCALE) jointScaleAccW++;
    else if(lastTreeMove == TM_NODEAGE){ parameterTree->recordNodeAgeMove(true); commitZoneBlock(); }
    else if(lastTreeMove == TM_CROWN) commitZoneBlock();
    if(lastMoveKind == MK_AZGIBBS){
        unresolvedFossils->updateForAcceptance();
        azAcc++;
        return;
    }
    if(lastMoveKind == MK_RATESHIFT){
        if(a.lambdaField != nullptr) a.lambdaField->commitProposed();
        else for(ParameterDouble* l : a.lambda) l->commitProposed();
        if(a.muField != nullptr) a.muField->commitProposed();
        else for(ParameterDouble* m : a.mu) m->commitProposed();
        if(unresolvedFossils != nullptr) unresolvedFossils->updateForAcceptance();
        rsAccW++;
        rsAcc++;
        return;
    }
    if(lastMoveKind == MK_RATEVEC){
        if(lastRateField != nullptr) lastRateField->commitProposed();
        else for(ParameterDouble* p : *lastRateVec) p->commitProposed();
        if(lastRateVecScale){ rvAccW++; rvAcc++; } else { seAccW++; seAcc++; }
        return;
    }
    if(lastMoveKind == MK_UPDOWN){
        for(ParameterDouble* l : a.lambda) l->commitProposed();
        for(ParameterDouble* m : a.mu) m->commitProposed();
        for(auto& pv : a.psi) for(ParameterDouble* p : pv) p->commitProposed();
        if(a.lambdaField != nullptr) a.lambdaField->commitProposed();
        if(a.muField != nullptr) a.muField->commitProposed();
        for(int tp = 0; tp < a.numPsiTypes; tp++) if(a.psiField[tp] != nullptr) a.psiField[tp]->commitProposed();
        if(originAge != nullptr) originAge->commitProposed();
        parameterTree->updateForAcceptance();
        if(unresolvedFossils != nullptr) unresolvedFossils->updateForAcceptance();
        upDownTotal++;
        upDownRecent.push_back(true);
        if(upDownRecent.size() > 1000) upDownRecent.pop_front();
    }else if(lastMoveKind == MK_JOINTSCALE){
        parameterTree->updateForAcceptance();
        unresolvedFossils->updateForAcceptance();
    }else{
        updatedParameter->updateForAcceptance();
    }
    if(isResolved && originAge != nullptr)
        parameterTree->getTree()->getRoot()->setTime(originAge->getValue());
}

void FBDTreeModel::updateForRejection(void){
    if(lastMoveKind == MK_PSIRJ){
        applyPsiModel(prevPsiModel);
        setPsiMainValues(psiRjSaved);
        return;
    }
    if(lastTreeMove != TM_NONE) tmAtt[lastTreeMove]++;
    if(lastTreeMove == TM_NODEAGE){ parameterTree->recordNodeAgeMove(false); restoreZoneBlock(); }
    if(lastTreeMove == TM_CROWN) restoreZoneBlock();
    if(lastMoveKind == MK_AZGIBBS){
        unresolvedFossils->updateForRejection();
        rebuildStalkIndex();
        return;
    }
    if(lastMoveKind == MK_RATESHIFT){
        if(a.lambdaField != nullptr) a.lambdaField->restoreProposed();
        else for(ParameterDouble* l : a.lambda) l->restoreProposed();
        if(a.muField != nullptr) a.muField->restoreProposed();
        else for(ParameterDouble* m : a.mu) m->restoreProposed();
        if(unresolvedFossils != nullptr){
            unresolvedFossils->updateForRejection();
            rebuildStalkIndex();
        }
        return;
    }
    if(lastMoveKind == MK_RATEVEC){
        if(lastRateField != nullptr) lastRateField->restoreProposed();
        else for(ParameterDouble* p : *lastRateVec) p->restoreProposed();
        return;
    }
    if(lastMoveKind == MK_UPDOWN){
        for(ParameterDouble* l : a.lambda) l->restoreProposed();
        for(ParameterDouble* m : a.mu) m->restoreProposed();
        for(auto& pv : a.psi) for(ParameterDouble* p : pv) p->restoreProposed();
        if(a.lambdaField != nullptr) a.lambdaField->restoreProposed();
        if(a.muField != nullptr) a.muField->restoreProposed();
        for(int tp = 0; tp < a.numPsiTypes; tp++) if(a.psiField[tp] != nullptr) a.psiField[tp]->restoreProposed();
        if(originAge != nullptr) originAge->restoreProposed();
        parameterTree->updateForRejection();
        if(unresolvedFossils != nullptr) unresolvedFossils->updateForRejection();
        upDownTotal++;
        upDownRecent.push_back(false);
        if(upDownRecent.size() > 1000) upDownRecent.pop_front();
    }else if(lastMoveKind == MK_JOINTSCALE){
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
    os << rateVecStep << ' ' << shrinkStep << ' ' << shiftStep << ' ' << saBatch << ' ' << saBatchF << ' ' << upDownStep << ' ' << upDownTotal << '\n';
    os << rvAccW << ' ' << rvAttW << ' ' << rvAtt << ' ' << seAccW << ' ' << seAttW << ' ' << seAtt << ' ' << rsAccW << ' ' << rsAttW << '\n';
    os << rsAcc << ' ' << rsTot << ' ' << rsAdapt << ' ' << azAcc << ' ' << azAtt << '\n';
    Serialize::writeBoolDeque(os, upDownRecent);
    for(int i = 0; i < TM_COUNT; i++) os << tmAcc[i] << ' ' << tmAtt[i] << ' ';
    os << '\n';
    os << subtreeStep << ' ' << subtreeAccW << ' ' << subtreeAttW << ' ' << subtreeAdapt << ' ' << jointScaleStep << ' ' << jointScaleAccW << ' ' << jointScaleAttW << ' ' << jointScaleAdapt << '\n';
    if(psiRjMain >= 0){
        os << "psi_rj " << psiModel << ' ' << rjAcc << ' ' << rjAtt;
        for(double x : getPsiRjState())
            os << ' ' << x;
        os << '\n';
    }
}

void FBDTreeModel::readState(std::istream& is){
    for(Parameter* p : parameters)
        p->readState(is);
    is >> rateVecStep >> shrinkStep >> shiftStep >> saBatch >> saBatchF >> upDownStep >> upDownTotal;
    is >> rvAccW >> rvAttW >> rvAtt >> seAccW >> seAttW >> seAtt >> rsAccW >> rsAttW;
    is >> rsAcc >> rsTot >> rsAdapt >> azAcc >> azAtt;
    Serialize::readBoolDeque(is, upDownRecent);
    for(int i = 0; i < TM_COUNT; i++) is >> tmAcc[i] >> tmAtt[i];
    is >> subtreeStep >> subtreeAccW >> subtreeAttW >> subtreeAdapt >> jointScaleStep >> jointScaleAccW >> jointScaleAttW >> jointScaleAdapt;
    if(psiRjMain >= 0){
        std::string tag;
        int m;
        is >> tag >> m >> rjAcc >> rjAtt;
        if(tag != "psi_rj")
            Msg::error("the checkpoint has no psi_rj line; it was not written by a psi_rj_main run.");
        std::vector<double> st(7 + 3 * psiRjShift.size());
        for(double& x : st)
            is >> x;
        setPsiRjState(st);
        applyPsiModel(m);
    }
    cacheInit = false;
    zoneInit = true;
}

double FBDTreeModel::calculateFBDProbability(void){
    Tree* tree = parameterTree->getTree();

    if(isResolved && originAge != nullptr)
        tree->getRoot()->setTime(originAge->getValue());

    int numInternalNodes = tree->getNumNodes() - tree->getNumBackbone();
    double crownAge = tree->getCrown()->getTime();
    std::vector<Node*>& dpseq = tree->getDownPassSequence();

    prepareIntervals();

    if(isResolved)
        return calculateResolvedFBD();

    bool useOrigin = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);
    double fbdProb = 0.0;
    static const bool fbdDbg = (std::getenv("FBD_DBG") != nullptr);
    double dbgT1 = 0.0, dbgT12 = 0.0;
    fbdDUses.assign(a.intervalStart.size(), 0);
    fbdDecUses.assign(a.intervalStart.size(), 0);
    std::vector<int>& lnDPrevUses = fbdDUses;
    std::vector<int>& decUses = fbdDecUses;
    Node* crownCond = useOrigin ? nullptr : tree->getCrown();

    //term 1: conditioning
    if(useOrigin){
        double x0 = originAge->getValue();
        if(x0 < crownAge)
            return -INFINITY;
        for(int i = 0; i < unresolvedFossils->getNumFossils(); i++){
            if(unresolvedFossils->getFossilAge(i) > x0)
                return -INFINITY;
            if(unresolvedFossils->getAttachAge(i) > x0)
                return -INFINITY;
        }
        int k = findIndex(x0);
        if(UserSettings::userSettings().getConditioningEvent() == ConditioningEvent::SURVIVAL)
            fbdProb += a.lnRPrev[k] + incDecLn(k, x0);
        else{
            fbdProb += lnDLocal(x0) - calculateLnConditioning(x0);
            lnDPrevUses[k]++;
        }
    }else{
        int k = findIndex(crownAge);
        fbdProb += a.lnRPrev[k] - std::log(rho) + incDecLn(k, crownAge) + lnDPlusDec(k, crownAge);
        lnDPrevUses[k]++;
        decUses[k]++;
    }

    dbgT1 = fbdProb;
    //term 2: main body
    int nDp = (int)dpseq.size();
    fbdTermNode.assign(nDp, 0.0);
    std::vector<double>& termNode = fbdTermNode;
    fbdNodeBin.resize(nDp);
    std::vector<int>& nodeBin = fbdNodeBin;
    ThreadPool::current().parallelFor(OP_FBD, nDp, [&](int lo, int hi){
        for(int idx = lo; idx < hi; idx++){
            nodeBin[idx] = -1;
            Node* n = dpseq[idx];
            if(n->getIsTip() || n == crownCond)
                continue;
            bool hasChild = false;
            for(Node* c : n->getNeighbors())
                if(c != n->getAncestor()){ hasChild = true; break; }
            if(hasChild){
                nodeBin[idx] = findIndex(n->getTime());
                termNode[idx] = std::log(a.lamVec[nodeBin[idx]]) + lnDLocal(n->getTime());
            }
        }
    });
    for(int idx = 0; idx < nDp; idx++){
        fbdProb += termNode[idx];
        if(nodeBin[idx] >= 0)
            lnDPrevUses[nodeBin[idx]]++;
    }
    fbdProb += numExtantTips * std::log(rho);
    dbgT12 = fbdProb;

    //term 3: fossil attachment
    int numFossils = unresolvedFossils->getNumFossils();
    updateGammaCache();
    int spineIdx = unresolvedFossils->getSpineIdx();
    std::vector<double> sig;
    sig.reserve((2 + a.numPsiTypes) * a.c1Vec.size());
    for(size_t q = 0; q < a.c1Vec.size(); q++){
        sig.push_back(a.lamVec[q]); sig.push_back(a.muVec[q]);
        for(int t = 0; t < a.numPsiTypes; t++)
            sig.push_back(a.lnPsiTypeVec[t][q]);
    }
    bool ratesSame = a.fbdMemoInit && sig == a.fbdRateSig;
    if(a.fbdMemoInit == false || (int)a.fbdTermFoss.size() != numFossils){
        a.fbdTermFoss.assign(numFossils, 0.0);
        a.fbdPrevGammaLn.assign(numFossils, std::numeric_limits<double>::quiet_NaN());
        a.fbdPrevFy.assign(numFossils, std::numeric_limits<double>::quiet_NaN());
        a.fbdPrevFz.assign(numFossils, std::numeric_limits<double>::quiet_NaN());
        a.fbdPrevKind.assign(numFossils, -1);
        a.fbdKy.assign(numFossils, 0);
        a.fbdKz.assign(numFossils, 0);
        ratesSame = false;
        a.fbdMemoInit = true;
    }
    a.fbdRateSig.swap(sig);

    std::vector<double>& termFoss = a.fbdTermFoss;
    ThreadPool::current().parallelFor(OP_FBD, numFossils, [&](int lo, int hi){
        for(int i = lo; i < hi; i++){
            double gy = unresolvedFossils->getFossilAge(i);
            double gz = unresolvedFossils->getAttachAge(i);
            int kind = (unresolvedFossils->isSA(i) ? 2 : 0) | (i == spineIdx ? 1 : 0);
            if(ratesSame && cachedGammaLn[i] == a.fbdPrevGammaLn[i] && gy == a.fbdPrevFy[i]
               && gz == a.fbdPrevFz[i] && kind == a.fbdPrevKind[i])
                continue;
            int ky = (gy == a.fbdPrevFy[i]) ? a.fbdKy[i] : findIndex(gy);
            int kz = (gz == a.fbdPrevFz[i]) ? a.fbdKz[i] : findIndex(gz);
            termFoss[i] = fossilTermLn(i, spineIdx, ky, kz, true);
            a.fbdPrevGammaLn[i] = cachedGammaLn[i];
            a.fbdPrevFy[i] = gy;
            a.fbdPrevFz[i] = gz;
            a.fbdKy[i] = ky;
            a.fbdKz[i] = kz;
            a.fbdPrevKind[i] = kind;
        }
    });
    for(int i = 0; i < numFossils; i++)
        fbdProb += termFoss[i];
    for(int i : ueFossils)
        if(i != spineIdx && unresolvedFossils->isSA(i) == false)
            lnDPrevUses[a.fbdKz[i]]++;
    if(spineIdx >= 0 && unresolvedFossils->isUE(spineIdx) == false && unresolvedFossils->isSA(spineIdx) == false)
        lnDPrevUses[a.fbdKy[spineIdx]]--;
    for(size_t k = 0; k < backboneFossils.size(); k++){
        double y = backboneFossils[k].tip->getTime();
        fbdProb += std::log(psiOfTypeAt(a.backboneType[k], findIndex(y))) + lnP0At(findIndex(y), y) - lnDLocal(y);
        lnDPrevUses[findIndex(y)]--;
    }
    // count x lnDInc per bin so lnDPrev offsets that cancel never gets added; absolue lnD terms lose every digits when c1*dt is large
    int olderD = 0, olderDec = 0;
    for(int j = (int)lnDPrevUses.size() - 2; j >= 0; j--){
        olderD += lnDPrevUses[j + 1];
        olderDec += decUses[j + 1];
        int pairs = (olderD > 0) ? std::min(olderD, olderDec) : 0;
        if(pairs > 0)
            fbdProb += pairs * a.lnDecPair[j];
        if(olderD != pairs)
            fbdProb += (olderD - pairs) * a.lnDInc[j];
        if(olderDec != pairs)
            fbdProb += (olderDec - pairs) * a.lnDecInc[j];
    }
    if(fbdDbg){
        double t2 = dbgT12 - dbgT1;
        double t3 = fbdProb - dbgT12;
        int amax = 0; double vmax = -1e300;
        for(int idx = 0; idx < nDp; idx++) if(termNode[idx] > vmax){ vmax = termNode[idx]; amax = idx; }
        double fsum = 0.0, fmin = 1e300, fmax = -1e300;
        for(int i = 0; i < numFossils; i++){ fsum += termFoss[i]; if(termFoss[i] < fmin) fmin = termFoss[i]; if(termFoss[i] > fmax) fmax = termFoss[i]; }
        double tn = dpseq[amax]->getTime();
        std::cout << "[FBD_DBG] total=" << fbdProb << " t1cond=" << dbgT1 << " t2spec=" << t2 << " t3foss=" << t3
                  << " | maxNodeTerm=" << vmax << " @t=" << tn << " bin=" << findIndex(tn) << " lam=" << lambdaAt(findIndex(tn)) << " lnD=" << lnD(tn)
                  << " | fossSum=" << fsum << " fossMin=" << fmin << " fossMax=" << fmax << "\n";
        int fam = 0; double fv = -1e300;
        for(int i = 0; i < numFossils; i++) if(termFoss[i] > fv){ fv = termFoss[i]; fam = i; }
        double yy = unresolvedFossils->getFossilAge(fam), zz = unresolvedFossils->getAttachAge(fam);
        int ty = a.fossilType[fam];
        std::cout << "[FBD_DBG_FOSS] i=" << fam << " term=" << fv << " isSA=" << unresolvedFossils->isSA(fam) << " isUE=" << unresolvedFossils->isUE(fam)
                  << " y=" << yy << " z=" << zz << " biny=" << findIndex(yy) << " binz=" << findIndex(zz)
                  << " logPsiY=" << std::log(psiOfTypeAt(ty, findIndex(yy))) << " log2LamZ=" << std::log(2 * lambdaAt(findIndex(zz)))
                  << " logP0Y=" << std::log(calculateP0(yy)) << " lnDz=" << lnD(zz) << " lnDy=" << lnD(yy) << " gammaLn=" << cachedGammaLn[fam] << "\n";
        int ni = (int)a.intervalStart.size();
        for(int i = 0; i < ni; i++){
            double p0 = calculateP0At(i, (i + 1 < ni) ? a.intervalStart[i + 1] : a.intervalStart[i] + 1.0);
            std::cout << "[FBD_DBG_IV] i=" << i << " start=" << a.intervalStart[i] << " lam=" << lambdaAt(i) << " mu=" << muAt(i) << " psi=" << psiTotalAt(i)
                      << " c1=" << a.c1Vec[i] << " P0end=" << p0 << " dPrev=" << a.dPrev[i] << "\n";
        }
        std::vector<double> divc(ni, 0.0), stalk(ni, 0.0), fy(ni, 0.0);
        double xtop = useOrigin ? originAge->getValue() : crownAge;
        for(int idx = 0; idx < nDp; idx++){
            Node* nn = dpseq[idx];
            if(nn->getIsTip()) continue;
            bool hc = false;
            for(Node* c : nn->getNeighbors()) if(c != nn->getAncestor()){ hc = true; break; }
            if(hc){ int bi = findIndex(nn->getTime()); if(bi < ni) divc[bi] += 1.0; }
        }
        for(int i = 0; i < numFossils; i++){
            double yv = unresolvedFossils->getFossilAge(i);
            int by = findIndex(yv); if(by < ni) fy[by] += 1.0;
            if(unresolvedFossils->isSA(i)) continue;
            double zv = unresolvedFossils->getAttachAge(i);
            int bz = findIndex(zv); if(bz < ni) divc[bz] += 1.0;
            for(int b = 0; b < ni; b++){
                double blo = a.intervalStart[b], bhi = (b + 1 < ni) ? a.intervalStart[b + 1] : xtop;
                double ov = std::min(zv, bhi) - std::max(yv, blo);
                if(ov > 0.0) stalk[b] += ov;
            }
        }
        for(int b = 0; b < ni; b++)
            std::cout << "[FBD_BAL] i=" << b << " [" << a.intervalStart[b] << "," << ((b + 1 < ni) ? a.intervalStart[b + 1] : xtop) << ") lam=" << lambdaAt(b)
                      << " nDiv=" << divc[b] << " stalkLen=" << stalk[b] << " nFossilY=" << fy[b] << "\n";
    }
    return fbdProb;
}

double FBDTreeModel::fossilTermLn(int i, int spineIdx, int ky, int kz, bool local){
    return fossilTermLn(unresolvedFossils, cachedGammaLn[i], i, spineIdx, ky, kz, local);
}

double FBDTreeModel::fossilTermLn(ParameterUnresolvedFossils* uf, double gLn, int i, int spineIdx, int ky, int kz, bool local){
    if(uf->isSA(i))
        return a.lnPsiTypeVec[a.fossilType[i]][ky >= 0 ? ky : findIndex(uf->getFossilAge(i))] + gLn;
    if(i == spineIdx && uf->isUE(i))
        return 0.0;
    if(uf->isUE(i)){
        double z = uf->getAttachAge(i);
        return std::log(rho) + a.ln2LambdaVec[findIndex(z)] + (local ? lnDLocal(z) : lnD(z)) + gLn;
    }
    if(i == spineIdx){
        double ys = uf->getFossilAge(i);
        return a.lnPsiTypeVec[a.fossilType[i]][findIndex(ys)] + lnP0At(findIndex(ys), ys) - (local ? lnDLocal(ys) : lnD(ys));
    }
    return fossilPqLn(uf->getFossilAge(i), uf->getAttachAge(i), a.fossilType[i], ky, kz) + gLn;
}

double FBDTreeModel::sweepTermLn(int i){
    double t = fossilTermLn(i, -1);
    if(pathPartner == nullptr)
        return t;
    return pathSelfW * t + pathPartnerW * pathPartner->fossilTermLn(unresolvedFossils, cachedGammaLn[i], i, -1, -1, -1, false);
}

double FBDTreeModel::term3Sum(void){
    int nf = unresolvedFossils->getNumFossils();
    int spineIdx = unresolvedFossils->getSpineIdx();
    bool useOrigin = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);
    double x0 = useOrigin ? originAge->getValue() : 0.0;
    double s = 0.0;
    for(int i = 0; i < nf; i++){
        if(useOrigin && (unresolvedFossils->getFossilAge(i) > x0 || unresolvedFossils->getAttachAge(i) > x0))
            return -INFINITY;
        s += fossilTermLn(i, spineIdx);
    }
    return s;
}

void FBDTreeModel::rebuildZoneStalks(int z, const std::vector<int>& fos){
    zoneStalks[z].y.clear();
    zoneStalks[z].z.clear();
    for(int i : fos){
        if(unresolvedFossils->isSA(i))
            continue;
        zoneStalks[z].y.push_back(unresolvedFossils->getFossilAge(i));
        zoneStalks[z].z.push_back(unresolvedFossils->getAttachAge(i));
    }
    std::sort(zoneStalks[z].y.begin(), zoneStalks[z].y.end());
    std::sort(zoneStalks[z].z.begin(), zoneStalks[z].z.end());
}

void FBDTreeModel::zoneRecomputeGamma(const std::vector<int>& fos){
    if(gammaValM.size() != cachedGammaLn.size())
        gammaValM.assign(cachedGammaLn.size(), 0.0);
    for(int i : fos){
        double g = computeGamma(unresolvedFossils->getAttachAge(i), i);
        cachedGammaLn[i] = (g > 0.0) ? std::log(g) : -INFINITY;
        gammaValM[i] = g;
    }
}

double FBDTreeModel::zoneTermSum(const std::vector<int>& fos, double x0, bool useOrigin){
    double s = 0.0;
    for(int i : fos){
        if(useOrigin && (unresolvedFossils->getFossilAge(i) > x0 || unresolvedFossils->getAttachAge(i) > x0))
            return -INFINITY;
        s += fossilTermLn(i, -1);
    }
    return s;
}

double FBDTreeModel::fossilSweep(void){
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);
    prepareIntervals();
    updateGammaCache();
    if(pathPartner != nullptr)
        pathPartner->prepareIntervals();
    double r = (unresolvedFossils->getSpineIdx() < 0) ? fossilSweepParallel() : fossilSweepSequential();
    RandomVariable::setActiveInstance(prevRng);
    return r;
}

double FBDTreeModel::fossilSweepSequential(void){
    int nf = unresolvedFossils->getNumFossils();
    std::vector<int> order(nf);
    for(int i = 0; i < nf; i++) order[i] = i;
    for(int i = nf - 1; i > 0; i--)
        std::swap(order[i], order[(int)(rng.uniformRv() * (i + 1))]);

    double t3 = (pathTarget != nullptr) ? pathTarget() : term3Sum();
    for(int idx : order){
        double ratio = unresolvedFossils->proposeOneFossil(idx);
        if(ratio == -INFINITY){
            unresolvedFossils->updateForRejection();
            continue;
        }
        updateGammaCache();
        double t3new = (pathTarget != nullptr) ? pathTarget() : term3Sum();
        if(std::log(rng.uniformRv()) < heat * (t3new - t3) + ratio){
            t3 = t3new;
            unresolvedFossils->updateForAcceptance();
        }else{
            unresolvedFossils->updateForRejection();
            updateGammaCache();
        }
    }
    return std::numeric_limits<double>::infinity();
}

double FBDTreeModel::fossilSweepParallel(void){
    int nf = unresolvedFossils->getNumFossils();
    int nsub = unresolvedFossils->subMoveCount();
    bool useOrigin = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);
    double x0 = useOrigin ? originAge->getValue() : 0.0;

    std::vector<std::vector<int> > zoneFos(numZones);
    for(int i = 0; i < nf; i++)
        zoneFos[unresolvedFossils->getAttachmentZone(i)].push_back(i);

    std::vector<unsigned int> seed(numZones);
    for(int z = 0; z < numZones; z++)
        seed[z] = (unsigned int)(rng.uniformRv() * 4294967295.0) + 1u;

    std::vector<std::vector<long> > zAtt(numZones, std::vector<long>(nsub, 0));
    std::vector<std::vector<long> > zAcc(numZones, std::vector<long>(nsub, 0));
    std::vector<long> zNAcc(numZones, 0), zNRej(numZones, 0);
    std::vector<double> zoneTerm(nf, 0.0);
    std::vector<double> gammaVal(nf, 0.0);
    int logMax = 2 * (nf + parameterTree->getTree()->getNumNodes()) + 8;
    std::vector<double> logHalf(logMax);
    for(int q = 0; q < logMax; q++)
        logHalf[q] = (q > 0) ? std::log(0.5 * (double)q) : -INFINITY;

    ThreadPool::current().parallelFor(OP_FBD, numZones, [&](int lo, int hi){
        for(int z = lo; z < hi; z++){
            std::vector<int>& fos = zoneFos[z];
            if(fos.empty())
                continue;
            RandomVariable zr(seed[z]);
            std::vector<int> ord = fos;
            for(int a = (int)ord.size() - 1; a > 0; a--)
                std::swap(ord[a], ord[(int)(zr.uniformRv() * (a + 1))]);
            rebuildZoneStalks(z, fos);
            zoneRecomputeGamma(fos);
            std::vector<std::pair<double,int> > zIdx;
            zIdx.reserve(fos.size());
            double zt = 0.0;
            for(int i : fos){
                gammaVal[i] = computeGamma(unresolvedFossils->getAttachAge(i), i);
                zoneTerm[i] = sweepTermLn(i);
                zt += zoneTerm[i];
                zIdx.push_back(std::make_pair(unresolvedFossils->getAttachAge(i), i));
            }
            std::sort(zIdx.begin(), zIdx.end());

            StalkBucket& sb = zoneStalks[z];
            std::vector<int> touched, keptIdx;
            std::vector<double> keptGamma, keptTerm, keptVal;
            bool halfFixSweep = (UserSettings::userSettings().getModel() == Model::UFBD);

            for(int idx : ord){
                double yO = unresolvedFossils->getFossilAge(idx);
                double zO = unresolvedFossils->getAttachAge(idx);
                bool saO = unresolvedFossils->isSA(idx);
                int sub;
                double ratio = unresolvedFossils->proposeFossilParallel(idx, zr, sub);
                zAtt[z][sub]++;
                if(ratio == -INFINITY){
                    unresolvedFossils->rejectFossil(idx);
                    zNRej[z]++;
                    continue;
                }
                double yN = unresolvedFossils->getFossilAge(idx);
                double zN = unresolvedFossils->getAttachAge(idx);
                bool saN = unresolvedFossils->isSA(idx);
                if(useOrigin && (yN > x0 || zN > x0)){
                    unresolvedFossils->rejectFossil(idx);
                    zNRej[z]++;
                    continue;
                }

                bool yKept = (saO == false && saN == false && yN == yO);
                if(saO == false){
                    if(yKept) eraseStalkZ(sb.z, zO);
                    else      eraseStalk(sb.y, sb.z, yO, zO);
                }
                if(saN == false){
                    if(yKept) insertStalkZ(sb.z, zN);
                    else      insertStalk(sb.y, sb.z, yN, zN);
                }
                eraseZIdx(zIdx, zO, idx);
                insertZIdx(zIdx, zN, idx);

                double b1lo, b1hi, b2lo, b2hi;
                stalkChangeBands(saO, yO, zO, saN, yN, zN, b1lo, b1hi, b2lo, b2hi);
                touched.clear();
                collectInBand(zIdx, b1lo, b1hi, touched);
                collectInBand(zIdx, b2lo, b2hi, touched);
                touched.push_back(idx);

                keptGamma.clear();
                keptTerm.clear();
                keptVal.clear();
                keptIdx.clear();
                double dz = 0.0;
                for(int a : touched){
                    double q = unresolvedFossils->getAttachAge(a);
                    double d = 0.0;
                    if(a != idx){
                        double inN = (saN == false && yN < q && q < zN) ? 1.0 : 0.0;
                        double inO = (saO == false && yO < q && q < zO) ? 1.0 : 0.0;
                        if(inN == inO)
                            continue;
                        d = ((halfFixSweep && unresolvedFossils->isSA(a) == false) ? 0.5 : 1.0) * (inN - inO);
                    }
                    keptGamma.push_back(cachedGammaLn[a]);
                    keptTerm.push_back(zoneTerm[a]);
                    keptVal.push_back(gammaVal[a]);
                    keptIdx.push_back(a);
                    if(a == idx){
                        gammaVal[a] = computeGamma(q, a);
                        cachedGammaLn[a] = (gammaVal[a] > 0.0) ? std::log(gammaVal[a]) : -INFINITY;
                        double t = sweepTermLn(a);
                        dz += t - zoneTerm[a];
                        zoneTerm[a] = t;
                    }else{
                        gammaVal[a] += d;
                        double ln = logHalf[(int)(2.0 * gammaVal[a] + 0.5)];
                        dz += ln - cachedGammaLn[a];
                        zoneTerm[a] += ln - cachedGammaLn[a];
                        cachedGammaLn[a] = ln;
                    }
                }

                if(std::log(zr.uniformRv()) < heat * dz + ratio){
                    zt += dz;
                    unresolvedFossils->acceptFossil(idx);
                    zAcc[z][sub]++;
                    zNAcc[z]++;
                }else{
                    unresolvedFossils->rejectFossil(idx);
                    zNRej[z]++;
                    if(saN == false){
                        if(yKept) eraseStalkZ(sb.z, zN);
                        else      eraseStalk(sb.y, sb.z, yN, zN);
                    }
                    if(saO == false){
                        if(yKept) insertStalkZ(sb.z, zO);
                        else      insertStalk(sb.y, sb.z, yO, zO);
                    }
                    eraseZIdx(zIdx, zN, idx);
                    insertZIdx(zIdx, zO, idx);
                    for(size_t q = keptIdx.size(); q-- > 0; ){
                        cachedGammaLn[keptIdx[q]] = keptGamma[q];
                        zoneTerm[keptIdx[q]] = keptTerm[q];
                        gammaVal[keptIdx[q]] = keptVal[q];
                    }
                }
            }
        }
    });

    std::vector<long> att(nsub, 0), acc(nsub, 0);
    long nAcc = 0, nRej = 0;
    for(int z = 0; z < numZones; z++){
        for(int s = 0; s < nsub; s++){ att[s] += zAtt[z][s]; acc[s] += zAcc[z][s]; }
        nAcc += zNAcc[z]; nRej += zNRej[z];
    }
    unresolvedFossils->mergeSubStats(att.data(), acc.data(), nAcc, nRej);
    stalkIndexFresh = true;
    updateGammaCache();
    return std::numeric_limits<double>::infinity();
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
            lnP += lnDDiff(n->getAncestor()->getTime(), n->getTime());
        if(n->getIsTip()){
            if(n->getIsFossil() == false)
                lnP += std::log(rho);
            else{
                int ft = 0;
                if(a.numPsiTypes > 1){ std::map<std::string,int>::iterator it = a.fossilTypeByName.find(n->getName()); if(it != a.fossilTypeByName.end()) ft = it->second; }
                if(n->getIsSA())
                    lnP += std::log(psiOfTypeAt(ft, findIndex(n->getTime())));
                else
                    lnP += std::log(psiOfTypeAt(ft, findIndex(n->getTime()))) + lnP0At(findIndex(n->getTime()), n->getTime());
            }
        }
        else if(n != root){
            bool fakeSplit = false;
            for(Node* c : n->getNeighbors())
                if(c != n->getAncestor() && c->getIsTip() && c->getIsFossil() && c->getIsSA()){ fakeSplit = true; break; }
            if(fakeSplit == false)
                lnP += std::log(2.0 * a.lamVec[findIndex(n->getTime())]);
        }
    }
    return lnP;
}

double FBDTreeModel::lnD(double t){
    return (t <= 0.0) ? 0.0 : a.lnDPrev[findIndex(t)] + std::log(4.0) - calculateLnQtAt(findIndex(t), t);
}

double FBDTreeModel::fossilPqLn(double y, double z, int type, int ky, int kz){
    if(ky < 0) ky = findIndex(y);
    if(kz < 0) kz = findIndex(z);
    return a.lnPsiTypeVec[type][ky] + a.ln2LambdaVec[kz] + lnP0At(ky, y)
           + calculateLnQtAt(ky, y) - calculateLnQtAt(kz, z) + lnDBetween(ky, kz);
}

double FBDTreeModel::lnDLocal(double t){
    return (t <= 0.0) ? 0.0 : std::log(4.0) - calculateLnQtAt(findIndex(t), t);
}

// mu - lambda - c1 = -(psi + c1 + beta) to keep lnD + survival decrement exact at large mu
double FBDTreeModel::lnDPlusDec(int i, double t){
    double linear;
    double lg = lnSurvDecParts(i, t, linear);
    if(linear == 0.0)
        return std::log(4.0) - calculateLnQtAt(i, t) + lg;
    return std::log(4.0) - 2.0 * lnDenomAt(i, t) + lg - (a.psiTotVec[i] + a.cpbVec[i]) * (t - a.intervalStart[i]);
}

double FBDTreeModel::lnQBracket(int i, double t){
    return 2.0 * lnDenomAt(i, t);
}

double FBDTreeModel::lnSurvDecrement(int i, double t){
    double linear;
    double lg = lnSurvDecParts(i, t, linear);
    return linear + lg;
}

double FBDTreeModel::lnSurvDecParts(int i, double t, double& linear){
    double lam = a.lamVec[i], mu = a.muVec[i];
    double tau = t - a.intervalStart[i];
    double c1 = std::abs(lam - mu);
    double u = (c1 > 0.0) ? -std::expm1(-c1 * tau) / c1 : tau;
    double s = std::exp(a.lnSPrevHat[i]);
    if(lam >= mu){
        linear = 0.0;
        return std::log(std::exp(-c1 * tau) + s * lam * u);
    }
    linear = c1 * tau;
    return std::log1p(s * lam * u);
}

double FBDTreeModel::incDecLn(int i, double t){
    double lam = a.lamVec[i], mu = a.muVec[i];
    if(lam >= mu)
        return std::log(4.0) - calculateLnQtAt(i, t) + lnSurvDecrement(i, t);
    double psi = a.psiTotVec[i], beta = lam - mu - psi;
    double delta = (beta >= 0.0 ? a.c1Vec[i] + beta : 4.0 * lam * psi / (a.c1Vec[i] - beta)) + psi;
    double tau = t - a.intervalStart[i], ch = mu - lam;
    double u = -std::expm1(-ch * tau) / ch;
    double s = std::exp(a.lnSPrevHat[i]);
    return -delta * tau + std::log(4.0) - lnQBracket(i, t) + std::log1p(s * lam * u);
}

double FBDTreeModel::lnDDiff(double zhi, double ylo){
    int kz = findIndex(zhi), ky = findIndex(ylo);
    return calculateLnQtAt(ky, ylo) - calculateLnQtAt(kz, zhi) + lnDBetween(ky, kz);
}

// summed per bin: lnDprev[kz] - lnDPrev[ky] loses every digit when c1*dt gets very large
double FBDTreeModel::lnDBetween(int ky, int kz){
    double v = 0.0;
    for(int j = ky; j < kz; j++)
        v += a.lnDInc[j];
    return v;
}

double FBDTreeModel::lnSurvivalWithin(int i, double t){
    double lam = a.lamVec[i], mu = a.muVec[i];
    double tau = t - a.intervalStart[i];
    double c1 = std::abs(lam - mu);
    double u = (c1 > 0.0) ? -std::expm1(-c1 * tau) / c1 : tau;
    double lnS = a.lnSPrevHat[i];
    double s = std::exp(lnS);
    if(lam >= mu)
        return lnS - std::log(std::exp(-c1 * tau) + s * lam * u);
    return lnS - c1 * tau - std::log1p(s * lam * u);
}

double FBDTreeModel::calculateLnSurvival(double t){
    return lnSurvivalWithin(findIndex(t), t);
}

double FBDTreeModel::calculateLnAnySample(double t){
    int k = findIndex(t);
    return std::log(sampledAt(k, t));
}

double FBDTreeModel::calculateLnConditioning(double t){
    switch(UserSettings::userSettings().getConditioningEvent()){
        case ConditioningEvent::SURVIVAL:
            return calculateLnSurvival(t);
        case ConditioningEvent::ANYSAMPLE:
            return calculateLnAnySample(t);
        case ConditioningEvent::EXTINCT: {
            int k = findIndex(t);
            double pExtinct = calculateSurvivalAt(k, t) - std::exp(lnSurvivalWithin(k, t));
            return (pExtinct > 0.0) ? std::log(pExtinct) : -INFINITY;
        }
    }
    return -INFINITY;
}

void FBDTreeModel::prepareIntervals(void){
    size_t n = a.intervalStart.size();
    a.c1Vec.assign(n, 0.0);
    a.cpbVec.resize(n);
    a.cmbVec.resize(n);
    a.dPrev.resize(n);
    a.lnDPrev.assign(n, 0.0);
    a.lnDInc.resize(n);
    a.lnDenEnd.resize(n);
    a.lamVec.resize(n);
    a.muVec.resize(n);
    a.psiTotVec.resize(n);
    for(size_t i = 0; i < n; i++){
        double li = lambdaAt(i);
        double mi = muAt(i);
        double pi = psiTotalAt(i);
        a.lamVec[i] = li;
        a.muVec[i] = mi;
        a.psiTotVec[i] = pi;
        double beta = li - mi - pi;
        a.c1Vec[i] = std::abs(std::sqrt(std::pow(beta, 2) + 4*li*pi));
        if(beta >= 0.0){
            a.cpbVec[i] = a.c1Vec[i] + beta;
            a.cmbVec[i] = 4*li*pi / a.cpbVec[i];
        }else{
            a.cmbVec[i] = a.c1Vec[i] - beta;
            a.cpbVec[i] = 4*li*pi / a.cmbVec[i];
        }
        if(i == 0){
            a.dPrev[0] = rho;
            a.lnDPrev[0] = 0.0;
        }else{
            int j = (int)i - 1;
            double x = a.c1Vec[j] * (a.intervalStart[i] - a.intervalStart[j]);
            double e = std::exp(-x), d = a.dPrev[j];
            double ome = (x < 0.7) ? -std::expm1(-x) : 1.0 - e;
            double den = a.cmbVec[j] + e * a.cpbVec[j] + 2.0 * d * a.lamVec[j] * ome;
            a.dPrev[i] = (2.0 * a.psiTotVec[j] * ome + d * (a.cpbVec[j] + e * a.cmbVec[j])) / den;
            a.lnDenEnd[j] = std::log(den / a.c1Vec[j]);
            double lnQ = x + 2.0 * a.lnDenEnd[j];
            a.lnDPrev[i] = a.lnDPrev[i-1] + std::log(4.0) - lnQ;
            a.lnDInc[j] = std::log(4.0) - lnQ;
        }
    }
    a.ln2LambdaVec.assign(n, 0.0);
    for(size_t i = 0; i < n; i++)
        a.ln2LambdaVec[i] = std::log(2 * a.lamVec[i]);
    a.lnPsiTypeVec.assign(a.numPsiTypes, std::vector<double>(n, 0.0));
    for(int t = 0; t < a.numPsiTypes; t++)
        for(size_t i = 0; i < n; i++)
            a.lnPsiTypeVec[t][i] = std::log(psiOfTypeAt(t, (int)i));

    a.lnSPrevHat.assign(n, 0.0);
    a.lnRPrev.assign(n, 0.0);
    a.lnDecInc.resize(n);
    a.lnDecPair.resize(n);
    a.lnSPrevHat[0] = std::log(rho);
    a.lnRPrev[0] = -std::log(rho);
    for(size_t i = 1; i < n; i++){
        int j = (int)i - 1;
        double linear;
        double lg = lnSurvDecParts(j, a.intervalStart[i], linear);
        double dec = linear + lg;
        double incdec = incDecLn(j, a.intervalStart[i]);
        a.lnDecInc[j] = dec;
        a.lnDecPair[j] = (linear == 0.0) ? a.lnDInc[j] + dec
                       : std::log(4.0) - 2.0 * a.lnDenEnd[j] + lg - (a.psiTotVec[j] + a.cpbVec[j]) * (a.intervalStart[i] - a.intervalStart[j]);
        a.lnRPrev[i] = a.lnRPrev[i-1] + incdec;
        a.lnSPrevHat[i] = a.lnSPrevHat[i-1] - dec;
    }
}

int FBDTreeModel::findIndex(double t){
    int i = (int)(std::upper_bound(a.intervalStart.begin(), a.intervalStart.end(), t) - a.intervalStart.begin()) - 1;
    return (i > 0) ? i : 0;
}

double FBDTreeModel::calculateLnQtAt(int i, double t){
    return a.c1Vec[i] * (t - a.intervalStart[i]) + 2.0 * lnDenomAt(i, t);
}

// d = 1 - E with nonnegative terms only; c2 cancels and rounds E to 1 when mu is very large
double FBDTreeModel::lnDenomAt(int i, double t){
    double e = std::exp(-a.c1Vec[i] * (t - a.intervalStart[i]));
    return std::log((a.cmbVec[i] + e * a.cpbVec[i] + 2.0 * a.dPrev[i] * a.lamVec[i] * (1.0 - e)) / a.c1Vec[i]);
}

double FBDTreeModel::sampledAt(int i, double t){
    double x = a.c1Vec[i] * (t - a.intervalStart[i]);
    double e = std::exp(-x), d = a.dPrev[i];
    double ome = (x < 0.7) ? -std::expm1(-x) : 1.0 - e;
    return (2.0 * a.psiTotVec[i] * ome + d * (a.cpbVec[i] + e * a.cmbVec[i]))
           / (a.cmbVec[i] + e * a.cpbVec[i] + 2.0 * d * a.lamVec[i] * ome);
}

double FBDTreeModel::calculateP0At(int i, double t){
    return 1.0 - sampledAt(i, t);
}

double FBDTreeModel::lnP0At(int i, double t){
    return std::log(1.0 - sampledAt(i, t));
}

double FBDTreeModel::calculateP0(double t){
    return calculateP0At(findIndex(t), t);
}

double FBDTreeModel::calculateSurvivalAt(int i, double t){
    return sampledAt(i, t);
}

double FBDTreeModel::calculateP0HatAt(int i, double t){
    return -std::expm1(lnSurvivalWithin(i, t));
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

void FBDTreeModel::buildZoneIndex(void){
    Tree* tree = parameterTree->getTree();
    if(eulerBuilt == false)
        buildEulerIndex();
    int nn = tree->getNumNodes();
    int nf = unresolvedFossils->getNumFossils();
    bool hasTrunk = (originAge != nullptr);

    std::vector<Node*> zCrown;
    std::vector<char>  zIsCrown, zIsStem;
    std::vector<int>   zoneOfFossil(nf, -1);
    zoneCladeName.clear();
    zoneType.clear();
    for(int i = 0; i < nf; i++){
        Node* c = unresolvedFossils->getCrownNode(i);
        char cr = unresolvedFossils->getIsCrown(i) ? 1 : 0;
        char st = unresolvedFossils->getIsStem(i) ? 1 : 0;
        for(int z = 0; z < (int)zCrown.size(); z++)
            if(zCrown[z] == c && zIsCrown[z] == cr && zIsStem[z] == st){
                zoneOfFossil[i] = z;
                break;
            }
        if(zoneOfFossil[i] < 0){
            zoneOfFossil[i] = (int)zCrown.size();
            zCrown.push_back(c);
            zIsCrown.push_back(cr);
            zIsStem.push_back(st);
            zoneCladeName.push_back(unresolvedFossils->getCladeName(i));
            zoneType.push_back(cr ? 'C' : (st ? 'S' : 'T'));
        }
    }
    numZones = (int)zCrown.size();

    std::vector<char> member((size_t)numZones * nn, 0);
    std::vector<char> trunkIn(numZones, 0);
    std::vector<int>  zSize(numZones, 0);
    for(int z = 0; z < numZones; z++){
        Node* c = zCrown[z];
        for(Node* n : tree->getDownPassSequence()){
            if(n == tree->getRoot())
                continue;
            bool in;
            if(zIsStem[z])
                in = (n == c);
            else
                in = inSub(n->getAncestor(), c) || (zIsCrown[z] == 0 && n == c);
            if(in){
                member[(size_t)z * nn + n->getOffset()] = 1;
                zSize[z]++;
            }
        }
        if(hasTrunk && c == tree->getCrown() && zIsCrown[z] == 0){
            trunkIn[z] = 1;
            zSize[z]++;
        }
    }

    std::vector<char> zoneSubsetMat((size_t)numZones * numZones, 0);
    for(int a = 0; a < numZones; a++)
        for(int b = 0; b < numZones; b++){
            bool sub = (trunkIn[a] == 0 || trunkIn[b] == 1);
            for(int o = 0; o < nn && sub; o++)
                if(member[(size_t)a * nn + o] && member[(size_t)b * nn + o] == 0)
                    sub = false;
            zoneSubsetMat[(size_t)a * numZones + b] = sub ? 1 : 0;
        }

    minZoneOfNode.assign(nn, -1);
    for(Node* n : tree->getDownPassSequence()){
        if(n == tree->getRoot())
            continue;
        int best = -1;
        for(int z = 0; z < numZones; z++){
            if(member[(size_t)z * nn + n->getOffset()] == 0)
                continue;
            if(best < 0 || zoneSubsetMat[(size_t)z * numZones + best])
                best = z;
        }
        minZoneOfNode[n->getOffset()] = best;
    }
    trunkMinZone = -1;
    if(hasTrunk)
        for(int z = 0; z < numZones; z++){
            if(trunkIn[z] == 0)
                continue;
            if(trunkMinZone < 0 || zoneSubsetMat[(size_t)z * numZones + trunkMinZone])
                trunkMinZone = z;
        }

    std::vector<char> occupied(numZones, 0);
    for(int o = 0; o < nn; o++)
        if(minZoneOfNode[o] >= 0)
            occupied[minZoneOfNode[o]] = 1;
    if(trunkMinZone >= 0)
        occupied[trunkMinZone] = 1;

    std::vector<std::vector<int> > domain(nf);
    for(int i = 0; i < nf; i++)
        for(int z = 0; z < numZones; z++)
            if(occupied[z] && zoneSubsetMat[(size_t)z * numZones + zoneOfFossil[i]])
                domain[i].push_back(z);

    azGibbsIdx.clear();
    for(int i = 0; i < nf; i++)
        if(domain[i].size() > 1)
            azGibbsIdx.push_back(i);

    zoneEdges.assign(numZones, ZoneEdges());
    zoneStalks.assign(numZones, StalkBucket());
    unresolvedFossils->setAttachmentZoneDomain(domain);
    zoneIndexBuilt = true;
}

void FBDTreeModel::rebuildStalkIndex(void){
    for(size_t k = 0; k < zoneStalks.size(); k++){
        zoneStalks[k].y.clear();
        zoneStalks[k].z.clear();
    }
    int nf = unresolvedFossils->getNumFossils();
    for(int i = 0; i < nf; i++){
        if(unresolvedFossils->isSA(i))
            continue;
        StalkBucket& b = zoneStalks[unresolvedFossils->getAttachmentZone(i)];
        b.y.push_back(unresolvedFossils->getFossilAge(i));
        b.z.push_back(unresolvedFossils->getAttachAge(i));
    }
    for(size_t k = 0; k < zoneStalks.size(); k++){
        std::sort(zoneStalks[k].y.begin(), zoneStalks[k].y.end());
        std::sort(zoneStalks[k].z.begin(), zoneStalks[k].z.end());
    }
}

void FBDTreeModel::moveStalk(int i, int from, int to){
    if(from == to || unresolvedFossils->isSA(i))
        return;
    double yi = unresolvedFossils->getFossilAge(i);
    double zi = unresolvedFossils->getAttachAge(i);
    eraseStalk(zoneStalks[from].y, zoneStalks[from].z, yi, zi);
    insertStalk(zoneStalks[to].y,  zoneStalks[to].z,  yi, zi);
}

double FBDTreeModel::validZoneSet(int i, int a, std::vector<std::pair<double,double> >& iv){
    Tree* tree = parameterTree->getTree();
    double lo = unresolvedFossils->getMinAttachAge(i);
    double hi = unresolvedFossils->getMaxAttachAge(i);
    std::vector<std::pair<double,double> > raw;
    for(Node* n : tree->getDownPassSequence()){
        if(n == tree->getRoot() || minZoneOfNode[n->getOffset()] != a)
            continue;
        raw.push_back(std::make_pair(n->getTime(), n->getAncestor()->getTime()));
    }
    if(a == trunkMinZone && tree->getNumBackbone() > 0)
        raw.push_back(std::make_pair(tree->getCrown()->getTime(), originAge->getValue()));
    int nf = unresolvedFossils->getNumFossils();
    for(int j = 0; j < nf; j++){
        if(j == i || unresolvedFossils->isSA(j) || unresolvedFossils->getAttachmentZone(j) != a)
            continue;
        raw.push_back(std::make_pair(unresolvedFossils->getFossilAge(j), unresolvedFossils->getAttachAge(j)));
    }
    std::vector<std::pair<double,double> > keep;
    for(std::pair<double,double>& p : raw){
        double a0 = (p.first  < lo) ? lo : p.first;
        double b0 = (p.second > hi) ? hi : p.second;
        if(b0 > a0)
            keep.push_back(std::make_pair(a0, b0));
    }
    std::sort(keep.begin(), keep.end());
    iv.clear();
    for(std::pair<double,double>& p : keep){
        if(iv.empty() == false && p.first <= iv.back().second){
            if(p.second > iv.back().second)
                iv.back().second = p.second;
        }else{
            iv.push_back(p);
        }
    }
    double measure = 0.0;
    for(std::pair<double,double>& p : iv)
        measure += p.second - p.first;
    return measure;
}

void FBDTreeModel::snapshotZoneEdges(ZoneEdgeSnap& s){
    Tree* tree = parameterTree->getTree();
    Node* root = tree->getRoot();
    s.e.assign(numZones, std::vector<std::pair<double,double> >());
    for(Node* nd : tree->getDownPassSequence()){
        if(nd == root)
            continue;
        int mz = minZoneOfNode[nd->getOffset()];
        if(mz < 0)
            continue;
        s.e[mz].push_back(std::make_pair(nd->getTime(), nd->getAncestor()->getTime()));
    }
    s.crownTime = tree->getCrown()->getTime();
    s.x0 = (originAge != nullptr) ? originAge->getValue() : 0.0;
    s.numBackbone = tree->getNumBackbone();
}

double FBDTreeModel::gammaWithZone(int i, int a, const ZoneEdgeSnap& s){
    double z = unresolvedFossils->getAttachAge(i);
    double count = 0.0;
    for(const std::pair<double,double>& e : s.e[a])
        if(e.first < z && z < e.second)
            count += 1.0;
    if(a == trunkMinZone){
        if(s.numBackbone == 0){
            if(z >= s.x0)
                count += 1.0;
        }else if(s.crownTime < z && z < s.x0){
            count += 1.0;
        }
    }
    bool halfFix = (UserSettings::userSettings().getModel() == Model::UFBD);
    bool focalIsTip = (unresolvedFossils->isSA(i) == false);
    double w = (halfFix && focalIsTip) ? 0.5 : 1.0;
    const StalkBucket& b = zoneStalks[a];
    count += w * (double)countStraddling(b.y, b.z, z);
    int sp = unresolvedFossils->getSpineIdx();
    if(halfFix && focalIsTip && sp >= 0 && sp != i && unresolvedFossils->isSA(sp) == false
       && unresolvedFossils->getAttachmentZone(sp) == a){
        double ys = unresolvedFossils->getFossilAge(sp), zs = unresolvedFossils->getAttachAge(sp);
        if(ys < z && z < zs)
            count += 0.5;
    }
    return count;
}

double FBDTreeModel::azBlockPass(const ZoneEdgeSnap& s, std::vector<int>& lab, bool draw){
    double lq = 0.0;
    std::vector<double> w;
    for(size_t k = 0; k < zbIdx.size(); k++){
        int i = zbIdx[k];
        std::vector<int>& dom = unresolvedFossils->getAttachmentZoneDomain(i);
        w.assign(dom.size(), 0.0);
        double sum = 0.0;
        for(size_t d = 0; d < dom.size(); d++){
            double g = gammaWithZone(i, dom[d], s);
            w[d] = (g > 0.0) ? g : 0.0;
            sum += w[d];
        }
        if(sum <= 0.0)
            return -INFINITY;
        if(draw){
            double u = rng.uniformRv() * sum;
            size_t pick = dom.size() - 1;
            for(size_t d = 0; d < dom.size(); d++){
                u -= w[d];
                if(u < 0.0){ pick = d; break; }
            }
            lab[k] = dom[pick];
            lq += std::log(w[pick] / sum);
        }else{
            double wl = 0.0;
            for(size_t d = 0; d < dom.size(); d++)
                if(dom[d] == lab[k])
                    wl = w[d];
            if(wl <= 0.0)
                return -INFINITY;
            lq += std::log(wl / sum);
        }
    }
    return lq;
}

double FBDTreeModel::relabelAcrossNode(Node* n, double oldAge){
    zbIdx.clear();
    if(isResolved || unresolvedFossils == nullptr || zoneIndexBuilt == false)
        return 0.0;
    double newAge = n->getTime();
    double lo = (oldAge < newAge) ? oldAge : newAge;
    double hi = (oldAge < newAge) ? newAge : oldAge;
    std::vector<char> aff(numZones, 0);
    int mz = minZoneOfNode[n->getOffset()];
    if(mz >= 0)
        aff[mz] = 1;
    for(Node* c : n->getDescendants()){
        int mc = minZoneOfNode[c->getOffset()];
        if(mc >= 0)
            aff[mc] = 1;
    }
    for(int i : azGibbsIdx){
        double t = unresolvedFossils->getAttachAge(i);
        if(t <= lo || t >= hi)
            continue;
        for(int a : unresolvedFossils->getAttachmentZoneDomain(i))
            if(aff[a]){ zbIdx.push_back(i); break; }
    }
    relAtt++;
    if(zbIdx.empty())
        return 0.0;
    relFire++;
    relSize += (long)zbIdx.size();

    zbOld.resize(zbIdx.size());
    for(size_t k = 0; k < zbIdx.size(); k++)
        zbOld[k] = unresolvedFossils->getAttachmentZone(zbIdx[k]);

    ZoneEdgeSnap snapNew, snapOld;
    snapshotZoneEdges(snapNew);
    n->setTime(oldAge);
    snapshotZoneEdges(snapOld);
    n->setTime(newAge);

    std::vector<int> lab(zbIdx.size(), -1);
    double lqf = azBlockPass(snapNew, lab, true);
    if(lqf == -INFINITY){
        zbIdx.clear();
        return -INFINITY;
    }
    unresolvedFossils->beginZoneBlock(zbIdx);
    if(inNodeAgeSweep){
        for(size_t k = 0; k < zbIdx.size(); k++){
            unresolvedFossils->setAttachmentZone(zbIdx[k], lab[k]);
            moveStalk(zbIdx[k], zbOld[k], lab[k]);
        }
    }else{
        for(size_t k = 0; k < zbIdx.size(); k++)
            unresolvedFossils->setAttachmentZone(zbIdx[k], lab[k]);
        rebuildStalkIndex();
    }

    std::vector<int> back = zbOld;
    double lqr = azBlockPass(snapOld, back, false);
    if(lqr == -INFINITY)
        return -INFINITY;
    return lqr - lqf;
}

void FBDTreeModel::commitZoneBlock(void){
    if(zbIdx.empty())
        return;
    for(size_t k = 0; k < zbIdx.size(); k++)
        if(unresolvedFossils->getAttachmentZone(zbIdx[k]) != zbOld[k]){
            relCross++;
            if(unresolvedFossils->isSA(zbIdx[k]))
                relCrossSA++;
        }
    unresolvedFossils->commitZoneBlock();
    zbIdx.clear();
}

void FBDTreeModel::restoreZoneBlock(void){
    if(zbIdx.empty())
        return;
    if(inNodeAgeSweep){
        for(size_t k = 0; k < zbIdx.size(); k++)
            moveStalk(zbIdx[k], unresolvedFossils->getAttachmentZone(zbIdx[k]), zbOld[k]);
        unresolvedFossils->restoreZoneBlock();
    }else{
        unresolvedFossils->restoreZoneBlock();
        rebuildStalkIndex();
    }
    zbIdx.clear();
}

double FBDTreeModel::doAttachmentZoneJump(int i){
    std::vector<int>& dom = unresolvedFossils->getAttachmentZoneDomain(i);
    int k = (int)dom.size();
    int cur = unresolvedFossils->getAttachmentZone(i);
    int pick = (int)(rng.uniformRv() * (k - 1));
    int cand = cur;
    for(int d = 0; d < k; d++){
        if(dom[d] == cur)
            continue;
        if(pick == 0){ cand = dom[d]; break; }
        pick--;
    }

    std::vector<std::pair<double,double> > ivNew, ivOld;
    double mNew = validZoneSet(i, cand, ivNew);
    double mOld = validZoneSet(i, cur, ivOld);
    unresolvedFossils->beginAttachmentZoneMove(i);
    if(mNew <= 0.0 || mOld <= 0.0)
        return -INFINITY;

    double u = rng.uniformRv() * mNew;
    double zNew = ivNew.back().second;
    for(std::pair<double,double>& p : ivNew){
        double w = p.second - p.first;
        if(u < w){ zNew = p.first + u; break; }
        u -= w;
    }

    unresolvedFossils->setAttachmentZone(i, cand);
    unresolvedFossils->setAttachAge(i, zNew);
    rebuildStalkIndex();
    return std::log(mNew) - std::log(mOld);
}

double FBDTreeModel::doAttachmentZoneGibbs(void){
    int i = azGibbsIdx[(int)(rng.uniformRv() * (int)azGibbsIdx.size())];
    if(unresolvedFossils->isSA(i) == false && i != unresolvedFossils->getSpineIdx()
       && rng.uniformRv() < 0.5)
        return doAttachmentZoneJump(i);
    std::vector<int>& dom = unresolvedFossils->getAttachmentZoneDomain(i);
    int k = (int)dom.size();
    int cur = unresolvedFossils->getAttachmentZone(i);
    int nf = unresolvedFossils->getNumFossils();

    std::vector<int> aff(1, i);                 // only i and the fossils its stalk can attachEdge change
    if(unresolvedFossils->isSA(i) == false){
        double yi = unresolvedFossils->getFossilAge(i);
        double zi = unresolvedFossils->getAttachAge(i);
        for(int j = 0; j < nf; j++){
            double tj = unresolvedFossils->getAttachAge(j);
            if(j != i && yi < tj && tj < zi)
                aff.push_back(j);
        }
    }

    std::vector<double> lnG(k, 0.0);
    int curIdx = 0;
    for(int d = 0; d < k; d++){
        if(dom[d] == cur)
            curIdx = d;
        unresolvedFossils->setAttachmentZone(i, dom[d]);
        rebuildStalkIndex();
        for(int j : aff){
            double g = computeGamma(unresolvedFossils->getAttachAge(j), j);
            lnG[d] += (g > 0.0) ? std::log(g) : -INFINITY;
        }
    }

    double mx = -INFINITY;
    for(int d = 0; d < k; d++)
        if(lnG[d] > mx) mx = lnG[d];
    std::vector<double> w(k);
    double sum = 0.0;
    for(int d = 0; d < k; d++){
        w[d] = std::exp(lnG[d] - mx);
        sum += w[d];
    }
    double u = rng.uniformRv() * sum;
    int pick = k - 1;
    for(int d = 0; d < k; d++){
        u -= w[d];
        if(u < 0.0){ pick = d; break; }
    }

    unresolvedFossils->beginAttachmentZoneMove(i);
    unresolvedFossils->setAttachmentZone(i, dom[pick]);
    rebuildStalkIndex();
    return lnG[curIdx] - lnG[pick];
}

double FBDTreeModel::computeGamma(double z, int i){
    int mz = unresolvedFossils->getAttachmentZone(i);
    double count = zoneBackboneEdges(mz, z);

    bool halfFix = (UserSettings::userSettings().getModel() == Model::UFBD);
    bool focalIsTip = (unresolvedFossils->isSA(i) == false);
    double w = (halfFix && focalIsTip) ? 0.5 : 1.0;

    const StalkBucket& b = zoneStalks[mz];
    count += w * (double)countStraddling(b.y, b.z, z);

    int sp = unresolvedFossils->getSpineIdx();
    if(halfFix && focalIsTip && sp >= 0 && sp != i && unresolvedFossils->isSA(sp) == false
       && unresolvedFossils->getAttachmentZone(sp) == mz){
        double ys = unresolvedFossils->getFossilAge(sp), zs = unresolvedFossils->getAttachAge(sp);
        if(ys < z && z < zs)
            count += 0.5;                              // spine acts as backbone: weight 1, not 1/2
    }
    return count;
}

namespace {
struct EdgeFlip {
    int    zone;
    double lo, hi;
    double sign;
    double gate;
    bool   gateIsUpper;
};
}

void FBDTreeModel::updateGammaCache(void){
    Tree* tree = parameterTree->getTree();
    int nf = unresolvedFossils->getNumFossils();
    std::vector<Node*>& dpseq = tree->getDownPassSequence();

    if(zoneIndexBuilt == false)
        buildZoneIndex();

    Node* root = tree->getRoot();
    int nEdge = 0;
    for(int z = 0; z < numZones; z++) nEdge += (int)zoneEdges[z].yng.size();
    int budget = 0;
    for(int x = nEdge; x > 1; x >>= 1) budget++;
    int edits = (cacheInit && nEdge > 0) ? 0 : -1;
    bool rootOrCrownMoved = false;
    std::vector<Node*> moved;
    if(edits == 0){
        for(Node* nd : dpseq)
            if(nd->getTime() != prevNodeAge[nd->getOffset()]){
                edits += (nd != root ? 1 : 0) + (int)nd->getDescendants().size();
                if(edits > budget){ edits = -1; break; }
                moved.push_back(nd);
                if(nd == root || nd == tree->getCrown())
                    rootOrCrownMoved = true;
            }
    }
    bool nodeDelta = (edits >= 0) && cacheInit && gammaValM.size() == (size_t)nf && rootOrCrownMoved == false;
    std::vector<EdgeFlip> flips;

    if(edits >= 0){
        for(Node* nd : moved){
            double oldA = prevNodeAge[nd->getOffset()];
            double newA = nd->getTime();
            if(nodeDelta){
                double blo = std::min(oldA, newA), bhi = std::max(oldA, newA);
                int mzf = minZoneOfNode[nd->getOffset()];
                if(mzf >= 0){
                    EdgeFlip f;
                    f.zone = mzf; f.lo = blo; f.hi = bhi;
                    f.sign = (newA < oldA) ? 1.0 : -1.0;
                    f.gate = nd->getAncestor()->getTime(); f.gateIsUpper = true;
                    flips.push_back(f);
                }
                for(Node* c : nd->getDescendants()){
                    int mcf = minZoneOfNode[c->getOffset()];
                    if(mcf < 0)
                        continue;
                    EdgeFlip f;
                    f.zone = mcf; f.lo = blo; f.hi = bhi;
                    f.sign = (newA > oldA) ? 1.0 : -1.0;
                    f.gate = prevNodeAge[c->getOffset()]; f.gateIsUpper = false;
                    flips.push_back(f);
                }
            }
            int mz = (nd != root) ? minZoneOfNode[nd->getOffset()] : -1;
            if(mz >= 0){
                std::vector<double>& v = zoneEdges[mz].yng;
                v.erase(std::lower_bound(v.begin(), v.end(), oldA));
                v.insert(std::lower_bound(v.begin(), v.end(), newA), newA);
            }
            for(Node* c : nd->getDescendants()){
                int mc = minZoneOfNode[c->getOffset()];
                if(mc < 0)
                    continue;
                std::vector<double>& v = zoneEdges[mc].old;
                v.erase(std::lower_bound(v.begin(), v.end(), oldA));
                v.insert(std::lower_bound(v.begin(), v.end(), newA), newA);
            }
        }
    }else{
        for(int z = 0; z < numZones; z++){
            zoneEdges[z].yng.clear();
            zoneEdges[z].old.clear();
        }
        for(Node* nd : dpseq){
            if(nd == root)
                continue;
            int mz = minZoneOfNode[nd->getOffset()];
            if(mz < 0)
                continue;
            zoneEdges[mz].yng.push_back(nd->getTime());
            zoneEdges[mz].old.push_back(nd->getAncestor()->getTime());
        }
        for(int z = 0; z < numZones; z++){
            std::sort(zoneEdges[z].yng.begin(), zoneEdges[z].yng.end());
            std::sort(zoneEdges[z].old.begin(), zoneEdges[z].old.end());
        }
    }

    if(cacheInit == false){
        cachedGammaLn.assign(nf, 0.0);
        gammaValM.assign(nf, 0.0);
        gammaStale.assign(nf, 1);
        prevY.assign(nf, -1.0);
        prevZ.assign(nf, -1.0);
        prevSA.assign(nf, -1);
        prevAttachmentZone.assign(nf, -1);
        prevNodeAge.assign(tree->getNumNodes(), -1.0);
        prevX0 = -1.0;
        cacheInit = true;
        if(zoneInit == false) [[unlikely]] {
            for(int i = 0; i < nf; i++){
                std::vector<int>& dom = unresolvedFossils->getAttachmentZoneDomain(i);
                double zi = unresolvedFossils->getAttachAge(i);
                int best = dom[0];
                double bestCount = -1.0;
                for(int d : dom){
                    double c = zoneBackboneEdges(d, zi);
                    if(c > bestCount){
                        bestCount = c;
                        best = d;
                    }
                }
                unresolvedFossils->initAttachmentZone(i, best);
            }
            zoneInit = true;
        }
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

    bool anyFossilMoved = false;
    std::vector<std::pair<double,int> > tSorted;
    for(int i = 0; i < nf; i++){
        double yi = unresolvedFossils->getFossilAge(i);
        double zi = unresolvedFossils->getAttachAge(i);
        int sai = unresolvedFossils->isSA(i) ? 1 : 0;
        int lzi = unresolvedFossils->getAttachmentZone(i);
        if(yi == prevY[i] && zi == prevZ[i] && sai == prevSA[i] && lzi == prevAttachmentZone[i])
            continue;
        anyFossilMoved = true;
        gammaStale[i] = 1;
        bool wasTerm = (prevSA[i] == 0);
        bool isTerm = (sai == 0);
        if((wasTerm || isTerm) && prevY[i] >= 0.0){
            bool pureZ = (yi == prevY[i] && wasTerm && isTerm && lzi == prevAttachmentZone[i]);
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
            if(tSorted.empty()){
                tSorted.reserve(nf);
                for(int j = 0; j < nf; j++)
                    tSorted.push_back(std::make_pair(unresolvedFossils->isSA(j)
                                                     ? unresolvedFossils->getFossilAge(j)
                                                     : unresolvedFossils->getAttachAge(j), j));
                std::sort(tSorted.begin(), tSorted.end());
            }
            const int lowI = std::numeric_limits<int>::min(), highI = std::numeric_limits<int>::max();
            std::vector<std::pair<double,int> >::const_iterator b = pureZ
                ? std::lower_bound(tSorted.begin(), tSorted.end(), std::make_pair(lo, lowI))
                : std::upper_bound(tSorted.begin(), tSorted.end(), std::make_pair(lo, highI));
            std::vector<std::pair<double,int> >::const_iterator e =
                std::lower_bound(tSorted.begin(), tSorted.end(), std::make_pair(hi, lowI));
            for(; b != e; ++b)
                gammaStale[b->second] = 1;
        }
        prevY[i] = yi;
        prevZ[i] = zi;
        prevSA[i] = sai;
        prevAttachmentZone[i] = lzi;
    }

    if(nodeDelta && flips.empty() == false){
        if((int)gammaLogTab.size() < 2 * (nf + tree->getNumNodes()) + 8){
            gammaLogTab.resize(2 * (nf + tree->getNumNodes()) + 8);
            for(size_t q = 0; q < gammaLogTab.size(); q++)
                gammaLogTab[q] = (q > 0) ? std::log(0.5 * (double)q) : -INFINITY;
        }
        flipsByZone.assign(numZones, std::vector<int>());
        for(size_t k = 0; k < flips.size(); k++)
            flipsByZone[flips[k].zone].push_back((int)k);
        for(int i = 0; i < nf; i++){
            if(gammaStale[i])
                continue;
            const std::vector<int>& fz = flipsByZone[unresolvedFossils->getAttachmentZone(i)];
            if(fz.empty())
                continue;
            double q = unresolvedFossils->getAttachAge(i);
            for(int k : fz){
                const EdgeFlip& f = flips[k];
                if(q <= f.lo || q >= f.hi)
                    continue;
                if(f.gateIsUpper ? (q >= f.gate) : (q <= f.gate))
                    continue;
                gammaValM[i] += f.sign;
                cachedGammaLn[i] = gammaLogTab[(int)(2.0 * gammaValM[i] + 0.5)];
            }
        }
    }else if(nodeDelta == false){
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
        for(int i = 0; i < nf && changedIntervals.empty() == false; i++){
            if(gammaStale[i]) continue;
            double ti = unresolvedFossils->isSA(i) ? unresolvedFossils->getFossilAge(i) : unresolvedFossils->getAttachAge(i);
            for(std::pair<double,double>& iv : changedIntervals)
                if(ti > iv.first && ti < iv.second){ gammaStale[i] = 1; break; }
        }
    }

    if(anyFossilMoved || stalkIndexBuilt == false){
        if(stalkIndexFresh == false)
            rebuildStalkIndex();
        stalkIndexBuilt = true;
    }
    stalkIndexFresh = false;

    std::vector<int> staleIdx;
    for(int i = 0; i < nf; i++)
        if(gammaStale[i])
            staleIdx.push_back(i);
    gammaStaleTot += (long long)staleIdx.size();
    gammaStaleCalls++;
    if(eulerBuilt == false)
        buildEulerIndex();
    ThreadPool::current().parallelFor(OP_GAMMA, (int)staleIdx.size(), [&](int a, int b){
        for(int k = a; k < b; k++){
            int i = staleIdx[k];
            double g = computeGamma(unresolvedFossils->getAttachAge(i), i);
            cachedGammaLn[i] = (g > 0.0) ? std::log(g) : -INFINITY;
            gammaValM[i] = g;
            gammaStale[i] = 0;
        }
    });

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
        if(unresolvedFossils->maxAttachIsNodeBounded(i) == false)
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
        std::vector<Node*> attachEdges;
        std::vector<double> los;
        std::vector<double> his;
        enumerateFossilAttachEdges(t, crown, origin, isCrown, y, attachEdges, los, his);
        int k = (int)(rv.uniformRv() * attachEdges.size());
        double z = los[k] + rv.uniformRv() * (his[k] - los[k]);
        t->insertFossilTip(attachEdges[k], f.getTaxon(), y, z);

        fossilName.push_back(f.getTaxon());
        fossilIsCrown.push_back(isCrown);
    }
}

void FBDTreeModel::enumerateFossilAttachEdges(Tree* t, Node* crown, Node* origin, bool isCrown, double y, std::vector<Node*>& attachEdges, std::vector<double>& los, std::vector<double>& his){
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
            attachEdges.push_back(n);
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
    Node* attachParent = split->getAncestor();
    Node* attachChild = nullptr;
    for(Node* nb : split->getNeighbors())
        if(nb != attachParent && nb != r)
            attachChild = nb;

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
    double oldRange = tree->isSATip(attachChild) ? 1.0 : std::min(ceilingS, attachParent->getTime()) - std::max(rAge, attachChild->getTime());

    attachParent->removeNeighbor(split);
    split->removeNeighbor(attachParent);
    attachChild->removeNeighbor(split);
    split->removeNeighbor(attachChild);
    attachParent->addNeighbor(attachChild);
    attachChild->addNeighbor(attachParent);
    attachChild->setAncestor(attachParent);
    tree->initializeDownPassSequence();

    std::vector<Node*> attachEdges;
    std::vector<double> los;
    std::vector<double> his;
    enumerateSubtreeAttachEdges(tree, crowns, isCrowns, origins, rAge, ceilingS, attachEdges, los, his);
    int k = (int)(rng.uniformRv() * (double)attachEdges.size());
    Node* newChild = attachEdges[k];
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
    if(f->getIsSA()){
        sp->setTime(y + rng.uniformRv() * range);
        f->setIsSA(false);
        return std::log(range);
    }
    if(y < sib->getTime())
        return -INFINITY;
    if(sib->getIsTip() && sib->getIsFossil() && sib->getIsSA())
        return -INFINITY;
    sp->setTime(y);
    f->setIsSA(true);
    return -std::log(range);
}

double FBDTreeModel::doUpDownScale(void){
    lastMoveKind = MK_UPDOWN;
    if(upDownTotal > 0 && upDownTotal % 100 == 0){
        double ar = 0.0;
        for(bool b : upDownRecent)
            if(b)
                ar++;
        if(upDownRecent.empty() == false)
            ar /= upDownRecent.size();
        upDownStep *= std::exp((1.0 / std::sqrt((double)(upDownTotal / 100))) * (ar - 0.3));
    }

    double mBact = 0.95;
    double sBact = std::sqrt(1.0 - mBact * mBact);
    double delta = mBact + Probability::Normal::rv(&rng) * sBact;
    if(Probability::Uniform::rv(&rng, 0.0, 1.0) < 0.5)
        delta = -delta;
    double c = std::exp(upDownStep * delta);
    double invc = 1.0 / c;
    double lnc = std::log(c);

    int nUp = 0;
    double h = 0.0;
    for(ParameterDouble* l : a.lambda) if(l->getProposalProbability() > 0.0){ l->scaleProposed(c); nUp++; }
    for(ParameterDouble* m : a.mu) if(m->getProposalProbability() > 0.0){ m->scaleProposed(c); nUp++; }
    for(auto& pv : a.psi) for(ParameterDouble* p : pv) if(p->getProposalProbability() > 0.0){ p->scaleProposed(c); nUp++; }
    if(a.lambdaField != nullptr) h += a.lambdaField->scaleAllProposed(c, true);
    if(a.muField != nullptr)     h += a.muField->scaleAllProposed(c, true);
    for(int tp = 0; tp < a.numPsiTypes; tp++)
        if(a.psiField[tp] != nullptr) h += a.psiField[tp]->scaleAllProposed(c, true);
    h += nUp * lnc;

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
    if(originAge != nullptr && originAge->getProposalProbability() > 0.0){
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
    parameterTree->setAgeFloors(floors);
    parameterTree->getTree()->setLastUpdateWasScale(false);
}

double FBDTreeModel::doJointScale(void){
    Tree* tree = parameterTree->getTree();
    jointScaleAttW++;
    if(jointScaleAttW >= 200){
        jointScaleAdapt++;
        double ar = (double)jointScaleAccW / (double)jointScaleAttW;
        double gain = 1.0 / std::sqrt((double)jointScaleAdapt);
        jointScaleStep *= std::exp(gain * (ar - 0.3));
        if(jointScaleStep < 1e-4) jointScaleStep = 1e-4;
        if(jointScaleStep > 5.0)  jointScaleStep = 5.0;
        jointScaleAccW = 0;
        jointScaleAttW = 0;
    }
    double m = std::exp(jointScaleStep * (rng.uniformRv() - 0.5));
    tree->setLastUpdateWasScale(true);
    double zJac = unresolvedFossils->scaleAllAttachAges(m);
    if(zJac == -INFINITY)
        return -INFINITY;
    int numScaled = tree->scaleInternalAges(m);
    for(Node* n : tree->getDownPassSequence())
        if(n != tree->getRoot() && tree->isSATip(n) == false && n->getTime() >= n->getAncestor()->getTime())
            return -INFINITY;
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
    subtreeAttW++;
    if(subtreeAttW >= 200){
        subtreeAdapt++;
        double ar = (double)subtreeAccW / (double)subtreeAttW;
        double gain = 1.0 / std::sqrt((double)subtreeAdapt);
        subtreeStep *= std::exp(gain * (ar - 0.3));
        if(subtreeStep < 1e-2) subtreeStep = 1e-2;
        if(subtreeStep > 2.0)  subtreeStep = 2.0;
        subtreeAccW = 0;
        subtreeAttW = 0;
    }
    double newAge = oldAge + subtreeStep * (parentAge - oldestTip) * (rng.uniformRv() - 0.5);
    while(newAge <= oldestTip || newAge >= parentAge){
        if(newAge <= oldestTip) newAge = 2.0 * oldestTip - newAge;
        if(newAge >= parentAge) newAge = 2.0 * parentAge - newAge;
    }
    double sf = newAge / oldAge;

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
    for(Node* n : tree->getDownPassSequence())
        if(n != tree->getRoot() && tree->isSATip(n) == false && n->getTime() >= n->getAncestor()->getTime())
            return -INFINITY;
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

void FBDTreeModel::enumerateSubtreeAttachEdges(Tree* t, std::vector<Node*>& crowns, std::vector<char>& isCrowns, std::vector<Node*>& origins, double rAge, double ceilingS, std::vector<Node*>& attachEdges, std::vector<double>& los, std::vector<double>& his){
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
            attachEdges.push_back(n);
            los.push_back(lo);
            his.push_back(hi);
        }
        if(n->getIsTip() && n->getIsFossil() && t->isSATip(n) == false && rAge < n->getTime() && n->getTime() < ceilingS){
            attachEdges.push_back(n);
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
