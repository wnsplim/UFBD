#include <map>
#include <functional>
#include <iomanip>
#include "Node.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "ParameterOUField.hpp"
#include "FBDTreeModel.hpp"
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

FBDTreeModel::FBDTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, unsigned int seed) :
    PhylogeneticModel(){

    lastMoveKind = MK_PARAM;
    lastRateVec = nullptr;
    rateVecStep = 0.2;
    shrinkStep = 0.2;
    rvAccW = rvAttW = rvAtt = seAccW = seAttW = seAtt = 0;
    lastTreeMove = TM_NONE;
    for(int i = 0; i < TM_COUNT; i++) tmAcc[i] = tmAtt[i] = 0;
    azAcc = 0;
    azAtt = 0;
    shiftStep = 0.005;
    saBatch = 0;
    saBatchF = 1.0;
    rsAccW = rsAttW = 0;
    rsAcc = rsTot = rsAdapt = 0;
    upDownStep = 0.1;
    upDownTotal = 0;
    cacheInit = false;
    zoneInit = false;
    lambdaField = nullptr;
    muField = nullptr;
    numPsiTypes = UserSettings::userSettings().getNumPsiTypes();
    rng.setSeed(seed);
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);

    parameterTree = new ParameterTree(1.0, this);
    isResolved = (UserSettings::userSettings().getModel() == Model::RFBD);

    const std::vector<std::string>& psiTypeNames = UserSettings::userSettings().getPsiTypeNames();
    std::map<std::string,int> typeIdx;
    for(size_t k = 0; k < psiTypeNames.size(); k++) typeIdx[psiTypeNames[k]] = (int)k;
    auto resolveType = [&](Fossil& f) -> int {
        if(psiTypeNames.empty() || f.getMaxAge() <= 0.0) return 0;
        if(f.getType().empty())
            Msg::error("fossil '" + f.getTaxon() + "' has no preservation type, but -psi_types declares multiple types.");
        std::map<std::string,int>::const_iterator it = typeIdx.find(f.getType());
        if(it == typeIdx.end())
            Msg::error("fossil '" + f.getTaxon() + "' has preservation type '" + f.getType() + "' not listed in -psi_types.");
        return it->second;
    };

    std::vector<Fossil> unrFossils;
    std::vector<std::string> bbClade;
    std::vector<Assignment> bbAsg;
    if(isResolved == false){
        for(Fossil& f : fossils){
            Node* tip = t->getTaxonNode(f.getTaxon());
            if(tip == nullptr || tip->getIsTip() == false){ unrFossils.push_back(f); continue; }
            tip->setIsFossil(true);
            tip->setTime(0.5 * (f.getMinAge() + f.getMaxAge()));
            tip->setFossilAgeRange(f.getMinAge(), f.getMaxAge());
            backboneFossils.push_back({tip, f.getMinAge(), f.getMaxAge(), resolveType(f)});
            bbClade.push_back(f.getClade());
            bbAsg.push_back(f.getAssignment());
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
                        if(g.getTaxon() == fTaxon) continue;
                        std::map<std::string, Node*>::iterator ci = cladeCrown.find(g.getClade());
                        if(ci == cladeCrown.end()) continue;
                        if(isUnder(fTip, ci->second) == false) continue;
                        if(g.getClade() == bbClade[k] &&
                           ((bbAsg[k] == Assignment::CROWN && g.getAssignment() == Assignment::STEM) ||
                            (bbAsg[k] == Assignment::STEM  && g.getAssignment() == Assignment::CROWN)))
                            continue;
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
            double hi = us.getConditionAgePriorP2();
            if(us.getConditionAgePrior() == Probability::PriorFamily::UNIFORM && x0init >= hi)
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
    
    UserSettings& rateUs = UserSettings::userSettings();
    std::vector<double> lambdaTimes, muTimes;
    lambdaTimes.push_back(0.0);
    for(double t : rateUs.getLambdaSkylineTimes()) lambdaTimes.push_back(t);
    muTimes.push_back(0.0);
    for(double t : rateUs.getMuSkylineTimes())     muTimes.push_back(t);
    std::vector<std::vector<double>> psiTimes(numPsiTypes);
    for(int tp = 0; tp < numPsiTypes; tp++){
        psiTimes[tp].push_back(0.0);
        for(double t : rateUs.getPsiSkylineTimes(tp)) psiTimes[tp].push_back(t);
    }
    intervalStart.push_back(0.0);
    for(double t : rateUs.getSkylineTimes())
        intervalStart.push_back(t);
    psiIdx.assign(numPsiTypes, std::vector<int>());
    for(double s : intervalStart){
        int li = 0, mi = 0;
        for(int k = 0; k < (int)lambdaTimes.size(); k++) if(lambdaTimes[k] <= s) li = k;
        for(int k = 0; k < (int)muTimes.size(); k++)     if(muTimes[k] <= s)     mi = k;
        lambdaIdx.push_back(li);
        muIdx.push_back(mi);
        for(int tp = 0; tp < numPsiTypes; tp++){
            int pi = 0;
            for(int k = 0; k < (int)psiTimes[tp].size(); k++) if(psiTimes[tp][k] <= s) pi = k;
            psiIdx[tp].push_back(pi);
        }
    }
    Probability::PriorSpec lp = rateUs.getLambdaPrior();
    Probability::PriorSpec mp = rateUs.getMuPrior();
    if(rateUs.getLambdaMode() == RateMode::OU && lp.set) Msg::error("-lambda_prior is only used under mode=indep; under mode=ou set the OU level with -lambda_ou_theta.");
    if(rateUs.getMuMode() == RateMode::OU && mp.set) Msg::error("-mu_prior is only used under mode=indep; under mode=ou set the OU level with -mu_ou_theta.");
    Probability::PriorSpec defRate{true, Probability::PriorFamily::EXPONENTIAL, 5.0, 1.0};
    if(!lp.set) lp = defRate;
    if(!mp.set) mp = defRate;
    int nLambda = (int)lambdaTimes.size();
    int nMu = (int)muTimes.size();
    double lam0 = Probability::priorMean(lp.family, lp.p1, lp.p2, lp.p3);
    if(!(lam0 > 0.0 && std::isfinite(lam0))) lam0 = 0.1;
    double mu0 = Probability::priorMean(mp.family, mp.p1, mp.p2, mp.p3);
    if(!(mu0 > 0.0 && std::isfinite(mu0))) mu0 = 0.1;
    mu0 = offsetFromLambda(mu0, lam0, mp);
    {
        std::vector<int> b2c = buildSkylineRates("lambda", "", nLambda, lambdaTimes, rateUs.getLambdaMode(), rateUs.getLambdaOU(), lp, lam0, rateUs.getLambdaGroups(), rateUs.getLambdaGroupPrior(), lambda, lambdaField, lambdaName);
        for(int& u : lambdaIdx) u = b2c[u];
        appendRateMap(lambdaTimes, b2c, lambdaName);
    }
    {
        std::vector<int> b2c = buildSkylineRates("mu", "", nMu, muTimes, rateUs.getMuMode(), rateUs.getMuOU(), mp, mu0, rateUs.getMuGroups(), rateUs.getMuGroupPrior(), mu, muField, muName);
        for(int& u : muIdx) u = b2c[u];
        appendRateMap(muTimes, b2c, muName);
    }
    psi.assign(numPsiTypes, std::vector<ParameterDouble*>());
    psiField.assign(numPsiTypes, nullptr);
    psiName.assign(numPsiTypes, std::vector<std::string>());
    for(int tp = 0; tp < numPsiTypes; tp++){
        Probability::PriorSpec pp = rateUs.getPsiPrior(tp);
        if(rateUs.getPsiMode(tp) == RateMode::OU && pp.set) Msg::error("-psi_prior is only used under mode=indep; under mode=ou set the OU level with -psi_ou_theta.");
        if(!pp.set) pp = defRate;
        double psi0 = Probability::priorMean(pp.family, pp.p1, pp.p2, pp.p3);
        if(!(psi0 > 0.0 && std::isfinite(psi0))) psi0 = 0.1;
        psi0 = offsetFromLambda(psi0, lam0, pp);
        int nPsi = (int)psiTimes[tp].size();
        std::string prefix = (numPsiTypes > 1) ? ("psi_" + rateUs.getPsiTypeNames()[tp]) : "psi";
        std::string sep = (numPsiTypes > 1) ? "_" : "";
        std::vector<int> b2c = buildSkylineRates(prefix, sep, nPsi, psiTimes[tp], rateUs.getPsiMode(tp), rateUs.getPsiOU(tp), pp, psi0, rateUs.getPsiGroups(tp), rateUs.getPsiGroupPrior(tp), psi[tp], psiField[tp], psiName[tp]);
        for(int& u : psiIdx[tp]) u = b2c[u];
        appendRateMap(psiTimes[tp], b2c, psiName[tp]);
    }
    rho = UserSettings::userSettings().getRho();
    fossilType.assign((int)unrFossils.size(), 0);
    for(size_t i = 0; i < unrFossils.size(); i++){
        fossilType[i] = resolveType(unrFossils[i]);
        fossilTypeByName[unrFossils[i].getTaxon()] = fossilType[i];
        unrFossilName.push_back(unrFossils[i].getTaxon());
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
    }

    if(isResolved){
        parameterTree->setProposalProbability(78.0);
        for(ParameterDouble* l : lambda) l->setProposalProbability(15.0);
        for(ParameterDouble* m : mu)     m->setProposalProbability(15.0);
        for(auto& pv : psi) for(ParameterDouble* p : pv) p->setProposalProbability(15.0);
    }
    double fieldBase = (isResolved ? 15.0 : 1.0);
    if(lambdaField) lambdaField->setProposalProbability(fieldBase * (double)lambdaField->getNumBins());
    if(muField)     muField->setProposalProbability(fieldBase * (double)muField->getNumBins());
    for(int tp = 0; tp < numPsiTypes; tp++) if(psiField[tp]) psiField[tp]->setProposalProbability(fieldBase * (double)psiField[tp]->getNumBins());

    double sum = 0.0;
    for(Parameter* p : parameters)
        sum += p->getProposalProbability();
    for(Parameter* p : parameters)
        p->setProposalProbability(p->getProposalProbability() / sum);

    RandomVariable::setActiveInstance(prevRng);
}

std::vector<int> FBDTreeModel::buildSkylineRates(const std::string& prefix, const std::string& sep, int nBins, const std::vector<double>& times, RateMode mode, const OUHyperSpec& ou, const Probability::PriorSpec& basePrior, double rate0, const std::vector<int>& groupIds, const std::map<int,Probability::PriorSpec>& groupPrior, std::vector<ParameterDouble*>& outVec, ParameterOUField*& outField, std::vector<std::string>& outNames){
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
    for(int c = 0; c < nChunks; c++)
        outNames.push_back(prefix + ((nChunks > 1) ? (sep + std::to_string(chunkMinBin[c])) : ""));

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
        double na = std::numeric_limits<double>::quiet_NaN();
        outField = new ParameterOUField(1.0, this, nChunks, loEdges, rate0, originAge,
            ou.thetaSet ? ou.thetaMedian : na, ou.thetaSet ? ou.thetaSd : na,
            ou.sdSet ? ou.sdShape : na, ou.sdSet ? ou.sdRate : na,
            ou.nuSet ? ou.nuShape : na, ou.nuSet ? ou.nuRate : na);
        parameters.push_back(outField);
        return binToChunk;
    }

    for(int c = 0; c < nChunks; c++){
        ParameterDouble* p = new ParameterDouble(1.0, this, outNames[c], 0.0, std::numeric_limits<double>::max());
        std::map<int,Probability::PriorSpec>::const_iterator pit = groupPrior.find(chunkGid[c]);
        bool ov = (pit != groupPrior.end() && pit->second.set);
        const Probability::PriorSpec& sp = ov ? pit->second : basePrior;
        p->setPrior(sp.family, sp.p1, sp.p2);
        p->setValue(ov ? Probability::priorMean(sp.family, sp.p1, sp.p2, sp.p3) : rate0);
        outVec.push_back(p);
        parameters.push_back(p);
    }
    return binToChunk;
}

void FBDTreeModel::appendRateMap(const std::vector<double>& times, const std::vector<int>& binToChunk, const std::vector<std::string>& names){
    int nBins = (int)binToChunk.size();
    for(int c = 0; c < (int)names.size(); c++){
        std::string iv;
        bool first = true;
        for(int i = 0; i < nBins; i++){
            if(binToChunk[i] != c) continue;
            char buf[64];
            if(i + 1 < nBins) std::snprintf(buf, sizeof(buf), "[%g, %g)", times[i], times[i + 1]);
            else              std::snprintf(buf, sizeof(buf), "[%g, inf)", times[i]);
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
    int j = lambdaIdx[i];
    return (lambdaField != nullptr) ? lambdaField->getRate(j) : lambda[j]->getValue();
}

double FBDTreeModel::muAt(int i){
    int j = muIdx[i];
    if(muField != nullptr)
        return muField->getRate(j);
    return mu[j]->getValue();
}

double FBDTreeModel::psiOfTypeAt(int type, int i){
    int j = psiIdx[type][i];
    return (psiField[type] != nullptr) ? psiField[type]->getRate(j) : psi[type][j]->getValue();
}

double FBDTreeModel::psiTotalAt(int i){
    double s = 0.0;
    for(int t = 0; t < numPsiTypes; t++)
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

std::vector<std::string> FBDTreeModel::getParameterNames(void){
    std::vector<std::string> names;
    if(lambdaField != nullptr)
        for(const std::string& n : lambdaName) names.push_back(n);
    else
        for(ParameterDouble* p : lambda) names.push_back(p->getName());
    if(muField != nullptr)
        for(const std::string& n : muName) names.push_back(n);
    else
        for(ParameterDouble* p : mu) names.push_back(p->getName());
    for(int tp = 0; tp < numPsiTypes; tp++){
        if(psiField[tp] != nullptr)
            for(const std::string& n : psiName[tp]) names.push_back(n);
        else
            for(ParameterDouble* p : psi[tp]) names.push_back(p->getName());
    }
    if(unresolvedFossils != nullptr || isResolved)
        names.push_back("nSA");
    if(originAge != nullptr)
        names.push_back("originAge");
    for(size_t i = 0; i < getAgeLogNodes().size(); i++)
        names.push_back("x" + std::to_string(i + 1));
    return names;
}

std::vector<bool> FBDTreeModel::getParameterFixedMask(void){
    std::vector<bool> m;
    if(lambdaField != nullptr) for(size_t i = 0; i < lambdaName.size(); i++) m.push_back(false);
    else for(ParameterDouble* p : lambda) m.push_back(p->getProposalProbability() == 0.0);
    if(muField != nullptr) for(size_t i = 0; i < muName.size(); i++) m.push_back(false);
    else for(ParameterDouble* p : mu) m.push_back(p->getProposalProbability() == 0.0);
    for(int tp = 0; tp < numPsiTypes; tp++){
        if(psiField[tp] != nullptr) for(size_t i = 0; i < psiName[tp].size(); i++) m.push_back(false);
        else for(ParameterDouble* p : psi[tp]) m.push_back(p->getProposalProbability() == 0.0);
    }
    if(unresolvedFossils != nullptr || isResolved) m.push_back(false);
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
    std::vector<double> vals;
    if(lambdaField != nullptr)
        for(int i = 0; i < lambdaField->getNumBins(); i++) vals.push_back(lambdaField->getRate(i));
    else
        for(ParameterDouble* p : lambda) vals.push_back(p->getValue());
    if(muField != nullptr)
        for(int i = 0; i < muField->getNumBins(); i++) vals.push_back(muField->getRate(i));
    else
        for(ParameterDouble* p : mu) vals.push_back(p->getValue());
    for(int tp = 0; tp < numPsiTypes; tp++){
        if(psiField[tp] != nullptr)
            for(int i = 0; i < psiField[tp]->getNumBins(); i++) vals.push_back(psiField[tp]->getRate(i));
        else
            for(ParameterDouble* p : psi[tp]) vals.push_back(p->getValue());
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
        lnP += Probability::priorLnPdf(us.getConditionAgePrior(), us.getConditionAgePriorP1(), us.getConditionAgePriorP2(), crownAge, 0.0, std::numeric_limits<double>::max(), us.getConditionAgePriorP3());
    }

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
        if(p != parameterTree)
            std::cout << p->getName() << " (A/R): " << p->getAcceptanceRatio() << "\t";
    }
    if(unresolvedFossils != nullptr)
        unresolvedFossils->print();
    static const char* tmName[TM_COUNT] = {"NE","WB","WE","treeScale","SARJ","upDown","jointScale","subtree","nodeAge","crown"};
    std::cout << "tree[";
    for(int i = 0; i < TM_COUNT; i++)
        if(tmAtt[i] > 0)
            std::cout << tmName[i] << ":" << (double)tmAcc[i] / (double)tmAtt[i] << " ";
    std::cout << "] treeScaleStep: " << parameterTree->getScaleLambda();
    if(rsTot > 0)
        std::cout << " rateShift[A/R:" << (double)rsAcc / (double)rsTot
                  << " step:" << std::setprecision(5) << shiftStep << std::setprecision(3)
                  << " batch:" << saBatch << "]";
    std::cout << "\n";
    os.flags(f);
    os.precision(pr);
}

std::vector<ParameterDouble*>* FBDTreeModel::pickIidRateVector(void){
    std::vector<std::vector<ParameterDouble*>*> cands;
    if(lambdaField == nullptr && lambda.size() >= 2) cands.push_back(&lambda);
    if(muField == nullptr && mu.size() >= 2)         cands.push_back(&mu);
    for(int tp = 0; tp < numPsiTypes; tp++) if(psiField[tp] == nullptr && psi[tp].size() >= 2) cands.push_back(&psi[tp]);
    if(cands.empty())
        return nullptr;
    return cands[(int)(rng.uniformRv() * cands.size())];
}

double FBDTreeModel::doRateVectorScale(void){
    lastMoveKind = MK_RATEVEC;
    lastRateVecScale = true;
    lastRateVec = pickIidRateVector();
    double mB = 0.95;
    double d = mB + Probability::Normal::rv(&rng) * std::sqrt(1.0 - mB * mB);
    if(rng.uniformRv() < 0.5) d = -d;
    double c = std::exp(rateVecStep * d);
    for(ParameterDouble* p : *lastRateVec)
        p->scaleProposed(c);
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
    return (double)lastRateVec->size() * std::log(c);
}

double FBDTreeModel::doRateShrinkExpand(void){
    lastMoveKind = MK_RATEVEC;
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
    seAtt++;
    if(seAttW >= 200){
        double gain = 1.0 / std::sqrt((double)(seAtt / 200));
        shrinkStep *= std::exp(gain * ((double)seAccW / seAttW - 0.3));
        if(shrinkStep < 1e-3) shrinkStep = 1e-3;
        if(shrinkStep > 10.0)  shrinkStep = 10.0;
        seAccW = 0;
        seAttW = 0;
    }
    return (double)(n - 1) * std::log(a);
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
    if(lambdaField != nullptr){
        double j = lambdaField->shiftRates(d);
        if(std::isinf(j)) return -std::numeric_limits<double>::infinity();
        lnH += j;
    }else{
        for(ParameterDouble* l : lambda)
            if(l->getValue() + d <= 0.0)
                return -std::numeric_limits<double>::infinity();
        for(ParameterDouble* l : lambda)
            l->shiftProposed(d);
    }
    if(muField != nullptr){
        double j = muField->shiftRates(d);
        if(std::isinf(j)) return -std::numeric_limits<double>::infinity();
        lnH += j;
    }else{
        for(ParameterDouble* mv : mu)
            if(mv->getValue() + d <= 0.0)
                return -std::numeric_limits<double>::infinity();
        for(ParameterDouble* mv : mu)
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
int countStraddling(const std::vector<double>& los, const std::vector<double>& his, double zq){
    int below = (int)(std::lower_bound(los.begin(), los.end(), zq) - los.begin());
    int above = (int)(std::upper_bound(his.begin(), his.end(), zq) - his.begin());
    return below - above;
}
// stalks as (z_attach,y) sorted by z_attach; y passed separately (sorted)
int countStraddling(const std::vector<double>& ySorted, const std::vector<std::pair<double,double>>& zy, double zq){
    int below = (int)(std::lower_bound(ySorted.begin(), ySorted.end(), zq) - ySorted.begin());
    int above = (int)(std::upper_bound(zy.begin(), zy.end(), zq,
                     [](double v, const std::pair<double,double>& p){ return v < p.first; }) - zy.begin());
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
    bool haveIid = (lambdaField == nullptr && lambda.size() >= 2)
                || (muField == nullptr && mu.size() >= 2);
    for(int tp = 0; tp < numPsiTypes && !haveIid; tp++)
        if(psiField[tp] == nullptr && psi[tp].size() >= 2) haveIid = true;
    lastTreeMove = TM_NONE;
    if(haveIid && rng.uniformRv() < 0.20){
        double r = (rng.uniformRv() < 0.5) ? doRateVectorScale() : doRateShrinkExpand();
        RandomVariable::setActiveInstance(prevRng);
        return r;
    }
    if((lambdaField != nullptr || muField != nullptr) && rng.uniformRv() < 0.10){
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
        double uMove = rng.uniformRv() * (slideAndCrown + 3.0 * blockUnit);
        if(uMove >= slideAndCrown + 2.0 * blockUnit){
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

    RandomVariable::setActiveInstance(prevRng);
    return ratio;
}

void FBDTreeModel::updateForAcceptance(void){
    if(lastTreeMove != TM_NONE){ tmAtt[lastTreeMove]++; tmAcc[lastTreeMove]++; }
    if(lastMoveKind == MK_AZGIBBS){
        unresolvedFossils->updateForAcceptance();
        azAcc++;
        return;
    }
    if(lastMoveKind == MK_RATESHIFT){
        if(lambdaField != nullptr) lambdaField->commitProposed();
        else for(ParameterDouble* l : lambda) l->commitProposed();
        if(muField != nullptr) muField->commitProposed();
        else for(ParameterDouble* m : mu) m->commitProposed();
        if(unresolvedFossils != nullptr) unresolvedFossils->updateForAcceptance();
        rsAccW++;
        rsAcc++;
        return;
    }
    if(lastMoveKind == MK_RATEVEC){
        for(ParameterDouble* p : *lastRateVec)
            p->commitProposed();
        if(lastRateVecScale) rvAccW++; else seAccW++;
        return;
    }
    if(lastMoveKind == MK_UPDOWN){
        for(ParameterDouble* l : lambda) l->commitProposed();
        for(ParameterDouble* m : mu) m->commitProposed();
        for(auto& pv : psi) for(ParameterDouble* p : pv) p->commitProposed();
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
    if(lastTreeMove != TM_NONE) tmAtt[lastTreeMove]++;
    if(lastMoveKind == MK_AZGIBBS){
        unresolvedFossils->updateForRejection();
        rebuildStalkIndex();
        return;
    }
    if(lastMoveKind == MK_RATESHIFT){
        if(lambdaField != nullptr) lambdaField->restoreProposed();
        else for(ParameterDouble* l : lambda) l->restoreProposed();
        if(muField != nullptr) muField->restoreProposed();
        else for(ParameterDouble* m : mu) m->restoreProposed();
        if(unresolvedFossils != nullptr){
            unresolvedFossils->updateForRejection();
            rebuildStalkIndex();
        }
        return;
    }
    if(lastMoveKind == MK_RATEVEC){
        for(ParameterDouble* p : *lastRateVec)
            p->restoreProposed();
        return;
    }
    if(lastMoveKind == MK_UPDOWN){
        for(ParameterDouble* l : lambda) l->restoreProposed();
        for(ParameterDouble* m : mu) m->restoreProposed();
        for(auto& pv : psi) for(ParameterDouble* p : pv) p->restoreProposed();
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
}

void FBDTreeModel::readState(std::istream& is){
    for(Parameter* p : parameters)
        p->readState(is);
    is >> rateVecStep >> shrinkStep >> shiftStep >> saBatch >> saBatchF >> upDownStep >> upDownTotal;
    is >> rvAccW >> rvAttW >> rvAtt >> seAccW >> seAttW >> seAtt >> rsAccW >> rsAttW;
    is >> rsAcc >> rsTot >> rsAdapt >> azAcc >> azAtt;
    Serialize::readBoolDeque(is, upDownRecent);
    for(int i = 0; i < TM_COUNT; i++) is >> tmAcc[i] >> tmAtt[i];
    cacheInit = false;
    zoneInit = true;
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
    static const bool fbdDbg = (std::getenv("FBD_DBG") != nullptr);
    double dbgT1 = 0.0, dbgT12 = 0.0;

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
        double lx0 = lambdaAt(findIndex(x0));
        fbdProb -= std::log(lx0);
        fbdProb -= calculateLnConditioning(x0);
        fbdProb += std::log(4 * lx0);
        fbdProb += lnD(x0) - std::log(4.0);
    }else{
        double lr = lambdaAt(findIndex(crownAge));
        fbdProb -= 2 * (std::log(lr) + calculateLnSurvival(crownAge));
        fbdProb += std::log(4 * lr);
        fbdProb += lnD(crownAge) - std::log(4.0);
    }

    dbgT1 = fbdProb;
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
                termNode[idx] = std::log(lambdaAt(findIndex(n->getTime()))) + lnD(n->getTime());
        }
    });
    for(int idx = 0; idx < nDp; idx++)
        fbdProb += termNode[idx];
    fbdProb += numExtantTips * std::log(rhoVal);
    dbgT12 = fbdProb;

    //term 3: fossil attachment
    int numFossils = unresolvedFossils->getNumFossils();
    updateGammaCache();
    int spineIdx = unresolvedFossils->getSpineIdx();
    std::vector<double> termFoss(numFossils, 0.0);
    ThreadPool::current().parallelFor(OP_FBD, numFossils, [&](int lo, int hi){
        for(int i = lo; i < hi; i++)
            termFoss[i] = fossilTermLn(i, spineIdx);
    });
    for(int i = 0; i < numFossils; i++)
        fbdProb += termFoss[i];
    for(const BackboneFossil& bf : backboneFossils){
        double y = bf.tip->getTime();
        fbdProb += std::log(psiOfTypeAt(bf.type, findIndex(y))) + std::log(calculateP0(y)) - lnD(y);
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
        int ty = fossilType[fam];
        std::cout << "[FBD_DBG_FOSS] i=" << fam << " term=" << fv << " isSA=" << unresolvedFossils->isSA(fam) << " isUE=" << unresolvedFossils->isUE(fam)
                  << " y=" << yy << " z=" << zz << " biny=" << findIndex(yy) << " binz=" << findIndex(zz)
                  << " logPsiY=" << std::log(psiOfTypeAt(ty, findIndex(yy))) << " log2LamZ=" << std::log(2 * lambdaAt(findIndex(zz)))
                  << " logP0Y=" << std::log(calculateP0(yy)) << " lnDz=" << lnD(zz) << " lnDy=" << lnD(yy) << " gammaLn=" << cachedGammaLn[fam] << "\n";
        int ni = (int)intervalStart.size();
        for(int i = 0; i < ni; i++){
            double p0 = calculateP0At(i, (i + 1 < ni) ? intervalStart[i + 1] : intervalStart[i] + 1.0);
            std::cout << "[FBD_DBG_IV] i=" << i << " start=" << intervalStart[i] << " lam=" << lambdaAt(i) << " mu=" << muAt(i) << " psi=" << psiTotalAt(i)
                      << " c1=" << c1Vec[i] << " c2=" << c2Vec[i] << " P0end=" << p0 << " ePrev=" << (i < (int)ePrev.size() ? ePrev[i] : 0.0) << "\n";
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
                double blo = intervalStart[b], bhi = (b + 1 < ni) ? intervalStart[b + 1] : xtop;
                double ov = std::min(zv, bhi) - std::max(yv, blo);
                if(ov > 0.0) stalk[b] += ov;
            }
        }
        for(int b = 0; b < ni; b++)
            std::cout << "[FBD_BAL] i=" << b << " [" << intervalStart[b] << "," << ((b + 1 < ni) ? intervalStart[b + 1] : xtop) << ") lam=" << lambdaAt(b)
                      << " nDiv=" << divc[b] << " stalkLen=" << stalk[b] << " nFossilY=" << fy[b] << "\n";
    }
    return fbdProb;
}

double FBDTreeModel::fossilTermLn(int i, int spineIdx){
    if(unresolvedFossils->isSA(i))
        return std::log(psiOfTypeAt(fossilType[i], findIndex(unresolvedFossils->getFossilAge(i)))) + cachedGammaLn[i];
    if(i == spineIdx && unresolvedFossils->isUE(i))
        return 0.0;
    if(unresolvedFossils->isUE(i))
        return uePqLn(unresolvedFossils->getAttachAge(i)) + cachedGammaLn[i];
    if(i == spineIdx){
        double ys = unresolvedFossils->getFossilAge(i);
        return std::log(psiOfTypeAt(fossilType[i], findIndex(ys))) + std::log(calculateP0(ys)) - lnD(ys);
    }
    return fossilPqLn(unresolvedFossils->getFossilAge(i), unresolvedFossils->getAttachAge(i), fossilType[i]) + cachedGammaLn[i];
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
    zoneStalks[z].zy.clear();
    for(int i : fos){
        if(unresolvedFossils->isSA(i))
            continue;
        double yi = unresolvedFossils->getFossilAge(i);
        zoneStalks[z].y.push_back(yi);
        zoneStalks[z].zy.push_back(std::make_pair(unresolvedFossils->getAttachAge(i), yi));
    }
    std::sort(zoneStalks[z].y.begin(), zoneStalks[z].y.end());
    std::sort(zoneStalks[z].zy.begin(), zoneStalks[z].zy.end());
}

void FBDTreeModel::zoneRecomputeGamma(const std::vector<int>& fos){
    for(int i : fos){
        double g = computeGamma(unresolvedFossils->getAttachAge(i), i);
        cachedGammaLn[i] = (g > 0.0) ? std::log(g) : -INFINITY;
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
    rhoVal = rho;
    prepareIntervals();
    updateGammaCache();
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

    double t3 = term3Sum();
    for(int idx : order){
        double ratio = unresolvedFossils->proposeOneFossil(idx);
        if(ratio == -INFINITY){
            unresolvedFossils->updateForRejection();
            continue;
        }
        updateGammaCache();
        double t3new = term3Sum();
        if(std::log(rng.uniformRv()) < (t3new - t3) + ratio){
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

    ThreadPool::current().parallelFor(OP_FBD, numZones, [&](int lo, int hi){
        for(int z = lo; z < hi; z++){
            std::vector<int>& fos = zoneFos[z];
            if(fos.empty())
                continue;
            RandomVariable zr;
            zr.setSeed(seed[z]);
            std::vector<int> ord = fos;
            for(int a = (int)ord.size() - 1; a > 0; a--)
                std::swap(ord[a], ord[(int)(zr.uniformRv() * (a + 1))]);
            double zt = zoneTermSum(fos, x0, useOrigin);
            for(int idx : ord){
                int sub;
                double ratio = unresolvedFossils->proposeFossilParallel(idx, zr, sub);
                zAtt[z][sub]++;
                if(ratio == -INFINITY){
                    unresolvedFossils->rejectFossil(idx);
                    zNRej[z]++;
                    continue;
                }
                rebuildZoneStalks(z, fos);
                zoneRecomputeGamma(fos);
                double ztnew = zoneTermSum(fos, x0, useOrigin);
                if(std::log(zr.uniformRv()) < (ztnew - zt) + ratio){
                    zt = ztnew;
                    unresolvedFossils->acceptFossil(idx);
                    zAcc[z][sub]++;
                    zNAcc[z]++;
                }else{
                    unresolvedFossils->rejectFossil(idx);
                    zNRej[z]++;
                    rebuildZoneStalks(z, fos);
                    zoneRecomputeGamma(fos);
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
            lnP += lnD(n->getAncestor()->getTime()) - lnD(n->getTime());
        if(n->getIsTip()){
            if(n->getIsFossil() == false)
                lnP += std::log(rhoVal);
            else{
                int ft = 0;
                if(numPsiTypes > 1){ std::map<std::string,int>::iterator it = fossilTypeByName.find(n->getName()); if(it != fossilTypeByName.end()) ft = it->second; }
                if(n->getIsSA())
                    lnP += std::log(psiOfTypeAt(ft, findIndex(n->getTime())));
                else
                    lnP += std::log(psiOfTypeAt(ft, findIndex(n->getTime()))) + std::log(calculateP0(n->getTime()));
            }
        }
        else if(n != root){
            bool fakeSplit = false;
            for(Node* c : n->getNeighbors())
                if(c != n->getAncestor() && c->getIsTip() && c->getIsFossil() && c->getIsSA()){ fakeSplit = true; break; }
            if(fakeSplit == false)
                lnP += std::log(2.0 * lambdaAt(findIndex(n->getTime())));
        }
    }
    return lnP;
}

double FBDTreeModel::lnD(double t){
    return (t <= 0.0) ? 0.0 : lnDPrev[findIndex(t)] + std::log(4.0) - calculateLnQtAt(findIndex(t), t);
}

double FBDTreeModel::fossilPqLn(double y, double z, int type){
    return std::log(psiOfTypeAt(type, findIndex(y))) + std::log(2*lambdaAt(findIndex(z))) + std::log(calculateP0(y)) + lnD(z) - lnD(y);
}

double FBDTreeModel::uePqLn(double z){
    return std::log(rhoVal) + std::log(2*lambdaAt(findIndex(z))) + lnD(z);
}

double FBDTreeModel::calculateLnSurvival(double t){
    int k = findIndex(t);
    double lam = lambdaAt(k), mu = muAt(k), c1 = c1HatVec[k], c2 = c2HatVec[k];
    double tau = t - intervalStart[k];
    if(lam < 1e-9 * c1){ // lambda->0: log(N)-log(2*lam) is a removable 0/0; exact lambda=0 limit
        double eb = (k == 0) ? (1.0 - rhoVal) : ePrevHat[k];
        return -mu * tau + std::log(1.0 - eb);
    }
    double E = std::exp(-c1 * tau);
    double N  = E * (1.0 - c2) * ((lam - mu) - c1) + (1.0 + c2) * ((lam - mu) + c1);
    double Dp = E * (1.0 - c2) + (1.0 + c2);
    return std::log(N) - std::log(2.0 * lam) - std::log(Dp);
}

double FBDTreeModel::calculateLnAnySample(double t){
    int k = findIndex(t);
    double lam = lambdaAt(k), mu = muAt(k), psi = psiTotalAt(k), c1 = c1Vec[k], c2 = c2Vec[k];
    double tau = t - intervalStart[k];
    if(lam < 1e-9 * c1){ // lambda->0: log(N)-log(2*lam) is a removable 0/0; exact lambda=0 limit
        double s = mu + psi;
        double eb = (k == 0) ? (1.0 - rhoVal) : ePrev[k];
        if(s <= 0.0)
            return std::log(1.0 - eb);
        double E = std::exp(-s * tau);
        return std::log((E * (1.0 - 2.0*eb) + (mu*E + psi*(2.0 - E)) / s) / 2.0);
    }
    double beta = lam - mu - psi, bpc, bmc;
    if(beta >= 0.0){ bpc = beta + c1; bmc = -4.0*lam*psi/bpc; }
    else           { bmc = beta - c1; bpc = -4.0*lam*psi/bmc; }
    double E = std::exp(-c1 * tau);
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
            int k = findIndex(t);
            double pExtinct = calculateP0HatAt(k, t) - calculateP0At(k, t);
            return (pExtinct > 0.0) ? std::log(pExtinct) : -INFINITY;
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
        double pi = psiTotalAt(i);
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
    double mi = muAt(i);
    double pi = psiTotalAt(i);
    if(li < 1e-7 * c1Vec[i]){ // lambda->0: tmp/(2*li) is an unstable 0/0; exact lambda=0 limit
        double s = mi + pi;
        double eb = (i == 0) ? (1.0 - rhoVal) : ePrev[i];
        if(s <= 0.0)
            return eb;
        double a = mi / s;
        return a + (eb - a) * std::exp(-s * tau);
    }
    double tmp = -li + mi + pi;
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
    double mi = muAt(i);
    if(li < 1e-7 * c1HatVec[i]){ // lambda->0: same 0/0 as P0At; exact lambda=0 limit (pure death)
        double eb = (i == 0) ? (1.0 - rhoVal) : ePrevHat[i];
        return 1.0 + (eb - 1.0) * std::exp(-mi * tau);
    }
    double tmp = -li + mi;
    tmp += c1HatVec[i] * (std::exp(-c1HatVec[i] * tau) * (1 - c2HatVec[i]) - (1+c2HatVec[i]) ) / ( std::exp(-c1HatVec[i] * tau) * (1 - c2HatVec[i]) + (1+c2HatVec[i])  );
    tmp /= 2*li;
    return 1 + tmp;
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
        zoneStalks[k].zy.clear();
    }
    int nf = unresolvedFossils->getNumFossils();
    for(int i = 0; i < nf; i++){
        if(unresolvedFossils->isSA(i))
            continue;
        StalkBucket& b = zoneStalks[unresolvedFossils->getAttachmentZone(i)];
        double yi = unresolvedFossils->getFossilAge(i);
        b.y.push_back(yi);
        b.zy.push_back(std::make_pair(unresolvedFossils->getAttachAge(i), yi));
    }
    for(size_t k = 0; k < zoneStalks.size(); k++){
        std::sort(zoneStalks[k].y.begin(),  zoneStalks[k].y.end());
        std::sort(zoneStalks[k].zy.begin(), zoneStalks[k].zy.end());
    }
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
    count += w * (double)countStraddling(b.y, b.zy, z);

    int sp = unresolvedFossils->getSpineIdx();
    if(halfFix && focalIsTip && sp >= 0 && sp != i && unresolvedFossils->isSA(sp) == false
       && unresolvedFossils->getAttachmentZone(sp) == mz){
        double ys = unresolvedFossils->getFossilAge(sp), zs = unresolvedFossils->getAttachAge(sp);
        if(ys < z && z < zs)
            count += 0.5;                              // spine acts as backbone: weight 1, not 1/2
    }
    return count;
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
        gammaStale.assign(nf, 1);
        prevY.assign(nf, -1.0);
        prevZ.assign(nf, -1.0);
        prevSA.assign(nf, -1);
        prevAttachmentZone.assign(nf, -1);
        prevNodeAge.assign(tree->getNumNodes(), -1.0);
        prevX0 = -1.0;
        cacheInit = true;
        if(zoneInit == false){
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

    for(int i = 0; i < nf; i++){
        double yi = unresolvedFossils->getFossilAge(i);
        double zi = unresolvedFossils->getAttachAge(i);
        int sai = unresolvedFossils->isSA(i) ? 1 : 0;
        int lzi = unresolvedFossils->getAttachmentZone(i);
        if(yi == prevY[i] && zi == prevZ[i] && sai == prevSA[i] && lzi == prevAttachmentZone[i])
            continue;
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
            for(int j = 0; j < nf; j++){
                if(gammaStale[j]) continue;
                double tj = unresolvedFossils->isSA(j) ? unresolvedFossils->getFossilAge(j) : unresolvedFossils->getAttachAge(j);
                if(pureZ ? (tj >= lo && tj < hi) : (tj > lo && tj < hi))
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

    rebuildStalkIndex();

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
        prevAttachmentZone[i] = unresolvedFossils->getAttachmentZone(i);
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
    for(ParameterDouble* l : lambda) if(l->getProposalProbability() > 0.0){ l->scaleProposed(c); nUp++; }
    for(ParameterDouble* m : mu) if(m->getProposalProbability() > 0.0){ m->scaleProposed(c); nUp++; }
    for(auto& pv : psi) for(ParameterDouble* p : pv) if(p->getProposalProbability() > 0.0){ p->scaleProposed(c); nUp++; }
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
    double m = std::exp(parameterTree->getScaleLambda() * (rng.uniformRv() - 0.5));
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
