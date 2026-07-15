#include "RelaxedClockTreeModel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include "ApproxBranchLengthLikelihood.hpp"
#include "FBDTreeModel.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "ParameterUnresolvedFossils.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Serialize.hpp"
#include "SequenceCTMCModel.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

void RelaxedClockTreeModel::buildClock(ClockModel clockModel, const double* rgeneParam, const double* sigma2Param){
    int nPart = (lik != nullptr) ? lik->getNumPartitions() : ctmc->getNumPartitions();
    std::vector<int> pgroup = UserSettings::userSettings().getClockGroups();
    if(pgroup.empty())
        pgroup = (lik != nullptr) ? lik->getPartitionGroups() : ctmc->getPartitionGroups();
    if(pgroup.empty() == false && (int)pgroup.size() != nPart)
        Msg::error("clock-partition assignment (" + std::to_string(pgroup.size()) + " entries) does not match the number of partitions (" + std::to_string(nPart) + ")");
    // CIR clock: halt — construction detached (ParameterBranchRatesCIR kept but never built)
    clock = new ParameterBranchRates(1.0, this, fbd->getTree(), nPart, pgroup, clockModel, rgeneParam, sigma2Param);
    clock->setUnresolvedFossils(fbd->getUnresolvedFossils());
    naSel.init(2);
}

void RelaxedClockTreeModel::crownInitScale(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils){
    UserSettings& us = UserSettings::userSettings();
    if(us.getConditionAgePriorSet() == false)
        return;
    Conditioning cond = us.getConditioning();
    if(cond != Conditioning::CROWN && cond != Conditioning::ORIGIN)
        return;
    double pm = Probability::priorMean(us.getConditionAgePrior(), us.getConditionAgePriorP1(), us.getConditionAgePriorP2());
    double crownAge = t->getCrown()->getTime();
    if(crownAge <= 0.0 || std::isfinite(pm) == false || pm <= 0.0)
        return;
    double target = (cond == Conditioning::ORIGIN) ? 0.9 * pm : pm;
    if(target > 0.0 && target != crownAge)
        t->scaleInternalAges(target / crownAge);

    std::map<Node*,double> floor;
    for(Fossil& f : fossils){
        Clade* clade = nullptr;
        for(Clade& c : clades)
            if(c.getName() == f.getClade()){ clade = &c; break; }
        Node* cr = (clade == nullptr || clade->getTaxa().empty()) ? t->getCrown() : t->getMRCA(clade->getTaxa());
        Node* node = (f.getAssignment() == Assignment::CROWN) ? cr : cr->getAncestor();
        if(node == nullptr)
            continue;
        std::map<Node*,double>::iterator it = floor.find(node);
        if(it == floor.end() || f.getMaxAge() > it->second)
            floor[node] = f.getMaxAge();
    }
    for(Node* n : t->getDownPassSequence()){
        if(n->getIsTip())
            continue;
        double a = n->getTime();
        std::map<Node*,double>::iterator it = floor.find(n);
        if(it != floor.end() && it->second * 1.05 > a)
            a = it->second * 1.05;
        double maxChild = 0.0;
        for(Node* c : n->getNeighbors())
            if(c != n->getAncestor() && c->getTime() > maxChild)
                maxChild = c->getTime();
        if(a <= maxChild)
            a = maxChild * 1.02 + 1e-9;
        n->setTime(a);
    }
}

double RelaxedClockTreeModel::nodeAgeJump2(void){
    std::vector<Node*> nodes = fbd->getTree()->getInternalAgeNodes();
    double s = 0.0;
    for(size_t i = 0; i < nodes.size() && i < naSnap.size(); i++){
        double d = std::log(nodes[i]->getTime()) - naSnap[i];
        s += d * d;
    }
    return s;
}

RelaxedClockTreeModel::RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& hessianFile, const std::string& mlTreeFile, int nStates, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param, unsigned int seed){
    rng.setSeed(seed);
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);
    crownInitScale(t, clades, fossils);
    fbd = new FBDTreeModel(t, clades, fossils, seed);
    ctmc = nullptr;
    std::vector<std::string> rogue;
    for(Fossil& f : fossils)
        rogue.push_back(f.getTaxon());
    lik = new ApproxBranchLengthLikelihood(hessianFile, mlTreeFile, rogue, nStates);
    buildClock(clockModel, rgeneParam, sigma2Param);
    parameters.push_back(fbd->getParameterTree());
    lastMoveType = 2;
    RandomVariable::setActiveInstance(prevRng);
}

RelaxedClockTreeModel::RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& sequenceFile, const std::string& partitionFile, int nStates, int numCats, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param, unsigned int seed){
    rng.setSeed(seed);
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);
    crownInitScale(t, clades, fossils);
    fbd = new FBDTreeModel(t, clades, fossils, seed);
    lik = nullptr;
    ctmc = new SequenceCTMCModel(this, sequenceFile, partitionFile, nStates, numCats);

    buildClock(clockModel, rgeneParam, sigma2Param);
    ctmc->buildParameters();

    parameters.push_back(fbd->getParameterTree());
    lastMoveType = 2;
    RandomVariable::setActiveInstance(prevRng);
}

void RelaxedClockTreeModel::invalidateLikelihoodCache(void){
    if(ctmc != nullptr)
        ctmc->invalidateCache();
}

double RelaxedClockTreeModel::lnLikelihood(void){
    if(ctmc == nullptr){
        return lik->computeLnL(fbd->getTree(), clock->getAbsoluteRates());
    }
    return ctmc->computeLnL(fbd->getTree(), clock->getAbsoluteRates(), clock->getBranchMGF());
}

void RelaxedClockTreeModel::invalidatePriorCache(void){
    fbd->invalidateGammaCache();
}

double RelaxedClockTreeModel::lnPriorProbability(void){
    double lnp = fbd->lnLikelihood() + fbd->lnPriorProbability() + clock->lnProbability();
    if(ctmc != nullptr)
        lnp += ctmc->lnPrior();
    return lnp;
}

double RelaxedClockTreeModel::nodeAgeSweep(void){
    RandomVariable& r = RandomVariable::randomVariableInstance();
    fbd->setupNodeAgeFloors();
    std::vector<Node*> nodes = fbd->getTree()->getInternalAgeNodes();
    for(int i = (int)nodes.size() - 1; i > 0; i--){
        int j = (int)(r.uniformRv() * (i + 1));
        std::swap(nodes[i], nodes[j]);
    }
    double curL = lnLikelihood();
    double curP = lnPriorProbability();
    for(Node* n : nodes){
        double ratio = fbd->getTree()->updateNodeAgeOnNode(n);
        double newL = lnLikelihood();
        double newP = lnPriorProbability();
        if(std::log(r.uniformRv()) < (newL - curL) + (newP - curP) + ratio){
            curL = newL;
            curP = newP;
            fbd->getParameterTree()->updateForAcceptance();
        }else{
            fbd->getParameterTree()->updateForRejection();
        }
    }
    return std::numeric_limits<double>::infinity(); // each node already MH-accepted above; force outer accept
}

double RelaxedClockTreeModel::update(void){
    RandomVariable& r = RandomVariable::randomVariableInstance();
    if(ctmc != nullptr && r.uniformRv() < 0.05){ lastMoveType = 6; return ctmc->update(); }
    double u = r.uniformRv();
    if(u < 0.20){ lastMoveType = 0; return clock->update(); }
    if(u < 0.516){
        fbd->setupNodeAgeFloors();
        if(nInternalAge == 0) nInternalAge = (int)fbd->getAgeLogNodes().size();
        double pCrown = (nInternalAge > 0) ? 1.0 / nInternalAge : 0.0;
        if(r.uniformRv() < pCrown){ lastMoveType = 8; return clock->simpleDistanceMove(); }
        if(static_cast<ParameterBranchRates*>(clock)->getClockModel() == ClockModel::UCLN && r.uniformRv() < pCrown){ lastMoveType = 10; return clock->smallPulleyMove(); }
        naSnap.clear();
        std::vector<Node*> nodes = fbd->getTree()->getInternalAgeNodes();
        for(Node* n : nodes)
            naSnap.push_back(std::log(n->getTime()));
        naOp = naSel.pick(r);
        if(naOp == 0){ lastMoveType = 1; return clock->constantDistanceMove(); }
        lastMoveType = 7;
        double h = nodeAgeSweep();
        naSel.record(1, nodeAgeJump2(), (double)naSnap.size());
        return h;
    }
    if(u < 0.572){
        lastMoveType = 4;
        return clock->rateAgeSubtreeMove();
    }
    if(u < 0.651){
        lastMoveType = 3;
        ageScaleAtt++;
        if(ageScaleAtt % 50 == 0){
            double ar = (double)ageScaleAcc / ageScaleAtt;
            double gain = 1.0 / std::sqrt((double)(ageScaleAtt / 50));
            ageScaleStep *= std::exp(gain * (ar - 0.3));
            if(ageScaleStep < 1e-4) ageScaleStep = 1e-4;
            if(ageScaleStep > 20.0)  ageScaleStep = 20.0;
        }
        double cc = std::exp(ageScaleStep * (r.uniformRv() - 0.5));
        int kAge = fbd->getTree()->scaleInternalAges(cc);
        clock->scaleAll(1.0 / cc);
        int nRate = clock->getNumClockPartitions() * (1 + clock->getNumBranchNodes());
        return ((double)kAge - (double)nRate) * std::log(cc);
    }
    if(fbd->getUnresolvedFossils() != nullptr && u < 0.751){
        lastMoveType = 9;
        return fbd->fossilSweep();
    }
    lastMoveType = 2; return fbd->update();
}

void RelaxedClockTreeModel::updateForAcceptance(void){
    if(lastMoveType == 7 || lastMoveType == 9) return;
    if(lastMoveType == 6)
        ctmc->updateForAcceptance();
    else if(lastMoveType == 0)
        clock->updateForAcceptance();
    else if(lastMoveType == 1){
        clock->updateForAcceptance();
        fbd->getParameterTree()->updateForAcceptance();
        naSel.record(0, nodeAgeJump2(), 1.0);
    }else if(lastMoveType == 3){
        ageScaleAcc++;
        clock->commitAll();
        fbd->getParameterTree()->updateForAcceptance();
    }else if(lastMoveType == 4){
        clock->updateForAcceptance();
        fbd->getParameterTree()->updateForAcceptance();
        if(fbd->getUnresolvedFossils() != nullptr)
            fbd->getUnresolvedFossils()->updateForAcceptance();
    }else if(lastMoveType == 8){
        clock->updateForAcceptance();
        fbd->getParameterTree()->updateForAcceptance();
    }else if(lastMoveType == 10){
        clock->updateForAcceptance();
    }else
        fbd->updateForAcceptance();
}

void RelaxedClockTreeModel::updateForRejection(void){
    if(lastMoveType == 7 || lastMoveType == 9) return;
    if(lastMoveType == 6)
        ctmc->updateForRejection();
    else if(lastMoveType == 0)
        clock->updateForRejection();
    else if(lastMoveType == 1){
        clock->updateForRejection();
        fbd->getParameterTree()->updateForRejection();
        naSel.record(0, 0.0, 1.0);
    }else if(lastMoveType == 3){
        clock->restoreAll();
        fbd->getParameterTree()->updateForRejection();
    }else if(lastMoveType == 4){
        clock->updateForRejection();
        fbd->getParameterTree()->updateForRejection();
        if(fbd->getUnresolvedFossils() != nullptr)
            fbd->getUnresolvedFossils()->updateForRejection();
    }else if(lastMoveType == 8){
        clock->updateForRejection();
        fbd->getParameterTree()->updateForRejection();
    }else if(lastMoveType == 10){
        clock->updateForRejection();
    }else
        fbd->updateForRejection();
}

void RelaxedClockTreeModel::writeState(std::ostream& os){
    fbd->getRng()->writeState(os);
    fbd->writeState(os);
    clock->writeState(os);
    if(ctmc != nullptr)
        ctmc->writeState(os);
    os << ageScaleStep << ' ' << ageScaleAtt << ' ' << ageScaleAcc << '\n';
    naSel.writeState(os);
}

void RelaxedClockTreeModel::readState(std::istream& is){
    fbd->getRng()->readState(is);
    fbd->readState(is);
    clock->readState(is);
    if(ctmc != nullptr)
        ctmc->readState(is);
    is >> ageScaleStep >> ageScaleAtt >> ageScaleAcc;
    naSel.readState(is);
}

void RelaxedClockTreeModel::collectNodeAges(std::vector<std::string>* names, std::vector<double>* vals){
    if(fbd->hasOrigin()){
        if(names) names->push_back("x0");
        if(vals)  vals->push_back(fbd->getOriginAgeValue());
    }
    std::vector<Node*> bb = fbd->getAgeLogNodes();
    for(size_t i = 0; i < bb.size(); i++){
        if(names) names->push_back("x" + std::to_string(i + 1));
        if(vals)  vals->push_back(bb[i]->getTime());
    }
}

static bool isNodeAgeName(const std::string& s){
    return s.size() >= 2 && s[0] == 'x' && s[1] >= '0' && s[1] <= '9';
}

std::vector<std::string> RelaxedClockTreeModel::getParameterNames(void){
    std::vector<std::string> n;
    collectNodeAges(&n, nullptr);
    for(const std::string& s : fbd->getParameterNames())
        if(s != "originAge" && isNodeAgeName(s) == false)
            n.push_back(s);
    for(int p = 0; p < clock->getNumClockPartitions(); p++){
        std::string suf = (clock->getNumClockPartitions() > 1) ? std::to_string(p) : "";
        n.push_back("clockMean" + suf);
        n.push_back("clockSigma2" + suf);
    }
    if(ctmc != nullptr)
        ctmc->appendParameterNames(n);
    return n;
}

std::vector<double> RelaxedClockTreeModel::getParameterString(void){
    std::vector<double> v;
    collectNodeAges(nullptr, &v);
    std::vector<std::string> fbdN = fbd->getParameterNames();
    std::vector<double> fbdV = fbd->getParameterString();
    for(size_t i = 0; i < fbdV.size(); i++)
        if(i >= fbdN.size() || (fbdN[i] != "originAge" && isNodeAgeName(fbdN[i]) == false))
            v.push_back(fbdV[i]);
    for(int p = 0; p < clock->getNumClockPartitions(); p++){
        v.push_back(clock->getClockPartitionRate(p));
        v.push_back(clock->getClockPartitionSigma2(p));
    }
    if(ctmc != nullptr)
        ctmc->appendParameterValues(v);
    return v;
}

bool RelaxedClockTreeModel::treeIncludesFossils(void){ return fbd->treeIncludesFossils(); }

std::vector<std::string> RelaxedClockTreeModel::getLatentNames(void){
    std::vector<std::string> n = fbd->getLatentNames();
    Tree* tr = fbd->getTree();
    std::map<int,std::string> lab;
    std::vector<Node*> bb = fbd->getAgeLogNodes();
    for(size_t i = 0; i < bb.size(); i++)
        lab[bb[i]->getOffset()] = "x" + std::to_string(i + 1);
    for(int p = 0; p < clock->getNumClockPartitions(); p++){
        std::string psuf = (clock->getNumClockPartitions() > 1) ? ("_p" + std::to_string(p)) : "";
        for(int i = 0; i < clock->getNumBranchNodes(); i++){
            int off = clock->getBranchNodeOffset(i);
            std::map<int,std::string>::iterator it = lab.find(off);
            n.push_back("rate_" + (it != lab.end() ? it->second : tr->getNodeByOffset(off)->getName()) + psuf);
        }
    }
    return n;
}

std::vector<double> RelaxedClockTreeModel::getLatentString(void){
    std::vector<double> v = fbd->getLatentString();
    std::vector<std::vector<double>> ar = clock->getAbsoluteRates();
    for(int p = 0; p < clock->getNumClockPartitions(); p++)
        for(int i = 0; i < clock->getNumBranchNodes(); i++)
            v.push_back(ar[p][clock->getBranchNodeOffset(i)]);
    return v;
}

void RelaxedClockTreeModel::print(void){
    fbd->print();
    clock->print();
}
