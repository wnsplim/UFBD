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
    int nLoci = (lik != nullptr) ? lik->getNumPartitions() : ctmc->getNumPartitions();
    // CIR clock: halt — construction detached (ParameterBranchRatesCIR kept but never built)
    clock = new ParameterBranchRates(1.0, this, fbd->getTree(), nLoci, clockModel, rgeneParam, sigma2Param);
    clock->setUnresolvedFossils(fbd->getUnresolvedFossils());
    naSel.init(2);
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
    fbd = new FBDTreeModel(t, clades, fossils, seed);
    ctmc = nullptr;
    std::vector<std::string> rogue;
    for(Fossil& f : fossils)
        rogue.push_back(f.getTaxon());
    lik = new ApproxBranchLengthLikelihood(hessianFile, mlTreeFile, rogue, 0, nStates);
    buildClock(clockModel, rgeneParam, sigma2Param);
    parameters.push_back(fbd->getParameterTree());
    lastMoveType = 2;
    RandomVariable::setActiveInstance(prevRng);
}

RelaxedClockTreeModel::RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& sequenceFile, const std::string& partitionFile, int nStates, int numCats, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param, unsigned int seed){
    rng.setSeed(seed);
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);
    fbd = new FBDTreeModel(t, clades, fossils, seed);
    lik = nullptr;
    ctmc = new SequenceCTMCModel(this, sequenceFile, partitionFile, nStates, numCats);

    buildClock(clockModel, rgeneParam, sigma2Param);
    ctmc->buildParameters();

    parameters.push_back(fbd->getParameterTree());
    lastMoveType = 2;
    RandomVariable::setActiveInstance(prevRng);
}

double RelaxedClockTreeModel::lnLikelihood(void){
    if(ctmc == nullptr)
        return lik->computeLnL(fbd->getTree(), clock->getAbsoluteRates());
    return ctmc->computeLnL(fbd->getTree(), clock->getAbsoluteRates(), clock->getBranchMGF());
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
    static const bool naSelOff = [](){ const char* e = std::getenv("FBD_NASEL"); return e != nullptr && e[0] == 'o'; }();
    RandomVariable& r = RandomVariable::randomVariableInstance();
    if(ctmc != nullptr && r.uniformRv() < 0.25){ lastMoveType = 6; return ctmc->update(); }
    double u = r.uniformRv();
    if(u < 0.25){ lastMoveType = 0; return clock->update(); }
    if(u < 0.65){
        naSnap.clear();
        std::vector<Node*> nodes = fbd->getTree()->getInternalAgeNodes();
        for(Node* n : nodes)
            naSnap.push_back(std::log(n->getTime()));
        naOp = naSelOff ? (u < 0.50 ? 0 : 1) : naSel.pick(r);
        if(naOp == 0){ lastMoveType = 1; return clock->constantDistanceMove(); }
        lastMoveType = 7;
        double h = nodeAgeSweep();
        naSel.record(1, nodeAgeJump2(), (double)naSnap.size());
        return h;
    }
    if(u < 0.72){
        lastMoveType = 4;
        return clock->rateAgeSubtreeMove();
    }
    if(u < 0.82){
        lastMoveType = 3;
        ageScaleAtt++;
        if(ageScaleAtt % 50 == 0){
            double ar = (double)ageScaleAcc / ageScaleAtt;
            double gain = 1.0 / std::sqrt((double)(ageScaleAtt / 50));
            ageScaleStep *= std::exp(gain * (ar - 0.44));
            if(ageScaleStep < 1e-4) ageScaleStep = 1e-4;
            if(ageScaleStep > 20.0)  ageScaleStep = 20.0;
        }
        double cc = std::exp(ageScaleStep * (r.uniformRv() - 0.5));
        int kAge = fbd->getTree()->scaleInternalAges(cc);
        clock->scaleAll(1.0 / cc);
        int nRate = clock->getNumLoci() * (1 + clock->getNumBranchNodes());
        return ((double)kAge - (double)nRate) * std::log(cc);
    }
    lastMoveType = 2; return fbd->update();
}

void RelaxedClockTreeModel::updateForAcceptance(void){
    if(lastMoveType == 7) return;
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
    }else
        fbd->updateForAcceptance();
}

void RelaxedClockTreeModel::updateForRejection(void){
    if(lastMoveType == 7) return;
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
    Node* crown = fbd->getTree()->getCrown();
    std::vector<Node*> stack(1, crown);
    int k = 1;
    while(stack.empty() == false){
        Node* n = stack.back();
        stack.pop_back();
        if(n->getIsTip())
            continue;
        if(names) names->push_back("x" + std::to_string(k));
        if(vals)  vals->push_back(n->getTime());
        k++;
        std::vector<Node*> kids;
        for(Node* c : n->getNeighbors())
            if(c != n->getAncestor())
                kids.push_back(c);
        std::sort(kids.begin(), kids.end(), [](Node* a, Node* b){ return a->getOffset() > b->getOffset(); });
        for(Node* c : kids)
            stack.push_back(c);
    }
}

std::vector<std::string> RelaxedClockTreeModel::getParameterNames(void){
    std::vector<std::string> n;
    collectNodeAges(&n, nullptr);
    for(const std::string& s : fbd->getParameterNames())
        if(s != "originAge")
            n.push_back(s);
    for(int p = 0; p < clock->getNumLoci(); p++){
        std::string suf = (clock->getNumLoci() > 1) ? std::to_string(p) : "";
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
        if(i >= fbdN.size() || fbdN[i] != "originAge")
            v.push_back(fbdV[i]);
    for(int p = 0; p < clock->getNumLoci(); p++){
        v.push_back(clock->getLocusRate(p));
        v.push_back(clock->getLocusSigma2(p));
    }
    if(ctmc != nullptr)
        ctmc->appendParameterValues(v);
    return v;
}

void RelaxedClockTreeModel::print(void){
    fbd->print();
    clock->print();
}
