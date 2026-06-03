#include "Node.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "FBDTreeModel.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"
#include "Probability.hpp"
#include "Msg.hpp"

#include <iostream>
#include <cmath>
#include <limits>

FBDTreeModel::FBDTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, unsigned int seed) :
    PhylogeneticModel(),
    c1(0.0),
    c2(0.0){

    lastWasJointScale = false;
    cacheInit = false;
    rng.setSeed(seed);
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);

    parameterTree = new ParameterTree(1.0, this);
    isFBD = (UserSettings::userSettings().getModel() == Model::FBD);
    if(isFBD){
        Tree* working = new Tree(*t);
        resolveFossils(working, clades, fossils);
        parameterTree->setTree(working);
        delete working;
    }else{
        parameterTree->setTree(t);
    }
    parameters.push_back(parameterTree);
    
    //instantiate FBD model parameters
    lambda = new ParameterDouble(1.0, this, "lambda", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(lambda);
    mu = new ParameterDouble(1.0, this, "mu", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(mu);
    psi = new ParameterDouble(1.0, this, "psi", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(psi);
    Probability::PriorSpec lp = UserSettings::userSettings().getLambdaPrior(); if(lp.set) lambda->setPrior(lp.family, lp.p1, lp.p2);
    Probability::PriorSpec mp = UserSettings::userSettings().getMuPrior();     if(mp.set) mu->setPrior(mp.family, mp.p1, mp.p2);
    Probability::PriorSpec pp = UserSettings::userSettings().getPsiPrior();    if(pp.set) psi->setPrior(pp.family, pp.p1, pp.p2);
    rho = UserSettings::userSettings().getRho();

    originAge = nullptr;
    if(UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN){
        double x0init = parameterTree->getTree()->getRoot()->getTime();
        for(Fossil& f : fossils)
            if(f.getMaxAge() > x0init)
                x0init = f.getMaxAge();
        originAge = new ParameterDouble(1.0, this, "originAge", 0.0, std::numeric_limits<double>::max());
        originAge->setValue(x0init * 1.05);
        parameters.push_back(originAge);
        UserSettings& us = UserSettings::userSettings();
        if(us.getConditionAgePriorSet())
            originAge->setPrior(us.getConditionAgePrior(), us.getConditionAgePriorP1(), us.getConditionAgePriorP2());
        else
            originAge->setPrior(Probability::PriorFamily::IMPROPER, 0.0, 0.0);
    }

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
            fossilOrigin.push_back(wt->getNodeByOffset(clade->getOrigin()->getOffset()));
        }
    }else{
        unresolvedFossils = new ParameterUnresolvedFossils(1.0, this, parameterTree->getTree(), clades, fossils, originAge);
        parameters.push_back(unresolvedFossils);
    }

    if(isFBD){
        parameterTree->setProposalProbability(78.0);
        lambda->setProposalProbability(15.0);
        mu->setProposalProbability(15.0);
        psi->setProposalProbability(15.0);
    }

    //normalize proposal probabilities
    double sum = 0.0;
    for(Parameter* p : parameters)
        sum += p->getProposalProbability();
    for(Parameter* p : parameters)
        p->setProposalProbability(p->getProposalProbability() / sum);

    RandomVariable::setActiveInstance(prevRng);
}

double FBDTreeModel::calculateFBDProbability(void){
    Tree* tree = parameterTree->getTree();
    
    int numInternalNodes = tree->getNumNodes() - tree->getNumTaxa();
    double rootAge = tree->getRoot()->getTime();
    std::vector<Node*> dpseq = tree->getDownPassSequence();
    
    lambdaVal = lambda->getValue();
    muVal = mu->getValue();
    rhoVal = rho;
    psiVal = psi->getValue();
    double log4LambdaRho = std::log(4*lambdaVal*rhoVal);
    
    calculateC1();
    calculateC2();

    if(isFBD)
        return calculateResolvedFBD();

    double qRoot = calculateQt(rootAge);

    //FBD probability in log-landia
    bool useOrigin = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);
    double fbdProb = 0.0;

    //term 1: conditioning
    if(useOrigin){
        double x0 = originAge->getValue();
        if(x0 < rootAge)
            return -INFINITY;
        fbdProb -= std::log(lambdaVal);
        fbdProb -= calculateLnSurvival(x0);
        fbdProb += log4LambdaRho;
        fbdProb -= std::log(calculateQt(x0));
    }else{
        fbdProb -= 2 * (std::log(lambdaVal) + calculateLnSurvival(rootAge));
        fbdProb += log4LambdaRho;
        fbdProb -= std::log(qRoot);
    }

    //term 2: main body
    fbdProb += numInternalNodes * log4LambdaRho;
    for(Node* n : dpseq)
        if(n->getIsTip() == false)
            fbdProb -= std::log(calculateQt(n->getTime()));

    //term 3: fossil attachment
    int numFossils = unresolvedFossils->getNumFossils();
    updateGammaCache();
    for(int i = 0; i < numFossils; i++){
        if(unresolvedFossils->isSA(i)){
            fbdProb += std::log(psiVal) + cachedGammaLn[i];
            continue;
        }
        fbdProb += fossilPqLn(unresolvedFossils->getFossilAge(i), unresolvedFossils->getAttachAge(i)) + cachedGammaLn[i];
    }
    return fbdProb;
}

double FBDTreeModel::lnD(double t){
    return (t <= 0.0) ? 0.0 : std::log(4.0) - std::log(calculateQt(t));
}

double FBDTreeModel::calculateResolvedFBD(void){
    Tree* tree = parameterTree->getTree();
    Node* root = tree->getRoot();
    double rootAge = root->getTime();
    bool useOrigin = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);
    double startAge = useOrigin ? originAge->getValue() : rootAge;
    if(useOrigin && startAge < rootAge)
        return -INFINITY;

    double logTwoLambda = std::log(2.0 * lambdaVal);
    int numInitialLineages = useOrigin ? 1 : 2;
    double lnP = -numInitialLineages * calculateLnSurvival(startAge);
    if(useOrigin)
        lnP += lnD(startAge) - lnD(rootAge);

    for(Node* n : tree->getDownPassSequence()){
        if(n != root)
            lnP += lnD(n->getAncestor()->getTime()) - lnD(n->getTime());
        if(n->getIsTip()){
            if(n->getIsFossil() == false)
                lnP += std::log(rhoVal);
            else if(n->getAncestor()->getTime() == n->getTime())
                lnP += std::log(psiVal);
            else
                lnP += std::log(psiVal) + std::log(calculatePo(n->getTime()));
        }
        else if(useOrigin || n != root){
            bool fakeSplit = false;
            for(Node* c : n->getNeighbors())
                if(c != n->getAncestor() && c->getIsTip() && c->getIsFossil() && c->getTime() == n->getTime()){ fakeSplit = true; break; }
            if(fakeSplit == false)
                lnP += logTwoLambda;
        }
    }
    return lnP;
}

static bool nodeInSubtree(Node* node, Node* subtreeRoot){
    for(Node* p = node; ; p = p->getAncestor()){
        if(p == subtreeRoot)
            return true;
        if(p->getAncestor() == p) // reached the root self-loop
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

double FBDTreeModel::computeGamma(double z, int i){
    Tree* tree = parameterTree->getTree();
    Node* crown = unresolvedFossils->getCrownNode(i);
    bool total = (unresolvedFossils->getIsCrown(i) == false);
    double count = 0.0;
    for(Node* n : tree->getDownPassSequence()){
        if(n == tree->getRoot())
            continue;
        Node* anc = n->getAncestor();
        if(n->getTime() < z && z < anc->getTime()){
            bool inZone = nodeInSubtree(anc, crown);
            if(inZone == false && total && n == crown)
                inZone = true;
            if(inZone)
                count++;
        }
    }
    if(total && crown == tree->getRoot() && originAge != nullptr){
        double x0 = originAge->getValue();
        if(tree->getRoot()->getTime() < z && z < x0)
            count++;
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
        if(nodeInSubtree(unresolvedFossils->getCrownNode(j), crown) == false)
            continue;
        if(total == false && zj > crown->getTime())
            continue;
        double w = 1.0;
        if(SymmetryCorrection && focalIsTip){
            Node* crownJ = unresolvedFossils->getCrownNode(j);
            bool jTotal = (unresolvedFossils->getIsCrown(j) == false);
            bool reciprocal = nodeInSubtree(crown, crownJ)
                              && (jTotal || z <= crownJ->getTime());
            if(reciprocal)
                w = 0.5;
        }
        count += w;
    }
    return count;
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
        Node* origin = t->getNodeByOffset(clade->getOrigin()->getOffset());
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
            if(n->getAncestor() != treeRoot && t->isFakeSplit(n->getAncestor()) == false)
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
    double oldRange = std::min(ceilingS, hostParent->getTime()) - std::max(rAge, hostChild->getTime());

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
    double newRange = his[k] - los[k];
    double z = los[k] + rng.uniformRv() * newRange;

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
    if(candidates.empty())
        return 0.0;

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

void FBDTreeModel::updateGammaCache(void){
    Tree* tree = parameterTree->getTree();
    int nf = unresolvedFossils->getNumFossils();
    std::vector<Node*>& dpseq = tree->getDownPassSequence();

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
            Node* root = tree->getRoot();
            for(int i = 0; i < nf; i++)
                if(unresolvedFossils->getCrownNode(i) == root && unresolvedFossils->getIsCrown(i) == false)
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

void FBDTreeModel::calculateC1(void){
    c1 =    std::abs(
                std::sqrt(
                    std::pow(lambdaVal - muVal - psiVal, 2) +
                    4*lambdaVal * psiVal
                )
            );
}

void FBDTreeModel::calculateC2(void){
    c2 = (-lambdaVal + muVal + 2*lambdaVal * rhoVal + psiVal);
    c2 /= c1;
}

double FBDTreeModel::calculateQt(double t){
    double tmp = 2 * (1 - std::pow(c2, 2));
    tmp += std::exp(-c1 * t) * std::pow(1-c2, 2);
    tmp += std::exp(c1 * t) * std::pow(1+c2, 2);
    return tmp;
}

double FBDTreeModel::calculatePo(double t){
    double tmp = -lambdaVal + muVal + psiVal;
    tmp += c1 * (std::exp(-c1 * t) * (1 - c2) - (1+c2) ) / ( std::exp(-c1 * t) * (1 - c2) + (1+c2)  );
    tmp /= 2*lambdaVal;
    return 1 + tmp;
}

double FBDTreeModel::calculatePoHat(double t){
    double tmp = rhoVal * (lambdaVal - muVal);
    tmp /= (lambdaVal * rhoVal + (lambdaVal*(1-rhoVal) - muVal)*std::exp(-1 * (lambdaVal - muVal) * t) );
    return 1 - tmp;
}

double FBDTreeModel::calculateLnSurvival(double t){
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

double FBDTreeModel::fossilPqLn(double y, double z){
    return std::log(psiVal) + std::log(2*lambdaVal) + std::log(calculatePo(y)) + std::log(calculateQt(y)) - std::log(calculateQt(z));
}

std::vector<std::string> FBDTreeModel::getParameterNames(void){
    std::vector<std::string> names;
    for(Parameter* p : parameters)
        if( p != parameterTree){ //by convention, exclude parameter tree from these getters
            ParameterUnresolvedFossils* uf = dynamic_cast<ParameterUnresolvedFossils*>(p);
            names.push_back(uf != nullptr ? "nSA" : p->getName());
        }
    return names;
}

std::vector<double> FBDTreeModel::getParameterString(void){
    std::vector<double> vals;
    for(Parameter* p : parameters)
        if( p != parameterTree){ //by convention, exclude parameter tree from these getters
            ParameterDouble* pd = dynamic_cast<ParameterDouble*>(p);
            if(pd != nullptr){
                vals.push_back(pd->getValue());
                continue;
            }
            ParameterUnresolvedFossils* uf = dynamic_cast<ParameterUnresolvedFossils*>(p);
            if(uf != nullptr)
                vals.push_back((double)uf->getNumSampledAncestors());
        }

    return vals;
}

double FBDTreeModel::lnLikelihood(void){
    return calculateFBDProbability();
}

double FBDTreeModel::lnPriorProbability(void){
    double lnP = 0.0;
    
    for(Parameter* p : parameters)
        if( p != parameterTree)
            lnP += p->lnProbability();

    UserSettings& us = UserSettings::userSettings();
    if(us.getConditionAgePriorSet() && us.getConditioning() == Conditioning::CROWN){
        double rootAge = parameterTree->getTree()->getRoot()->getTime();
        lnP += Probability::priorLnPdf(us.getConditionAgePrior(), us.getConditionAgePriorP1(), us.getConditionAgePriorP2(), rootAge, 0.0, std::numeric_limits<double>::max());
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

double FBDTreeModel::update(void){
    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&rng);

    double u = rng.uniformRv();
    double sum = 0.0;
    for(Parameter* p  : parameters){
        sum += p->getProposalProbability();
        if(u < sum){
            updatedParameter = p;
            break;
        }
    }
    lastWasJointScale = false;
    if(updatedParameter == parameterTree && isFBD == false){
        Tree* t = parameterTree->getTree();
        int numSlideable = 0;
        for(Node* n : t->getDownPassSequence())
            if(n != t->getRoot() && n->getIsTip() == false)
                numSlideable++;
        double fixedWeight = 3.0;
        double slideAndRoot = numSlideable + fixedWeight;
        double uMove = rng.uniformRv() * (slideAndRoot + 2.0 * fixedWeight);
        if(uMove >= slideAndRoot){
            lastWasJointScale = true;
            double r = (uMove < slideAndRoot + fixedWeight) ? doJointScale() : doSubtreeScale();
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
        double wAge = 40.0;
        double uMove = rng.uniformRv() * (wNE + wWB + wWE + wTreeScale + wSA + wAge);
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
    }
    double ratio = updatedParameter->update();

    RandomVariable::setActiveInstance(prevRng);
    return ratio;
}

void FBDTreeModel::updateForAcceptance(void){
    if(lastWasJointScale){
        parameterTree->updateForAcceptance();
        unresolvedFossils->updateForAcceptance();
    }else{
        updatedParameter->updateForAcceptance();
    }
}

void FBDTreeModel::updateForRejection(void){
    if(lastWasJointScale){
        parameterTree->updateForRejection();
        unresolvedFossils->updateForRejection();
    }else{
        updatedParameter->updateForRejection();
    }
}
