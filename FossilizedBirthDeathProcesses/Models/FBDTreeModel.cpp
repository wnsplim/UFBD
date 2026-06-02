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
    if(isFBD)
        resolveFossils(t, clades, fossils);
    parameterTree->setTree(t);
    parameters.push_back(parameterTree);
    
    //instantiate FBD model parameters
    lambda = new ParameterDouble(1.0, this, "lambda", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(lambda);
    mu = new ParameterDouble(1.0, this, "mu", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(mu);
    psi = new ParameterDouble(1.0, this, "psi", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(psi);
    rho = UserSettings::userSettings().getRho();

    originAge = nullptr;
    if(UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN && isFBD == false){
        double x0init = parameterTree->getTree()->getRoot()->getTime();
        for(Fossil& f : fossils)
            if(f.getMaxAge() > x0init)
                x0init = f.getMaxAge();
        originAge = new ParameterDouble(1.0, this, "originAge", 0.0, std::numeric_limits<double>::max());
        originAge->setValue(x0init * 1.05);
        parameters.push_back(originAge);
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
        }
    }else{
        unresolvedFossils = new ParameterUnresolvedFossils(1.0, this, parameterTree->getTree(), clades, fossils, originAge);
        parameters.push_back(unresolvedFossils);
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
    double qRoot = calculateQt(rootAge);
    
    //FBD probability in log-landia
    bool useOrigin = (UserSettings::userSettings().getConditioning() == Conditioning::ORIGIN);
    double fbdProb = 0.0;

    //term 1: conditioning
    if(useOrigin){
        double x0 = (originAge != nullptr) ? originAge->getValue() : rootAge;
        if(x0 < rootAge)
            return -INFINITY;
        fbdProb -= std::log(lambdaVal);
        fbdProb -= calculateLnSurvival(x0);
        if(originAge != nullptr){
            fbdProb += log4LambdaRho;
            fbdProb -= std::log(calculateQt(x0));
        }
    }else{
        fbdProb -= 2 * (std::log(lambdaVal) + calculateLnSurvival(rootAge));
        fbdProb += log4LambdaRho;
        fbdProb -= std::log(qRoot);
    }

    //term 2: main body
    if(isFBD){//this chunk is for resolved (vanilla) FBD
        for(Node* n : dpseq){
            if(n->getIsTip())
                continue;
            bool isSplit = false;
            for(Node* nb : n->getNeighbors())
                if(nb != n->getAncestor() && nb->getIsTip() && nb->getIsFossil()){
                    isSplit = true;
                    break;
                }
            if(isSplit)
                continue;
            fbdProb += log4LambdaRho;
            fbdProb -= std::log(calculateQt(n->getTime()));
        }
        for(Node* n : dpseq)
            if(n->getIsTip() && n->getIsFossil())
                fbdProb += fossilPqLn(n->getTime(), n->getAncestor()->getTime());
        return fbdProb;
    }

    fbdProb += numInternalNodes * log4LambdaRho;
    for(Node* n : dpseq)
        if(n->getIsTip() == false)
            fbdProb -= std::log(calculateQt(n->getTime()));

    //term 3: fossil attachment
    int numFossils = unresolvedFossils->getNumFossils();
    updateGammaCache();
    for(int i = 0; i < numFossils; i++){
        if(unresolvedFossils->isSampledAncestor(i)){
            fbdProb += std::log(psiVal);
            continue;
        }
        fbdProb += fossilPqLn(unresolvedFossils->getFossilAge(i), unresolvedFossils->getAttachAge(i)) + cachedGammaLn[i];
    }
    return fbdProb;
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
    bool focalIsTip = (unresolvedFossils->isSampledAncestor(i) == false);
    int numFossils = unresolvedFossils->getNumFossils();
    for(int j = 0; j < numFossils; j++){
        if(j == i)
            continue;
        if(unresolvedFossils->isSampledAncestor(j))
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
        Node* origin = clade->getOrigin();
        bool isCrown = (f.getAssignment() == Assignment::CROWN);
        if(crown->getAncestor() == crown)
            isCrown = true;
        double y = 0.5 * (f.getMinAge() + f.getMaxAge());
        double ceiling = isCrown ? crown->getTime() : origin->getTime();

        std::vector<Node*> hostChildren;
        std::vector<double> hostLo;
        std::vector<double> hostHi;
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
                hostChildren.push_back(n);
                hostLo.push_back(lo);
                hostHi.push_back(hi);
            }
        }
        int k = (int)(rv.uniformRv() * hostChildren.size());
        double z = hostLo[k] + rv.uniformRv() * (hostHi[k] - hostLo[k]);
        t->insertFossilTip(hostChildren[k], f.getTaxon(), y, z);

        fossilName.push_back(f.getTaxon());
        fossilIsCrown.push_back(isCrown);
        fossilY.push_back(y);
    }
}

void FBDTreeModel::computeAgeFloors(std::map<Node*,double>& floors){
    int numFossils = unresolvedFossils->getNumFossils();
    for(int i = 0; i < numFossils; i++){
        if(unresolvedFossils->isSampledAncestor(i))
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
        if(unresolvedFossils->isSampledAncestor(i))
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
        prevSa.assign(nf, -1);
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
        int sai = unresolvedFossils->isSampledAncestor(i) ? 1 : 0;
        if(yi == prevY[i] && zi == prevZ[i] && sai == prevSa[i])
            continue;
        gammaStale[i] = 1;
        bool wasTerm = (prevSa[i] == 0);
        bool isTerm = (sai == 0);
        if((wasTerm || isTerm) && prevY[i] >= 0.0){
            double lo = std::min(prevY[i], yi);
            double hi = -INFINITY;
            if(wasTerm) hi = std::max(hi, prevZ[i]);
            if(isTerm)  hi = std::max(hi, zi);
            for(int j = 0; j < nf; j++){
                if(gammaStale[j]) continue;
                double tj = unresolvedFossils->isSampledAncestor(j) ? unresolvedFossils->getFossilAge(j) : unresolvedFossils->getAttachAge(j);
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
            double ti = unresolvedFossils->isSampledAncestor(i) ? unresolvedFossils->getFossilAge(i) : unresolvedFossils->getAttachAge(i);
            for(std::pair<double,double>& iv : changedIntervals)
                if(ti > iv.first && ti < iv.second){ gammaStale[i] = 1; break; }
        }
    }

    for(int i = 0; i < nf; i++){
        if(gammaStale[i] == 0) continue;
        if(unresolvedFossils->isSampledAncestor(i) == false){
            double g = computeGamma(unresolvedFossils->getAttachAge(i), i);
            cachedGammaLn[i] = (g > 0.0) ? std::log(g) : -INFINITY;
        }
        gammaStale[i] = 0;
    }

    for(int i = 0; i < nf; i++){
        prevY[i] = unresolvedFossils->getFossilAge(i);
        prevZ[i] = unresolvedFossils->getAttachAge(i);
        prevSa[i] = unresolvedFossils->isSampledAncestor(i) ? 1 : 0;
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
        if( p != parameterTree) //by convention, exclude parameter tree from these getters
            names.push_back(p->getName());
    return names;
}

std::vector<double> FBDTreeModel::getParameterString(void){
    std::vector<double> vals;
    for(Parameter* p : parameters)
        if( p != parameterTree){ //by convention, exclude parameter tree from these getters
            ParameterDouble* pd = dynamic_cast<ParameterDouble*>(p);
            if(pd != nullptr)
                vals.push_back(pd->getValue());
        }

    return vals;
}

double FBDTreeModel::lnLikelihood(void){
    return calculateFBDProbability();
}

double FBDTreeModel::lnPriorProbability(void){
    double lnP = 0.0;
    
    //calcualte lnP on FBD parameters
    for(Parameter* p : parameters){
        //exclude parameter tree (prior accounted for in FBD likelihood) and originAge (prior is the conditioning-age prior below)
        if( p != parameterTree && p != originAge)
            lnP += p->lnProbability();
    }

    UserSettings& us = UserSettings::userSettings();
    if(us.getConditionAgePriorSet()){
        double x = (originAge != nullptr) ? originAge->getValue() : parameterTree->getTree()->getRoot()->getTime();
        double p1 = us.getConditionAgePriorP1();
        double p2 = us.getConditionAgePriorP2();
        switch(us.getConditionAgePrior()){
            case ConditionAgePriorFamily::EXP:       lnP += Probability::Exponential::lnPdf(p1, x); break;
            case ConditionAgePriorFamily::GAMMA:     lnP += Probability::Gamma::lnPdf(p1, p2, x); break;
            case ConditionAgePriorFamily::LOGNORMAL: lnP += Probability::Normal::lnPdf(p1, p2 * p2, std::log(x)) - std::log(x); break;
            case ConditionAgePriorFamily::UNIFORM:   lnP += (x < p1 || x > p2) ? -INFINITY : Probability::Uniform::lnPdf(p1, p2, x); break;
        }
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
