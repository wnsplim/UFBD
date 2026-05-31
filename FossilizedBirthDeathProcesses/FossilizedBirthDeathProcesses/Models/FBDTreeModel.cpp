#include "Node.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "FBDTreeModel.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"

#include <iostream>
#include <cmath>
#include <limits>

FBDTreeModel::FBDTreeModel(Tree* t) :
    PhylogeneticModel(),
    c1(0.0),
    c2(0.0){
    
    parameterTree = new ParameterTree(0.0, this); // proposal probability = 0
    parameterTree->setTree(t);
    parameterTree->getTree()->initializeTimes();
    parameters.push_back(parameterTree);
    
    //instantiate FBD model parameters
    lambda = new ParameterDouble(1.0, this, "lambda", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(lambda);
    mu = new ParameterDouble(1.0, this, "mu", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(mu);
    psi = new ParameterDouble(1.0, this, "psi", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(psi);
    rho = new ParameterDouble(1.0, this, "rho", 0.0, 1.0);
    parameters.push_back(rho);
    
    //normalize proposal probabilities
    double sum = 0.0;
    for(Parameter* p : parameters)
        sum += p->getProposalProbability();
    for(Parameter* p : parameters)
        p->setProposalProbability(p->getProposalProbability() / sum);
}

double FBDTreeModel::calculateFBDProbability(void){
    Tree* tree = parameterTree->getTree();
    
    int numInternalNodes = tree->getNumNodes() - tree->getNumTaxa();
    double rootAge = tree->getRoot()->getTime();
    std::vector<Node*> dpseq = tree->getDownPassSequence();
    
    lambdaVal = lambda->getValue();
    muVal = mu->getValue();
    rhoVal = rho->getValue();
    psiVal = psi->getValue();
    double log4LambdaRho = std::log(4*lambdaVal*rhoVal);
    double log2Lambda = std::log(2*lambdaVal);
    
    calculateC1();
    calculateC2();
    double qRoot = calculateQt(rootAge);
    double poHatRoot = calculatePoHat(rootAge);
    
    //FBD probability in log land
    //term 1
    double fbdProb = 0.0;
    fbdProb -= 2 * std::log(lambdaVal * (1 - poHatRoot));

    //term 2
    fbdProb += log4LambdaRho;
    fbdProb -= std::log(qRoot);

    //term 3
    fbdProb += numInternalNodes * log4LambdaRho;
    for(Node* n : dpseq)
        if(n->getIsTip() == false)
            fbdProb -= std::log(calculateQt(n->getTime()));

    //term 4
        //by convention: ancestral fossils are still flagged as tips
    for(Node* n : dpseq)
        if(n->getIsFossil())
            fbdProb += std::log(psiVal * tree->getNumLineagesAtTime(n->getTime()));

    for(Node* n : dpseq)
        if(n->getIsFossil() && (n->getDescendants().size() == 0)){
            double nTime = n->getTime();
            double nAttachTime = nTime + tree->getBranchLength(n, n->getAncestor());
            //fossil is not ancestral
            fbdProb += log2Lambda;
            fbdProb += std::log(calculatePo(nTime)) + std::log(calculateQt(nTime));
            fbdProb -= std::log(calculateQt(nAttachTime));
        }
    return fbdProb;
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
        //exclude parameter tree; we account for its prior probability in FBD likelihood
        if( p != parameterTree)
            lnP += p->lnProbability();
    }
    
    return lnP;
}

void FBDTreeModel::print(void){
    for(Parameter* p : parameters){
        if(p != parameterTree && p->getParmPrintConsole() == true)
            std::cout << p->getName() << " (A/R): " << p->getAcceptanceRatio() << "\t";
    }
    std::cout << "\n";
}

double FBDTreeModel::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    //update parameters
    double u = rng.uniformRv();
    double sum = 0.0;
    for(Parameter* p  : parameters){
        sum += p->getProposalProbability();
        if(u < sum){
            updatedParameter = p;
            break;
        }
    }
    return updatedParameter->update();
}

void FBDTreeModel::updateForAcceptance(void){
   updatedParameter->updateForAcceptance();
}

void FBDTreeModel::updateForRejection(void){
    updatedParameter->updateForRejection();
}
