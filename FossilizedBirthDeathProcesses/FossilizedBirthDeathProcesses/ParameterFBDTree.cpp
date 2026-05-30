#include "Node.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"
#include "ParameterFBDTree.hpp"
#include "RandomVariable.hpp"

ParameterFBDTree::ParameterFBDTree(double prob, PhylogeneticModel* m, Tree* t) :
    ParameterTree(prob, m),
    c1(0.0),
    c2(0.0){
    
    ParameterTree::trees[0] = new Tree(*t); //use explicit namespacing when accessing base class member objs.
    ParameterTree::trees[1] = new Tree(*t);
    
    //instantiate FBD model parameters
    lambda = new ParameterDouble(1.0, m, "lambda", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(lambda);
    mu = new ParameterDouble(1.0, m, "mu", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(mu);
    psi = new ParameterDouble(1.0, m, "psi", 0.0, std::numeric_limits<double>::max());
    parameters.push_back(psi);
    rho = new ParameterDouble(1.0, m, "rho", 0.0, 1.0);
    parameters.push_back(rho);
}

double ParameterFBDTree::calculateFBDProbability(void){
    double rootAge = ParameterTree::getTree()->getRoot()->getTime();
    lambdaVal = lambda->getValue();
    muVal = mu->getValue();
    rhoVal = rho->getValue();
    psiVal = psi->getValue();
    
    //expects all calculations to be in log land
    calculateC1();
    calculateC2();
    double qRoot = calculateQt(rootAge);
    double poHatRoot = calculatePoHat(rootAge);
    return 0.0;
}

void ParameterFBDTree::calculateC1(void){
    c1 =    std::abs(
                std::sqrt(
                    std::pow(lambda->getValue() - mu->getValue() - psi->getValue(), 2) +
                    4*lambda->getValue() * psi->getValue()
                )
            );
}

void ParameterFBDTree::calculateC2(void){
    c2 = (-lambda->getValue() + mu->getValue() + 2*lambda->getValue() * rho->getValue() + psi->getValue());
    c2 /= c1;
}

double ParameterFBDTree::calculateQt(double t){
    double tmp = 2 * (1 - std::pow(c2, 2));
    tmp += std::exp(-c1 * t) * std::pow(1-c2, 2);
    tmp += std::exp(c1 * t) * std::pow(1+c2, 2);
    return tmp;
}

double ParameterFBDTree::calculatePo(double t){
    double tmp = -lambda->getValue() + mu->getValue() + psi->getValue();
    tmp += c1 * (std::exp(-c1 * t) * (1 - c2) - (1+c2) ) / ( std::exp(-c1 * t) * (1 - c2) + (1+c2)  );
    tmp /= 2*lambda->getValue();
    return 1 + tmp;
}

double ParameterFBDTree::calculatePoHat(double t){
    double tmp = rho->getValue() * (lambda->getValue() - mu->getValue());
    tmp /= (lambda->getValue() * rho->getValue() + (lambda->getValue()*(1-rho->getValue()) - mu->getValue())*std::exp(-1 * (lambda->getValue() - mu->getValue()) * t) );
    return 1 + tmp;
}

bool ParameterFBDTree::getAdaptiveProposalActive(void){
    for(Parameter* p : parameters)
        if(p->getAdaptiveProposalActive() == true){
            return true;
        }
    return false;
}

double ParameterFBDTree::lnProbability(void){
    cachedLnP = 0.0;
    
    //calcualte lnP on FBD parameters
    for(Parameter* p : parameters)
        cachedLnP += p->lnProbability();
        
    //calcualte FBD probabity
    cachedLnP += calcualteFBDProbability();
   
    return 0.0;
}

void ParameterFBDTree::print(void){
    
}

double ParameterFBDTree::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    ParameterTree::useCachedLnP = false;
    
    if(rng.uniformRv() < 0.0){
        //update topology/branch length
        updatedParameter = nullptr;
        return ParameterTree::update();
    }else{
        //update FBD model parameters
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
}

void ParameterFBDTree::updateForAcceptance(void){
    if(updatedParameter == nullptr){
        ParameterTree::updateForAcceptance();
    }else{
        updatedParameter->updateForAcceptance();
    }
    ParameterTree::useCachedLnP = false;
    ParameterTree::numAcceptances++;
}

void ParameterFBDTree::updateForRejection(void){
    if(updatedParameter == nullptr){
        ParameterTree::updateForRejection();
    }else{
        updatedParameter->updateForRejection();
    }
    ParameterTree::useCachedLnP = false;
    ParameterTree::numRejections++;
}
