#include "Msg.hpp"
#include "ParameterDouble.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

ParameterDouble::ParameterDouble(double prob, PhylogeneticModel* m, std::string n, double lb, double ub) :
    Parameter(prob, m, n),
    lowerBound(lb),
    upperBound(ub),
    numRejections(0), numAcceptances(0),
    numAdaptive(10000),
    targetAr(0.3){
    
    adaptiveProposalActive = true;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    double u = Probability::TruncatedNormal::rv(&rng, 0.0, 1.0, lowerBound, upperBound);
    
    value.push_back(u);
    value.push_back(u);
    windowSize = 1;
    priorFamily = Probability::PriorFamily::TRUNCATED_NORMAL;
    priorP1 = 0.0;
    priorP2 = 1.0;
}

double ParameterDouble::lnProbability(void){
    return Probability::priorLnPdf(priorFamily, priorP1, priorP2, value[0], lowerBound, upperBound);
}

void ParameterDouble::setPrior(Probability::PriorFamily f, double p1, double p2){
    priorFamily = f;
    priorP1 = p1;
    priorP2 = p2;
    if(f == Probability::PriorFamily::FIXED){
        setValue(p1);
        setProposalProbability(0.0);
    } else if(f == Probability::PriorFamily::UNIFORM && (value[0] < p1 || value[0] > p2))
        setValue(0.5 * (p1 + p2));
}

void ParameterDouble::print(void){

}

double ParameterDouble::update(void) {
    return updateBactrianScale();
}

double ParameterDouble::updateSlidingWindow(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    double acceptRej = 0.0;
    for(bool b : recentAcceptRej)
        if(b == true)
            acceptRej++;
    acceptRej /= recentAcceptRej.size();
    
    if((numRejections + numAcceptances) % 100 ==0 && ((numRejections + numAcceptances) < numAdaptive)){
            if(acceptRej < targetAr - 0.2)
                windowSize /= 1.1;
            else if (acceptRej > targetAr + 0.2)
                windowSize *= 1.1;
            
    }else if ((numRejections + numAcceptances) == numAdaptive){
        adaptiveProposalActive = false;
    }

    double u = Probability::Uniform::rv(&rng, value[1] - windowSize, value[1] + windowSize);
    if(u <= lowerBound){ //bounce
        u = lowerBound + (lowerBound - u);
    }else if ( u >= upperBound){
        u = upperBound - ( u - upperBound);
    }
    value[0] = u;
    return 0.0;
}

double ParameterDouble::updateBactrianScale(){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double acceptRej = 0.0;
    for(bool b : recentAcceptRej)
        if(b == true)
            acceptRej++;
    acceptRej /= recentAcceptRej.size();

    int total = numRejections + numAcceptances;
    if(total > 0 && total % 100 == 0){
        double gain = 1.0 / std::sqrt((double)(total / 100));
        windowSize *= std::exp(gain * (acceptRej - targetAr));
    }

    double m = 0.95;
    double s = std::sqrt(1.0 - m * m);
    double delta = m + Probability::Normal::rv(&rng) * s;
    if(Probability::Uniform::rv(&rng, 0.0, 1.0) < 0.5)
        delta = -delta;
    double logScale = windowSize * delta;
    value[0] = value[1] * std::exp(logScale);
    return logScale;
}

void ParameterDouble::updateForAcceptance(void) {
    numAcceptances++;
    value[1] = value[0];
    recentAcceptRej.push_back(true);
    if(recentAcceptRej.size() > 1000)
            recentAcceptRej.pop_front();
}

void ParameterDouble::updateForRejection(void) {
    numRejections++;
    value[0] = value[1];
    recentAcceptRej.push_back(false);
    if(recentAcceptRej.size() > 1000)
            recentAcceptRej.pop_front();
}
