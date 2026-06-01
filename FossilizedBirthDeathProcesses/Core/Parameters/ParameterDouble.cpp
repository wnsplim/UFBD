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
    targetAr(0.43){
    
    adaptiveProposalActive = true;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    double u = Probability::TruncatedNormal::rv(&rng, 0.0, 1.0, lowerBound, upperBound);
    
    value.push_back(u);
    value.push_back(u);
    windowSize = 1;
}

double ParameterDouble::lnProbability(void){
    return Probability::TruncatedNormal::lnPdf(value[0], 0.0, 1.0, lowerBound, upperBound); //this is kinda sloppy but most flexible for chaning bounds
}

void ParameterDouble::print(void){

}

double ParameterDouble::update(void) {
    return updateSlidingWindow(); //done this way so that if we want, we can implement other proposal mechanisms etc fairly painlessly
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
        std::cout << parmName << " done adaptive sampling | final acceptRej: " << acceptRej << std::endl;
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
