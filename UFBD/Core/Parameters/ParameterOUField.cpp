#include <cmath>
#include <limits>
#include <vector>

#include "ParameterOUField.hpp"

#include "PhylogeneticModel.hpp"
#include "Serialize.hpp"

ParameterOUField::ParameterOUField(double prob, PhylogeneticModel* m, int nB, const std::vector<double>& gridSpacing, double init0) : Parameter(prob, m, "rateField"){
    nBins = nB;
    grid = gridSpacing;
    double present0 = (init0 > 0.0) ? init0 : 1.0;
    rateVal[0].assign(nBins, present0);
    rateVal[1].assign(nBins, present0);
}

double ParameterOUField::shiftRates(double d){
    double rBar = 0.0;
    for(int k = 0; k < nBins; k++)
        rBar += rateVal[1][k];
    rBar /= (double)nBins;
    if(rBar + d <= 0.0)
        return -INFINITY;
    double c = (rBar + d) / rBar;
    for(int k = 0; k < nBins; k++)
        rateVal[0][k] = rateVal[1][k] * c;
    return std::log(rBar) - std::log(rBar + d);
}

void ParameterOUField::writeState(std::ostream& os){
    Serialize::writeVec(os, rateVal[1]);
}

void ParameterOUField::readState(std::istream& is){
    Serialize::readVec(is, rateVal[1]);
    rateVal[0] = rateVal[1];
}
