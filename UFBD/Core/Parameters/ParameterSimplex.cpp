#include "ParameterSimplex.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Serialize.hpp"

#include <cmath>
#include <limits>

ParameterSimplex::ParameterSimplex(double prob, PhylogeneticModel* m, std::string n, int dimension, double concentration, double tuning) :
    Parameter(prob, m, n),
    tuning(tuning),
    numAcceptances(0),
    numRejections(0)
{
    value.assign(2, std::vector<double>(dimension, 1.0 / dimension));
    priorAlpha.assign(dimension, concentration);
}

double ParameterSimplex::lnProbability(void){
    return Probability::Dirichlet::lnPdf(priorAlpha, value[0]);
}

void ParameterSimplex::print(void){
}

double ParameterSimplex::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    int k = (int)value[1].size();
    value[0] = value[1];

    double acceptRej = 0.0;
    for(bool b : recentAcceptRej)
        if(b)
            acceptRej++;
    if(recentAcceptRej.empty() == false)
        acceptRej /= recentAcceptRej.size();
    int total = numAcceptances + numRejections;
    if(total > 0 && total % 100 == 0){
        double gain = 1.0 / std::sqrt((double)(total / 100));
        tuning *= std::exp(gain * (acceptRej - 0.3));
    }

    int i = (int)(rng.uniformRv() * k);
    int j = (int)(rng.uniformRv() * (k - 1));
    if(j >= i)
        j++;
    double m = 0.95;
    double s = std::sqrt(1.0 - m * m);
    double d = m + Probability::Normal::rv(&rng) * s;
    if(Probability::Uniform::rv(&rng, 0.0, 1.0) < 0.5)
        d = -d;
    double t = tuning * d;
    value[0][i] = value[1][i] - t;
    value[0][j] = value[1][j] + t;
    if(value[0][i] <= 0.0 || value[0][i] >= 1.0 || value[0][j] <= 0.0 || value[0][j] >= 1.0)
        return -std::numeric_limits<double>::infinity();
    return 0.0;
}

void ParameterSimplex::updateForAcceptance(void){
    value[1] = value[0];
    numAcceptances++;
    recentAcceptRej.push_back(true);
    if(recentAcceptRej.size() > 1000)
        recentAcceptRej.pop_front();
}

void ParameterSimplex::updateForRejection(void){
    value[0] = value[1];
    numRejections++;
    recentAcceptRej.push_back(false);
    if(recentAcceptRej.size() > 1000)
        recentAcceptRej.pop_front();
}

void ParameterSimplex::writeState(std::ostream& os){
    Serialize::writeVec(os, value[1]);
    os << tuning << ' ' << numAcceptances << ' ' << numRejections << '\n';
    Serialize::writeBoolDeque(os, recentAcceptRej);
}

void ParameterSimplex::readState(std::istream& is){
    Serialize::readVec(is, value[1]);
    value[0] = value[1];
    is >> tuning >> numAcceptances >> numRejections;
    Serialize::readBoolDeque(is, recentAcceptRej);
}
