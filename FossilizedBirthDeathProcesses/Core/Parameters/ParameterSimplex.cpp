#include "ParameterSimplex.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"

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

double ParameterSimplex::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    int k = (int)value[0].size();
    std::vector<double> a(k);
    for(int i = 0; i < k; i++)
        a[i] = tuning * value[1][i];
    Probability::Dirichlet::rv(&rng, a, value[0]);
    std::vector<double> aProposed(k);
    for(int i = 0; i < k; i++)
        aProposed[i] = tuning * value[0][i];
    double forward = Probability::Dirichlet::lnPdf(a, value[0]);
    double reverse = Probability::Dirichlet::lnPdf(aProposed, value[1]);
    return reverse - forward;
}

void ParameterSimplex::updateForAcceptance(void){
    value[1] = value[0];
    numAcceptances++;
}

void ParameterSimplex::updateForRejection(void){
    value[0] = value[1];
    numRejections++;
}

void ParameterSimplex::print(void){
}
