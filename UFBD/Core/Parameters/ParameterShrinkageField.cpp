#include <cmath>
#include <limits>
#include <vector>

#include "ParameterShrinkageField.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Serialize.hpp"

ParameterShrinkageField::ParameterShrinkageField(double prob, PhylogeneticModel* m, int nB, Probability::PriorSpec ap, double priorNShifts, double shiftSize, double init0) : Parameter(prob, m, "shrinkageField"){
    nBins = nB;
    nDelta = nB - 1;
    anchorPrior = ap;
    zeta = calibrateGlobalScale(nB, priorNShifts, shiftSize);
    lastMove = -1;
    for(int i = 0; i < 4; i++){
        acc[i] = 0;
        rej[i] = 0;
    }
    step[0] = 0.2;
    step[1] = 0.05;
    step[2] = 0.4;
    step[3] = 0.4;
    moveWeight[0] = 0.15;
    moveWeight[1] = 0.15;
    moveWeight[2] = 0.30;
    moveWeight[3] = 0.40;
    double present0 = (init0 > 0.0) ? init0 : (ap.set ? Probability::priorMean(ap.family, ap.p1, ap.p2) : 1.0);
    if(present0 <= 0.0)
        present0 = 1.0;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    anchor[0] = anchor[1] = std::log(present0) + Probability::Normal::rv(&rng) * zeta;
    gamma[0] = gamma[1] = zeta;
    delta[0].assign(nDelta, 0.0);
    sigma[0].assign(nDelta, 1.0);
    rateVal[0].assign(nBins, present0);
    rateVal[1].assign(nBins, present0);
    for(int k = 0; k < nDelta; k++)
        delta[0][k] = Probability::Normal::rv(&rng) * zeta;
    recomputeRates();
    delta[1] = delta[0];
    sigma[1] = sigma[0];
    rateVal[1] = rateVal[0];
}

double ParameterShrinkageField::calibrateGlobalScale(int nBins, double priorNShifts, double shiftSize){
    static std::vector<int>    cN;
    static std::vector<double> cPn, cSs, cZ;
    for(size_t i = 0; i < cN.size(); i++)
        if(cN[i] == nBins && cPn[i] == priorNShifts && cSs[i] == shiftSize)
            return cZ[i];
    const int G = 2000;
    const double sqrt2 = std::sqrt(2.0);
    double shift = std::log(shiftSize);
    std::vector<double> bs(G), bg(G);
    for(int i = 0; i < G; i++){
        double q = 1e-10 + (1.0 - 2e-10) * i / (G - 1);
        q = 1.0 - (1.0 - q) / 2.0;
        bs[i] = std::tan(PI * (q - 0.5));
    }
    for(int i = 0; i < G; i++){
        double q = 1e-4 + (1.0 - 2e-4) * i / (G - 1);
        q = 1.0 - (1.0 - q) / 2.0;
        bg[i] = std::tan(PI * (q - 0.5));
    }
    double lo = 1e-12, hi = 1.0;
    for(int it = 0; it < 60; it++){
        double mid = 0.5 * (lo + hi);
        double E = 0.0;
        for(int g = 0; g < G; g++){
            double gam = mid * bg[g];
            double pt = 0.0;
            for(int i = 0; i < G; i++)
                pt += 0.5 * std::erfc(shift / (gam * bs[i] * sqrt2));
            E += (pt / (double)G / 0.5) * (double)(nBins - 1);
        }
        E /= (double)G;
        if(E < priorNShifts)
            lo = mid;
        else
            hi = mid;
    }
    double zeta = 0.5 * (lo + hi);
    cN.push_back(nBins);
    cPn.push_back(priorNShifts);
    cSs.push_back(shiftSize);
    cZ.push_back(zeta);
    return zeta;
}

double ParameterShrinkageField::halfCauchyLnP(double x, double scale){
    if(x <= 0.0)
        return -INFINITY;
    double z = x / scale;
    return std::log(2.0 / PI) - std::log(scale) - std::log(1.0 + z * z);
}

double ParameterShrinkageField::lnProbability(void){
    if(gamma[0] <= 0.0)
        return -INFINITY;
    double present = std::exp(anchor[0]);
    double lnp = Probability::priorLnPdf(anchorPrior.family, anchorPrior.p1, anchorPrior.p2, present, 0.0, std::numeric_limits<double>::max());
    lnp += anchor[0];
    lnp += halfCauchyLnP(gamma[0], zeta);
    double g2 = gamma[0] * gamma[0];
    for(int k = 0; k < nDelta; k++){
        double s = sigma[0][k];
        if(s <= 0.0)
            return -INFINITY;
        lnp += halfCauchyLnP(s, 1.0);
        lnp += Probability::Normal::lnPdf(0.0, g2 * s * s, delta[0][k]);
    }
    return lnp;
}

void ParameterShrinkageField::recomputeRates(void){
    double f = anchor[0];
    rateVal[0][0] = std::exp(f);
    for(int k = 0; k < nDelta; k++){
        f += delta[0][k];
        rateVal[0][k + 1] = std::exp(f);
    }
}

double ParameterShrinkageField::bactrianStep(int mt){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double ar = 0.0;
    for(bool b : recentAR[mt])
        if(b)
            ar++;
    if(recentAR[mt].empty() == false)
        ar /= recentAR[mt].size();
    int total = acc[mt] + rej[mt];
    if(total > 0 && total % 100 == 0){
        double gain = 1.0 / std::sqrt((double)(total / 100));
        step[mt] *= std::exp(gain * (ar - 0.3));
    }
    double m = 0.95;
    double sd = std::sqrt(1.0 - m * m);
    double d = m + Probability::Normal::rv(&rng) * sd;
    if(Probability::Uniform::rv(&rng, 0.0, 1.0) < 0.5)
        d = -d;
    return step[mt] * d;
}

double ParameterShrinkageField::gibbsScales(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double before = lnProbability();
    double g2 = gamma[0] * gamma[0];
    for(int k = 0; k < nDelta; k++){
        double s2 = sigma[0][k] * sigma[0][k];
        double sigScale = Probability::Gamma::rv(&rng, 1.0, 1.0 + 1.0 / s2) + delta[0][k] * delta[0][k] / (2.0 * g2);
        sigma[0][k] = std::sqrt(1.0 / Probability::Gamma::rv(&rng, 1.0, sigScale));
    }
    double gamScale = Probability::Gamma::rv(&rng, 1.0, 1.0 / (zeta * zeta) + 1.0 / g2);
    for(int k = 0; k < nDelta; k++)
        gamScale += 0.5 * delta[0][k] * delta[0][k] / (sigma[0][k] * sigma[0][k]);
    gamma[0] = std::sqrt(1.0 / Probability::Gamma::rv(&rng, (nDelta + 1.0) / 2.0, gamScale));
    double after = lnProbability();
    return before - after;
}

double ParameterShrinkageField::ellipticalSliceDelta(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double lnL0 = model->lnLikelihood();
    std::vector<double> f0 = delta[0];
    std::vector<double> nu(nDelta);
    for(int k = 0; k < nDelta; k++)
        nu[k] = gamma[0] * sigma[0][k] * Probability::Normal::rv(&rng);
    double logy = lnL0 + std::log(rng.uniformRv());
    double theta = 2.0 * PI * rng.uniformRv();
    double tmin = theta - 2.0 * PI;
    double tmax = theta;
    while(true){
        double c = std::cos(theta);
        double s = std::sin(theta);
        for(int k = 0; k < nDelta; k++)
            delta[0][k] = f0[k] * c + nu[k] * s;
        recomputeRates();
        if(model->lnLikelihood() > logy)
            break;
        if(theta < 0.0)
            tmin = theta;
        else
            tmax = theta;
        theta = tmin + (tmax - tmin) * rng.uniformRv();
    }
    return std::numeric_limits<double>::infinity();
}

double ParameterShrinkageField::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double u = rng.uniformRv();
    double cum = 0.0;
    int mv = 3;
    for(int i = 0; i < 4; i++){
        cum += moveWeight[i];
        if(u < cum){
            mv = i;
            break;
        }
    }
    lastMove = mv;
    if(mv == 0){
        anchor[0] += bactrianStep(0);
        recomputeRates();
        return 0.0;
    }
    if(mv == 1){
        int k = (int)(rng.uniformRv() * nDelta);
        delta[0][k] += bactrianStep(1);
        recomputeRates();
        return 0.0;
    }
    if(mv == 2)
        return gibbsScales();
    return ellipticalSliceDelta();
}

double ParameterShrinkageField::getAcceptanceRatio(void){
    int a = 0, r = 0;
    for(int i = 0; i < 4; i++){
        a += acc[i];
        r += rej[i];
    }
    return ((double)a) / ((double)a + (double)r);
}

void ParameterShrinkageField::updateForAcceptance(void){
    acc[lastMove]++;
    recentAR[lastMove].push_back(true);
    if(recentAR[lastMove].size() > 1000)
        recentAR[lastMove].pop_front();
    anchor[1] = anchor[0];
    gamma[1] = gamma[0];
    delta[1] = delta[0];
    sigma[1] = sigma[0];
    rateVal[1] = rateVal[0];
}

void ParameterShrinkageField::updateForRejection(void){
    rej[lastMove]++;
    recentAR[lastMove].push_back(false);
    if(recentAR[lastMove].size() > 1000)
        recentAR[lastMove].pop_front();
    anchor[0] = anchor[1];
    gamma[0] = gamma[1];
    delta[0] = delta[1];
    sigma[0] = sigma[1];
    rateVal[0] = rateVal[1];
}

void ParameterShrinkageField::writeState(std::ostream& os){
    os << anchor[1] << ' ' << gamma[1] << '\n';
    Serialize::writeVec(os, delta[1]);
    Serialize::writeVec(os, sigma[1]);
    Serialize::writeVec(os, rateVal[1]);
    for(int k = 0; k < 4; k++) os << step[k] << ' ';
    os << '\n';
    for(int k = 0; k < 4; k++) os << acc[k] << ' ' << rej[k] << ' ';
    os << '\n';
    for(int k = 0; k < 4; k++) Serialize::writeBoolDeque(os, recentAR[k]);
}

void ParameterShrinkageField::readState(std::istream& is){
    is >> anchor[1] >> gamma[1];
    anchor[0] = anchor[1];
    gamma[0] = gamma[1];
    Serialize::readVec(is, delta[1]);
    delta[0] = delta[1];
    Serialize::readVec(is, sigma[1]);
    sigma[0] = sigma[1];
    Serialize::readVec(is, rateVal[1]);
    rateVal[0] = rateVal[1];
    for(int k = 0; k < 4; k++) is >> step[k];
    for(int k = 0; k < 4; k++) is >> acc[k] >> rej[k];
    for(int k = 0; k < 4; k++) Serialize::readBoolDeque(is, recentAR[k]);
}
