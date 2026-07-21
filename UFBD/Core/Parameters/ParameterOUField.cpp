#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

#include "ParameterOUField.hpp"

#include "Node.hpp"
#include "ParameterDouble.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Serialize.hpp"
#include "Tree.hpp"

enum { OU_BIN = 0, OU_THETA = 1, OU_SDEQ = 2, OU_NU = 3 };

ParameterOUField::ParameterOUField(double prob, PhylogeneticModel* m, int nB, const std::vector<double>& lo, double init0, ParameterDouble* orig, double thetaMedianOv, double thetaSdOv, double sdShapeOv, double sdRateOv, double nuShapeOv, double nuRateOv) : Parameter(prob, m, "rateField"){
    nBins = nB;
    loEdges = lo;
    originAge = orig;

    double thMed = std::isnan(thetaMedianOv) ? 0.2 : thetaMedianOv;
    thetaMean = std::log(thMed);
    thetaSd = std::isnan(thetaSdOv) ? 0.5 : thetaSdOv;
    theta[0] = theta[1] = thetaMean;

    double r0 = (init0 > 0.0) ? init0 : thMed;
    rateVal[0].assign(nBins, r0);
    rateVal[1].assign(nBins, r0);

    if(std::isnan(sdShapeOv)){ sdShape = 6.0; sdRate = sdShape / 1.2; }
    else { sdShape = sdShapeOv; sdRate = sdRateOv; }
    sdEq[0] = sdEq[1] = sdShape / sdRate;

    std::vector<double> dts;
    for(int c = 1; c < nBins - 1; c++){
        double midPrev = 0.5 * (loEdges[c - 1] + loEdges[c]);
        double midCur  = 0.5 * (loEdges[c] + loEdges[c + 1]);
        dts.push_back(midCur - midPrev);
    }
    double medDt = 5.0;
    if(dts.empty() == false){
        std::sort(dts.begin(), dts.end());
        medDt = dts[dts.size() / 2];
    }
    if(std::isnan(nuShapeOv)){ nuShape = 4.0; nuRate = nuShape / (-std::log(0.7) / medDt); }
    else { nuShape = nuShapeOv; nuRate = nuRateOv; }
    nu[0] = nu[1] = nuShape / nuRate;

    for(int k = 0; k < 4; k++){ step[k] = 0.3; attW[k] = 0; accW[k] = 0; adaptN[k] = 0; moveAtt[k] = 0; moveAcc[k] = 0; }
    step[OU_BIN] = 0.5;
    lastMove = OU_BIN;
    lastBin = -1;
    binAtt.assign(nBins, 0);
    binAcc.assign(nBins, 0);
    numAcc = 0;
    numRej = 0;

    perBinStep = (getenv("FBD_OU_PERBIN_STEP") != nullptr);
    stepBin.assign(nBins, step[OU_BIN]);
    attWBin.assign(nBins, 0);
    accWBin.assign(nBins, 0);
    adaptNBin.assign(nBins, 0);
}

double ParameterOUField::topAge(void){
    if(originAge != nullptr)
        return originAge->getValue();
    return model->getTree()->getCrown()->getTime();
}

double ParameterOUField::lnProbability(void){
    double th = theta[0], sd = sdEq[0], nuv = nu[0];
    if(sd <= 0.0 || nuv <= 0.0)
        return -std::numeric_limits<double>::infinity();

    double lp = Probability::Normal::lnPdf(thetaMean, thetaSd * thetaSd, th);
    lp += Probability::Gamma::lnPdf(sdShape, sdRate, sd);
    lp += Probability::Gamma::lnPdf(nuShape, nuRate, nuv);

    double var0 = sd * sd;
    double x = std::log(rateVal[0][0]);
    lp += Probability::Normal::lnPdf(th, var0, x) - x;

    double top = topAge();
    double prevMid = 0.5 * (loEdges[0] + loEdges[1]);
    double prevX = x;
    for(int c = 1; c < nBins; c++){
        double mid = (c < nBins - 1) ? 0.5 * (loEdges[c] + loEdges[c + 1]) : 0.5 * (loEdges[c] + top);
        double dt = mid - prevMid;
        double rho = std::exp(-nuv * dt);
        double mean = th + (prevX - th) * rho;
        double var = var0 * (1.0 - rho * rho);
        double xc = std::log(rateVal[0][c]);
        lp += Probability::Normal::lnPdf(mean, var, xc) - xc;
        prevMid = mid;
        prevX = xc;
    }
    return lp;
}

double ParameterOUField::bactrianDelta(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double m = 0.95;
    double delta = m + Probability::Normal::rv(&rng) * std::sqrt(1.0 - m * m);
    if(rng.uniformRv() < 0.5)
        delta = -delta;
    return delta;
}

double ParameterOUField::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double pHyper = 3.0 / ((double)nBins + 3.0);
    if(rng.uniformRv() >= pHyper){
        lastMove = OU_BIN;
        int k = (int)(rng.uniformRv() * (double)nBins);
        if(k >= nBins) k = nBins - 1;
        lastBin = k;
        double s = perBinStep ? stepBin[k] : step[OU_BIN];
        double c = std::exp(s * bactrianDelta());
        rateVal[0][k] = rateVal[1][k] * c;
        return std::log(c);
    }
    double v = rng.uniformRv();
    lastBin = -1;
    if(v < 1.0 / 3.0){
        lastMove = OU_THETA;
        theta[0] = theta[1] + step[OU_THETA] * bactrianDelta();
        return 0.0;
    }
    if(v < 2.0 / 3.0){
        lastMove = OU_SDEQ;
        double c = std::exp(step[OU_SDEQ] * bactrianDelta());
        sdEq[0] = sdEq[1] * c;
        return std::log(c);
    }
    lastMove = OU_NU;
    double c = std::exp(step[OU_NU] * bactrianDelta());
    nu[0] = nu[1] * c;
    return std::log(c);
}

void ParameterOUField::adaptStep(int m, bool accepted){
    attW[m]++;
    if(accepted) accW[m]++;
    if(attW[m] >= 100){
        adaptN[m]++;
        double gain = 1.0 / std::sqrt((double)adaptN[m]);
        step[m] *= std::exp(gain * ((double)accW[m] / (double)attW[m] - 0.3));
        if(step[m] < 1e-4) step[m] = 1e-4;
        if(step[m] > 10.0)  step[m] = 10.0;
        attW[m] = 0;
        accW[m] = 0;
    }
}

void ParameterOUField::adaptStepBin(int k, bool accepted){
    attWBin[k]++;
    if(accepted) accWBin[k]++;
    if(attWBin[k] >= 100){
        adaptNBin[k]++;
        double gain = 1.0 / std::sqrt((double)adaptNBin[k]);
        stepBin[k] *= std::exp(gain * ((double)accWBin[k] / (double)attWBin[k] - 0.3));
        if(stepBin[k] < 1e-4) stepBin[k] = 1e-4;
        if(stepBin[k] > 10.0)  stepBin[k] = 10.0;
        attWBin[k] = 0;
        accWBin[k] = 0;
    }
}

void ParameterOUField::commitProposed(void){
    rateVal[1] = rateVal[0];
    theta[1] = theta[0];
    sdEq[1] = sdEq[0];
    nu[1] = nu[0];
}

void ParameterOUField::restoreProposed(void){
    rateVal[0] = rateVal[1];
    theta[0] = theta[1];
    sdEq[0] = sdEq[1];
    nu[0] = nu[1];
}

void ParameterOUField::updateForAcceptance(void){
    commitProposed();
    numAcc++;
    moveAtt[lastMove]++;
    moveAcc[lastMove]++;
    if(perBinStep && lastMove == OU_BIN && lastBin >= 0)
        adaptStepBin(lastBin, true);
    else
        adaptStep(lastMove, true);
    if(lastMove == OU_BIN && lastBin >= 0){ binAtt[lastBin]++; binAcc[lastBin]++; }
}

void ParameterOUField::updateForRejection(void){
    restoreProposed();
    numRej++;
    moveAtt[lastMove]++;
    if(perBinStep && lastMove == OU_BIN && lastBin >= 0)
        adaptStepBin(lastBin, false);
    else
        adaptStep(lastMove, false);
    if(lastMove == OU_BIN && lastBin >= 0) binAtt[lastBin]++;
}

void ParameterOUField::printPerBinAccept(std::ostream& os, const char* label) const {
    os << " " << label << "_binAR[";
    for(int k = 0; k < nBins; k++)
        os << k << ":" << (binAtt[k] > 0 ? (double)binAcc[k] / (double)binAtt[k] : 0.0)
           << "(" << binAcc[k] << "/" << binAtt[k] << ",step " << (perBinStep ? stepBin[k] : step[OU_BIN]) << ") ";
    os << "]";
    static const char* hyName[3] = {"theta", "sdEq", "nu"};
    static const int hyMove[3] = {OU_THETA, OU_SDEQ, OU_NU};
    os << " " << label << "_hyper[";
    for(int h = 0; h < 3; h++){
        int m = hyMove[h];
        os << hyName[h] << ":" << (moveAtt[m] > 0 ? (double)moveAcc[m] / (double)moveAtt[m] : 0.0)
           << "(" << moveAcc[m] << "/" << moveAtt[m] << ",step " << step[m] << ") ";
    }
    os << "]";
}

double ParameterOUField::scaleAllProposed(double c){
    for(int k = 0; k < nBins; k++)
        rateVal[0][k] = rateVal[1][k] * c;
    return (double)nBins * std::log(c);
}

double ParameterOUField::shrinkExpandProposed(double a){
    double logMean = 0.0;
    for(int k = 0; k < nBins; k++)
        logMean += std::log(rateVal[1][k]);
    logMean /= (double)nBins;
    for(int k = 0; k < nBins; k++)
        rateVal[0][k] = std::exp(logMean + a * (std::log(rateVal[1][k]) - logMean));
    return (double)(nBins - 1) * std::log(a);
}

double ParameterOUField::shiftRates(double d){
    double rBar = 0.0;
    for(int k = 0; k < nBins; k++)
        rBar += rateVal[1][k];
    rBar /= (double)nBins;
    if(rBar + d <= 0.0)
        return -std::numeric_limits<double>::infinity();
    double c = (rBar + d) / rBar;
    for(int k = 0; k < nBins; k++)
        rateVal[0][k] = rateVal[1][k] * c;
    theta[0] = theta[1] + std::log(c);
    return std::log(rBar) - std::log(rBar + d);
}

void ParameterOUField::writeState(std::ostream& os){
    Serialize::writeVec(os, rateVal[1]);
    os << theta[1] << ' ' << sdEq[1] << ' ' << nu[1] << ' '
       << step[OU_BIN] << ' ' << step[OU_THETA] << ' ' << step[OU_SDEQ] << ' ' << step[OU_NU] << '\n';
}

void ParameterOUField::readState(std::istream& is){
    Serialize::readVec(is, rateVal[1]);
    rateVal[0] = rateVal[1];
    is >> theta[1] >> sdEq[1] >> nu[1]
       >> step[OU_BIN] >> step[OU_THETA] >> step[OU_SDEQ] >> step[OU_NU];
    theta[0] = theta[1];
    sdEq[0] = sdEq[1];
    nu[0] = nu[1];
}
