#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "ParameterOUField.hpp"

#include "Node.hpp"
#include "ParameterDouble.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Serialize.hpp"
#include "Tree.hpp"

ParameterOUField::ParameterOUField(double prob, PhylogeneticModel* m, int nB, const std::vector<double>& lo, double init0, ParameterDouble* orig, const OUHyperSpec& ou) : Parameter(prob, m, "rateField"){
    nBins = nB;
    loEdges = lo;
    originAge = orig;

    thetaPrior = ou.thetaPrior.set ? ou.thetaPrior : Probability::PriorSpec{true, Probability::PriorFamily::NORMAL, std::log(0.2), 1.6, 0.0};
    theta[0] = theta[1] = Probability::priorMean(thetaPrior.family, thetaPrior.p1, thetaPrior.p2, thetaPrior.p3);

    double r0 = (init0 > 0.0) ? init0 : std::exp(theta[0]);
    rateVal[0].assign(nBins, r0);
    rateVal[1].assign(nBins, r0);

    sdPrior = ou.sdPrior.set ? ou.sdPrior : Probability::PriorSpec{true, Probability::PriorFamily::EXPONENTIAL, 1.0 / 0.3, 1.0, 0.0};
    nuPrior = ou.nuPrior.set ? ou.nuPrior : Probability::PriorSpec{true, Probability::PriorFamily::EXPONENTIAL, 2.0, 1.0, 0.0};
    sdEq[0] = sdEq[1] = Probability::priorMean(sdPrior.family, sdPrior.p1, sdPrior.p2, sdPrior.p3);
    nu[0] = nu[1] = Probability::priorMean(nuPrior.family, nuPrior.p1, nuPrior.p2, nuPrior.p3);
    if(isThetaFixed() == false) hyperMoves.push_back(OU_THETA);
    if(isSdFixed() == false) hyperMoves.push_back(OU_SDEQ);
    if(isNuFixed() == false) hyperMoves.push_back(OU_NU);
    if(isSdFixed() == false) hyperMoves.push_back(OU_JOINTSD);
    if(isNuFixed() == false) hyperMoves.push_back(OU_NUNC);

    for(int k = 0; k < OU_NMOVES; k++){ step[k] = 0.3; attW[k] = 0; accW[k] = 0; adaptN[k] = 0; moveAtt[k] = 0; moveAcc[k] = 0; }
    step[OU_BIN] = 0.5;
    lastMove = OU_BIN;
    lastBin = -1;
    binAtt.assign(nBins, 0);
    binAcc.assign(nBins, 0);
    numAcc = 0;
    numRej = 0;

    stepBin.assign(nBins, step[OU_BIN]);
    attWBin.assign(nBins, 0);
    accWBin.assign(nBins, 0);
    adaptNBin.assign(nBins, 0);
}

double ParameterOUField::lnProbability(void){
    double th = theta[0], sd = sdEq[0], nuv = nu[0];
    if(sd <= 0.0 || nuv <= 0.0)
        return -std::numeric_limits<double>::infinity();

    double lp = Probability::priorLnPdf(thetaPrior.family, thetaPrior.p1, thetaPrior.p2, th, -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), thetaPrior.p3);
    lp += Probability::priorLnPdf(sdPrior.family, sdPrior.p1, sdPrior.p2, sd, 0.0, std::numeric_limits<double>::infinity(), sdPrior.p3);
    lp += Probability::priorLnPdf(nuPrior.family, nuPrior.p1, nuPrior.p2, nuv, 0.0, std::numeric_limits<double>::infinity(), nuPrior.p3);

    double var0 = sd * sd;
    double prevX = std::log(rateVal[0][0]);
    lp += Probability::Normal::lnPdf(th, var0, prevX) - prevX;
    for(int c = 1; c < nBins; c++){
        double rho = std::exp(-nuv * (binMid(c) - binMid(c - 1)));
        double mean = th + (prevX - th) * rho;
        double var = var0 * (1.0 - rho * rho);
        double xc = std::log(rateVal[0][c]);
        lp += Probability::Normal::lnPdf(mean, var, xc) - xc;
        prevX = xc;
    }
    return lp;
}

double ParameterOUField::binMid(int c){
    if(c == 0){
        double w0 = loEdges[1] - loEdges[0];
        if(nBins > 2)
            w0 = std::min(w0, loEdges[2] - loEdges[1]);
        return loEdges[0] + 0.5 * w0;
    }
    if(c < nBins - 1)
        return 0.5 * (loEdges[c] + loEdges[c + 1]);
    return loEdges[c] + 0.5 * (loEdges[c] - loEdges[c - 1]);
}

void ParameterOUField::drawFromPrior(RandomVariable* rng){
    double th = Probability::priorRv(thetaPrior.family, thetaPrior.p1, thetaPrior.p2, rng, thetaPrior.p3);
    double sd = Probability::priorRv(sdPrior.family, sdPrior.p1, sdPrior.p2, rng, sdPrior.p3);
    double nuv = Probability::priorRv(nuPrior.family, nuPrior.p1, nuPrior.p2, rng, nuPrior.p3);
    std::vector<double> r(nBins);
    double x = Probability::Normal::rv(rng, th, sd);
    r[0] = std::exp(x);
    for(int c = 1; c < nBins; c++){
        double rho = std::exp(-nuv * (binMid(c) - binMid(c - 1)));
        x = Probability::Normal::rv(rng, th + (x - th) * rho, sd * std::sqrt(1.0 - rho * rho));
        r[c] = std::exp(x);
    }
    setState(r, th, sd, nuv);
}

void ParameterOUField::setState(const std::vector<double>& r, double th, double sd, double n){
    rateVal[0] = r;
    theta[0] = th;
    sdEq[0] = sd;
    nu[0] = n;
    commitProposed();
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
    double nHyper = (double)hyperMoves.size();
    double pHyper = nHyper / ((double)nBins + nHyper);
    if(rng.uniformRv() >= pHyper){
        lastMove = OU_BIN;
        int k = (int)(rng.uniformRv() * (double)nBins);
        if(k >= nBins) k = nBins - 1;
        lastBin = k;
        double s = stepBin[k];
        double c = std::exp(s * bactrianDelta());
        rateVal[0][k] = rateVal[1][k] * c;
        return std::log(c);
    }
    lastBin = -1;
    lastMove = hyperMoves[(int)(rng.uniformRv() * nHyper)];
    if(lastMove == OU_THETA){
        theta[0] = theta[1] + step[OU_THETA] * bactrianDelta();
//         if(partner != nullptr)
//             partner->theta[0] = theta[0];
        return 0.0;
    }
    if(lastMove == OU_SDEQ){
        double c = std::exp(step[OU_SDEQ] * bactrianDelta());
        sdEq[0] = sdEq[1] * c;
//         if(partner != nullptr)
//             partner->sdEq[0] = sdEq[0];
        return std::log(c);
    }
    if(lastMove == OU_JOINTSD)
        return jointSdScale();
    if(lastMove == OU_NUNC)
        return nuNonCentered();
    double c = std::exp(step[OU_NU] * bactrianDelta());
    nu[0] = nu[1] * c;
//     if(partner != nullptr)
//         partner->nu[0] = nu[0];
    return std::log(c);
}

// double ParameterOUField::standardizedDeviations(std::vector<double>& u){
//     return standardizedDeviations(rateVal[0], theta[0], sdEq[0], nu[0], u);
// }
// 
// double ParameterOUField::standardizedDeviations(const std::vector<double>& rates, double th, double sd, double nuv, std::vector<double>& u){
//     u.resize(nBins);
//     double xPrev = std::log(rates[0]);
//     u[0] = (xPrev - th) / sd;
//     double lnJac = -xPrev - (double)nBins * std::log(sd);
//     for(int i = 1; i < nBins; i++){
//         double dt = binMid(i) - binMid(i - 1);
//         double q = -std::expm1(-2.0 * nuv * dt);
//         double x = std::log(rates[i]);
//         u[i] = (x - th - std::exp(-nuv * dt) * (xPrev - th)) / (sd * std::sqrt(q));
//         lnJac -= x + 0.5 * std::log(q);
//         xPrev = x;
//     }
//     return lnJac;
// }
// 
// void ParameterOUField::setStandardizedDeviations(const std::vector<double>& u){
//     double th = theta[0], sd = sdEq[0];
//     double x = th + sd * u[0];
//     rateVal[0][0] = std::exp(x);
//     for(int i = 1; i < nBins; i++){
//         double dt = binMid(i) - binMid(i - 1);
//         x = th + std::exp(-nu[0] * dt) * (x - th) + sd * std::sqrt(-std::expm1(-2.0 * nu[0] * dt)) * u[i];
//         rateVal[0][i] = std::exp(x);
//     }
//     commitProposed();
// }
// 
// void ParameterOUField::shareHyperparameters(ParameterOUField* p){
//     p->copyHyper(*this);
//     p->hyperMoves.clear();
//     p->partner = this;
//     partner = p;
// }
// 
double ParameterOUField::nuNonCentered(void){
    double lnc = step[OU_NUNC] * bactrianDelta();
    double nuOld = nu[1], nuNew = nuOld * std::exp(lnc);
    nu[0] = nuNew;
//     double h = nuNonCenteredRates(nuOld, nuNew, lnc);
//     if(partner != nullptr){
//         partner->nu[0] = nuNew;
//         h = partner->nuNonCenteredRates(nuOld, nuNew, h);
//     }
//     return h;
    return nuNonCenteredRates(nuOld, nuNew, lnc);
}

double ParameterOUField::nuNonCenteredRates(double nuOld, double nuNew, double h){
    double th = theta[1], sd = sdEq[1];
    double xPrev = std::log(rateVal[1][0]), xPrevNew = xPrev;
    for(int i = 1; i < nBins; i++){
        double dt = binMid(i) - binMid(i - 1);
        double q = -std::expm1(-2.0 * nuOld * dt), qNew = -std::expm1(-2.0 * nuNew * dt);
        double x = std::log(rateVal[1][i]);
        double u = (x - th - std::exp(-nuOld * dt) * (xPrev - th)) / (sd * std::sqrt(q));
        double xn = th + std::exp(-nuNew * dt) * (xPrevNew - th) + sd * std::sqrt(qNew) * u;
        h += (xn - x) + 0.5 * std::log(qNew / q);
        rateVal[0][i] = std::exp(xn);
        xPrev = x;
        xPrevNew = xn;
    }
    return h;
}

double ParameterOUField::jointSdScale(void){
    double c = std::exp(step[OU_JOINTSD] * bactrianDelta());
    sdEq[0] = sdEq[1] * c;
//     double h = jointSdScaleRates(c, (double)(nBins + 1) * std::log(c));
//     if(partner != nullptr){
//         partner->sdEq[0] = sdEq[0];
//         h = partner->jointSdScaleRates(c, h + (double)partner->nBins * std::log(c));
//     }
//     return h;
    return jointSdScaleRates(c, (double)(nBins + 1) * std::log(c));
}

double ParameterOUField::jointSdScaleRates(double c, double h){
    for(int k = 0; k < nBins; k++){
        double x = std::log(rateVal[1][k]);
        double xn = theta[1] + c * (x - theta[1]);
        rateVal[0][k] = std::exp(xn);
        h += xn - x;
    }
    return h;
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
//     if(partner != nullptr){
//         partner->rateVal[1] = partner->rateVal[0];
//         partner->copyHyper(*this);
//     }
}

void ParameterOUField::restoreProposed(void){
    rateVal[0] = rateVal[1];
    theta[0] = theta[1];
    sdEq[0] = sdEq[1];
    nu[0] = nu[1];
//     if(partner != nullptr){
//         partner->rateVal[0] = partner->rateVal[1];
//         partner->copyHyper(*this);
//     }
}

void ParameterOUField::updateForAcceptance(void){
    commitProposed();
    numAcc++;
    moveAtt[lastMove]++;
    moveAcc[lastMove]++;
    if(lastMove == OU_BIN && lastBin >= 0)
        adaptStepBin(lastBin, true);
    else
        adaptStep(lastMove, true);
    if(lastMove == OU_BIN && lastBin >= 0){ binAtt[lastBin]++; binAcc[lastBin]++; }
}

void ParameterOUField::updateForRejection(void){
    restoreProposed();
    numRej++;
    moveAtt[lastMove]++;
    if(lastMove == OU_BIN && lastBin >= 0)
        adaptStepBin(lastBin, false);
    else
        adaptStep(lastMove, false);
    if(lastMove == OU_BIN && lastBin >= 0) binAtt[lastBin]++;
}

void ParameterOUField::printPerBinAccept(std::ostream& os, const char* label) const {
    os << " " << label << "_binAR[";
    for(int k = 0; k < nBins; k++)
        os << k << ":" << (binAtt[k] > 0 ? (double)binAcc[k] / (double)binAtt[k] : 0.0)
           << "(" << binAcc[k] << "/" << binAtt[k] << ",step " << stepBin[k] << ") ";
    os << "]";
    static const char* hyName[5] = {"theta", "sdEq", "nu", "jointSd", "nuNC"};
    static const int hyMove[5] = {OU_THETA, OU_SDEQ, OU_NU, OU_JOINTSD, OU_NUNC};
    os << " " << label << "_hyper[";
    for(int h = 0; h < 5; h++){
        int m = hyMove[h];
        os << hyName[h] << ":" << (moveAtt[m] > 0 ? (double)moveAcc[m] / (double)moveAtt[m] : 0.0)
           << "(" << moveAcc[m] << "/" << moveAtt[m] << ",step " << step[m] << ") ";
    }
    os << "]";
}

double ParameterOUField::scaleAllProposed(double c, bool withTheta){
    for(int k = 0; k < nBins; k++)
        rateVal[0][k] = rateVal[1][k] * c;
//     if(withTheta && isThetaFixed() == false){
    if(withTheta && isThetaFixed() == false)
        theta[0] = theta[1] + std::log(c);
//         if(partner != nullptr)
//             partner->theta[0] = theta[0];
//     }
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
//     if(isThetaFixed() == false){
    if(isThetaFixed() == false)
        theta[0] = theta[1] + std::log(c);
//         if(partner != nullptr)
//             partner->theta[0] = theta[0];
//     }
    return (double)(nBins - 1) * std::log(c);
}

void ParameterOUField::writeState(std::ostream& os){
    Serialize::writeVec(os, rateVal[1]);
    os << theta[1] << ' ' << sdEq[1] << ' ' << nu[1] << ' '
       << step[OU_BIN] << ' ' << step[OU_THETA] << ' ' << step[OU_SDEQ] << ' ' << step[OU_NU] << '\n';
    os << "perbin\n";
    Serialize::writeVec(os, stepBin);
    Serialize::writeLVec(os, attWBin);
    Serialize::writeLVec(os, accWBin);
    Serialize::writeLVec(os, adaptNBin);
    for(int m = 0; m < 4; m++) os << attW[m] << ' ' << accW[m] << ' ' << adaptN[m] << ' ';
    os << '\n';
    os << "jointsd " << step[OU_JOINTSD] << ' ' << attW[OU_JOINTSD] << ' ' << accW[OU_JOINTSD] << ' ' << adaptN[OU_JOINTSD] << '\n';
    os << "nunc " << step[OU_NUNC] << ' ' << attW[OU_NUNC] << ' ' << accW[OU_NUNC] << ' ' << adaptN[OU_NUNC] << '\n';
}

bool ParameterOUField::readUntaggedPerBin(std::istream& is, size_t n){
    std::vector<double> sb(n);
    std::vector<long> aw(n), cw(n), an(n);
    long long c = 0;
    bool ok = (bool)(is >> c) && c == (long long)n;
    for(size_t i = 0; ok && i < n; i++)
        ok = (bool)(is >> sb[i]);
    std::vector<long>* dst[3] = { &aw, &cw, &an };
    for(int v = 0; v < 3 && ok; v++){
        ok = (bool)(is >> c) && c == (long long)n;
        for(size_t i = 0; ok && i < n; i++)
            ok = (bool)(is >> (*dst[v])[i]);
    }
    if(ok == false){
        is.clear();
        return false;
    }
    stepBin = sb;
    attWBin = aw;
    accWBin = cw;
    adaptNBin = an;
    return true;
}

void ParameterOUField::readState(std::istream& is){
    Serialize::readVec(is, rateVal[1]);
    rateVal[0] = rateVal[1];
    is >> theta[1] >> sdEq[1] >> nu[1]
       >> step[OU_BIN] >> step[OU_THETA] >> step[OU_SDEQ] >> step[OU_NU];
    size_t n = rateVal[1].size();
    std::streampos here = is.tellg();
    std::string tag;
    is >> tag;
    bool block = (tag == "perbin");
    if(block){
        Serialize::readVec(is, stepBin);
        Serialize::readLVec(is, attWBin);
        Serialize::readLVec(is, accWBin);
        Serialize::readLVec(is, adaptNBin);
    }else{
        is.clear();
        is.seekg(here);
        block = readUntaggedPerBin(is, n);
        if(block == false)
            is.seekg(here);
    }
    if(block){
        for(int m = 0; m < 4; m++) is >> attW[m] >> accW[m] >> adaptN[m];
    }else{
        stepBin.assign(n, step[OU_BIN]);
        attWBin.assign(n, 0);
        accWBin.assign(n, 0);
        adaptNBin.assign(n, 0);
        for(int m = 0; m < 4; m++){ attW[m] = 0; accW[m] = 0; adaptN[m] = 0; }
    }
    for(;;){
        std::streampos at = is.tellg();
        int m = -1;
        if(is >> tag)
            m = (tag == "jointsd") ? OU_JOINTSD : (tag == "nunc" ? OU_NUNC : -1);
        if(m < 0){
            is.clear();
            is.seekg(at);
            break;
        }
        is >> step[m] >> attW[m] >> accW[m] >> adaptN[m];
    }
    theta[0] = theta[1];
    sdEq[0] = sdEq[1];
    nu[0] = nu[1];
}
