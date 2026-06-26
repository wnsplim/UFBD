#ifndef ParameterBranchRates_hpp
#define ParameterBranchRates_hpp

#include <cmath>
#include <deque>
#include <iosfwd>
#include <vector>

#include "Parameter.hpp"

class Tree;
class ParameterUnresolvedFossils;
class RandomVariable;

class AdaptiveMixSelector {
    public:
        void init(int nOps);
        int  pick(RandomVariable& rng);
        void record(int op, double jump2, double cpu);
        void writeState(std::ostream& os) const;
        void readState(std::istream& is);
    private:
        std::vector<double> cumJ2;
        std::vector<double> cumCpu;
        std::vector<long>   tries;
        int                 nOps;
};

enum class ClockModel { UCLN, WN, GBM, CIR, GBMC }; // CIR clock: halt — detached, not selectable

struct BranchMGF {
    int    kind;
    double rho, rhoUp, Ln, sigmaPB, theta, muH;
    double alpha, beta;
};

// CIR clock: halt — detached dead code (kept, not wired)
inline double cirLogBesselI(double nu, double u){
    if(nu < 0.0 || u < 0.0)
        return -INFINITY;
    if(u < 700.0)
        return std::log(std::cyl_bessel_i(nu, u));
    double z = 4.0 * nu * nu;
    double corr = 1.0 - (z - 1.0) / (8.0 * u) + (z - 1.0) * (z - 9.0) / (128.0 * u * u);
    return u - 0.5 * std::log(2.0 * M_PI * u) + std::log(corr);
}

inline double cirBridgeLogG(double eta, double r0, double rt, double t, double s2, double b){
    double bbar = std::sqrt(b * b - 2.0 * eta * s2);
    double decay = std::exp(-bbar * t);
    double om = 1.0 - decay;
    double c = 2.0 * bbar / (s2 * om);
    double nu = 2.0 * b / s2 - 1.0;
    double nu0 = b / s2 - 0.5;
    double u = 2.0 * c * std::sqrt(r0 * rt * decay);
    double A = -(b * t / s2) * (bbar - b) + ((b - bbar) / s2) * r0 - ((b + bbar) / s2) * rt - c * (r0 + rt) * decay;
    return std::log(c) + A + nu0 * std::log(rt / (r0 * decay)) + cirLogBesselI(nu, u);
}

inline double cirBridgeMGF(double eta, double r0, double rt, double t, double s2, double b){
    if(t * b < 0.0001)
        return std::exp(eta * (r0 + rt) * t * 0.5);
    return std::exp(cirBridgeLogG(eta, r0, rt, t, s2, b) - cirBridgeLogG(0.0, r0, rt, t, s2, b));
}

inline double normalCdf(double x){ return 0.5 * std::erfc(-x * M_SQRT1_2); }

inline void gbmBridgeMoments(double T, double A, double B, double u, double* mean, double* var){
    double uuT = u * u * T;
    double s = std::sqrt(uuT);
    double lz = std::log(B / A);
    double pa1 = lz / s + 0.5 * s;
    double pa2 = pa1 - s;
    if((pa1 > 2.0 && pa2 > 2.0) || (pa1 < -2.0 && pa2 < -2.0)){
        *mean = ((B - A) / lz + u * u * T * ((A + B) / (2.0 * lz * lz) - (B - A) / (lz * lz * lz)));
        *var = 0.0;
        return;
    }
    double m = (A / (u * u)) * std::sqrt(2.0 * M_PI * uuT) * std::exp((uuT / 2.0 + lz) * (uuT / 2.0 + lz) / (2.0 * uuT)) * (normalCdf(pa1) - normalCdf(pa2)) / T;
    *mean = m;
    double U = uuT, z = B / A, sU = std::sqrt(U);
    double ez2 = 2.0 * std::sqrt(2.0 * M_PI * U) * std::exp((U + lz) * (U + lz) / (2.0 * U)) * (normalCdf(lz / sU + sU) - normalCdf(lz / sU - sU))
               - (1.0 + z) * 2.0 * std::sqrt(2.0 * M_PI * U) * std::exp((U / 2.0 + lz) * (U / 2.0 + lz) / (2.0 * U)) * (normalCdf(lz / sU + 0.5 * sU) - normalCdf(lz / sU - 0.5 * sU));
    ez2 = ez2 * A * A / (u * u * u * u) / (T * T);
    double v = ez2 - m * m;
    *var = (v < 0.0) ? 0.0 : v;
}

class BranchRateModel : public Parameter {

    public:
                                    BranchRateModel(double prob, PhylogeneticModel* m, Tree* tree, int numLoci, const double* rgeneParam, const double* sigma2Param);
        double                      getAcceptanceRatio(void);
        int                         getNumLoci(void) { return numLoci; }
        double                      getLocusRate(int p) { return mu[0][p]; }
        double                      getLocusSigma2(int p) { return sigma2[0][p]; }
        int                         getNumBranchNodes(void) { return (int)branchNodes.size(); }
        void                        scaleAll(double sf);
        void                        commitAll(void);
        void                        restoreAll(void);
        double                      constantDistanceMove(void);
        double                      rateAgeSubtreeMove(void);
        void                        setUnresolvedFossils(ParameterUnresolvedFossils* u) { uf = u; }
        void                        print(void) {}
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);
        virtual std::vector<std::vector<double>> getAbsoluteRates(void) = 0;
        virtual std::vector<std::vector<BranchMGF>> getBranchMGF(void){ return std::vector<std::vector<BranchMGF>>(numLoci, std::vector<BranchMGF>(numNodes, BranchMGF{0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0})); }

    protected:
        double                      scaleLocusRate(int p);
        double                      scaleLocusSigma2(int p);
        double                      scaleBranchRate(int p, int b);
        double                      globalRateBranchRatesScale(int p);
        double                      bactrianMultiplier(int moveType);
        double                      gammaDirichletLnP(const std::vector<double>& v, const double* param);
        double                      gammaLnPdf(double a, double b, double x);
        Tree*                       tree;
        int                         numLoci;
        int                         numNodes;
        std::vector<int>            branchNodes;
        double                      rgeneParam[3];
        double                      sigma2Param[3];
        std::vector<double>         mu[2];
        std::vector<double>         sigma2[2];
        std::vector<std::vector<double>> rate[2];
        double                      step[4];
        std::deque<bool>            recentAR[4];
        int                         acc[4];
        int                         rej[4];
        int                         lastMove;
        int                         lastLocus;
        int                         lastNode;
        ParameterUnresolvedFossils* uf;
        double                      cdStep;
        long                        cdAccW;
        long                        cdAttW;
        std::vector<int>            cdNodes;
        double                      ncStep;
        long                        ncAccW;
        long                        ncAttW;
        AdaptiveMixSelector         sigSel;
        double                      sigPreLog;
        std::vector<std::vector<double>> sigTauL;
        std::vector<std::vector<double>> sigEllB;
        long                        sigRefresh;
};

class ParameterBranchRates : public BranchRateModel {

    public:
                                    ParameterBranchRates(void) = delete;
                                    ParameterBranchRates(double prob, PhylogeneticModel* m, Tree* tree, int numLoci, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param);
        std::vector<std::vector<double>> getAbsoluteRates(void);
        std::vector<std::vector<BranchMGF>> getBranchMGF(void);
        double                      lnProbability(void);
        double                      update(void);

    private:
        double                      lognormalLnP(double r, double s2, double m);
        double                      whiteNoiseLnP(double r, double s2, double t, double m);
        double                      gbmLnP(void);
        double                      gbmContinuousLnP(void);
        double                      sigmaNonCenteredMove(int p);
        double                      sigmaNonCenteredMoveGBMC(int p);
        double                      sigmaNonCenteredMoveGBM(int p);
        double                      sigmaNonCenteredMoveWN(int p);
        void                        branchLikePrecision(int p, std::vector<double>& tauL, std::vector<double>& ellB);
        ClockModel                  clockModel;
};

// CIR clock: halt — detached dead code (kept, never constructed)
class ParameterBranchRatesCIR : public BranchRateModel {

    public:
                                    ParameterBranchRatesCIR(void) = delete;
                                    ParameterBranchRatesCIR(double prob, PhylogeneticModel* m, Tree* tree, int numLoci, const double* rgeneParam, const double* sigma2Param);
        std::vector<std::vector<double>> getAbsoluteRates(void);
        std::vector<std::vector<BranchMGF>> getBranchMGF(void);
        double                      lnProbability(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);

    private:
        double                      scaleLocusTheta(int p);
        double                      scaleMuRates(int p);
        double                      sigmaNonCenteredMoveCIR(int p);
        double                      cirLnP(void);
        double                      getMeanTau(double rho, double rhoUp, double t, double sigma, double theta);
        double                      besselIRatio(double nu, double x);
        double                      thetaParam[3];
        std::vector<double>         theta[2];
};

#endif
