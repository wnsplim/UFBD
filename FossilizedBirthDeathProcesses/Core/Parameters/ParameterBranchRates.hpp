#ifndef ParameterBranchRates_hpp
#define ParameterBranchRates_hpp

#include <cmath>
#include <deque>
#include <vector>

#include "Parameter.hpp"

class Tree;
class ParameterUnresolvedFossils;

enum class ClockModel { UCLN, WN, GBM, CIR };

struct CirBranch {
    double rho, rhoUp, Ln, sigmaPB, theta, muH;
    bool   active;
};

inline double cirLogBesselI(double nu, double u){
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
        virtual std::vector<std::vector<double>> getAbsoluteRates(void) = 0;
        virtual std::vector<std::vector<CirBranch>> getBranchCir(void){ return std::vector<std::vector<CirBranch>>(numLoci, std::vector<CirBranch>(numNodes, CirBranch{0.0,0.0,0.0,0.0,0.0,0.0,false})); }

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
};

class ParameterBranchRates : public BranchRateModel {

    public:
                                    ParameterBranchRates(void) = delete;
                                    ParameterBranchRates(double prob, PhylogeneticModel* m, Tree* tree, int numLoci, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param);
        std::vector<std::vector<double>> getAbsoluteRates(void);
        double                      lnProbability(void);
        double                      update(void);

    private:
        double                      lognormalLnP(double r, double s2, double m);
        double                      whiteNoiseLnP(double r, double s2, double t, double m);
        double                      gbmLnP(void);
        double                      sigmaNonCenteredMove(int p);
        ClockModel                  clockModel;
};

class ParameterBranchRatesCIR : public BranchRateModel {

    public:
                                    ParameterBranchRatesCIR(void) = delete;
                                    ParameterBranchRatesCIR(double prob, PhylogeneticModel* m, Tree* tree, int numLoci, const double* rgeneParam, const double* sigma2Param);
        std::vector<std::vector<double>> getAbsoluteRates(void);
        std::vector<std::vector<CirBranch>> getBranchCir(void);
        double                      lnProbability(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);

    private:
        double                      scaleLocusTheta(int p);
        double                      scaleMuRates(int p);
        double                      cirLnP(void);
        double                      getMeanTau(double rho, double rhoUp, double t, double sigma, double theta);
        double                      besselIRatio(double nu, double x);
        double                      thetaParam[3];
        std::vector<double>         theta[2];
};

#endif
