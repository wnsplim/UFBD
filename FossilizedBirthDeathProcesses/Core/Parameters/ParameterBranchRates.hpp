#ifndef ParameterBranchRates_hpp
#define ParameterBranchRates_hpp

#include <deque>
#include <vector>

#include "Parameter.hpp"

class Tree;

enum class ClockModel { UCLN, WN, GBM, CIR };

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
        void                        print(void) {}
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        virtual std::vector<std::vector<double>> getAbsoluteRates(void) = 0;

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
        double                      step[3];
        std::deque<bool>            recentAR[3];
        int                         acc[3];
        int                         rej[3];
        int                         lastMove;
        int                         lastLocus;
        int                         lastNode;
        double                      cdStep;
        long                        cdAccW;
        long                        cdAttW;
        std::vector<int>            cdNodes;
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
        ClockModel                  clockModel;
};

class ParameterBranchRatesCIR : public BranchRateModel {

    public:
                                    ParameterBranchRatesCIR(void) = delete;
                                    ParameterBranchRatesCIR(double prob, PhylogeneticModel* m, Tree* tree, int numLoci, const double* rgeneParam, const double* sigma2Param);
        std::vector<std::vector<double>> getAbsoluteRates(void);
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
