#ifndef ParameterBranchRates_hpp
#define ParameterBranchRates_hpp

#include <deque>
#include <vector>

#include "Parameter.hpp"

class Tree;

enum class ClockModel { UCLN, WN, GBM, CIR };

class BranchRateModel : public Parameter {

    public:
                                    BranchRateModel(double prob, PhylogeneticModel* m, std::string n) : Parameter(prob, m, n) {}
        virtual std::vector<std::vector<double>> getAbsoluteRates(void) = 0;
        virtual int                 getNumLoci(void) = 0;
        virtual double              getLocusRate(int p) = 0;
        virtual double              getLocusSigma2(int p) = 0;
        virtual int                 getNumBranchNodes(void) = 0;
        virtual void                scaleAll(double sf) = 0;
        virtual void                commitAll(void) = 0;
        virtual void                restoreAll(void) = 0;
        virtual double              constantDistanceMove(void) = 0;
};

class ParameterBranchRates : public BranchRateModel {

    public:
                                    ParameterBranchRates(void) = delete;
                                    ParameterBranchRates(double prob, PhylogeneticModel* m, Tree* tree, int numLoci, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param);
        double                      getAcceptanceRatio(void);
        std::vector<std::vector<double>> getAbsoluteRates(void);
        int                         getNumLoci(void) { return numLoci; }
        double                      getLocusRate(int p) { return mu[0][p]; }
        double                      getLocusSigma2(int p) { return sigma2[0][p]; }
        double                      lnProbability(void);
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        double                      constantDistanceMove(void);
        int                         getNumBranchNodes(void) { return (int)branchNodes.size(); }
        void                        scaleAll(double sf);
        void                        commitAll(void);
        void                        restoreAll(void);

    private:
        double                      scaleLocusRate(int p);
        double                      scaleLocusSigma2(int p);
        double                      scaleLocusTheta(int p);
        double                      scaleBranchRate(int p, int b);
        double                      globalRateBranchRatesScale(int p);
        double                      vectorScale(int p);
        double                      bactrianMultiplier(int moveType);
        double                      gammaDirichletLnP(const std::vector<double>& v, const double* param);
        double                      gammaLnPdf(double a, double b, double x);
        double                      lognormalLnP(double r, double s2, double m);
        double                      whiteNoiseLnP(double r, double s2, double t, double m);
        double                      gbmLnP(void);
        double                      cirLnP(void);
        double                      getMeanTau(double rho, double rhoUp, double t, double sigma, double theta);
        double                      besselIRatio(double nu, double x);
        Tree*                       tree;
        int                         numLoci;
        ClockModel                  clockModel;
        int                         numNodes;
        std::vector<int>            branchNodes;
        double                      rgeneParam[3];
        double                      sigma2Param[3];
        double                      thetaParam[3];
        std::vector<double>         mu[2];
        std::vector<double>         sigma2[2];
        std::vector<double>         theta[2];
        std::vector<std::vector<double>> rate[2];
        double                      step[4];
        std::deque<bool>            recentAR[4];
        int                         acc[4];
        int                         rej[4];
        int                         lastMove;
        int                         lastLocus;
        int                         lastNode;
        double                      cdStep;
        long                        cdAccW;
        long                        cdAttW;
        std::vector<int>            cdNodes;
};

#endif
