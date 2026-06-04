#ifndef ParameterBranchRates_hpp
#define ParameterBranchRates_hpp

#include <deque>
#include <vector>

#include "Parameter.hpp"

class Tree;

enum class ClockModel { UCLN, WN, GBM, CIR };

class ParameterBranchRates : public Parameter {

    public:
                                    ParameterBranchRates(void) = delete;
                                    ParameterBranchRates(double prob, PhylogeneticModel* m, Tree* tree, int numLoci, ClockModel clockModel, const double* rgenePara, const double* sigma2Para);
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

    private:
        double                      scaleLocusRate(int p);
        double                      scaleLocusSigma2(int p);
        double                      scaleLocusTheta(int p);
        double                      scaleBranchRate(int p, int b);
        double                      bactrianMultiplier(int moveType);
        double                      gammaDirichletLnP(const std::vector<double>& v, const double* para);
        double                      gammaLnPdf(double a, double b, double x);
        double                      lognormalLnP(double r, double s2);
        double                      whiteNoiseLnP(double r, double s2, double t);
        double                      gbmLnP(void);
        double                      cirLnP(void);
        double                      getMeanTau(double rho, double rhoUp, double t, double sigma, double theta);
        double                      besselIRatio(double nu, double x);
        Tree*                       tree;
        int                         numLoci;
        ClockModel                  clockModel;
        int                         numNodes;
        std::vector<int>            branchNodes;
        double                      rgenePara[3];
        double                      sigma2Para[3];
        double                      thetaPara[3];
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
};

#endif
