#ifndef ParameterShrinkageField_hpp
#define ParameterShrinkageField_hpp

#include <deque>
#include <vector>

#include "Parameter.hpp"
#include "Probability.hpp"

class ParameterShrinkageField : public Parameter {

    public:
                                    ParameterShrinkageField(void) = delete;
                                    ParameterShrinkageField(double prob, PhylogeneticModel* m, int nBins, Probability::PriorSpec anchorPrior, double priorNShifts, double shiftSize);
        double                      getRate(int i) { return rateVal[0][i]; }
        int                         getNumBins(void) { return nBins; }
        double                      getGlobalScaleHyperprior(void) { return zeta; }
        double                      getAcceptanceRatio(void);
        double                      lnProbability(void);
        void                        print(void) {}
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        static double               calibrateGlobalScale(int nBins, double priorNShifts, double shiftSize);

    private:
        double                      halfCauchyLnP(double x, double scale);
        double                      bactrianStep(int moveType);
        double                      gibbsScales(void);
        double                      ellipticalSliceDelta(void);
        void                        recomputeRates(void);
        int                         nBins;
        int                         nDelta;
        double                      zeta;
        Probability::PriorSpec      anchorPrior;
        double                      anchor[2];
        double                      gamma[2];
        std::vector<double>         delta[2];
        std::vector<double>         sigma[2];
        std::vector<double>         rateVal[2];
        double                      step[4];
        std::deque<bool>            recentAR[4];
        int                         acc[4];
        int                         rej[4];
        int                         lastMove;
};

#endif
