#ifndef ParameterShrinkageField_hpp
#define ParameterShrinkageField_hpp

#include <deque>
#include <vector>

#include "AdaptiveOperatorSampler.hpp"
#include "Parameter.hpp"
#include "Probability.hpp"

class ParameterShrinkageField : public Parameter {

    public:
                                    ParameterShrinkageField(void) = delete;
                                    ParameterShrinkageField(double prob, PhylogeneticModel* m, int nBins, Probability::PriorSpec anchorPrior, double priorNShifts, double shiftSize);
        double                      getRate(int i) { return rateVal[0][i]; }
        int                         getNumBins(void) { return nBins; }
        double                      getGlobalScaleHyperprior(void) { return zeta; }
        double                      getGlobalScale(void) { return gamma[0]; }
        void                        setInterweave(bool b) { interweave = b; sampler.setActive(4, b); sampler.setWeight(0,0.15); sampler.setWeight(1,0.15); sampler.setWeight(2,0.30); sampler.setWeight(3,b?0.20:0.40); sampler.setWeight(4,b?0.20:0.0); }
        void                        setAdaptive(bool b) { adaptive = b; sampler.setAdapting(b); }
        double                      getMoveWeight(int m) { return sampler.weightOf(m); }
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
        double                      interweaveGlobalScale(void);
        void                        recordReward(bool accepted);
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
        double                      step[5];
        std::deque<bool>            recentAR[5];
        int                         acc[5];
        int                         rej[5];
        int                         lastMove;
        bool                        interweave;
        bool                        adaptive;
        AdaptiveOperatorSampler     sampler;
        std::vector<double>         snapLogR;
        double                      snapLogGamma;
        std::vector<double>         logrM;
        std::vector<double>         logrV;
        double                      lgM;
        double                      lgV;
};

#endif
