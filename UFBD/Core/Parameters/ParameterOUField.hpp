#ifndef ParameterOUField_hpp
#define ParameterOUField_hpp

#include <cmath>
#include <vector>

#include "Parameter.hpp"
#include "UserSettings.hpp"

class ParameterDouble;
class RandomVariable;

class ParameterOUField : public Parameter {

    public:
                                    ParameterOUField(void) = delete;
                                    ParameterOUField(double prob, PhylogeneticModel* m, int nBins, const std::vector<double>& loEdges, double init0, ParameterDouble* originAge, const OUHyperSpec& ou);
        void                        drawFromPrior(RandomVariable* rng);
        void                        setState(const std::vector<double>& r, double th, double sd, double n);
        double                      getRate(int i) { return rateVal[0][i]; }
        double                      shiftRates(double d);
        double                      scaleAllProposed(double c, bool withTheta = false);
        double                      shrinkExpandProposed(double a);
        void                        commitProposed(void);
        void                        restoreProposed(void);
//         void                        copyState(const ParameterOUField& o) { rateVal[0] = o.rateVal[0]; rateVal[1] = o.rateVal[1]; copyHyper(o); }
        int                         getNumBins(void) { return nBins; }
        double                      getTheta(void) { return theta[0]; }
        double                      getSdEq(void) { return sdEq[0]; }
        double                      getNu(void) { return nu[0]; }
        bool                        isThetaFixed(void) { return thetaPrior.family == Probability::PriorFamily::FIXED; }
        bool                        isSdFixed(void) { return sdPrior.family == Probability::PriorFamily::FIXED; }
        bool                        isNuFixed(void) { return nuPrior.family == Probability::PriorFamily::FIXED; }
//         std::pair<double,double>    sdSupport(void) { return Probability::priorSupport(sdPrior.family, sdPrior.p1, sdPrior.p2, 0.0, INFINITY, sdPrior.p3); }
//         std::pair<double,double>    nuSupport(void) { return Probability::priorSupport(nuPrior.family, nuPrior.p1, nuPrior.p2, 0.0, INFINITY, nuPrior.p3); }
//         double                      standardizedDeviations(std::vector<double>& u);
//         double                      standardizedDeviations(const std::vector<double>& rates, double th, double sd, double nuv, std::vector<double>& u);
//         void                        setStandardizedDeviations(const std::vector<double>& u);
//         void                        shareHyperparameters(ParameterOUField* p);
//         bool                        sharesHyperparameters(void) { return partner != nullptr; }
        double                      getAcceptanceRatio(void) { return (numAcc + numRej > 0) ? (double)numAcc / (double)(numAcc + numRej) : 0.0; }
        double                      lnProbability(void);
        void                        print(void) {}
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);
        bool                        readUntaggedPerBin(std::istream& is, size_t n);
        void                        printPerBinAccept(std::ostream& os, const char* label) const;

    private:
        enum OuMove { OU_BIN = 0, OU_THETA = 1, OU_SDEQ = 2, OU_NU = 3, OU_JOINTSD = 4, OU_NUNC = 5, OU_NMOVES = 6 };

        double                      binMid(int c);
        double                      bactrianDelta(void);
        double                      jointSdScale(void);
        double                      jointSdScaleRates(double c, double h);
        double                      nuNonCentered(void);
        double                      nuNonCenteredRates(double nuOld, double nuNew, double h);
//         void                        copyHyper(const ParameterOUField& o) { for(int k = 0; k < 2; k++){ theta[k] = o.theta[k]; sdEq[k] = o.sdEq[k]; nu[k] = o.nu[k]; } }
        void                        adaptStep(int m, bool accepted);
        void                        adaptStepBin(int k, bool accepted);

        int                         nBins;
        std::vector<double>         loEdges;
        std::vector<double>         rateVal[2];
        double                      theta[2];
        double                      sdEq[2];
        double                      nu[2];
//         ParameterOUField*           partner = nullptr;
        ParameterDouble*            originAge;

        Probability::PriorSpec      thetaPrior, sdPrior, nuPrior;
        std::vector<int>            hyperMoves;

        int                         lastMove;
        double                      step[OU_NMOVES];
        long                        attW[OU_NMOVES], accW[OU_NMOVES], adaptN[OU_NMOVES];
        long                        moveAtt[OU_NMOVES], moveAcc[OU_NMOVES];
        long                        numAcc, numRej;
        int                         lastBin;
        std::vector<long>           binAtt, binAcc;

        std::vector<double>         stepBin;
        std::vector<long>           attWBin, accWBin, adaptNBin;
};

#endif
