#ifndef ParameterOUField_hpp
#define ParameterOUField_hpp

#include <vector>

#include "Parameter.hpp"

class ParameterDouble;

class ParameterOUField : public Parameter {

    public:
                                    ParameterOUField(void) = delete;
                                    ParameterOUField(double prob, PhylogeneticModel* m, int nBins, const std::vector<double>& loEdges, ParameterDouble* originAge, double thetaMedianOv, double thetaSdOv, double sdShapeOv, double sdRateOv, double nuShapeOv, double nuRateOv);
        double                      getRate(int i) { return rateVal[0][i]; }
        double                      shiftRates(double d);
        void                        commitProposed(void);
        void                        restoreProposed(void);
        int                         getNumBins(void) { return nBins; }
        double                      getAcceptanceRatio(void) { return (numAcc + numRej > 0) ? (double)numAcc / (double)(numAcc + numRej) : 0.0; }
        double                      lnProbability(void);
        void                        print(void) {}
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);

    private:
        double                      bactrianDelta(void);
        double                      topAge(void);
        void                        adaptStep(int m, bool accepted);

        int                         nBins;
        std::vector<double>         loEdges;
        std::vector<double>         rateVal[2];
        double                      theta[2];
        double                      sdEq[2];
        double                      nu[2];
        ParameterDouble*            originAge;

        double                      thetaMean, thetaSd;
        double                      sdShape, sdRate;
        double                      nuShape, nuRate;

        int                         lastMove;
        double                      step[4];
        long                        attW[4], accW[4], adaptN[4];
        long                        numAcc, numRej;
};

#endif
