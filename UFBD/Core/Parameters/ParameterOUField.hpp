#ifndef ParameterOUField_hpp
#define ParameterOUField_hpp

#include <vector>

#include "Parameter.hpp"

class ParameterOUField : public Parameter {

    public:
                                    ParameterOUField(void) = delete;
                                    ParameterOUField(double prob, PhylogeneticModel* m, int nBins, const std::vector<double>& gridSpacing, double init0);
        double                      getRate(int i) { return rateVal[0][i]; }
        double                      shiftRates(double d);
        void                        commitProposed(void) { rateVal[1] = rateVal[0]; }
        void                        restoreProposed(void) { rateVal[0] = rateVal[1]; }
        int                         getNumBins(void) { return nBins; }
        double                      getAcceptanceRatio(void) { return 0.0; }
        double                      lnProbability(void) { return 0.0; }
        void                        print(void) {}
        double                      update(void) { return 0.0; }
        void                        updateForAcceptance(void) { rateVal[1] = rateVal[0]; }
        void                        updateForRejection(void) { rateVal[0] = rateVal[1]; }
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);

    private:
        int                         nBins;
        std::vector<double>         grid;
        std::vector<double>         rateVal[2];
};

#endif
