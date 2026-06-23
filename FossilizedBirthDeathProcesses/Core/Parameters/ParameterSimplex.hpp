#ifndef ParameterSimplex_hpp
#define ParameterSimplex_hpp

#include "Parameter.hpp"
#include <string>
#include <vector>

class ParameterSimplex : public Parameter {

    public:
                                    ParameterSimplex(void) = delete;
                                    ParameterSimplex(double prob, PhylogeneticModel* m, std::string n, int dimension, double concentration, double tuning);
        double                      getAcceptanceRatio(void) { return ((double)numAcceptances) / ((double)(numAcceptances + numRejections)); }
        const std::vector<double>&  getValue(void) { return value[0]; }
        double                      lnProbability(void);
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);

    private:
        std::vector<std::vector<double>> value;
        std::vector<double>              priorAlpha;
        double                           tuning;
        int                              numAcceptances;
        int                              numRejections;
};

#endif
