#ifndef ParameterSimplex_hpp
#define ParameterSimplex_hpp

#include "Parameter.hpp"
#include <iosfwd>
#include <string>
#include <vector>

class ParameterSimplex : public Parameter {

    public:
                                    ParameterSimplex(void) = delete;
                                    ParameterSimplex(double prob, PhylogeneticModel* m, std::string n, int dimension, double concentration, double tuning);
        double                      getAcceptanceRatio(void) { return ((double)numAcceptances) / ((double)(numAcceptances + numRejections)); }
        double                      getStep(void) { return tuning; }
        const std::vector<double>&  getValue(void) { return value[0]; }
        void                        setValue(const std::vector<double>& v) { value[0] = v; value[1] = v; }
        double                      lnProbability(void);
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);

    private:
        std::vector<std::vector<double>> value;
        std::vector<double>              priorAlpha;
        double                           tuning;
        int                              numAcceptances;
        int                              numRejections;
};

#endif
