#ifndef GtrRateModel_hpp
#define GtrRateModel_hpp

#include <vector>

class GtrRateModel {

    public:
                            GtrRateModel(int numStates);
        void                setParameters(const std::vector<double>& exchangeability, const std::vector<double>& frequency);
        void                transitionProbabilities(double t, double* P) const;
        int                 getNumStates(void) const { return numStates; }

    private:
        int                 numStates;
        std::vector<double> eigenvalue;
        std::vector<double> cijk;
};

#endif
