#ifndef GTRrateModel_hpp
#define GTRrateModel_hpp

#include <vector>

class GTRrateModel {

    public:
                            GTRrateModel(int numStates);
        void                setParameters(const std::vector<double>& exchangeability, const std::vector<double>& frequency);
        void                transitionProbabilities(double bl, double cat, double relVar, double* P) const;
        int                 getNumStates(void) const { return numStates; }

    private:
        int                 numStates;
        std::vector<double> eigenvalue;
        std::vector<double> cijk;
};

#endif
