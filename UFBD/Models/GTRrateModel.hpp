#ifndef GTRrateModel_hpp
#define GTRrateModel_hpp

#include <vector>

struct BranchMGF;

class GTRrateModel {

    public:
                            GTRrateModel(int numStates);
        void                setParameters(const std::vector<double>& exchangeability, const std::vector<double>& frequency);
        void                transitionProbabilities(double bl, double cat, const BranchMGF& cb, double* P) const;

    private:
        int                 numStates;
        std::vector<double> eigenvalue;
        std::vector<double> cijk;
};

#endif
