#ifndef RelaxedClockTreeModel_hpp
#define RelaxedClockTreeModel_hpp

#include "FBDInput.hpp"
#include "ParameterBranchRates.hpp"
#include "PhylogeneticModel.hpp"

#include <string>
#include <vector>

class ApproxBranchLengthLikelihood;
class FBDTreeModel;
class ParameterTree;
class Tree;

class RelaxedClockTreeModel : public PhylogeneticModel {

    public:
                                    RelaxedClockTreeModel(void) = delete;
                                    RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& hessianFile, const std::string& mlTreeFile, int nStates, ClockModel clockModel, const double* rgenePara, const double* sigma2Para, unsigned int seed);
        std::vector<std::string>    getParameterNames(void);
        std::vector<double>         getParameterString(void);
        double                      lnLikelihood(void);
        double                      lnPriorProbability(void);
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);

    private:
        double                      bdPrior(void);
        FBDTreeModel*               fbd;
        ParameterTree*              parameterTree;
        ParameterBranchRates*       clock;
        ApproxBranchLengthLikelihood* lik;
        double                      clockWeight;
        int                         lastMoveType;
        bool                        fixedRoot;
        bool                        isOrigin;
        double                      fixAge;
};

#endif
