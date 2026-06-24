#ifndef RelaxedClockTreeModel_hpp
#define RelaxedClockTreeModel_hpp

#include "FBDInput.hpp"
#include "ParameterBranchRates.hpp"
#include "PhylogeneticModel.hpp"

#include <string>
#include <vector>

class ApproxBranchLengthLikelihood;
class FBDTreeModel;
class ParameterDouble;
class ParameterSimplex;
class ParameterTree;
class SequenceCTMCModel;
class SequenceLikelihood;
class Tree;

class RelaxedClockTreeModel : public PhylogeneticModel {

    public:
                                    RelaxedClockTreeModel(void) = delete;
                                    RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& hessianFile, const std::string& mlTreeFile, int nStates, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param, unsigned int seed);
                                    RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& sequenceFile, const std::string& partitionFile, int nStates, int numCats, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param, unsigned int seed);
        std::vector<std::string>    getParameterNames(void);
        std::vector<double>         getParameterString(void);
        double                      lnLikelihood(void);
        double                      lnPriorProbability(void);
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);

    private:
        void                        buildClock(ClockModel clockModel, const double* rgeneParam, const double* sigma2Param);
        double                      nodeAgeSweep(void);
        FBDTreeModel*               fbd;
        BranchRateModel*            clock;
        ApproxBranchLengthLikelihood* lik;
        SequenceCTMCModel*          ctmc;
        int                         lastMoveType;
        double                      ageScaleStep = 0.5;
        int                         ageScaleAtt = 0;
        int                         ageScaleAcc = 0;
};

#endif
