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
class Tree;

class RelaxedClockTreeModel : public PhylogeneticModel {

    public:
                                    RelaxedClockTreeModel(void) = delete;
                                    RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& hessianFile, const std::string& mlTreeFile, int nStates, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param, unsigned int seed);
                                    RelaxedClockTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, const std::string& sequenceFile, const std::string& partitionFile, int nStates, int numCats, ClockModel clockModel, const double* rgeneParam, const double* sigma2Param, unsigned int seed);
        std::vector<std::string>    getParameterNames(void);
        std::vector<double>         getParameterString(void);
        std::vector<std::string>    getLatentNames(void);
        std::vector<double>         getLatentString(void);
        std::vector<std::string>    getZoneLegend(void);
        bool                        treeIncludesFossils(void);
        std::string                 getRateMap(void);
        double                      lnLikelihood(void);
        void                        invalidateLikelihoodCache(void);
        void                        invalidatePriorCache(void);
        double                      lnPriorProbability(void);
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);
        void                        freezePncpTuning(void) { clock->freezePncp(); }
        void                        setChainLabel(int c) { clock->setChainLabel(c); }

    private:
        void                        buildClock(ClockModel clockModel, const double* rgeneParam, const double* sigma2Param);
        void                        crownInitScale(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils);
        void                        collectNodeAges(std::vector<std::string>* names, std::vector<double>* vals);
        double                      nodeAgeSweep(void);
        double                      nodeAgeJump2(void);
        FBDTreeModel*               fbd;
        BranchRateModel*            clock;
        ApproxBranchLengthLikelihood* lik;
        SequenceCTMCModel*          ctmc;
        int                         lastMoveType;
    public:
        int                         getLastMoveType(void){ return lastMoveType; }
    private:
        double                      ageScaleStep = 0.5;
        int                         ageScaleAtt = 0;
        int                         ageScaleAcc = 0;
        long                        naSweepAtt = 0;
        long                        naSweepAcc = 0;
        AdaptiveMixSelector         naSel;
        std::vector<double>         naSnap;
        int                         naOp = 0;
        int                         nInternalAge = 0;
};

#endif
