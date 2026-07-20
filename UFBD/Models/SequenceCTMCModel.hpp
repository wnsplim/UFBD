#ifndef SequenceCTMCModel_hpp
#define SequenceCTMCModel_hpp

#include <iosfwd>
#include <string>
#include <vector>

class Parameter;
class ParameterDouble;
class ParameterSimplex;
class PhylogeneticModel;
class SequenceLikelihood;
class Tree;
struct BranchMGF;

class SequenceCTMCModel {

    public:
                                        SequenceCTMCModel(PhylogeneticModel* owner, const std::string& sequenceFile, const std::string& partitionFile, int nStates, int numCats);
        void                            buildParameters(void);
        double                          computeLnL(Tree* tree, const std::vector<std::vector<double> >& branchRates, const std::vector<std::vector<BranchMGF> >& branchMGF);
        void                            invalidateCache(void);
        double                          lnPrior(void);
        double                          update(void);
        void                            updateForAcceptance(void);
        void                            updateForRejection(void);
        void                            print(void);
        void                            appendParameterNames(std::vector<std::string>& names);
        void                            appendParameterValues(std::vector<double>& values);
        void                            writeState(std::ostream& os);
        void                            readState(std::istream& is);
        int                             getNumPartitions(void) const;
        std::vector<int>                getPartitionGroups(void) const { return partitionGroup; }

    private:
        PhylogeneticModel*              owner;
        SequenceLikelihood*             seqLik;
        std::vector<int>                partitionGroup;
        std::vector<ParameterSimplex*>  exch;
        std::vector<ParameterSimplex*>  freq;
        std::vector<ParameterDouble*>   alpha;
        std::vector<ParameterDouble*>   pinv;
        std::vector<std::vector<double> > observedFreq;
        Parameter*                      lastSubstParm;
        int                             nStates;
        int                             numCats;
        bool                            empirical;
        bool                            freqEstimated;
        bool                            useGammaHet;
        bool                            useInvariant;
};

#endif
