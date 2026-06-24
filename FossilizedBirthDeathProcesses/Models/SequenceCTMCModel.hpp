#ifndef SequenceCTMCModel_hpp
#define SequenceCTMCModel_hpp

#include <string>
#include <vector>

class Parameter;
class ParameterDouble;
class ParameterSimplex;
class PhylogeneticModel;
class SequenceLikelihood;
class Tree;

class SequenceCTMCModel {

    public:
                                        SequenceCTMCModel(PhylogeneticModel* owner, const std::string& sequenceFile, const std::string& partitionFile, int nStates, int numCats);
        void                            buildParameters(void);
        double                          computeLnL(Tree* tree, const std::vector<std::vector<double> >& branchRates);
        double                          lnPrior(void);
        double                          update(void);
        void                            updateForAcceptance(void);
        void                            updateForRejection(void);
        void                            appendParameterNames(std::vector<std::string>& names);
        void                            appendParameterValues(std::vector<double>& values);
        int                             getNumPartitions(void) const;

    private:
        PhylogeneticModel*              owner;
        SequenceLikelihood*             seqLik;
        std::vector<ParameterSimplex*>  exch;
        std::vector<ParameterSimplex*>  freq;
        std::vector<ParameterDouble*>   alpha;
        std::vector<ParameterDouble*>   pinv;
        Parameter*                      lastSubstParm;
        int                             nStates;
        int                             numCats;
};

#endif
