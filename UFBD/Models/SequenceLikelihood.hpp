#ifndef SequenceLikelihood_hpp
#define SequenceLikelihood_hpp

#include <string>
#include <vector>

#include "GTRrateModel.hpp"
#include "ParameterBranchRates.hpp"

class Tree;
class Node;

class SequenceLikelihood {

    public:
                            SequenceLikelihood(int numStates, int numCats);
        void                addPartition(const std::vector<std::string>& taxa, const std::vector<std::vector<int>>& patterns, const std::vector<int>& weight);
        int                 getNumPartitions(void) const { return numPartitions; }
        double              computeLnL(Tree* tree,
                                       const std::vector<std::vector<double>>& branchRates,
                                       const std::vector<std::vector<double>>& exchangeability,
                                       const std::vector<std::vector<double>>& frequency,
                                       const std::vector<double>& alpha,
                                       const std::vector<double>& proportionInvariant,
                                       const std::vector<std::vector<BranchMGF>>& branchMGF);

    private:
        void                mapTaxaToNodes(Tree* tree);
        void                buildRepeats(Tree* tree);
        double              computePartitionLnL(int p, Tree* tree,
                                       const std::vector<std::vector<double>>& branchRates,
                                       const std::vector<std::vector<double>>& exchangeability,
                                       const std::vector<std::vector<double>>& frequency,
                                       const std::vector<double>& alpha,
                                       const std::vector<double>& proportionInvariant,
                                       const std::vector<std::vector<BranchMGF>>& branchMGF,
                                       bool parallelPatterns);

        int                                          numStates;
        int                                          numCats;
        int                                          numPartitions;
        std::vector<std::vector<std::string>>        taxonNames;
        std::vector<std::vector<std::vector<int>>>   patternState;
        std::vector<std::vector<int>>                patternWeight;
        std::vector<std::vector<int>>                constantState;
        std::vector<GTRrateModel>                    rateModel;

        Node*                                        mappedRoot;
        std::vector<std::vector<std::vector<int>>>   tipStateByOffset;
        std::vector<std::vector<char>>               tipMissing;
        std::vector<std::vector<std::vector<int>>>   clsId;
        std::vector<std::vector<std::vector<int>>>   clsRep;

        bool                                            cacheValid;
    public:
        void                                            invalidateCache(void) { cacheValid = false; }
    private:
        std::vector<std::vector<std::vector<double>>>   conP;
        std::vector<std::vector<std::vector<double>>>   cumScale;
        std::vector<std::vector<double>>                lastBl;
        std::vector<std::vector<BranchMGF>>             lastMGF;
        std::vector<std::vector<double>>                lastExch;
        std::vector<std::vector<double>>                lastFreq;
        std::vector<double>                             lastAlpha;
        std::vector<double>                             lastPinv;
};

#endif
