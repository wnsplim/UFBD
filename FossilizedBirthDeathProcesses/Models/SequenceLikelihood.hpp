#ifndef SequenceLikelihood_hpp
#define SequenceLikelihood_hpp

#include <string>
#include <vector>

#include "GTRrateModel.hpp"

class Tree;
class Node;

class SequenceLikelihood {

    public:
                            SequenceLikelihood(int numStates, int numCats);
        void                addPartition(const std::vector<std::string>& taxa, const std::vector<std::vector<int>>& patterns, const std::vector<int>& weight);
        int                 getNumStates(void) const { return numStates; }
        int                 getNumPartitions(void) const { return numPartitions; }
        double              computeLnL(Tree* tree,
                                       const std::vector<std::vector<double>>& branchRates,
                                       const std::vector<std::vector<double>>& exchangeability,
                                       const std::vector<std::vector<double>>& frequency,
                                       const std::vector<double>& alpha,
                                       const std::vector<double>& proportionInvariant);

    private:
        void                mapTaxaToNodes(Tree* tree);

        int                                          numStates;
        int                                          numCats;
        int                                          numPartitions;
        std::vector<std::vector<std::string>>        taxonNames;
        std::vector<std::vector<std::vector<int>>>   patternState;
        std::vector<std::vector<int>>                patternWeight;
        std::vector<std::vector<int>>                constantState;
        std::vector<GTRrateModel>                    rateModel;

        Node*                                        mappedCrown;
        std::vector<std::vector<std::vector<int>>>   tipStateByOffset;

        bool                                            cacheValid;
        std::vector<std::vector<std::vector<double>>>   conP;
        std::vector<std::vector<double>>                lastBl;
        std::vector<std::vector<double>>                lastExch;
        std::vector<std::vector<double>>                lastFreq;
        std::vector<double>                             lastAlpha;
        std::vector<double>                             lastPinv;
};

#endif
