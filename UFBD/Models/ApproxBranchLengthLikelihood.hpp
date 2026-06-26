#ifndef ApproxBranchLengthLikelihood_hpp
#define ApproxBranchLengthLikelihood_hpp

#include <set>
#include <string>
#include <vector>

class Tree;
class Node;

class ApproxBranchLengthLikelihood {

    public:
                                    ApproxBranchLengthLikelihood(void) = delete;
                                    ApproxBranchLengthLikelihood(const std::string& hessianFile, const std::string& mlTreeFile, const std::vector<std::string>& rogueTaxa, int partitionIndex, int nStates);
        double                      computeLnL(Tree* tree, const std::vector<std::vector<double>>& branchRates);
        int                         getNumPartitions(void) { return nPartitions; }

    private:
        void                        readHessianFile(const std::string& fn);
        void                        applyArcsinTransform(int p);
        void                        buildBranchOrder(const std::string& backboneNewick);
        void                        newickCanonBiparts(const std::string& nwk, std::set<std::set<std::string>>& out);
        Node*                       findNodeByBipartition(const std::set<std::string>& bp, Tree* tree);
        std::set<std::string>       backboneTipsBelow(Node* n, Tree* tree);
        std::set<std::string>       canonicalize(const std::set<std::string>& clade);

        int                         nb;
        int                         nPartitions;
        int                         crownBranchIdx;
        int                         partitionIndex;
        double                      cJc;
        std::vector<std::vector<double>>    blMle;
        std::vector<std::vector<double>>    gradient;
        std::vector<std::vector<double>>    hessian;
        std::vector<std::set<std::string>>  bipartitions;
        std::string                 hessianNewick;
        std::set<std::string>       backboneTaxa;
        std::set<std::string>       rogueTaxa;
        std::vector<int>            branchNodeIdx;
        Node*                       cachedRoot;
};

#endif
