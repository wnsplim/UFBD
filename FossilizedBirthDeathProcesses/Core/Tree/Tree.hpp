#ifndef Tree_hpp
#define Tree_hpp

#include <map>
#include <sstream>
#include <string>
#include <vector>

class Node;
class RandomVariable;

typedef std::map<std::pair<Node*,Node*>,double> BranchLengths;


class Tree {

    public:
                                            Tree(void) = delete;
                                            Tree(const Tree& t); //copy constructor
#if 0
                                            Tree(std::vector<std::string> taxonNames, double lambda); //random tree w/ taxon names + branches from exp dist
#endif
                                            Tree(std::string newick);
                                           ~Tree(void);
        Tree&                               operator=(const Tree& t);
#if 0
        void                                calculateTreeHeight(void);
#endif
        void                                checkBranchLengthsNeg(void);
#if 0
        void                                dropTip(std::string tip);
        void                                dropTip(int idx);
        void                                forceBinary(void);
#endif
        std::vector<Node*>                  getAllDescendants(Node* n);
        double                              getBranchLength(Node* e1, Node* e2);
        std::vector<Node*>&                 getDownPassSequence(void) { return downPassSequence; }
        std::string                         getNewickString(void);
        int                                 getNumLineagesAtTime(double t);
        bool                                isSATip(Node* n);
        bool                                isFakeSplit(Node* n);
        int                                 getNumNodes(void) { return (int)nodes.size(); }
        Node*                               getNodeByOffset(int o) { return nodes[o]; }
        int                                 getNumTaxa(void);
        Node*                               getCrown(void) { return crown; }
        bool                                getLastUpdateWasScale(void) { return lastUpdateWasScale; }
        bool                                isBinary(void);
        void                                validateBackbone(void);
        static bool                         isValidNewick(const std::string& s);
        Node*                               getTaxonNode(std::string name);
        Node*                               getMRCA(const std::vector<std::string>& taxonNames);
        void                                initializeDownPassSequence(void);
        void                                initializeTimes(void); //starting node ages from topology; extant tips at 0
#if 0
        void                                keepTips(std::vector<std::string> t);
#endif
        double                              update(double scaleLambda);
        int                                 scaleInternalAges(double m);
        int                                 scaleSubtreeAges(Node* subtreeCrown, double m);
        void                                setLastUpdateWasScale(bool b) { lastUpdateWasScale = b; }
        void                                assignStartingAges(const std::map<Node*,double>& minAges, double unit);
        Node*                               insertFossilTip(Node* hostChild, std::string name, double y, double z);
        void                                setAgeFloors(const std::map<Node*,double>& f) { ageFloors = f; }
        void                                print(void);
        void                                print(std::string header);
#if 0
        std::pair<Node*,Node*>              randomlyChooseBranch(void);
#endif
        void                                reindexNodes(void);
        
    private:
        Node*                               addNode(void);
        void                                clone(const Tree& t);
#if 0
        void                                collapseNode(Node* n);
#endif
        void                                deleteNodes(void);
#if 0
        double                              branchLengthFromMap(Node* e1, Node* e2);
#endif
        void                                initializeBranchLengthKey(std::pair<Node*,Node*>& key, Node* e1, Node* e2);
        std::vector<std::string>            parseNewickString(std::string);
        void                                passDown(Node* p, Node* from);
#if 0
        void                                removeAllBranches(void);
        void                                removeBranch(Node* e1, Node* e2);
        void                                rSPR(void);
        void                                rSPR(std::string s); // for a specifc tip
        void                                reroot(Node* r);
#endif
#if 0
        double                              roundDecimal(double value, int n);
#endif
        void                                setBranch(Node* e1, Node* e2, double x);
        void                                showNode(Node* p, int indent);
        double                              updateNodeAge(void);
        double                              updateCrownAge(double scaleLambda);
        void                                writeTree(Node* p, std::stringstream& strm);
        BranchLengths                       branchLengths;
        std::map<Node*,double>              ageFloors;
        std::vector<Node*>                  downPassSequence;
        std::vector<Node*>                  nodes;
        std::vector<Node*>                  tips;
        Node*                               freeNode;
        Node*                               crown;
        double                              treeHeight;
        int                                 numTaxa;
        int                                 numInternalNodes;
        bool                                lastUpdateWasScale = false;
};

#endif
