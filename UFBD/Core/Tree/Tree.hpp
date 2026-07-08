#ifndef Tree_hpp
#define Tree_hpp

#include <iosfwd>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

class Node;
class RandomVariable;

class Tree {

    public:
                                            Tree(void) = delete;
                                            Tree(const Tree& t); //copy constructor
                                            Tree(std::string newick);
                                           ~Tree(void);
        Tree&                               operator=(const Tree& t);
        std::vector<Node*>                  getAllDescendants(Node* n);
        double                              getBranchLength(Node* e1, Node* e2);
        std::vector<Node*>&                 getDownPassSequence(void) { return downPassSequence; }
        std::string                         getNewickString(void);
        std::string                         getBackboneNewickString(bool keepFossils = false);
        std::vector<Node*>                  getBackboneAgeNodes(void);
        std::vector<Node*>                  getAllAgeNodes(void);
        std::string                         getSummaryNewickString(const std::map<Node*,double>& age, const std::map<Node*,std::pair<double,double>>& hpd, bool keepFossils = false, double y0 = 0.0);
        const std::vector<Node*>&           getBackboneRateNodes(void);
        void                                ensureBackboneCache(void);
        Node*                               getBackboneRoot(void);
        Node*                               getBackboneParent(Node* n);
        const std::vector<Node*>&           getBackboneChildren(Node* n);
        bool                                isBackboneNode(Node* n);
        int                                 getNumLineagesAtTime(double t);
        bool                                isSATip(Node* n);
        bool                                isFakeSplit(Node* n);
        int                                 getNumNodes(void) { return (int)nodes.size(); }
        Node*                               getNodeByOffset(int o) { return nodes[o]; }
        int                                 getNumBackbone(void);
        int                                 getNumExtant(void);
        void                                liftInternalAgesAboveChildren(void);
        Node*                               getCrown(void) { return crown; }
        Node*                               getRoot(void) { return origin != nullptr ? origin : crown; }
        Node*                               getOrigin(void) { return origin; }
        Node*                               addOriginNode(double x0);
        bool                                getLastUpdateWasScale(void) { return lastUpdateWasScale; }
        bool                                isBinary(void);
        void                                validateBackbone(void);
        static bool                         isValidNewick(const std::string& s);
        Node*                               getTaxonNode(std::string name);
        Node*                               getMRCA(const std::vector<std::string>& taxonNames);
        void                                initializeDownPassSequence(void);
        void                                initializeTimes(void); //starting node ages from topology; extant tips at 0
        double                              update(double scaleLambda);
        int                                 scaleInternalAges(double m);
        int                                 scaleSubtreeAges(Node* subtreeCrown, double m);
        double                              updateNodeAge(void);
        double                              updateNodeAgeOnNode(Node* n);
        std::vector<Node*>                  getInternalAgeNodes(void);
        std::vector<Node*>                  getFossilTipAgeNodes(void);
        double                              updateFossilTipAge(Node* n);
        void                                setLastUpdateWasScale(bool b) { lastUpdateWasScale = b; }
        void                                assignStartingAges(const std::map<Node*,double>& minAges, double unit);
        Node*                               insertFossilTip(Node* hostChild, std::string name, double y, double z);
        void                                setAgeFloors(const std::map<Node*,double>& f) { ageFloors = f; }
        void                                print(void);
        void                                print(std::string header);
        void                                reindexNodes(void);
        void                                writeState(std::ostream& os);
        void                                readState(std::istream& is);
        
    private:
        Node*                               addNode(void);
        void                                buildBackboneCache(void);
        void                                clone(const Tree& t);
        void                                deleteNodes(void);
        std::vector<std::string>            parseNewickString(std::string);
        void                                passDown(Node* p, Node* from);
        void                                showNode(Node* p, int indent);
        double                              updateCrownAge(double scaleLambda);
        void                                writeTree(Node* p, std::stringstream& strm);
        std::map<Node*,double>              ageFloors;
        std::vector<Node*>                  downPassSequence;
        std::vector<Node*>                  nodes;
        Node*                               crown;
        Node*                               origin = nullptr;
        int                                 numTaxa;
        bool                                lastUpdateWasScale = false;
        bool                                bbCacheValid = false;
        Node*                               bbRootNode = nullptr;
        std::vector<Node*>                  bbRateNodes;
        std::vector<Node*>                  bbParentByOffset;
        std::vector<std::vector<Node*>>     bbChildrenByOffset;
};

bool writeSummaryTree(Tree* tree, const std::vector<std::string>& names, const std::vector<std::vector<double>>& cols, const std::vector<std::string>& latentNames, const std::vector<std::vector<double>>& latentCols, double burninFrac, const std::string& path, bool keepFossils = false);

#endif
