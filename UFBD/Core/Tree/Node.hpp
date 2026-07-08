#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
#include <vector>
class RandomVariable;


class Node {

    public:
                            Node(void);
        void                addNeighbor(Node* p);
        Node*               getAncestor(void) { return ancestor; }
        std::vector<Node*>& getDescendants(void);
        bool                getFlag(void) { return flag; }
        int                 getIndex(void) { return index; }
        bool                getIsFossil(void) { return fossil; }
        bool                getIsTip(void) { return isTip; }
        std::string         getName(void) { return name; }
        std::set<Node*>&    getNeighbors(void) { return neighbors; }
        int                 getOffset(void) { return offset; }
        double              getTime(void) { return time; }
        double              getFossilYMin(void) { return fossilYMin; }
        double              getFossilYMax(void) { return fossilYMax; }
        void                setFossilAgeRange(double lo, double hi) { fossilYMin = lo; fossilYMax = hi; }
        void                removeNeighbor(Node* p);
        void                removeAllNeighbors(void) { neighbors.clear(); }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setFlag(bool tf) { flag = tf; }
        void                setIndex(int x) { index = x; }
        void                setIsFossil(bool x) { fossil = x; }
        void                setIsTip(bool tf) { isTip = tf; }
        void                setName(std::string s) { name = s; }
        void                setOffset(int x) { offset = x; }
        void                setTime(double x) { time = x; }
    
    private:
        std::vector<Node*>  descendantsVector;
        std::set<Node*>     neighbors;
        Node*               ancestor;
        std::string         name;
        double              time;
        double              fossilYMin;
        double              fossilYMax;
        bool                isTip;
        bool                flag;
        bool                fossil;
        int                 index;
        int                 offset;
};

#endif
