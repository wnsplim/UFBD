#ifndef ParameterTree_hpp
#define ParameterTree_hpp

#include <deque>
#include <vector>

#include "Parameter.hpp"
#include "Tree.hpp"

class Node;

class ParameterTree : public Parameter {

    public:
                                    ParameterTree(void) = delete;
                                    ParameterTree(double prob, PhylogeneticModel* m);
        double                      getAcceptanceRatio(void) { return ((double) numAcceptances) /( (double)numAcceptances + (double)numRejections ) ;}
        double                      getScaleLambda(void) { return scaleLambda; }
        double                      getNodeAgeStep(void) { return nodeAgeStep; }
        void                        recordNodeAgeMove(bool accepted);
        Tree*                       getTree(void) { return trees[0]; }
        void                        setAgeFloors(const std::map<Node*,double>& f);
        double                      lnProbability(void);
        void                        print(void);
        void                        setTree(Tree* t);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);

    protected:
        void                        tuneScale(bool accepted);
        //ordered by menmory footprint
        Tree*                       trees[2];
        int                         numAcceptances;
        int                         numRejections;
        double                      scaleLambda;
        std::deque<bool>            recentScaleAcceptRej;
        int                         numScaleMoves;
        double                      nodeAgeStep = 0.5;
        long                        naAccW = 0;
        long                        naAttW = 0;
        long                        naAdapt = 0;
};

#endif
