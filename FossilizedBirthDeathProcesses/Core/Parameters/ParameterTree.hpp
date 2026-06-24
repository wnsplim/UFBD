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
        virtual bool                getAdaptiveProposalActive(void) { return false; }
        Tree*                       getTree(void) { return trees[0]; }
        double                      lnProbability(void);
        void                        print(void);
        void                        setTree(Tree* t);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        
    protected:
        void                        tuneScale(bool accepted);
        //ordered by menmory footprint
        Tree*                       trees[2];
        double                      cachedLnP;
        bool                        useCachedLnP;
        int                         numAcceptances;
        int                         numRejections;
        double                      scaleLambda;
        std::deque<bool>            recentScaleAcceptRej;
        int                         numScaleMoves;
};

#endif
