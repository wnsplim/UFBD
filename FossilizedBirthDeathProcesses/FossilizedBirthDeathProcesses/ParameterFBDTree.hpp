#ifndef ParameterFBDTree_hpp
#define ParameterFBDTree_hpp

#include "Parameter.hpp"
#include "ParameterTree.hpp"
#include "Tree.hpp"

#include <vector>

class Parameter;

class ParameterFBDTree : protected ParameterTree {

    public:
                                    ParameterFBDTree(void) = delete;
                                    ParameterFBDTree(double prob, PhylogeneticModel* m, Tree* t);
        bool                        getAdaptiveProposalActive(void) override;
        double                      lnProbability(void) override;
        void                        print(void) override;
        double                      update(void) override;
        void                        updateForAcceptance(void) override;
        void                        updateForRejection(void) override;
        
    private:
        double                      calculateFBDProbability(void);
        void                        calculateCs(void);
        void                        calculateC1(void);
        void                        calculateC2(void);
        double                      calculateQt(double t);
        double                      calculatePo(double t);
        double                      calculatePoHat(double t);
        //ordered by menmory footprint
        std::vector<Parameter*>     parameters;
        Parameter*                  updatedParameter; //by convention, nullptr indicates topology/branch length move (i.e., ParameterTree::update())
        ParameterDouble*            lambda;
        ParameterDouble*            mu;
        ParameterDouble*            psi;
        ParameterDouble*            rho;
        double                      c1;
        double                      c2;
        double                      lambdaVal;
        double                      muVal;
        double                      rhoVal;
        double                      psiVal;
};


#endif
