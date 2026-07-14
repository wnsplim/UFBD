#ifndef PhylogeneticModel_hpp
#define PhylogeneticModel_hpp

#include <iosfwd>
#include <string>
#include <vector>

#include "RandomVariable.hpp"

class Parameter;
class Tree;

class PhylogeneticModel {

    public:
                                                PhylogeneticModel(void);
                                               ~PhylogeneticModel(void);
        virtual std::vector<std::string>        getParameterNames(void) = 0;
        virtual std::vector<double>             getParameterString(void) = 0;
        virtual std::vector<std::string>        getLatentNames(void) { return {}; }
        virtual std::vector<double>             getLatentString(void) { return {}; }
        virtual bool                            treeIncludesFossils(void) { return false; }
        Tree*                                   getTree(void);
        RandomVariable*                         getRng(void) { return &rng; }
        const Parameter*                        getUpdatedParameter(void) { return updatedParameter; };
        virtual double                          lnLikelihood(void) = 0;
        virtual void                            invalidateLikelihoodCache(void) {}
        virtual void                            invalidatePriorCache(void) {}
        virtual int                             getLastMoveType(void) { return -1; }
        virtual double                          lnPriorProbability(void) = 0;
        virtual void                            print(void) = 0;
        virtual double                          update(void) = 0;
        virtual void                            updateForAcceptance(void) = 0;
        virtual void                            updateForRejection(void) = 0;
        virtual void                            writeState(std::ostream& os) {}
        virtual void                            readState(std::istream& is) {}
        std::vector<Parameter*>                 parameters;
        Parameter*                              updatedParameter;
        RandomVariable                          rng;
};

#endif
