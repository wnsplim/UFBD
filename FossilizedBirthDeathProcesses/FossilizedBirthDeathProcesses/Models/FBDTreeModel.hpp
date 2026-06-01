#ifndef FBDTreeModel_hpp
#define FBDTreeModel_hpp

#include "Parameter.hpp"
#include "ParameterTree.hpp"
#include "ParameterUnresolvedFossils.hpp"
#include "PhylogeneticModel.hpp"
#include "Tree.hpp"

#include <map>
#include <vector>

class Parameter;
class ParameterDouble;

class FBDTreeModel : public PhylogeneticModel {

    public:
                                    FBDTreeModel(void) = delete;
                                    FBDTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, unsigned int seed);
        std::vector<std::string>    getParameterNames(void);
        std::vector<double>         getParameterString(void);
        double                      lnLikelihood(void);
        double                      lnPriorProbability(void);
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
    private:
        double                      calculateFBDProbability(void);
        double                      computeGamma(double z, int i);
        void                        computeAgeFloors(std::map<Node*,double>& floors);
        double                      doJointScale(void);
        void                        calculateCs(void);
        void                        calculateC1(void);
        void                        calculateC2(void);
        double                      calculateQt(double t);
        double                      calculatePo(double t);
        double                      calculatePoHat(double t);
        //ordered by memory footprint
        ParameterDouble*            lambda;
        ParameterDouble*            mu;
        ParameterDouble*            psi;
        ParameterTree*              parameterTree;
        ParameterUnresolvedFossils* unresolvedFossils;
        double                      rho;
        double                      c1;
        double                      c2;
        double                      lambdaVal;
        double                      muVal;
        double                      rhoVal;
        double                      psiVal;
        bool                        lastWasJointScale;
};


#endif
