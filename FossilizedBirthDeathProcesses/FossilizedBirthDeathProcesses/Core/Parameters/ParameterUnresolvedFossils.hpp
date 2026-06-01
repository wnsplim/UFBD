#ifndef ParameterUnresolvedFossils_hpp
#define ParameterUnresolvedFossils_hpp

#include <string>
#include <vector>

#include "FBDInput.hpp"
#include "Parameter.hpp"

class Node;
class Tree;

class ParameterUnresolvedFossils : public Parameter {

    public:
                                    ParameterUnresolvedFossils(void) = delete;
                                    ParameterUnresolvedFossils(double prob, PhylogeneticModel* m, Tree* backbone, std::vector<Clade>& clades, std::vector<Fossil>& fossils);
        double                      getAcceptanceRatio(void) { return ((double)numAcceptances) / ((double)numAcceptances + (double)numRejections); }
        int                         getNumFossils(void) { return numFossils; }
        double                      getFossilAge(int i) { return y[0][i]; }
        double                      getAttachAge(int i) { return z[0][i]; }
        bool                        isSampledAncestor(int i) { return sa[0][i] != 0; }
        Node*                       getCrownNode(int i) { return crownNode[i]; }
        bool                        getIsCrown(int i) { return isCrown[i]; }
        double                      getMinAttachAge(int i);
        double                      getMaxAttachAge(int i);
        double                      lnProbability(void) { return 0.0; } // density is in the model's lnLikelihood, not a prior here
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);

    private:
        double                      updateFossilAge(int i);
        double                      updateAttachAge(int i);
        Tree*                       backbone;
        int                         numFossils;
        std::vector<double>         yMin;
        std::vector<double>         yMax;
        std::vector<Node*>          crownNode;
        std::vector<Node*>          originNode;
        std::vector<bool>           isCrown;
        std::vector<double>         y[2]; // 0 = working (scored), 1 = last accepted
        std::vector<double>         z[2];
        std::vector<int>            sa[2];
        int                         lastFossil;
        bool                        lastWasAttach;
        int                         numAcceptances;
        int                         numRejections;
};

#endif
