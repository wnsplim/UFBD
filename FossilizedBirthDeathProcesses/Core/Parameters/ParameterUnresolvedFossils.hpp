#ifndef ParameterUnresolvedFossils_hpp
#define ParameterUnresolvedFossils_hpp

#include <string>
#include <vector>

#include "FBDInput.hpp"
#include "Parameter.hpp"

class Node;
class ParameterDouble;
class Tree;

class ParameterUnresolvedFossils : public Parameter {

    public:
                                    ParameterUnresolvedFossils(void) = delete;
                                    ParameterUnresolvedFossils(double prob, PhylogeneticModel* m, Tree* backbone, std::vector<Clade>& clades, std::vector<Fossil>& fossils, ParameterDouble* originAge);
        double                      getAcceptanceRatio(void) { return ((double)numAcceptances) / ((double)numAcceptances + (double)numRejections); }
        int                         getNumFossils(void) { return numFossils; }
        double                      getFossilAge(int i) { return y[0][i]; }
        double                      getAttachAge(int i) { return z[0][i]; }
        void                        setAttachAge(int i, double zv) { z[0][i] = zv; z[1][i] = zv; }
        bool                        isUE(int i) { return ue[i]; }
        bool                        isSA(int i) { return ue[i] == false && z[0][i] == y[0][i]; }
        int                         getNumSampledAncestors(void) { int n = 0; for(int i = 0; i < numFossils; i++) if(ue[i] == false && z[0][i] == y[0][i]) n++; return n; }
        int                         getSpineIdx(void) { return spineIdx; }
        void                        syncSpine(double x0) { if(spineIdx >= 0){ z[0][spineIdx] = x0; z[1][spineIdx] = x0; } }
        Node*                       getCrownNode(int i) { return crownNode[i]; }
        Node*                       getMaxAttachNode(int i) { return isCrown[i] ? crownNode[i] : originNode[i]; }
        bool                        getIsCrown(int i) { return isCrown[i]; }
        double                      getMinAttachAge(int i);
        double                      getMaxAttachAge(int i);
        double                      lnProbability(void) { return 0.0; } // density is in the model's lnLikelihood, not a prior here
        void                        print(void);
        double                      update(void);
        double                      scaleAllAttachAges(double m);
        double                      scaleAttachAges(const std::vector<int>& indices, double m);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);

    private:
        double                      updateFossilAge(int i);
        double                      updateAttachAge(int i);
        double                      updateSampledAncestor(int i);
        Tree*                       backbone;
        ParameterDouble*            originAge;
        int                         numFossils;
        std::vector<double>         yMin;
        std::vector<double>         yMax;
        std::vector<Node*>          crownNode;
        std::vector<Node*>          originNode;
        std::vector<bool>           isCrown;
        std::vector<char>           ue;
        std::vector<double>         y[2]; // 0 = working (scored), 1 = last accepted
        std::vector<double>         z[2];
        int                         lastFossil;
        int                         spineIdx;
        bool                        lastWasBulk;
        int                         numAcceptances;
        int                         numRejections;
};

#endif
