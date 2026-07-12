#ifndef ParameterUnresolvedFossils_hpp
#define ParameterUnresolvedFossils_hpp

#include <string>
#include <vector>

#include "FBDInput.hpp"
#include "Parameter.hpp"
#include "ParameterDouble.hpp"

class Node;
class Tree;

class ParameterUnresolvedFossils : public Parameter {

    public:
                                    ParameterUnresolvedFossils(void) = delete;
                                    ParameterUnresolvedFossils(double prob, PhylogeneticModel* m, Tree* backbone, std::vector<Clade>& clades, std::vector<Fossil>& fossils, ParameterDouble* originAge);
        double                      getAcceptanceRatio(void) { return ((double)numAcceptances) / ((double)numAcceptances + (double)numRejections); }
        int                         getNumFossils(void) { return numFossils; }
        double                      getFossilAge(int i) { return y[0][i]; }
        double                      getAttachAge(int i) { return (i == spineIdx) ? originAge->getValue() : z[0][i]; }
        double                      getYMin(int i) { return yMin[i]; }
        double                      getYMax(int i) { return yMax[i]; }
        bool                        isUE(int i) { return ue[i]; }
        bool                        isSA(int i) { return ue[i] == false && getAttachAge(i) == y[0][i]; }
        int                         getNumSampledAncestors(void) { int n = 0; for(int i = 0; i < numFossils; i++) if(isSA(i)) n++; return n; }
        int                         getSpineIdx(void) { return spineIdx; }
        int                         getLandingZone(int i) { return lz[0][i]; }
        std::vector<int>&           getLandingZoneDomain(int i) { return lzDomain[i]; }
        void                        initLandingZone(int i, int v) { lz[0][i] = v; lz[1][i] = v; }
        void                        setLandingZone(int i, int v) { lz[0][i] = v; }
        void                        setAttachAge(int i, double v) { z[0][i] = v; }
        void                        beginLandingZoneMove(int i) { lastFossil = i; lastMove = SINGLE; }
        void                        setLandingZoneDomain(std::vector<std::vector<int> >& domain);
        Node*                       getCrownNode(int i) { return crownNode[i]; }
        Node*                       getMaxAttachNode(int i) { return isCrown[i] ? crownNode[i] : originNode[i]; }
        bool                        getIsCrown(int i) { return isCrown[i]; }
        bool                        getIsStem(int i) { return isStem[i]; }
        double                      getMinAttachAge(int i);
        double                      getMaxAttachAge(int i);
        double                      lnProbability(void) { return 0.0; } // density is in the model's lnLikelihood, not a prior here
        void                        print(void);
        double                      update(void);
        double                      scaleAllAttachAges(double m);
        double                      scaleAttachAges(const std::vector<int>& indices, double m);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);

    private:
        enum MoveKind { SINGLE, FLIP, BULK };
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
        std::vector<bool>           isStem;
        std::vector<char>           ue;
        std::vector<double>         y[2]; // 0 = working (scored), 1 = last accepted
        std::vector<double>         z[2];
        std::vector<int>            lz[2];
        std::vector<std::vector<int> > lzDomain;
        int                         lastFossil;
        int                         spineIdx;
        MoveKind                    lastMove;
        int                         flipS;
        int                         flipT;
        double                      flipSy;
        double                      flipSz;
        double                      flipTy;
        double                      flipTz;
        int                         numAcceptances;
        int                         numRejections;
};

#endif
