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
        bool                        isSA(int i) { return sa[0][i] != 0; }
        int                         getNumSampledAncestors(void) { int n = 0; for(int i = 0; i < numFossils; i++) if(isSA(i)) n++; return n; }
        int                         getSpineIdx(void) { return spineIdx; }
        int                         getAttachmentZone(int i) { return az[0][i]; }
        std::vector<int>&           getAttachmentZoneDomain(int i) { return azDomain[i]; }
        void                        initAttachmentZone(int i, int v) { az[0][i] = v; az[1][i] = v; }
        void                        setAttachmentZone(int i, int v) { az[0][i] = v; }
        void                        setAttachAge(int i, double v) { z[0][i] = v; }
        void                        beginAttachmentZoneMove(int i) { lastFossil = i; lastMove = SINGLE; }
        void                        setAttachmentZoneDomain(std::vector<std::vector<int> >& domain);
        Node*                       getCrownNode(int i) { return crownNode[i]; }
        Node*                       getMaxAttachNode(int i) { return isCrown[i] ? crownNode[i] : originNode[i]; }
        bool                        getIsCrown(int i) { return isCrown[i]; }
        bool                        getIsStem(int i) { return isStem[i]; }
        double                      getMinAttachAge(int i);
        double                      getMaxAttachAge(int i);
        double                      lnProbability(void) { return 0.0; } // density is in the model's lnLikelihood, not a prior here
        void                        print(void);
        double                      update(void);
        double                      proposeOneFossil(int i);
        bool                        saEligible(int i) { return ue[i] == false && i != spineIdx && getMaxAttachAge(i) > y[0][i]; }
        void                        beginBatchMove(void) { lastMove = BULK; }
        double                      flipSA(int i, bool toSA);
        double                      scaleAllAttachAges(double m);
        double                      scaleAttachAges(const std::vector<int>& indices, double m);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);

    private:
        enum MoveKind { SINGLE, FLIP, BULK };
        enum SubMove { SUB_NOOP, SUB_SA, SUB_Y, SUB_Z, SUB_COUNT };
        SubMove                     lastSub;
        long                        subAcc[SUB_COUNT];
        long                        subAtt[SUB_COUNT];
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
        std::vector<char>           sa[2];
        std::vector<int>            az[2];
        std::vector<std::vector<int> > azDomain;
        int                         lastFossil;
        int                         spineIdx;
        MoveKind                    lastMove;
        int                         flipS;
        int                         flipT;
        double                      flipSy;
        double                      flipSz;
        double                      flipTy;
        double                      flipTz;
        char                        flipSsa;
        char                        flipTsa;
        int                         numAcceptances;
        int                         numRejections;
};

#endif
