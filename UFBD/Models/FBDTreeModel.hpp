#ifndef FBDTreeModel_hpp
#define FBDTreeModel_hpp

#include "Parameter.hpp"
#include "ParameterBranchRates.hpp"
#include "ParameterTree.hpp"
#include "ParameterUnresolvedFossils.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "Tree.hpp"

#include <deque>
#include <map>
#include <vector>

class Parameter;
class ParameterDouble;
class ParameterShrinkageField;

class FBDTreeModel : public PhylogeneticModel {

    public:
                                    FBDTreeModel(void) = delete;
                                    FBDTreeModel(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils, unsigned int seed);
        std::vector<std::string>    getParameterNames(void);
        std::vector<double>         getParameterString(void);
        std::vector<bool>           getParameterFixedMask(void);
        std::vector<std::string>    getLatentNames(void);
        std::vector<double>         getLatentString(void);
        bool                        treeIncludesFossils(void) { return isResolved == false && backboneFossils.empty() == false; }
        std::vector<Node*>          getAgeLogNodes(void) { return isResolved ? parameterTree->getTree()->getBackboneAgeNodes() : parameterTree->getTree()->getAllAgeNodes(); }
        double                      lnLikelihood(void);
        double                      lnPriorProbability(void);
        void                        print(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);
        ParameterTree*              getParameterTree(void) { return parameterTree; }
        ParameterUnresolvedFossils* getUnresolvedFossils(void) { return unresolvedFossils; }
        std::string                 getRateMap(void);
        bool                        hasOrigin(void) { return originAge != nullptr; }
        double                      getOriginAgeValue(void);
        void                        setupNodeAgeFloors(void);
        double                      fossilSweep(void);
    private:
        double                      fossilTermLn(int i, int spineIdx);
        double                      term3Sum(void);
        double                      fossilSweepSequential(void);
        double                      fossilSweepParallel(void);
        void                        rebuildZoneStalks(int z, const std::vector<int>& fos);
        void                        zoneRecomputeGamma(const std::vector<int>& fos);
        double                      zoneTermSum(const std::vector<int>& fos, double x0, bool useOrigin);
        double                      calculateFBDProbability(void);
        double                      calculateResolvedFBD(void);
        int                         countResolvedSA(void);
        double                      lnD(double t);
        double                      fossilPqLn(double y, double z, int type);
        double                      uePqLn(double z);
        double                      calculateLnSurvival(double t);
        double                      calculateLnAnySample(double t);
        double                      calculateLnConditioning(double t);
        void                        prepareIntervals(void);
        int                         findIndex(double t);
        double                      calculateLnQtAt(int i, double t);
        double                      calculateP0At(int i, double t);
        double                      calculateP0(double t);
        double                      calculateP0HatAt(int i, double t);
        double                      computeGamma(double z, int i);
        static double               lnChoose(int n, int k);
        void                        buildEulerIndex(void);
        bool                        inSub(Node* node, Node* subtreeCrown);
        void                        buildZoneIndex(void);
        double                      zoneBackboneEdges(int mz, double z);
        void                        rebuildStalkIndex(void);
        double                      doAttachmentZoneGibbs(void);
        double                      doAttachmentZoneJump(int i);
        double                      validZoneSet(int i, int a, std::vector<std::pair<double,double> >& iv);
        void                        updateGammaCache(void);
        void                        computeAgeFloors(std::map<Node*,double>& floors);
        void                        resolveFossils(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils);
        void                        enumerateFossilAttachEdges(Tree* t, Node* crown, Node* origin, bool isCrown, double y, std::vector<Node*>& attachEdges, std::vector<double>& los, std::vector<double>& his);
        int                         fossilIndexByName(const std::string& nm);
        double                      doWilsonBalding(void);
        double                      doNarrowExchange(void);
        double                      doWideExchange(void);
        double                      doTreeScale(void);
        double                      doSARJMCMC(void);
        double                      doUpDownScale(void);
        double                      doJointScale(void);
        double                      doSubtreeScale(void);
        double                      doRateVectorScale(void);
        double                      doRateShrinkExpand(void);
        double                      doRateShift(void);
        std::vector<ParameterDouble*>* pickIidRateVector(void);
        void                        enumeratePrunableRoots(Tree* t, std::vector<Node*>& roots);
        void                        enumerateSubtreeAttachEdges(Tree* t, std::vector<Node*>& crowns, std::vector<char>& isCrowns, std::vector<Node*>& origins, double rAge, double ceilingS, std::vector<Node*>& attachEdges, std::vector<double>& los, std::vector<double>& his);
        bool                        subtreeFossilsValidAt(Tree* t, Node* s, Node* g);
        bool                        subtreeAllFossil(Node* n);
        double                      lambdaAt(int i);
        double                      muAt(int i);
        double                      psiTotalAt(int i);
        double                      psiOfTypeAt(int type, int i);
        std::vector<int>            buildSkylineRates(const std::string& prefix, const std::string& sep, int nBins, const std::vector<double>& times, double topAge, const Probability::PriorSpec& basePrior, double rate0, bool smooth, bool gmrf, const std::vector<int>& groupIds, const std::map<int,Probability::PriorSpec>& groupPrior, double nShifts, double shiftSize, std::vector<ParameterDouble*>& outVec, ParameterShrinkageField*& outField, std::vector<std::string>& outNames);
        void                        appendRateMap(const std::vector<double>& times, const std::vector<int>& binToChunk, const std::vector<std::string>& names);
        //ordered by memory footprint
        std::vector<ParameterDouble*> lambda;
        std::vector<ParameterDouble*> mu;
        std::vector<std::vector<ParameterDouble*>> psi;
        std::vector<std::string>    lambdaName, muName;
        std::vector<std::vector<std::string>> psiName;
        std::vector<std::pair<std::string,std::string>> rateMapRows;
        ParameterShrinkageField*    lambdaField;
        ParameterShrinkageField*    muField;
        std::vector<ParameterShrinkageField*> psiField;
        int                         numPsiTypes;
        int                         numExtantTips;
        std::vector<int>            fossilType;
        std::map<std::string,int>   fossilTypeByName;
        struct BackboneFossil { Node* tip; double yMin; double yMax; int type; };
        std::vector<BackboneFossil> backboneFossils;
        std::vector<std::string>    unrFossilName;
        std::vector<double>           intervalStart;
        std::vector<int>            lambdaIdx;
        std::vector<int>            muIdx;
        std::vector<std::vector<int>> psiIdx;
        ParameterDouble*            originAge;
        ParameterTree*              parameterTree;
        ParameterUnresolvedFossils* unresolvedFossils;
        std::vector<std::string>    fossilName;
        std::vector<Node*>          fossilCrown;
        std::vector<Node*>          fossilOrigin;
        std::vector<bool>           fossilIsCrown;
        double                      rho;
        std::vector<double>         c1Vec;
        std::vector<double>         c2Vec;
        std::vector<double>         ePrev;
        std::vector<double>         lnDPrev;
        std::vector<double>         c1HatVec;
        std::vector<double>         c2HatVec;
        std::vector<double>         ePrevHat;
        double                      rhoVal;
        enum MoveKind { MK_PARAM, MK_AZGIBBS, MK_RATESHIFT, MK_RATEVEC, MK_UPDOWN, MK_JOINTSCALE };
        MoveKind                    lastMoveKind;
        double                      shiftStep;
        int                         saBatch;
        double                      saBatchF;
        long                        rsAccW;
        long                        rsAttW;
        long                        rsAcc;
        long                        rsTot;
        long                        rsAdapt;
        bool                        lastRateVecScale;
        std::vector<ParameterDouble*>* lastRateVec;
        double                      rateVecStep;
        double                      shrinkStep;
        long                        rvAccW;
        long                        rvAttW;
        long                        rvAtt;
        long                        seAccW;
        long                        seAttW;
        long                        seAtt;
        enum TreeMove { TM_NONE = -1, TM_NE, TM_WB, TM_WE, TM_TREESCALE, TM_SARJ, TM_UPDOWN, TM_JOINTSCALE, TM_SUBTREE, TM_NODEAGE, TM_CROWN, TM_COUNT };
        int                         lastTreeMove;
        long                        tmAcc[TM_COUNT];
        long                        tmAtt[TM_COUNT];
        double                      upDownStep;
        int                         upDownTotal;
        std::deque<bool>            upDownRecent;
        std::vector<double>         cachedGammaLn;
        std::vector<char>           gammaStale;
        std::vector<double>         prevY;
        std::vector<double>         prevZ;
        std::vector<int>            prevSA;
        std::vector<double>         prevNodeAge;
        double                      prevX0;
        bool                        cacheInit;
    public:
        void                        invalidateGammaCache(void){ cacheInit = false; }
    private:
        bool                        zoneInit;
        bool                        isResolved;
        std::vector<int>            subPre;
        std::vector<int>            subSize;
        std::vector<Node*>          nodesByPre;
        bool                        eulerBuilt = false;
        int                         numZones = 0;
        std::vector<int>            minZoneOfNode;
        int                         trunkMinZone = -1;
        struct ZoneEdges {
            std::vector<double> yng, old;
        };
        std::vector<ZoneEdges>      zoneEdges;
        struct StalkBucket {
            std::vector<double> y;
            std::vector<std::pair<double,double>> zy;
        };
        std::vector<StalkBucket>    zoneStalks;
        bool                        zoneIndexBuilt = false;
        std::vector<int>            prevAttachmentZone;
        std::vector<int>            azGibbsIdx;
        long                        azAcc;
        long                        azAtt;
};


#endif
