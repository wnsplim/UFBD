#ifndef FBDTreeModel_hpp
#define FBDTreeModel_hpp

#include "Parameter.hpp"
#include "ParameterBranchRates.hpp"
#include "ParameterTree.hpp"
#include "ParameterUnresolvedFossils.hpp"
#include "PhylogeneticModel.hpp"
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
        bool                        hasOrigin(void) { return originAge != nullptr; }
        double                      getOriginAgeValue(void);
        void                        setupNodeAgeFloors(void);
    private:
        double                      calculateFBDProbability(void);
        double                      calculateResolvedFBD(void);
        int                         countResolvedSA(void);
        double                      lnD(double t);
        double                      fossilPqLn(double y, double z);
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
        double                      calculateP0Hat(double t);
        double                      computeGamma(double z, int i);
        void                        buildEulerIndex(void);
        bool                        inSub(Node* node, Node* subtreeCrown);
        void                        updateGammaCache(void);
        void                        computeAgeFloors(std::map<Node*,double>& floors);
        void                        resolveFossils(Tree* t, std::vector<Clade>& clades, std::vector<Fossil>& fossils);
        void                        enumerateFossilHosts(Tree* t, Node* crown, Node* origin, bool isCrown, double y, std::vector<Node*>& hosts, std::vector<double>& los, std::vector<double>& his);
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
        double                      doTurnoverMove(void);
        double                      cladeBackboneLineages(int i, double z);
        std::vector<ParameterDouble*>* pickIidRateVector(void);
        void                        enumeratePrunableRoots(Tree* t, std::vector<Node*>& roots);
        void                        enumerateSubtreeHosts(Tree* t, std::vector<Node*>& crowns, std::vector<char>& isCrowns, std::vector<Node*>& origins, double rAge, double ceilingS, std::vector<Node*>& hosts, std::vector<double>& los, std::vector<double>& his);
        bool                        subtreeFossilsValidAt(Tree* t, Node* s, Node* g);
        bool                        subtreeAllFossil(Node* n);
        double                      lambdaAt(int i);
        double                      muAt(int i);
        double                      psiAt(int i);
        //ordered by memory footprint
        std::vector<ParameterDouble*> lambda;
        std::vector<ParameterDouble*> mu;
        std::vector<ParameterDouble*> psi;
        ParameterShrinkageField*    lambdaField;
        ParameterShrinkageField*    muField;
        ParameterShrinkageField*    psiField;
        std::vector<double>           intervalStart;
        std::vector<int>            lambdaIdx;
        std::vector<int>            muIdx;
        std::vector<int>            psiIdx;
        ParameterDouble*            originAge;
        ParameterTree*              parameterTree;
        ParameterUnresolvedFossils* unresolvedFossils;
        std::vector<std::string>    fossilName;
        std::vector<Node*>          fossilCrown;
        std::vector<Node*>          fossilOrigin;
        std::vector<bool>           fossilIsCrown;
        std::vector<double>         fossilY;
        double                      rho;
        std::vector<double>         c1Vec;
        std::vector<double>         c2Vec;
        std::vector<double>         ePrev;
        std::vector<double>         lnDPrev;
        std::vector<double>         c1HatVec;
        std::vector<double>         c2HatVec;
        std::vector<double>         ePrevHat;
        double                      rhoVal;
        bool                        lastWasJointScale;
        bool                        lastWasUpDown;
        bool                        lastWasRateVec;
        bool                        lastRateVecScale;
        std::vector<ParameterDouble*>* lastRateVec;
        double                      rateVecStep;
        double                      shrinkStep;
        long                        rvAccW;
        long                        rvAttW;
        long                        seAccW;
        long                        seAttW;
        bool                        lastWasFbdRate;
        enum TreeMove { TM_NONE = -1, TM_NE, TM_WB, TM_WE, TM_TREESCALE, TM_SARJ, TM_UPDOWN, TM_JOINTSCALE, TM_SUBTREE, TM_NODEAGE, TM_CROWN, TM_COUNT };
        int                         lastTreeMove;
        long                        tmAcc[TM_COUNT];
        long                        tmAtt[TM_COUNT];
        double                      turnoverStep;
        long                        frAccW;
        long                        frAttW;
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
        bool                        isResolved;
        std::vector<int>            subPre;
        std::vector<int>            subSize;
        std::vector<Node*>          nodesByPre;
        bool                        eulerBuilt = false;
        std::vector<double>         sortedYounger;
        std::vector<double>         sortedOlder;
        std::vector<double>         sortedFossilY;
        std::vector<double>         sortedFossilZ;
        std::vector<std::pair<double,double>> sortedZY;
        int                         wholeTreeTotalFast = -1;
        int                         fastIsCrown = 0;
        struct CladeGammaIndex {
            std::vector<double> subY;
            std::vector<std::pair<double,double>> subZY;
            std::vector<double> totY, crY;
            std::vector<std::pair<double,double>> totZY, crZY;
        };
        std::map<Node*,CladeGammaIndex> cladeGamma;
        std::vector<Node*>          activeClades;
        int                         multiCladeFast = -1;
};


#endif
