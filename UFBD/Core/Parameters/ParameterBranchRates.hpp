#ifndef ParameterBranchRates_hpp
#define ParameterBranchRates_hpp

#include <cmath>
#include <deque>
#include <iosfwd>
#include <vector>

#include "Parameter.hpp"
#include "UserSettings.hpp"

class Tree;
class Node;
class ParameterUnresolvedFossils;
class RandomVariable;

class AdaptiveMixSelector {
    public:
        void init(int nOps);
        int  pick(RandomVariable& rng);
        void record(int op, double jump2, double cpu);
        void writeState(std::ostream& os) const;
        void readState(std::istream& is);
    private:
        std::vector<double> cumJ2;
        std::vector<double> cumCpu;
        std::vector<long>   tries;
        int                 nOps;
};

enum class ClockModel { UCLN, WN, GBM, CIR, GBMC }; // WN + CIR + GBMC: halt — detached, not selectable

class BranchRateModel : public Parameter {

    public:
                                    BranchRateModel(double prob, PhylogeneticModel* m, Tree* tree, int numPartitions, const std::vector<int>& partitionGroup, const std::vector<Sigma2Param>& sigma2ParamPerGroup, const double* rgeneParam, const double* sigma2Param);
        double                      getAcceptanceRatio(void);
        int                         getNumPartitions(void) { return numPartitions; }
        int                         getNumClockGroups(void) { return numClockGroups; }
        int                         getFirstPartitionOfGroup(int g) { for(int p = 0; p < numPartitions; p++) if(groupOf[p] == g) return p; return 0; }
        double                      getClockGroupRate(int g) { return mu[0][g]; }
        double                      getClockGroupSigma2(int g) { return sigma2[0][g]; }
        int                         getNumBranchNodes(void) { return (int)branchNodes.size(); }
        int                         getBranchNodeOffset(int i) { return branchNodes[i]; }
        void                        scaleAll(double sf);
        void                        commitAll(void);
        void                        restoreAll(void);
        double                      constantDistanceMove(void);
        Node*                       getCdNode(void) { return cdMovedNode; }
        double                      getCdOldAge(void) { return cdOldAge; }
        double                      rateAgeSubtreeMove(void);
        double                      simpleDistanceMove(void);
        double                      smallPulleyMove(void);
        void                        setUnresolvedFossils(ParameterUnresolvedFossils* u) { uf = u; }
        void                        print(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);
        void                        writeState(std::ostream& os);
        void                        readState(std::istream& is);
        void                        freezePncp(void);
        void                        setChainLabel(int c) { chainLabel = c; }
        void                        setPncpReporter(bool b) { pncpReporter = b; }
        virtual std::vector<std::vector<double>> getAbsoluteRates(void) = 0;

    protected:
        double                      scalePartitionRate(int p);
        double                      scalePartitionSigma2(int p);
        double                      scaleBranchRate(int p, int b);
        double                      globalRateBranchRatesScale(int p);
        double                      bactrianMultiplier(int moveType, int p);
        void                        adaptStep(int moveType, int p, bool accepted);
        void                        adaptNcStep(int p, bool accepted);
        double                      groupedGammaDirichletLnP(const std::vector<double>& v, const double* param);
        double                      gammaLnPdf(double a, double b, double x);
        Tree*                       tree;
        int                         numPartitions;
        int                         numClockGroups;
        std::vector<int>            groupOf;
        std::vector<Sigma2Param>    sigma2ParamByGroup;
        int                         numNodes;
        std::vector<int>            branchNodes;
        double                      rgeneParam[3];
        double                      sigma2Param[3];
        std::vector<double>         mu[2];
        std::vector<double>         sigma2[2];
        std::vector<std::vector<double>> rate[2];
        std::vector<double>         step[4];
        std::vector<int>            acc[4];
        std::vector<int>            rej[4];
        int                         lastMove;
        int                         lastPartition;
        int                         lastNode;
        ParameterUnresolvedFossils* uf;
        std::vector<double>         cdStepNode;
        std::vector<long>           cdAccNode;
        std::vector<long>           cdAttNode;
        std::vector<long>           cdTotNode;
        int                         lastCdNode;
        Node*                       cdMovedNode = nullptr;
        double                      cdOldAge = 0.0;
        std::vector<int>            cdNodes;
        double                      sdStep;
        long                        sdAccW;
        long                        sdAttW;
        double                      spStep;
        long                        spAccW;
        long                        spAttW;
        long                        sdAcc;
        long                        sdAtt;
        long                        spAcc;
        long                        spAtt;
        std::vector<double>         ncStep;
        std::vector<long>           ncAtt;
        std::vector<long>           ncAcc;
        long                        cdAccTot = 0;
        long                        cdAttTot = 0;
        long                        rasAcc = 0;
        long                        rasAtt = 0;
        double                      rasStep = 0.5;
        long                        rasAccW = 0;
        long                        rasAttW = 0;
        std::vector<std::vector<double>> sigTauL;
        std::vector<std::vector<double>> sigEllB;
        std::vector<long>           sigRefresh;
        std::vector<long>           sigCount;
        std::vector<double>         centeredness;
        std::vector<double>         pncpMeanTau;
        bool                        pncpFrozen = false;
        bool                        pncpReporter = true;
        int                         chainLabel = 0;
};

class ParameterBranchRates : public BranchRateModel {

    public:
                                    ParameterBranchRates(void) = delete;
                                    ParameterBranchRates(double prob, PhylogeneticModel* m, Tree* tree, int numPartitions, const std::vector<int>& partitionGroup, const std::vector<ClockModel>& clockModel, const std::vector<Sigma2Param>& sigma2ParamPerGroup, const double* rgeneParam, const double* sigma2Param);
        std::vector<std::vector<double>> getAbsoluteRates(void);
        double                      lnProbability(void);
        double                      update(void);
        double                      updatePncpPartition(int p);
        int                         getLastPncpPartition(void) { return lastMove == 8 ? lastPartition : -1; }
        ClockModel                  getClockModel(int g) { return clockModelOf[g]; }
        bool                        allGroupsAre(ClockModel c) { for(ClockModel m : clockModelOf) if(m != c) return false; return true; }

    private:
        double                      lognormalLnP(double r, double s2, double m);
        double                      whiteNoiseLnP(double r, double s2, double t, double m);
        double                      gbmLnP(void);
        double                      gbmContinuousLnP(void);
        double                      sigmaPncpMove(int p);
        double                      sigmaPncpMoveGBMC(int p);
        double                      sigmaPncpMoveGBM(int p);
        double                      sigmaPncpMoveWN(int p);
        void                        branchLikePrecision(int p, std::vector<double>& tauL, std::vector<double>& ellB);
        std::vector<ClockModel>     clockModelOf;
};

// CIR clock: halt — detached dead code (kept, never constructed)
class ParameterBranchRatesCIR : public BranchRateModel {

    public:
                                    ParameterBranchRatesCIR(void) = delete;
                                    ParameterBranchRatesCIR(double prob, PhylogeneticModel* m, Tree* tree, int numPartitions, const std::vector<int>& partitionGroup, const std::vector<Sigma2Param>& sigma2ParamPerGroup, const double* rgeneParam, const double* sigma2Param);
        std::vector<std::vector<double>> getAbsoluteRates(void);
        double                      lnProbability(void);
        double                      update(void);
        void                        updateForAcceptance(void);
        void                        updateForRejection(void);

    private:
        double                      scaleClockPartitionTheta(int p);
        double                      scaleMuRates(int p);
        double                      sigmaPncpMoveCIR(int p);
        double                      cirLnP(void);
        double                      getMeanTau(double rho, double rhoUp, double t, double sigma, double theta);
        double                      besselIRatio(double nu, double x);
        double                      thetaParam[3];
        std::vector<double>         theta[2];
};

#endif
