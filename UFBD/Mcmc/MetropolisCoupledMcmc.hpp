#ifndef MetropolisCoupledMcmc_hpp
#define MetropolisCoupledMcmc_hpp

#include <fstream>
#include <string>
#include <vector>

#include "ChainRunner.hpp"
#include "RandomVariable.hpp"
#include "ThreadPool.hpp"
#include "WriteTSV.hpp"

class PhylogeneticModel;

class MetropolisCoupledMcmc : public ChainRunner {
    public:
                                            MetropolisCoupledMcmc(void) = delete;
                                            MetropolisCoupledMcmc(unsigned long ng, int thin, std::vector<PhylogeneticModel*> m, unsigned int masterSeed);
                                            ~MetropolisCoupledMcmc(void);
        void                                init(void);
        void                                advance(unsigned long nGens);
        void                                finalize(void);
        bool                                treeHasFossils(void);
        void                                printMoveDiagnostics(int rep);
        void                                writeCheckpoint(void);
        bool                                loadCheckpoint(void);
        void                                setVerbose(bool b) { verbose = b; }
        void                                setLabel(int i) { runLabel = i; }
        Tree*                               getTree(void);

    protected:
        RandomVariable*                     resumeRng(void);
        std::vector<std::string>            resumeParameterNames(void);
        std::vector<std::string>            resumeLatentNames(void);

    private:
        //functions
        double                              calcHeating(int idx);
        void                                sample(unsigned long n);
        void                                writeColdSample(unsigned long n);
        void                                adaptDeltaT(bool swapAccepted);
        //objects ordered by memory footprint
        std::deque<bool>                    recentAcceptRej;
        std::vector<PhylogeneticModel*>     models;
        std::vector<double>                 currLnL;
        std::vector<double>                 newLnL;
        std::vector<double>                 currLnP;
        std::vector<double>                 newLnP;
        std::vector<double>                 lnProposalRatio;
        std::vector<double>                 lnLikelihoodRatio;
        std::vector<double>                 lnPriorRatio;
        std::vector<double>                 lnAcceptanceProbabilities;
        std::vector<int>                    indices;
        ThreadPool                          threadPool;
        std::vector<ThreadPool*>            chainPools;
        long                                chainCalibGen;
        double                              chainSeqT;
        double                              chainParT;
        int                                 chainDecision;
        RandomVariable                      swapRng;
        bool                                verbose = false;
        bool                                tuning = false;
        int                                 runLabel = 0;
        unsigned long                       numCycles;
        double                              deltaT;
        unsigned long                       resampleEvery;
        long                                swapAdaptCount;
        int                                 swapAdaptAcc;
        int                                 swapAdaptAtt;
        int                                 coldModelIdx;
        int                                 numModels;
        int                                 thinning;
};

#endif
