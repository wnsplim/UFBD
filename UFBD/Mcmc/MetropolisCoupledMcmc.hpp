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
                                            MetropolisCoupledMcmc(unsigned long ng, int pf, int sf, std::vector<PhylogeneticModel*> m, unsigned int masterSeed);
                                            ~MetropolisCoupledMcmc(void);
        void                                run(void);
        void                                init(void);
        void                                advance(unsigned long nGens);
        void                                finalize(void);
        void                                setOutputPaths(const std::string& po, const std::string& to) { paramOut = po; treeOut = to; }
        const std::vector<std::vector<double>>& traceColumns(void) { return traceCols; }
        const std::vector<std::string>&     traceNames(void) { return traceNms; }
        void                                writeCheckpoint(void);
        bool                                loadCheckpoint(void);
        void                                resumeOutputs(void);
        unsigned long                       currentGen(void) { return gen; }

    private:
        //functions
        double                              calcHeating(int idx);
        void                                sample(unsigned long n);
        void                                updateDeltaT(void);
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
        WriteTSV                            params;
        WriteTSV                            trees;
        std::string                         treeOut;
        std::string                         paramOut;
        bool                                writeTrees;
        unsigned long                       numCycles;
        double                              deltaT;
        int                                 coldModelIdx;
        int                                 numModels;
        int                                 printFrequency;
        int                                 sampleFrequency;
        unsigned long                       gen;
        std::vector<std::vector<double>>    traceCols;
        std::vector<std::string>            traceNms;
};

#endif
