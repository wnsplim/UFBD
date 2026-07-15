#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <string>
#include <vector>

#include "ChainRunner.hpp"
#include "WriteTSV.hpp"
class PhylogeneticModel;
class RandomVariable;
class Tree;

class Mcmc : public ChainRunner {
    public:
                                Mcmc(void) = delete;
                                Mcmc(int ng, int thin, PhylogeneticModel* m);
        void                    init(void);
        void                    advance(unsigned long nGens);
        void                    finalize(void);
        bool                    treeHasFossils(void);
        void                    printMoveDiagnostics(int rep);
        void                    writeCheckpoint(void);
        bool                    loadCheckpoint(void);
        Tree*                   getTree(void);
        void                    setVerbose(bool b) { verbose = b; }
        void                    setLabel(int i) { runLabel = i; }

    protected:
        RandomVariable*         resumeRng(void);
        std::vector<std::string> resumeParameterNames(void);
        std::vector<std::string> resumeLatentNames(void);

    private:
        void                    sample(unsigned long n, double lnL, double lnP);
        void                    checkCheckpointRoundTrip(const std::string& path);
        PhylogeneticModel*      model;
        bool                    verbose = false;
        int                     runLabel = 0;
        int                     numCycles;
        int                     thinning;
        bool                    tuning = false;
        double                  curLnL;
        double                  curLnP;
};

#endif
