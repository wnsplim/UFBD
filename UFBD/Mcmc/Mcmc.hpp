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
        const std::vector<std::vector<double>>& latentColumns(void) { return latentCols; }
        const std::vector<std::string>&         latentNames(void) { return latentNms; }
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

    private:
        void                    sample(unsigned long n, double lnL, double lnP);
        PhylogeneticModel*      model;
        WriteTSV                latent;
        std::string             latentOut;
        bool                    writeLatent = false;
        bool                    verbose = false;
        int                     runLabel = 0;
        int                     numCycles;
        int                     thinning;
        double                  curLnL;
        double                  curLnP;
        std::vector<std::vector<double>> latentCols;
        std::vector<std::string>         latentNms;
};

#endif
