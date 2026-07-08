#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <string>
#include <vector>

#include "ChainRunner.hpp"
#include "WriteTSV.hpp"
class PhylogeneticModel;
class Tree;

class Mcmc : public ChainRunner {
    public:
                                Mcmc(void) = delete;
                                Mcmc(int ng, int thin, PhylogeneticModel* m);
        void                    init(void);
        void                    advance(unsigned long nGens);
        void                    finalize(void);
        void                    setOutputPaths(const std::string& po, const std::string& to) { paramOut = po; treeOut = to; }
        const std::vector<std::vector<double>>& traceColumns(void) { return traceCols; }
        const std::vector<std::string>&         traceNames(void) { return traceNms; }
        const std::vector<std::vector<double>>& latentColumns(void) { return latentCols; }
        const std::vector<std::string>&         latentNames(void) { return latentNms; }
        bool                    treeHasFossils(void);
        void                    writeCheckpoint(void);
        bool                    loadCheckpoint(void);
        void                    resumeOutputs(void);
        unsigned long           currentGen(void) { return gen; }
        Tree*                   getTree(void);
        void                    setVerbose(bool b) { verbose = b; }
        void                    setLabel(int i) { runLabel = i; }

    private:
        void                    sample(unsigned long n, double lnL, double lnP);
        PhylogeneticModel*      model;
        WriteTSV                params;
        WriteTSV                trees;
        WriteTSV                latent;
        std::string             treeOut;
        std::string             paramOut;
        std::string             latentOut;
        bool                    writeTrees;
        bool                    writeLatent = false;
        bool                    verbose = false;
        int                     runLabel = 0;
        int                     numCycles;
        int                     thinning;
        unsigned long           gen;
        double                  curLnL;
        double                  curLnP;
        std::vector<std::vector<double>> traceCols;
        std::vector<std::string>         traceNms;
        std::vector<std::vector<double>> latentCols;
        std::vector<std::string>         latentNms;
};

#endif
