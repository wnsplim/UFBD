#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <string>
#include <vector>

#include "ChainRunner.hpp"
#include "WriteTSV.hpp"
class PhylogeneticModel;

class Mcmc : public ChainRunner {
    public:
                                Mcmc(void) = delete;
                                Mcmc(int ng, int pf, int sf, PhylogeneticModel* m);
        void                    run(void);
        void                    init(void);
        void                    advance(unsigned long nGens);
        void                    finalize(void);
        void                    setOutputPaths(const std::string& po, const std::string& to) { paramOut = po; treeOut = to; }
        const std::vector<std::vector<double>>& traceColumns(void) { return traceCols; }
        const std::vector<std::string>&         traceNames(void) { return traceNms; }

    private:
        void                    sample(unsigned long n, double lnL, double lnP);
        PhylogeneticModel*      model;
        WriteTSV                params;
        WriteTSV                trees;
        std::string             treeOut;
        std::string             paramOut;
        int                     numCycles;
        int                     printFrequency;
        int                     sampleFrequency;
        unsigned long           gen;
        double                  curLnL;
        double                  curLnP;
        std::vector<std::vector<double>> traceCols;
        std::vector<std::string>         traceNms;
};

#endif
