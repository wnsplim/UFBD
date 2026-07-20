#ifndef ConvergenceRunner_hpp
#define ConvergenceRunner_hpp

#include <string>
#include <vector>

class ChainRunner;

class ConvergenceRunner {

    public:
                            ConvergenceRunner(std::vector<ChainRunner*> reps, const std::string& paramOut, const std::string& treeOut);
        bool                run(void);
        double              getMaxRhat(void) const { return lastMaxRhat; }
        double              getMinBulkEss(void) const { return lastMinBulkEss; }
        double              getMinTailEss(void) const { return lastMinTailEss; }
        void                setEmitSummaryTree(bool b) { emitSummary = b; }
        void                setVerbose(bool b) { verbose = b; }

    private:
        bool                report(unsigned long gen, bool finalPass);
        void                writeMerged(void);
        std::vector<ChainRunner*> replicates;
        std::string         paramBase;
        std::string         treeBase;
        std::vector<std::string> repParamFiles;
        std::vector<std::string> repTreeFiles;
        double              lastMaxRhat = 1.0;
        double              lastMinBulkEss = -1.0;
        double              lastMinTailEss = -1.0;
        long                lastReportedGen = -1;
        bool                emitSummary = false;
        bool                verbose = false;
};

#endif
