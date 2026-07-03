#ifndef ConvergenceRunner_hpp
#define ConvergenceRunner_hpp

#include <string>
#include <vector>

class ChainRunner;

class ConvergenceRunner {

    public:
                            ConvergenceRunner(std::vector<ChainRunner*> reps, const std::string& paramOut, const std::string& treeOut);
        bool                run(void);

    private:
        bool                report(unsigned long gen, bool finalPass);
        void                writeMerged(void);
        std::vector<ChainRunner*> replicates;
        std::string         paramBase;
        std::string         treeBase;
        std::vector<std::string> repParamFiles;
        std::vector<std::string> repTreeFiles;
};

#endif
