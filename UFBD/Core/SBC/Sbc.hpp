#ifndef Sbc_hpp
#define Sbc_hpp

#include "ForwardSimulator.hpp"
#include "Probability.hpp"

#include <string>
#include <vector>

class RandomVariable;

struct SbcConfig {
    int                     numReps;
    bool                    simulateOnly;
    bool                    emitFiles;
    bool                    originConditioning;
    ConditioningEvent       condEvent;
    double                  rho;
    double                  bb;
    std::vector<double>     intervalStart;
    Probability::PriorSpec  lambdaPrior;
    Probability::PriorSpec  muPrior;
    Probability::PriorSpec  psiPrior;
    Probability::PriorSpec  startAgePrior;
    long                    mcmcGen;
    double                  burninFraction;
    int                     mcmcThin;
    int                     rankBins;
    std::string             dumpPrefix;
};

class Sbc {

    public:
                        Sbc(const SbcConfig& c, RandomVariable* r) : cfg(c), rng(r) {}
        void            run(void);

    private:
        SimParams       drawParams(void);
        void            runSimulateOnly(void);
        void            runInference(void);
        void            runEmit(void);
        SbcConfig       cfg;
        RandomVariable* rng;
};

#endif
