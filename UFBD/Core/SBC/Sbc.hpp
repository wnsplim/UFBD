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
    std::vector<double>     lambdaTimes;
    std::vector<double>     muTimes;
    int                     numPsiTypes;
    std::vector<std::string>            psiTypeNames;
    std::vector<std::vector<double>>    psiTimes;
    std::vector<Probability::PriorSpec> psiPriors;
    Probability::PriorSpec  lambdaPrior;
    Probability::PriorSpec  muPrior;
    Probability::PriorSpec  startAgePrior;
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
