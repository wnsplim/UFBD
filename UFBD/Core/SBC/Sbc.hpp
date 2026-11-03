#ifndef Sbc_hpp
#define Sbc_hpp

#include "ForwardSimulator.hpp"
#include "Probability.hpp"

#include <map>
#include <string>
#include <vector>

class Clade;
class Fossil;
class RandomVariable;
class Tree;

struct SbcConfig {
    int                     numReps;
    bool                    simulateOnly;
    bool                    emitFiles;
    bool                    originConditioning;
    ConditioningEvent       condEvent;
    double                  rho;
    double                  ue;
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
    std::vector<int>        lambdaGroups;
    std::vector<int>        muGroups;
    std::vector<std::vector<int>> psiGroups;
    std::map<int,Probability::PriorSpec> lambdaGroupPrior;
    std::map<int,Probability::PriorSpec> muGroupPrior;
    std::vector<std::map<int,Probability::PriorSpec>> psiGroupPriors;
    RateMode                lambdaMode;
    RateMode                muMode;
    std::vector<RateMode>   psiModes;
    OUHyperSpec             lambdaOU;
    OUHyperSpec             muOU;
    std::vector<OUHyperSpec> psiOU;
    std::string             dumpPrefix;
    double                  ageRange = 0.0;
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
        void            buildInputs(const SimResult& r, Tree* tree, std::vector<Clade>& clades, std::vector<Fossil>& fossils);
        SbcConfig       cfg;
        RandomVariable* rng;
};

#endif
