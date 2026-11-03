#ifndef ForwardSimulator_hpp
#define ForwardSimulator_hpp

#include <string>
#include <vector>

#include "UserSettings.hpp"

class RandomVariable;

struct SimParams {
    std::vector<double> intervalStart;
    std::vector<int>    lambdaIdx;
    std::vector<int>    muIdx;
    std::vector<std::vector<int>> psiIdx;
    std::vector<double> lambda;
    std::vector<double> mu;
    std::vector<std::vector<double>> psi;
    std::vector<double> lambdaOU;
    std::vector<double> muOU;
    std::vector<std::vector<double>> psiOU;
    double              rho = 1.0;
    double              ue = 0.0;
    double              startAge = 0.0;
    bool                originConditioning = false;
    ConditioningEvent   condEvent = ConditioningEvent::SURVIVAL;
};

struct SimResult {
    std::string         backboneNewick;
    std::vector<double> fossilAges;
    std::vector<int>    fossilTypes;
    int                 numExtantSampled = 0;
    int                 numBackbone = 0;
    int                 numUE = 0;
    int                 numFossils = 0;
    std::vector<double> fossilAttachAge;
    std::vector<char>   fossilIsSA;
    std::vector<std::vector<int>> fossilEdgeTaxa;
    std::vector<double> ueAttachAge;
    std::vector<double> backboneNodeAges;
    long                attempts = 0;
};

class ForwardSimulator {

    public:
                        ForwardSimulator(RandomVariable* r) : rng(r), maxActiveEdges(100000) {}
        SimResult       simulate(const SimParams& p);

    private:
        RandomVariable* rng;
        int             maxActiveEdges;
};

#endif
