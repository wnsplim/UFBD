#ifndef ForwardSimulator_hpp
#define ForwardSimulator_hpp

#include <string>
#include <vector>

class RandomVariable;

struct SimParams {
    std::vector<double> intervalStart;
    std::vector<double> lambda;
    std::vector<double> mu;
    std::vector<double> psi;
    double              rho = 1.0;
    double              bb = 1.0;
    double              startAge = 0.0;
    bool                originConditioning = false;
};

struct SimResult {
    std::string         backboneNewick;
    std::vector<double> fossilAges;
    int                 numExtantSampled = 0;
    int                 numBackbone = 0;
    int                 numUE = 0;
    int                 numFossils = 0;
    double              crownAge = 0.0;
    double              originAge = 0.0;
};

class ForwardSimulator {

    public:
                        ForwardSimulator(RandomVariable* r) : rng(r), maxLineages(100000) {}
        SimResult       simulate(const SimParams& p);

    private:
        RandomVariable* rng;
        int             maxLineages;
};

#endif
