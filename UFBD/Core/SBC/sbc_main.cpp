#ifdef BUILD_SBC

#include "Msg.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Sbc.hpp"
#include "UserSettings.hpp"

#include <random>
#include <string>
#include <vector>

int main(int argc, const char* argv[]){

    std::string mode, out;
    bool repsSet = false, genSet = false, burninSet = false, thinSet = false, binsSet = false;
    int reps = 0, thin = 0, bins = 0;
    long gen = 0;
    bool autoGen = false;
    double burnin = 0.0;
    double bb = 1.0;
    std::vector<const char*> rest;
    rest.push_back(argv[0]);
    for(int i = 1; i < argc; i++){
        std::string a = argv[i];
        if(a == "-sbc_mode" && i + 1 < argc)        mode = argv[++i];
        else if(a == "-sbc_reps" && i + 1 < argc){  reps = std::stoi(argv[++i]); repsSet = true; }
        else if(a == "-sbc_gen" && i + 1 < argc){   std::string g = argv[++i]; if(g == "auto") autoGen = true; else gen = std::stol(g); genSet = true; }
        else if(a == "-sbc_burnin" && i + 1 < argc){ burnin = std::stod(argv[++i]); burninSet = true; }
        else if(a == "-sbc_thin" && i + 1 < argc){  thin = std::stoi(argv[++i]); thinSet = true; }
        else if(a == "-sbc_bins" && i + 1 < argc){  bins = std::stoi(argv[++i]); binsSet = true; }
        else if(a == "-sbc_bb" && i + 1 < argc)     bb = std::stod(argv[++i]);
        else if(a == "-sbc_out" && i + 1 < argc)    out = argv[++i];
        else                                        rest.push_back(argv[i]);
    }

    if(mode != "simulate" && mode != "infer" && mode != "emit")
        Msg::error("SBC: -sbc_mode must be supplied as 'simulate', 'infer', or 'emit'.");
    if(repsSet == false)
        Msg::error("SBC: -sbc_reps is required.");
    if(mode == "emit" && out.empty())
        Msg::error("SBC emit: -sbc_out (output file prefix) is required.");

    UserSettings& s = UserSettings::userSettings();
    s.initializeSettings((int)rest.size(), rest.data(), true);

    Probability::PriorSpec lambdaPrior = s.getLambdaPrior();
    Probability::PriorSpec muPrior     = s.getMuPrior();
    Probability::PriorSpec psiPrior    = s.getPsiPrior();
    if(lambdaPrior.set == false || muPrior.set == false || psiPrior.set == false)
        Msg::error("SBC: -lambda-prior, -mu-prior and -psi-prior are all required.");
    if(s.getConditionAgePriorSet() == false)
        Msg::error("SBC: -cond must supply a distribution to draw the origin age from.");

    SbcConfig cfg;
    cfg.numReps = reps;
    cfg.simulateOnly = (mode == "simulate");
    cfg.emitFiles = (mode == "emit");
    cfg.originConditioning = (s.getConditioning() == Conditioning::ORIGIN);
    cfg.condEvent = s.getConditioningEvent();
    cfg.rho = s.getRho();
    cfg.bb = bb;
    if(cfg.originConditioning == false && cfg.bb < 1.0)
        Msg::error("SBC: crown conditioning requires bb=1. Use origin conditioning for bb<1.");
    cfg.intervalStart.push_back(0.0);
    for(double t : s.getSkylineTimes())
        cfg.intervalStart.push_back(t);
    cfg.lambdaPrior = lambdaPrior;
    cfg.muPrior = muPrior;
    cfg.psiPrior = psiPrior;
    cfg.startAgePrior = { true, s.getConditionAgePrior(), s.getConditionAgePriorP1(), s.getConditionAgePriorP2() };
    cfg.dumpPrefix = out;
    (void)gen; (void)autoGen; (void)burnin; (void)thin; (void)bins;
    (void)genSet; (void)burninSet; (void)thinSet; (void)binsSet;

    unsigned int seed = s.getSeedSet() ? s.getSeed() : std::random_device{}();
    RandomVariable rng(seed);
    RandomVariable::setActiveInstance(&rng);

    Sbc sbc(cfg, &rng);
    sbc.run();
    return 0;
}

#endif
