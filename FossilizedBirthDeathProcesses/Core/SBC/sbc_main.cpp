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
    double burnin = 0.0;
    double bb = 1.0;
    std::string postFossils;
    std::vector<const char*> rest;
    rest.push_back(argv[0]);
    for(int i = 1; i < argc; i++){
        std::string a = argv[i];
        if(a == "-sbc-mode" && i + 1 < argc)        mode = argv[++i];
        else if(a == "-sbc-reps" && i + 1 < argc){  reps = std::stoi(argv[++i]); repsSet = true; }
        else if(a == "-sbc-gen" && i + 1 < argc){   gen = std::stol(argv[++i]); genSet = true; }
        else if(a == "-sbc-burnin" && i + 1 < argc){ burnin = std::stod(argv[++i]); burninSet = true; }
        else if(a == "-sbc-thin" && i + 1 < argc){  thin = std::stoi(argv[++i]); thinSet = true; }
        else if(a == "-sbc-bins" && i + 1 < argc){  bins = std::stoi(argv[++i]); binsSet = true; }
        else if(a == "-sbc-bb" && i + 1 < argc)     bb = std::stod(argv[++i]);
        else if(a == "-sbc-out" && i + 1 < argc)    out = argv[++i];
        else if(a == "-post-fossils" && i + 1 < argc) postFossils = argv[++i];
        else                                        rest.push_back(argv[i]);
    }

    if(mode != "simulate" && mode != "infer" && mode != "posterior")
        Msg::error("SBC: -sbc-mode must be supplied as 'simulate', 'infer', or 'posterior'.");
    if(repsSet == false)
        Msg::error("SBC: -sbc-reps is required.");

    UserSettings& s = UserSettings::userSettings();
    s.initializeSettings((int)rest.size(), rest.data(), true);

    Probability::PriorSpec lambdaPrior = s.getLambdaPrior();
    Probability::PriorSpec muPrior     = s.getMuPrior();
    Probability::PriorSpec psiPrior    = s.getPsiPrior();
    if(lambdaPrior.set == false || muPrior.set == false || psiPrior.set == false)
        Msg::error("SBC: -lambda-prior, -mu-prior and -psi-prior are all required (these ARE the inference priors).");
    if(s.getConditionAgePriorSet() == false)
        Msg::error("SBC: -cond must supply a distribution to draw the start age x from (e.g. -cond crown \"lognormal(4.094,0.2)\").");

    SbcConfig cfg;
    cfg.numReps = reps;
    cfg.simulateOnly = (mode == "simulate");
    cfg.originConditioning = (s.getConditioning() == Conditioning::ORIGIN);
    cfg.condEvent = s.getConditioningEvent();
    cfg.rho = s.getRho();
    cfg.bb = bb;
    cfg.intervalStart.push_back(0.0);
    for(double t : s.getSkylineTimes())
        cfg.intervalStart.push_back(t);
    cfg.lambdaPrior = lambdaPrior;
    cfg.muPrior = muPrior;
    cfg.psiPrior = psiPrior;
    cfg.startAgePrior = { true, s.getConditionAgePrior(), s.getConditionAgePriorP1(), s.getConditionAgePriorP2() };
    cfg.dumpPrefix = out;
    cfg.posteriorMode = (mode == "posterior");
    if(cfg.posteriorMode){
        size_t p = 0;
        while(p < postFossils.size()){
            size_t comma = postFossils.find(',', p);
            std::string tok = postFossils.substr(p, (comma == std::string::npos) ? std::string::npos : comma - p);
            size_t colon = tok.find(':');
            cfg.postFossilMin.push_back(std::stod(tok.substr(0, colon)));
            cfg.postFossilMax.push_back(std::stod(tok.substr(colon + 1)));
            if(comma == std::string::npos) break;
            p = comma + 1;
        }
        if(cfg.postFossilMin.empty())
            Msg::error("SBC posterior: -post-fossils 'min:max,min:max,...' is required.");
    }

    if(cfg.simulateOnly == false){
        if(genSet == false)    Msg::error("SBC infer: -sbc-gen (MCMC generations per replicate) is required.");
        if(burninSet == false) Msg::error("SBC infer: -sbc-burnin (burn-in fraction in [0,1)) is required.");
        if(thinSet == false)   Msg::error("SBC infer: -sbc-thin (sample thinning interval) is required.");
        if(binsSet == false)   Msg::error("SBC infer: -sbc-bins (rank-histogram bin count) is required.");
        cfg.mcmcGen = gen;
        cfg.burninFraction = burnin;
        cfg.mcmcThin = thin;
        cfg.rankBins = bins;
    }else{
        cfg.mcmcGen = 0;
        cfg.burninFraction = 0.0;
        cfg.mcmcThin = 1;
        cfg.rankBins = 1;
    }

    unsigned int seed = s.getSeedSet() ? s.getSeed() : std::random_device{}();
    RandomVariable rng(seed);
    RandomVariable::setActiveInstance(&rng);

    Sbc sbc(cfg, &rng);
    sbc.run();
    return 0;
}

#endif
