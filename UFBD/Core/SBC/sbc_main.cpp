#ifdef BUILD_SBC

#include "Msg.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Sbc.hpp"
#include "UserSettings.hpp"

#include <cstdio>
#include <random>
#include <string>
#include <vector>

int main(int argc, const char* argv[]){

    std::string invocation;
    for(int i = 0; i < argc; i++){
        std::string a = argv[i];
        bool q = a.find_first_of(" (),") != std::string::npos;
        invocation += (i ? " " : "");
        invocation += q ? "\"" + a + "\"" : a;
    }
    printf("%s\n", invocation.c_str());

    std::string mode, out;
    bool repsSet = false;
    int reps = 0;
    double bb = 1.0;
    std::vector<const char*> rest;
    rest.push_back(argv[0]);
    for(int i = 1; i < argc; i++){
        std::string a = argv[i];
        if(a == "-sbc_mode" && i + 1 < argc)        mode = argv[++i];
        else if(a == "-sbc_reps" && i + 1 < argc){  reps = std::stoi(argv[++i]); repsSet = true; }
        else if(a == "-sbc_bb" && i + 1 < argc)     bb = std::stod(argv[++i]);
        else if(a == "-sbc_out" && i + 1 < argc)    out = argv[++i];
        else                                        rest.push_back(argv[i]);
    }

    if(mode != "simulate" && mode != "infer" && mode != "emit")
        Msg::error("-sbc_mode must be supplied as 'simulate', 'infer', or 'emit'.");
    if(repsSet == false)
        Msg::error("-sbc_reps is required.");
    if(mode == "emit" && out.empty())
        Msg::error("-sbc_out (output file prefix) is required.");

    UserSettings& s = UserSettings::userSettings();
    s.initializeSettings((int)rest.size(), rest.data(), true);

    Probability::PriorSpec lambdaPrior = s.getLambdaPrior();
    Probability::PriorSpec muPrior     = s.getMuPrior();
    if(lambdaPrior.set == false || muPrior.set == false)
        Msg::error("-lambda_prior, -mu_prior, and -psi_prior are all required.");
    for(int t = 0; t < s.getNumPsiTypes(); t++)
        if(s.getPsiPrior(t).set == false)
            Msg::error("-lambda_prior, -mu_prior, and -psi_prior are all required.");
    if(s.getConditionAgePriorSet() == false)
        Msg::error("-conditioning must supply a distribution to draw the root age from.");

    SbcConfig cfg;
    cfg.numReps = reps;
    cfg.simulateOnly = (mode == "simulate");
    cfg.emitFiles = (mode == "emit");
    cfg.originConditioning = (s.getConditioning() == Conditioning::ORIGIN);
    cfg.condEvent = s.getConditioningEvent();
    cfg.rho = s.getRho();
    cfg.bb = bb;
    if(cfg.originConditioning == false && cfg.bb < 1.0)
        Msg::error("crown conditioning requires bb=1. Use origin conditioning for bb<1.");
    cfg.intervalStart.push_back(0.0);
    for(double t : s.getSkylineTimes())
        cfg.intervalStart.push_back(t);
    cfg.lambdaTimes.push_back(0.0);
    for(double t : s.getLambdaSkylineTimes()) cfg.lambdaTimes.push_back(t);
    cfg.muTimes.push_back(0.0);
    for(double t : s.getMuSkylineTimes())     cfg.muTimes.push_back(t);
    cfg.numPsiTypes = s.getNumPsiTypes();
    cfg.psiTypeNames = s.getPsiTypeNames();
    for(int t = 0; t < cfg.numPsiTypes; t++){
        std::vector<double> times;
        times.push_back(0.0);
        for(double tv : s.getPsiSkylineTimes(t)) times.push_back(tv);
        cfg.psiTimes.push_back(times);
        cfg.psiPriors.push_back(s.getPsiPrior(t));
        cfg.psiGroups.push_back(s.getPsiGroups(t));
        cfg.psiGroupPriors.push_back(s.getPsiGroupPrior(t));
    }
    cfg.lambdaPrior = lambdaPrior;
    cfg.muPrior = muPrior;
    cfg.lambdaGroups = s.getLambdaGroups();
    cfg.muGroups = s.getMuGroups();
    cfg.lambdaGroupPrior = s.getLambdaGroupPrior();
    cfg.muGroupPrior = s.getMuGroupPrior();
    cfg.startAgePrior = { true, s.getConditionAgePrior(), s.getConditionAgePriorP1(), s.getConditionAgePriorP2() };
    cfg.dumpPrefix = out;

    unsigned int seed = s.getSeed();
    RandomVariable rng(seed);
    RandomVariable::setActiveInstance(&rng);

    Sbc sbc(cfg, &rng);
    sbc.run();
    return 0;
}

#endif
