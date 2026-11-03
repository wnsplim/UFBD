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

static void printSbcHelp(){
    printf("%s", R"SBC(USAGE

  sbc -sbc_mode <emit|simulate|infer> -sbc_reps <N> [flags]

MODES
  emit       write N simulated datasets to disk
  simulate   simulate N datasets and print summary statistics only
  infer      draw -> simulate -> infer -> rank; outputs <prefix>_ranks.tsv

SBC FLAGS
  -sbc_mode <emit|simulate|infer>
                          engine mode
  -sbc_reps <N>           number of replicate datasets
  -sbc_ue <prob>          probability that each extant taxon is designated an
                          unsampled extant (UE) (default 0; a value > 0 is
                          incompatible with crown conditioning)
  -sbc_out <prefix>       output prefix (required for -sbc_mode emit)
  -sbc_age_range <width>  give each fossil an age range of this width placed
                          uniformly at random around its simulated age and
                          truncated at the present, so fossil ages are sampled
                          (default 0 = fixed ages)
  -fbd_model <UFBD|RFBD|HEA14>
                          model variant

DRAW PRIORS AND CHAIN CONTROL
  These reuse the main ufbd flags: -lambda_prior, -mu_prior, -psi_prior
  (required for every rate under mode=indep), -*_prior_mode, -*_ou_theta,
  -*_ou_sd, -*_ou_nu, -conditioning <event> <prior> (required), -rho,
  -*_skyline_times, -psi_types, and, for -sbc_mode infer, the MCMC/convergence flags
  -chain_length, -min_ess, -max_gen, -thinning, -burn_in, -parallel_chains,
  -cores, -seed. Run 'ufbd -help' for their descriptions.
)SBC");
}

int main(int argc, const char* argv[]){

    for(int i = 1; i < argc; i++)
        if(std::string(argv[i]) == "-help" || std::string(argv[i]) == "-h"){ printSbcHelp(); return 0; }

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
    double ue = 0.0, ageRange = 0.0;
    std::vector<const char*> rest;
    rest.push_back(argv[0]);
    for(int i = 1; i < argc; i++){
        std::string a = argv[i];
        if(a == "-sbc_mode" && i + 1 < argc)        mode = argv[++i];
        else if(a == "-sbc_reps" && i + 1 < argc){  reps = std::stoi(argv[++i]); repsSet = true; }
        else if(a == "-sbc_ue" && i + 1 < argc)     ue = std::stod(argv[++i]);
        else if(a == "-sbc_age_range" && i + 1 < argc) ageRange = std::stod(argv[++i]);
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
    bool missing = (s.getLambdaMode() == RateMode::INDEP && lambdaPrior.set == false)
                || (s.getMuMode() == RateMode::INDEP && muPrior.set == false);
    for(int t = 0; t < s.getNumPsiTypes(); t++)
        if(s.getPsiMode(t) == RateMode::INDEP && s.getPsiPrior(t).set == false)
            missing = true;
    if(missing)
        Msg::error("-lambda_prior, -mu_prior, and -psi_prior are required for every rate under mode=indep.");
    if(s.getConditionAgePriorSet() == false)
        Msg::error("-conditioning must supply a distribution to draw the root age from.");

    SbcConfig cfg;
    cfg.numReps = reps;
    cfg.simulateOnly = (mode == "simulate");
    cfg.emitFiles = (mode == "emit");
    cfg.originConditioning = (s.getConditioning() == Conditioning::ORIGIN);
    cfg.condEvent = s.getConditioningEvent();
    cfg.rho = s.getRho();
    cfg.ue = ue;
    if(cfg.originConditioning == false && cfg.ue > 0.0)
        Msg::error("crown conditioning requires -sbc_ue 0.");
    if(ageRange < 0.0)
        Msg::error("-sbc_age_range must be >= 0.");
    cfg.ageRange = ageRange;
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
        cfg.psiModes.push_back(s.getPsiMode(t));
        cfg.psiOU.push_back(s.getPsiOU(t));
    }
    cfg.lambdaPrior = lambdaPrior;
    cfg.muPrior = muPrior;
    cfg.lambdaMode = s.getLambdaMode();
    cfg.muMode = s.getMuMode();
    cfg.lambdaOU = s.getLambdaOU();
    cfg.muOU = s.getMuOU();
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
