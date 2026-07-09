#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <random>

#include "Msg.hpp"
#include "Probability.hpp"
#include "UserSettings.hpp"


static void parseSkylineTimes(const std::string& flag, const std::string& val, std::vector<double>& out){
    std::stringstream ss(val);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty())
            continue;
        double t;
        try { t = std::stod(tok); }
        catch (...) { Msg::error("flag \"" + flag + "\" expects comma-separated numbers, but got \"" + tok + "\"."); }
        if (t <= 0.0)
            Msg::error("flag \"" + flag + "\" values must be positive, but got \"" + tok + "\".");
        if (out.empty() == false && t <= out.back())
            Msg::error("flag \"" + flag + "\" must be strictly increasing.");
        out.push_back(t);
    }
}

static bool isPriorFamilyKeyword(const std::string& s){
    std::string k = s;
    for(char& ch : k) ch = (char)std::tolower((unsigned char)ch);
    return k == "exp" || k == "exponential" || k == "gamma" || k == "lognormal"
        || k == "unif" || k == "uniform" || k == "truncnormal" || k == "truncnorm" || k == "normal"
        || k == "improper";
}

std::vector<double> UserSettings::getSkylineTimes(void){
    std::set<double> u;
    for(double t : lambdaSkylineTimes) u.insert(t);
    for(double t : muSkylineTimes)     u.insert(t);
    for(double t : psiSkylineTimes)    u.insert(t);
    for(std::map<std::string, std::vector<double>>::iterator it = psiTimesByName.begin(); it != psiTimesByName.end(); ++it)
        for(double t : it->second) u.insert(t);
    return std::vector<double>(u.begin(), u.end());
}

Probability::PriorSpec UserSettings::getPsiPrior(int t){
    if(psiTypeNames.empty() == false){
        std::map<std::string, Probability::PriorSpec>::iterator it = psiPriorByName.find(psiTypeNames[t]);
        if(it != psiPriorByName.end()) return it->second;
    }
    return psiPrior;
}

std::vector<double> UserSettings::getPsiSkylineTimes(int t){
    if(psiTypeNames.empty() == false){
        std::map<std::string, std::vector<double>>::iterator it = psiTimesByName.find(psiTypeNames[t]);
        if(it != psiTimesByName.end()) return it->second;
        return std::vector<double>();
    }
    return psiSkylineTimes;
}

RateMode UserSettings::getPsiMode(int t){
    if(psiTypeNames.empty() == false){
        std::map<std::string, RateMode>::iterator it = psiModeByName.find(psiTypeNames[t]);
        if(it != psiModeByName.end()) return it->second;
        return RateMode::IID;
    }
    return psiMode;
}

std::vector<int> UserSettings::getPsiGroups(int t){
    if(psiTypeNames.empty() == false){
        std::map<std::string, std::vector<int>>::iterator it = psiGroupsByName.find(psiTypeNames[t]);
        if(it != psiGroupsByName.end()) return it->second;
        return std::vector<int>();
    }
    return psiGroups;
}

std::map<int,Probability::PriorSpec> UserSettings::getPsiGroupPrior(int t){
    if(psiTypeNames.empty() == false){
        std::map<std::string, std::map<int,Probability::PriorSpec>>::iterator it = psiGroupPriorByName.find(psiTypeNames[t]);
        if(it != psiGroupPriorByName.end()) return it->second;
        return std::map<int,Probability::PriorSpec>();
    }
    return psiGroupPrior;
}

void UserSettings::initializeSettings(int argc, const char* argv[], bool sbcMode) {
    // Defaults
    treeOut         = "";
    parametersOut   = "";
    treeFile        = "";
    cladesFile      = "";
    fossilFile      = "";
    conditioningSet = false;
    conditioningEvent = ConditioningEvent::SURVIVAL;
    conditionAgePriorSet = false;
    conditionAgePrior    = Probability::PriorFamily::IMPROPER;
    conditionAgePriorP1  = 0.0;
    conditionAgePriorP2  = 0.0;
    model           = Model::UFBD;
    rho             = 1.0;
    seed            = 0;
    seedSet         = false;
    chainLength     = 1000000;
    numCoupledChains       = 1;
    numRuns         = 4;
    resume          = false;
    arLog           = false;
    autoChainLength = false;
    burninFraction  = 0.25;
    rhatThreshold   = 1.01;
    essThreshold    = 200.0;
    maxGen          = 1000000000;
    numCores        = 1;
    thinning        = 1000;
    hessianFile     = "";
    clockModelName  = "ucln";
    nStates         = 4;
    sequenceFile    = "";
    partitionFile   = "";
    numCats         = 4;
    seqDataType     = "nt";
    substModel      = "gtr";
    freqMode        = "model";
    useInvariant    = false;
    nstatesProvided = false;
    datatypeProvided = false;
    coresProvided   = false;
    rgeneGamma[0]   = 2.0;  rgeneGamma[1]  = 2000.0; rgeneGamma[2]  = 1.0;
    sigma2Gamma[0]  = 1.0;  sigma2Gamma[1] = 10.0;   sigma2Gamma[2] = 1.0;
    lambdaMode      = RateMode::IID;
    muMode          = RateMode::IID;
    psiMode         = RateMode::IID;
    hsmrfShifts     = std::log(2.0);
    hsmrfShiftSize  = 2.0;
    cpuTime         = false;

    std::vector<std::string> arguments;
    for (int i = 0; i < argc; i++)
        arguments.push_back(std::string(argv[i]));

    std::set<std::string> knownFlags = {
        "-tree_output", "-log_output", "-backbone_tree", "-clade_def", "-fossils", "-conditioning", "-rho", "-seed", "-chain_length", "-thinning", "-coupled_chains", "-cores", "-help", "-h",
        "-lambda_prior", "-mu_prior", "-psi_prior", "-psi_types",
        "-lambda_skyline_times", "-mu_skyline_times", "-psi_skyline_times", "-clock_partitions",
        "-lambda_prior_mode", "-mu_prior_mode", "-psi_prior_mode", "-lambda_groups", "-mu_groups", "-psi_groups", "-lambda_group_prior", "-mu_group_prior", "-psi_group_prior", "-hsmrf_shifts", "-hsmrf_shift_size", "-cpu_time",
        "-hessian", "-clock_model", "-n_states", "-rgene_gamma", "-sigma2_gamma",
        "-sequence", "-partition", "-ctmc_gamma_cat", "-datatype", "-ctmc_model", "-ctmc_inv", "-ctmc_freq",
        "-parallel_chains", "-burn_in", "-rhat", "-min_ess", "-max_gen", "-resume", "-ar_log"
    };
    std::set<std::string> valueFlags = {
        "-tree_output", "-log_output", "-backbone_tree", "-clade_def", "-fossils", "-conditioning", "-rho", "-seed", "-chain_length", "-thinning", "-coupled_chains", "-cores",
        "-lambda_prior", "-mu_prior", "-psi_prior", "-psi_types",
        "-lambda_skyline_times", "-mu_skyline_times", "-psi_skyline_times", "-clock_partitions",
        "-lambda_prior_mode", "-mu_prior_mode", "-psi_prior_mode", "-lambda_groups", "-mu_groups", "-psi_groups", "-lambda_group_prior", "-mu_group_prior", "-psi_group_prior", "-hsmrf_shifts", "-hsmrf_shift_size", "-cpu_time",
        "-hessian", "-clock_model", "-n_states", "-rgene_gamma", "-sigma2_gamma",
        "-sequence", "-partition", "-ctmc_gamma_cat", "-datatype", "-ctmc_model", "-ctmc_inv", "-ctmc_freq",
        "-parallel_chains", "-burn_in", "-rhat", "-min_ess", "-max_gen"
    };
    if (sbcMode) {
        knownFlags.insert("-fbd_model");
        valueFlags.insert("-fbd_model");
    }

    bool hsmrfProvided = false;
    bool coupledChainsProvided = false;
    for (int i = 1; i < (int)arguments.size(); i++) {
        std::string arg = arguments[i];

        // Check it looks like a flag
        if (arg.empty())
            Msg::error("empty argument at position " + std::to_string(i));

        if (knownFlags.find(arg) == knownFlags.end())
            Msg::error("unknown flag \"" + arg + "\".");

        // Help flag (no value)
        if (arg == "-help" || arg == "-h") {
            printHelp();
            return;
        }

        if (arg == "-resume") {
            resume = true;
            continue;
        }

        if (arg == "-ar_log") {
            arLog = true;
            continue;
        }


        // All remaining flags require a value — check it exists
        if (valueFlags.count(arg)) {
            if (i + 1 >= (int)arguments.size())
                Msg::error("flag \"" + arg + "\" requires a value but none was provided.");

            std::string val = arguments[++i];

            // Catch accidentally passing another flag as a value
            if (knownFlags.count(val))
                Msg::error("flag \"" + arg + "\" expects a value, but got another flag \"" + val + "\".");

            if (arg == "-tree_output") {
                treeOut = val + ".trees";
            } else if (arg == "-log_output") {
                parametersOut = val + ".log";
            } else if (arg == "-backbone_tree") {
                treeFile = val;
            } else if (arg == "-clade_def") {
                cladesFile = val;
            } else if (arg == "-fossils") {
                fossilFile = val;
            } else if (arg == "-conditioning") {
                std::string v = val;
                for (char& ch : v) ch = std::toupper((unsigned char)ch);
                if (v == "CROWN")          { conditioning = Conditioning::CROWN;  conditioningEvent = ConditioningEvent::SURVIVAL; }
                else if (v == "ORIGIN")    { conditioning = Conditioning::ORIGIN; conditioningEvent = ConditioningEvent::SURVIVAL; }
                else if (v == "ANYSAMPLE") { conditioning = Conditioning::ORIGIN; conditioningEvent = ConditioningEvent::ANYSAMPLE; }
                else if (v == "EXTINCT")   { conditioning = Conditioning::ORIGIN; conditioningEvent = ConditioningEvent::EXTINCT; }
                else Msg::error("flag \"-conditioning\" expects crown, origin, anysample, or extinct, but got \"" + val + "\".");
                conditioningSet = true;
                if (i + 1 < (int)arguments.size() && arguments[i + 1].empty() == false && arguments[i + 1][0] != '-') {
                    parsePriorInto(arguments[++i], conditionAgePrior, conditionAgePriorP1, conditionAgePriorP2);
                    conditionAgePriorSet = true;
                }
            } else if (arg == "-fbd_model") {
                std::string v = val;
                for (char& ch : v) ch = std::toupper((unsigned char)ch);
                if (v == "RFBD")        model = Model::RFBD;
                else if (v == "HEA14")  model = Model::HEA14;
                else if (v == "UFBD")   model = Model::UFBD;
                else Msg::error("flag \"-fbd_model\" expects RFBD, HEA14, or UFBD, but got \"" + val + "\".");
            } else if (arg == "-rho") {
                try {
                    rho = std::stod(val);
                } catch (...) {
                    Msg::error("flag \"-rho\" expects a number, but got \"" + val + "\".");
                }
                if (rho <= 0.0 || rho > 1.0)
                    Msg::error("flag \"-rho\" must be in (0, 1].");
            } else if (arg == "-seed") {
                try {
                    seed = (unsigned int)std::stoul(val);
                } catch (...) {
                    Msg::error("flag \"-seed\" expects a non-negative integer, but got \"" + val + "\".");
                }
                seedSet = true;
            } else if (arg == "-lambda_prior") {
                parsePriorInto(val, lambdaPrior.family, lambdaPrior.p1, lambdaPrior.p2); lambdaPrior.set = true;
            } else if (arg == "-mu_prior") {
                parsePriorInto(val, muPrior.family, muPrior.p1, muPrior.p2); muPrior.set = true;
            } else if (arg == "-psi_prior") {
                size_t c = val.find(':');
                bool typed = (c != std::string::npos) && (isPriorFamilyKeyword(val.substr(0, c)) == false);
                if (typed == false) { parsePriorInto(val, psiPrior.family, psiPrior.p1, psiPrior.p2); psiPrior.set = true; }
                else { Probability::PriorSpec ps; parsePriorInto(val.substr(c + 1), ps.family, ps.p1, ps.p2); ps.set = true; psiPriorByName[val.substr(0, c)] = ps; }
            } else if (arg == "-lambda_groups" || arg == "-mu_groups" || arg == "-psi_groups") {
                std::string nm, lst = val;
                if (arg == "-psi_groups") { size_t c = val.find(':'); if (c != std::string::npos) { nm = val.substr(0, c); lst = val.substr(c + 1); } }
                std::vector<int> g;
                std::stringstream ss(lst); std::string tok;
                while (std::getline(ss, tok, ',')) if (tok.empty() == false) {
                    try { g.push_back(std::stoi(tok)); }
                    catch (...) { Msg::error("flag \"" + arg + "\" expects comma-separated group IDs, but got \"" + val + "\"."); }
                }
                if (arg == "-lambda_groups") lambdaGroups = g;
                else if (arg == "-mu_groups") muGroups = g;
                else if (nm.empty()) psiGroups = g; else psiGroupsByName[nm] = g;
            } else if (arg == "-lambda_group_prior" || arg == "-mu_group_prior" || arg == "-psi_group_prior") {
                size_t c1 = val.find(':');
                if (c1 == std::string::npos) Msg::error("flag \"" + arg + "\" expects [type:]<group_id>:<prior>.");
                std::string first = val.substr(0, c1), nm, gidStr, priorSpec;
                bool firstIsInt = first.empty() == false && std::all_of(first.begin(), first.end(), ::isdigit);
                if (firstIsInt) {
                    gidStr = first; priorSpec = val.substr(c1 + 1);
                } else {
                    nm = first;
                    size_t c2 = val.find(':', c1 + 1);
                    if (c2 == std::string::npos) Msg::error("flag \"" + arg + "\" expects [type:]<group_id>:<prior>.");
                    gidStr = val.substr(c1 + 1, c2 - c1 - 1); priorSpec = val.substr(c2 + 1);
                }
                int gid = 0;
                try { gid = std::stoi(gidStr); }
                catch (...) { Msg::error("flag \"" + arg + "\" expects [type:]<group_id>:<prior>."); }
                Probability::PriorSpec ps; parsePriorInto(priorSpec, ps.family, ps.p1, ps.p2); ps.set = true;
                if (arg == "-lambda_group_prior") lambdaGroupPrior[gid] = ps;
                else if (arg == "-mu_group_prior") muGroupPrior[gid] = ps;
                else if (nm.empty()) psiGroupPrior[gid] = ps; else psiGroupPriorByName[nm][gid] = ps;
            } else if (arg == "-psi_types") {
                std::stringstream ss(val); std::string tok;
                while (std::getline(ss, tok, ',')) if (tok.empty() == false) psiTypeNames.push_back(tok);
            } else if (arg == "-lambda_prior_mode" || arg == "-mu_prior_mode") {
                std::string v = val;
                for (char& ch : v) ch = std::tolower((unsigned char)ch);
                RateMode rm = RateMode::IID;
                if (v == "iid") rm = RateMode::IID;
                else if (v == "smooth" || v == "hsmrf") rm = RateMode::SMOOTH;
                else Msg::error("flag \"" + arg + "\" expects iid or smooth, but got \"" + val + "\".");
                if (arg == "-lambda_prior_mode") lambdaMode = rm;
                else muMode = rm;
            } else if (arg == "-psi_prior_mode") {
                std::string nm, ms = val;
                size_t c = val.find(':');
                if (c != std::string::npos) { nm = val.substr(0, c); ms = val.substr(c + 1); }
                for (char& ch : ms) ch = std::tolower((unsigned char)ch);
                RateMode rm = RateMode::IID;
                if (ms == "iid") rm = RateMode::IID;
                else if (ms == "smooth" || ms == "hsmrf") rm = RateMode::SMOOTH;
                else Msg::error("flag \"-psi_prior_mode\" expects iid or smooth, but got \"" + ms + "\".");
                if (nm.empty()) psiMode = rm; else psiModeByName[nm] = rm;
            } else if (arg == "-hsmrf_shifts") {
                try { hsmrfShifts = std::stod(val); } catch (...) { Msg::error("flag \"-hsmrf_shifts\" expects a number, but got \"" + val + "\"."); }
                hsmrfProvided = true;
            } else if (arg == "-hsmrf_shift_size") {
                try { hsmrfShiftSize = std::stod(val); } catch (...) { Msg::error("flag \"-hsmrf_shift_size\" expects a number, but got \"" + val + "\"."); }
                hsmrfProvided = true;
            } else if (arg == "-cpu_time") {
                std::string v = val;
                for (char& ch : v) ch = std::tolower((unsigned char)ch);
                if (v == "on" || v == "true" || v == "1") cpuTime = true;
                else if (v == "off" || v == "false" || v == "0") cpuTime = false;
                else Msg::error("flag \"" + arg + "\" expects on or off, but got \"" + val + "\".");
            } else if (arg == "-lambda_skyline_times") {
                parseSkylineTimes(arg, val, lambdaSkylineTimes);
            } else if (arg == "-mu_skyline_times") {
                parseSkylineTimes(arg, val, muSkylineTimes);
            } else if (arg == "-psi_skyline_times") {
                size_t c = val.find(':');
                if (c == std::string::npos) parseSkylineTimes(arg, val, psiSkylineTimes);
                else { std::vector<double> tv; parseSkylineTimes(arg, val.substr(c + 1), tv); psiTimesByName[val.substr(0, c)] = tv; }
            } else if (arg == "-hessian") {
                hessianFile = val;
            } else if (arg == "-sequence") {
                sequenceFile = val;
            } else if (arg == "-partition") {
                partitionFile = val;
            } else if (arg == "-clock_partitions") {
                clockGroups.clear();
                std::stringstream cgs(val);
                std::string ct;
                while (std::getline(cgs, ct, ',')) {
                    try { clockGroups.push_back(std::stoi(ct)); }
                    catch (...) { Msg::error("flag \"-clock_partitions\" expects comma-separated integers, but got \"" + val + "\"."); }
                }
                std::vector<int> seen;
                for (int& g : clockGroups) {
                    int idx = -1;
                    for (size_t k = 0; k < seen.size(); k++) if (seen[k] == g) { idx = (int)k; break; }
                    if (idx < 0) { idx = (int)seen.size(); seen.push_back(g); }
                    g = idx;
                }
            } else if (arg == "-clock_model") {
                clockModelName = val;
                for (char& ch : clockModelName) ch = std::tolower((unsigned char)ch);
            } else if (arg == "-rgene_gamma") {
                std::stringstream ss(val); std::string tok; int k = 0;
                while (std::getline(ss, tok, ',') && k < 3) if (tok.empty() == false) rgeneGamma[k++] = std::stod(tok);
            } else if (arg == "-sigma2_gamma") {
                std::stringstream ss(val); std::string tok; int k = 0;
                while (std::getline(ss, tok, ',') && k < 3) if (tok.empty() == false) sigma2Gamma[k++] = std::stod(tok);
            } else if (arg == "-datatype") {
                datatypeProvided = true;
                std::string v = val;
                for (char& ch : v) ch = std::tolower((unsigned char)ch);
                if (v == "nt" || v == "dna" || v == "nucleotide") seqDataType = "nt";
                else if (v == "aa" || v == "protein" || v == "aminoacid") seqDataType = "aa";
                else Msg::error("flag \"-datatype\" expects nt or aa, but got \"" + val + "\".");
            } else if (arg == "-ctmc_model") {
                substModel = val;
            } else if (arg == "-ctmc_inv") {
                std::string v = val;
                for (char& ch : v) ch = std::tolower((unsigned char)ch);
                if (v == "on" || v == "true" || v == "1") useInvariant = true;
                else if (v == "off" || v == "false" || v == "0") useInvariant = false;
                else Msg::error("flag \"-ctmc_inv\" expects on or off, but got \"" + val + "\".");
            } else if (arg == "-ctmc_freq") {
                std::string v = val;
                for (char& ch : v) ch = std::tolower((unsigned char)ch);
                if (v == "model") freqMode = "model";
                else if (v == "empirical") freqMode = "empirical";
                else if (v == "estimated") freqMode = "estimate";
                else Msg::error("flag \"-ctmc_freq\" expects model, empirical, or estimated, but got \"" + val + "\".");
            } else if (arg == "-burn_in") {
                burninFraction = std::stod(val);
                if(burninFraction < 0.0 || burninFraction >= 1.0)
                    Msg::error("flag \"-burn_in\" expects a fraction in [0, 1).");
            } else if (arg == "-rhat") {
                rhatThreshold = std::stod(val);
                rhatThresholdSet = true;
            } else if (arg == "-min_ess") {
                essThreshold = std::stod(val);
                essThresholdSet = true;
            } else if (arg == "-max_gen") {
                maxGen = std::stoull(val);
            } else if (arg == "-chain_length" && val == "auto") {
                autoChainLength = true;
            }else {
                // Integer-valued flags
                // Check all characters are digits (allowing leading minus for negative detection)
                bool isNegative = (val[0] == '-');
                std::string digits = isNegative ? val.substr(1) : val;
                bool isInt = !digits.empty() && std::all_of(digits.begin(), digits.end(), ::isdigit);

                if (!isInt)
                    Msg::error("flag \"" + arg + "\" expects an integer, but got \"" + val + "\".");

                int intVal = std::stoi(val);

                if (arg == "-chain_length")        chainLength     = intVal;
                else if (arg == "-thinning") thinning = intVal;
                else if (arg == "-coupled_chains") { numCoupledChains = intVal; coupledChainsProvided = true; }
                else if (arg == "-cores") { numCores = intVal; coresProvided = true; }
                else if (arg == "-n_states") { nStates = intVal; nstatesProvided = true; }
                else if (arg == "-ctmc_gamma_cat")    numCats     = intVal;
                else if (arg == "-parallel_chains")    numRuns     = intVal;
            }
        }
    }

    // ── Post-parse validation  ──────────────────────────────

    if (conditioningSet == false)
        Msg::error("flag \"-conditioning\" is required.");

    if (sbcMode == false && fossilFile.empty())
        Msg::error("flag \"-fossils\" is required.");

    if (sbcMode == false && (sequenceFile.empty() == false || hessianFile.empty() == false) && treeFile.empty())
        Msg::error("when -sequence or -hessian is set, a backbone tree is required (-backbone_tree).");

    if (seedSet == false)
        seed = std::random_device{}();

    if (conditioningEvent == ConditioningEvent::EXTINCT && rho != 1.0) {
        Msg::warning("extinct conditioning has no extant sampling; forcing rho = 1 (was " + std::to_string(rho) + ").");
        rho = 1.0;
    }

    int maxNumThreads = (int)std::thread::hardware_concurrency();
    if (maxNumThreads <= 0) {
        Msg::warning("could not determine hardware thread count; using single thread.");
        maxNumThreads = 1;
    }

    if (numCoupledChains < 1) {
        Msg::warning("chains must be >= 1; setting to 1.");
        numCoupledChains = 1;
    }

    if (coupledChainsProvided && numCoupledChains == 1)
        Msg::warning("-coupled_chains 1 runs plain MCMC, not Metropolis-coupled MCMC; use -coupled_chains > 1 to enable MC3.");

    if (chainLength < 1)
        Msg::error("flag \"-chain_length\" must be a positive integer (or \"auto\").");

    if (sbcMode == false && autoChainLength == false && (essThresholdSet || rhatThresholdSet))
        Msg::warning("-min_ess / -rhat apply only under -chain_length auto; ignored for a fixed chain length.");

    if (thinning < 1)
        Msg::error("flag \"-thinning\" must be a positive integer.");

    if (coresProvided == false && numCoupledChains > 1)
        numCores = numCoupledChains;

    if (numCores > maxNumThreads) {
        if (coresProvided)
            Msg::warning("requested " + std::to_string(numCores) +
                         " cores, but only " + std::to_string(maxNumThreads) +
                         " available; capping at " + std::to_string(maxNumThreads) + ".");
        numCores = maxNumThreads;
    }
    if (numCores < 1)
        numCores = 1;

    bool empiricalModel = (substModel != "gtr");
    if (empiricalModel && seqDataType == "nt")
        Msg::error("empirical rate matrix (-ctmc_model " + substModel + ") is for amino-acid; use -datatype aa.");
    if (seqDataType == "aa" && empiricalModel == false)
        Msg::warning("amino-acid data (-datatype aa) with -ctmc_model gtr estimates 190 exchangeability parameters.");
    if (sequenceFile.empty() == false && nstatesProvided)
        Msg::warning("-n_states applies to the approximate-dating (-hessian) path only; ignoring it.");

    bool anySmooth = (lambdaMode == RateMode::SMOOTH || muMode == RateMode::SMOOTH || psiMode == RateMode::SMOOTH);
    for (std::map<std::string, RateMode>::iterator it = psiModeByName.begin(); it != psiModeByName.end(); ++it)
        if (it->second == RateMode::SMOOTH) anySmooth = true;
    if (hsmrfProvided && anySmooth == false)
        Msg::warning("-hsmrf_shifts / -hsmrf_shift_size have no effect: no rate uses HSMRF smoothing.");

    if (sbcMode == false) {
        if (treeOut.empty())
            Msg::warning("no tree output file specified (-tree_output).");
        if (parametersOut.empty())
            Msg::warning("no parameter output file specified (-log_output).");
    }
}

void UserSettings::parsePriorInto(const std::string& spec, Probability::PriorFamily& family, double& p1, double& p2) {
    const std::string help = "prior must be improper, a fixed number, or a distribution written exp:rate, gamma:shape,rate, lognormal:mu,sigma, unif:a,b, or truncnormal:mean,sd; got \"" + spec + "\".";
    size_t cp = spec.find(':');
    if (cp == std::string::npos) {
        std::string s = spec;
        for (char& ch : s) ch = std::tolower((unsigned char)ch);
        if (s == "improper") { family = Probability::PriorFamily::IMPROPER; return; }
        try { p1 = std::stod(spec); family = Probability::PriorFamily::FIXED; }
        catch (...) { Msg::error(help); }
        return;
    }
    std::string fam = spec.substr(0, cp);
    for (char& ch : fam) ch = std::tolower((unsigned char)ch);
    std::vector<double> ps;
    std::stringstream ss(spec.substr(cp + 1));
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        try { ps.push_back(std::stod(tok)); }
        catch (...) { Msg::error("prior parameter \"" + tok + "\" is not a number in \"" + spec + "\"."); }
    }
    if (fam == "exp" || fam == "exponential") {
        if (ps.size() != 1 || ps[0] <= 0.0) Msg::error("exp prior needs one positive rate: exp:rate.");
        family = Probability::PriorFamily::EXPONENTIAL; p1 = ps[0];
    } else if (fam == "gamma") {
        if (ps.size() != 2 || ps[0] <= 0.0 || ps[1] <= 0.0) Msg::error("gamma prior needs positive shape,rate: gamma:shape,rate.");
        family = Probability::PriorFamily::GAMMA; p1 = ps[0]; p2 = ps[1];
    } else if (fam == "lognormal") {
        if (ps.size() != 2 || ps[1] <= 0.0) Msg::error("lognormal prior needs mu,sigma with sigma>0: lognormal:mu,sigma.");
        family = Probability::PriorFamily::LOGNORMAL; p1 = ps[0]; p2 = ps[1];
    } else if (fam == "unif" || fam == "uniform") {
        if (ps.size() != 2 || ps[0] >= ps[1]) Msg::error("unif prior needs a<b: unif:a,b.");
        family = Probability::PriorFamily::UNIFORM; p1 = ps[0]; p2 = ps[1];
    } else if (fam == "truncnormal" || fam == "truncnorm" || fam == "normal") {
        if (ps.size() != 2 || ps[1] <= 0.0) Msg::error("truncnormal prior needs mean,sd with sd>0: truncnormal:mean,sd.");
        family = Probability::PriorFamily::TRUNCATED_NORMAL; p1 = ps[0]; p2 = ps[1];
    } else if (fam == "improper") {
        family = Probability::PriorFamily::IMPROPER;
    } else {
        Msg::error(help);
    }
}

void UserSettings::print(void) {

    if (treeFile.empty() == false)
        std::cout << "Tree input file:              " << treeFile << std::endl;
    if (cladesFile.empty() == false)
        std::cout << "Clades input file:            " << cladesFile << std::endl;
    if (fossilFile.empty() == false)
        std::cout << "Fossils input file:           " << fossilFile << std::endl;
    if (sequenceFile.empty() == false)
        std::cout << "Sequence input file:          " << sequenceFile << std::endl;
    if (hessianFile.empty() == false)
        std::cout << "Hessian input file:           " << hessianFile << std::endl;
    if (parametersOut.empty() == false)
        std::cout << "MCMC log output file:         " << parametersOut << std::endl;
    if (treeOut.empty() == false)
        std::cout << "Tree output file:             " << treeOut << std::endl;
    if (parametersOut.empty() == false) {
        std::string base = parametersOut;
        size_t d = base.rfind(".log");
        if (d != std::string::npos) base = base.substr(0, d);
        std::cout << "Mean tree output file:        " << base << ".tree" << std::endl;
    }
    std::cout << "Random number seed:           " << seed << std::endl;
    if (autoChainLength)
        std::cout << "Maximum chain length:         " << maxGen << std::endl;
    else
        std::cout << "Chain length:                 " << chainLength << std::endl;
    if (numCoupledChains > 1)
        std::cout << "Coupled chains per run (MC3): " << numCoupledChains << std::endl;
    std::cout << "Number of parallel chains:    " << numRuns << std::endl;
    std::cout << "Thinning:                     " << thinning << std::endl;
    std::cout << "Number of cores:              " << numCores << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;

}

void UserSettings::printHelp(void){
    std::cout << ".\n";
}
