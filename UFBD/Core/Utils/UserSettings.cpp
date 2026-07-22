#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <set>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <random>

#include "ConfigReader.hpp"
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
        return RateMode::INDEP;
    }
    return psiMode;
}

OUHyperSpec UserSettings::getPsiOU(int t){
    if(psiTypeNames.empty() == false){
        std::map<std::string, OUHyperSpec>::iterator it = psiOUByName.find(psiTypeNames[t]);
        if(it != psiOUByName.end()) return it->second;
        return OUHyperSpec();
    }
    return psiOU;
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
    conditionAgePriorP3  = 0.0;
    model           = Model::UFBD;
    rho             = 1.0;
    seed            = 0;
    seedSet         = false;
    chainLength     = 1000000;
    numCoupledChains       = 1;
    numRuns         = 4;
    resume          = false;
    arLog           = false;
    noLatentLog     = false;
    autoChainLength = false;
    burninFraction  = 0.25;
    rhatThreshold   = 1.01;
    essThreshold    = 0.0;
    essThresholdSet = false;
    rhatThresholdSet = false;
    maxGen          = 1000000000;
    numCores        = 1;
    thinning        = 1000;
    hessianFile     = "";
    clockModelName  = "ucln";
    sigma2Param     = Sigma2Param::PNCP;
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
    lambdaMode      = RateMode::INDEP;
    muMode          = RateMode::INDEP;
    psiMode         = RateMode::INDEP;
    ageOffset       = 0.0;
    cpuTime         = false;

    std::vector<std::string> arguments;
    for (int i = 0; i < argc; i++)
        arguments.push_back(std::string(argv[i]));

    if (sbcMode == false) {
        std::vector<std::string> cfgTokens;
        for (size_t i = 1; i < arguments.size(); ) {
            if (arguments[i] == "-config") {
                if (i + 1 >= arguments.size()) Msg::error("flag \"-config\" expects a file path.");
                configFilePath = arguments[i + 1];
                cfgTokens = ConfigReader::translate(configFilePath);
                arguments.erase(arguments.begin() + i, arguments.begin() + i + 2);
            } else {
                i++;
            }
        }
        for (size_t i = 0; i < cfgTokens.size(); i++)
            arguments.push_back(cfgTokens[i]);
    }

    for (size_t i = 1; i < arguments.size(); i++)
        if (arguments[i] == "-help" || arguments[i] == "-h") { printHelp(); std::exit(0); }

    for (size_t i = 0; i < arguments.size(); i++) {
        bool q = arguments[i].find_first_of(" (),") != std::string::npos;
        invocation += (i ? " " : "");
        invocation += q ? "\"" + arguments[i] + "\"" : arguments[i];
    }
    if (sbcMode == false) {
        std::cout << invocation << "\n";
        if (configFilePath.empty() == false)
            std::cout << "Config file:                  " << configFilePath << "\n";
    }

    std::set<std::string> knownFlags = {
        "-tree_output", "-log_output", "-backbone_tree", "-clade_def", "-fossils", "-conditioning", "-rho", "-seed", "-chain_length", "-thinning", "-coupled_chains", "-cores", "-help", "-h",
        "-lambda_prior", "-mu_prior", "-psi_prior", "-psi_types",
        "-lambda_skyline_times", "-mu_skyline_times", "-psi_skyline_times", "-clock_partitions",
        "-lambda_prior_mode", "-mu_prior_mode", "-psi_prior_mode",
        "-lambda_ou_theta", "-mu_ou_theta", "-psi_ou_theta", "-lambda_ou_sd", "-mu_ou_sd", "-psi_ou_sd", "-lambda_ou_nu", "-mu_ou_nu", "-psi_ou_nu",
        "-lambda_groups", "-mu_groups", "-psi_groups", "-cpu_time", "-age_offset",
        "-hessian", "-clock_model", "-n_states", "-rgene_gamma", "-sigma2_gamma", "-sigma2_param", "-pncp_tuning",
        "-sequence", "-partition", "-ctmc_gamma_cat", "-datatype", "-ctmc_model", "-ctmc_inv", "-ctmc_freq",
        "-parallel_chains", "-burn_in", "-rhat", "-min_ess", "-max_gen", "-delta_temperature", "-swap_interval", "-resume", "-ar_log", "-no_latent_log"
    };
    std::set<std::string> valueFlags = {
        "-tree_output", "-log_output", "-backbone_tree", "-clade_def", "-fossils", "-conditioning", "-rho", "-seed", "-chain_length", "-thinning", "-coupled_chains", "-cores",
        "-lambda_prior", "-mu_prior", "-psi_prior", "-psi_types",
        "-lambda_skyline_times", "-mu_skyline_times", "-psi_skyline_times", "-clock_partitions",
        "-lambda_prior_mode", "-mu_prior_mode", "-psi_prior_mode",
        "-lambda_ou_theta", "-mu_ou_theta", "-psi_ou_theta", "-lambda_ou_sd", "-mu_ou_sd", "-psi_ou_sd", "-lambda_ou_nu", "-mu_ou_nu", "-psi_ou_nu",
        "-lambda_groups", "-mu_groups", "-psi_groups", "-cpu_time", "-age_offset",
        "-hessian", "-clock_model", "-n_states", "-rgene_gamma", "-sigma2_gamma", "-sigma2_param", "-pncp_tuning",
        "-sequence", "-partition", "-ctmc_gamma_cat", "-datatype", "-ctmc_model", "-ctmc_inv", "-ctmc_freq",
        "-parallel_chains", "-burn_in", "-rhat", "-min_ess", "-max_gen", "-delta_temperature", "-swap_interval"
    };
    if (sbcMode) {
        knownFlags.insert("-fbd_model");
        valueFlags.insert("-fbd_model");
    }

    bool coupledChainsProvided = false;
    for (int i = 1; i < (int)arguments.size(); i++) {
        std::string arg = arguments[i];

        // Check it looks like a flag
        if (arg.empty())
            Msg::error("empty argument at position " + std::to_string(i));

        if (knownFlags.find(arg) == knownFlags.end())
            Msg::error("unknown flag \"" + arg + "\".");

        if (arg == "-resume") {
            resume = true;
            continue;
        }

        if (arg == "-ar_log") {
            arLog = true;
            continue;
        }

        if (arg == "-no_latent_log") {
            noLatentLog = true;
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
                    parsePriorInto(arguments[++i], conditionAgePrior, conditionAgePriorP1, conditionAgePriorP2, conditionAgePriorP3);
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
            } else if (arg == "-lambda_prior" || arg == "-mu_prior") {
                // <prior> (all bins) OR <group_id>:<prior> (per-chunk override)
                size_t c = val.find(':');
                std::string first = (c == std::string::npos) ? val : val.substr(0, c);
                bool isGid = (c != std::string::npos) && first.empty() == false && std::all_of(first.begin(), first.end(), ::isdigit);
                Probability::PriorSpec& base = (arg == "-lambda_prior") ? lambdaPrior : muPrior;
                std::map<int, Probability::PriorSpec>& gp = (arg == "-lambda_prior") ? lambdaGroupPrior : muGroupPrior;
                if (isGid) {
                    Probability::PriorSpec ps; parsePriorInto(val.substr(c + 1), ps.family, ps.p1, ps.p2, ps.p3); ps.set = true;
                    gp[std::stoi(first)] = ps;
                } else {
                    parsePriorInto(val, base.family, base.p1, base.p2, base.p3); base.set = true;
                }
            } else if (arg == "-psi_prior") {
                // <prior> | <type>:<prior> | <group_id>:<prior> | <type>:<group_id>:<prior>
                size_t c = val.find(':');
                std::string first = (c == std::string::npos) ? val : val.substr(0, c);
                if (c == std::string::npos || isPriorFamilyKeyword(first)) {
                    parsePriorInto(val, psiPrior.family, psiPrior.p1, psiPrior.p2, psiPrior.p3); psiPrior.set = true;
                } else if (first.empty() == false && std::all_of(first.begin(), first.end(), ::isdigit)) {
                    Probability::PriorSpec ps; parsePriorInto(val.substr(c + 1), ps.family, ps.p1, ps.p2, ps.p3); ps.set = true;
                    psiGroupPrior[std::stoi(first)] = ps;
                } else {
                    std::string rest = val.substr(c + 1);
                    size_t c2 = rest.find(':');
                    std::string r0 = (c2 == std::string::npos) ? rest : rest.substr(0, c2);
                    Probability::PriorSpec ps;
                    if (c2 != std::string::npos && r0.empty() == false && std::all_of(r0.begin(), r0.end(), ::isdigit)) {
                        parsePriorInto(rest.substr(c2 + 1), ps.family, ps.p1, ps.p2, ps.p3); ps.set = true;
                        psiGroupPriorByName[first][std::stoi(r0)] = ps;
                    } else {
                        parsePriorInto(rest, ps.family, ps.p1, ps.p2, ps.p3); ps.set = true;
                        psiPriorByName[first] = ps;
                    }
                }
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
            } else if (arg == "-psi_types") {
                std::stringstream ss(val); std::string tok;
                while (std::getline(ss, tok, ',')) if (tok.empty() == false) psiTypeNames.push_back(tok);
            } else if (arg == "-lambda_prior_mode" || arg == "-mu_prior_mode") {
                std::string v = val;
                for (char& ch : v) ch = std::tolower((unsigned char)ch);
                if (v != "indep" && v != "ou") Msg::error("flag \"" + arg + "\" expects indep|ou, but got \"" + val + "\".");
                RateMode m = (v == "ou") ? RateMode::OU : RateMode::INDEP;
                if (arg == "-lambda_prior_mode") lambdaMode = m;
                else muMode = m;
            } else if (arg == "-psi_prior_mode") {
                std::string nm, ms = val;
                size_t c = val.find(':');
                if (c != std::string::npos) { nm = val.substr(0, c); ms = val.substr(c + 1); }
                for (char& ch : ms) ch = std::tolower((unsigned char)ch);
                if (ms != "indep" && ms != "ou") Msg::error("flag \"-psi_prior_mode\" expects indep|ou, but got \"" + ms + "\".");
                RateMode m = (ms == "ou") ? RateMode::OU : RateMode::INDEP;
                if (nm.empty()) psiMode = m; else psiModeByName[nm] = m;
            } else if (arg == "-lambda_ou_theta" || arg == "-mu_ou_theta" || arg == "-psi_ou_theta") {
                bool isPsi = (arg.rfind("-psi_", 0) == 0);
                std::string nm, body = val;
                if (isPsi) {
                    size_t c = val.find(':');
                    if (c != std::string::npos) { nm = val.substr(0, c); body = val.substr(c + 1); }
                }
                OUHyperSpec* target;
                if (arg.rfind("-lambda_", 0) == 0)   target = &lambdaOU;
                else if (arg.rfind("-mu_", 0) == 0)   target = &muOU;
                else                                  target = nm.empty() ? &psiOU : &psiOUByName[nm];
                double med, sd;
                size_t comma = body.find(',');
                if (comma == std::string::npos) Msg::error("flag \"" + arg + "\" expects median,sd, but got \"" + val + "\".");
                try { med = std::stod(body.substr(0, comma)); sd = std::stod(body.substr(comma + 1)); }
                catch (...) { Msg::error("flag \"" + arg + "\" expects median,sd, but got \"" + val + "\"."); }
                if (med <= 0.0 || sd <= 0.0) Msg::error("flag \"" + arg + "\" median and sd must be > 0.");
                target->thetaSet = true; target->thetaMedian = med; target->thetaSd = sd;
            } else if (arg == "-lambda_ou_sd" || arg == "-mu_ou_sd" || arg == "-psi_ou_sd"
                    || arg == "-lambda_ou_nu" || arg == "-mu_ou_nu" || arg == "-psi_ou_nu") {
                bool isPsi = (arg.rfind("-psi_", 0) == 0);
                std::string nm, body = val;
                if (isPsi) {
                    size_t c = val.find(':');
                    if (c != std::string::npos) { nm = val.substr(0, c); body = val.substr(c + 1); }
                }
                OUHyperSpec* target;
                if (arg.rfind("-lambda_", 0) == 0)   target = &lambdaOU;
                else if (arg.rfind("-mu_", 0) == 0)   target = &muOU;
                else                                  target = nm.empty() ? &psiOU : &psiOUByName[nm];
                double shape, rate;
                size_t comma = body.find(',');
                if (comma == std::string::npos) Msg::error("flag \"" + arg + "\" expects shape,rate, but got \"" + val + "\".");
                try { shape = std::stod(body.substr(0, comma)); rate = std::stod(body.substr(comma + 1)); }
                catch (...) { Msg::error("flag \"" + arg + "\" expects shape,rate, but got \"" + val + "\"."); }
                if (shape <= 0.0 || rate <= 0.0) Msg::error("flag \"" + arg + "\" shape and rate must be > 0.");
                if (arg.find("_ou_sd") != std::string::npos) { target->sdSet = true; target->sdShape = shape; target->sdRate = rate; }
                else                                          { target->nuSet = true; target->nuShape = shape; target->nuRate = rate; }
            } else if (arg == "-age_offset") {
                try { ageOffset = std::stod(val); } catch (...) { Msg::error("flag \"-age_offset\" expects a number, but got \"" + val + "\"."); }
                if (ageOffset < 0.0) Msg::error("flag \"-age_offset\" must be >= 0.");
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
            } else if (arg == "-sigma2_param") {
                std::string v = val;
                for (char& ch : v) ch = std::tolower((unsigned char)ch);
                if (v == "pncp") sigma2Param = Sigma2Param::PNCP;
                else if (v == "c") sigma2Param = Sigma2Param::C;
                else if (v == "nc") sigma2Param = Sigma2Param::NC;
                else Msg::error("flag \"-sigma2_param\" expects pncp, c or nc, but got \"" + val + "\".");
            } else if (arg == "-pncp_tuning") {
                pncpTuningGens = std::stoul(val);
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
            } else if (arg == "-delta_temperature") {
                deltaTemperature = std::stod(val);
                if(deltaTemperature <= 0.0)
                    Msg::error("flag \"-delta_temperature\" must be positive.");
                deltaTemperatureProvided = true;
            } else if (arg == "-swap_interval") {
                resampleEvery = std::stoul(val);
                if(resampleEvery == 0)
                    Msg::error("flag \"-swap_interval\" must be a positive integer.");
                resampleEveryProvided = true;
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

    if ((sequenceFile.empty() == false || hessianFile.empty() == false) && sigma2Param == Sigma2Param::PNCP && pncpTuningGens == 0)
        Msg::error("sigma2_param=pncp requires pncp_tuning > 0");

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

    if (essThresholdSet == false)
        essThreshold = 100.0 * (double)numRuns;

    if (numCoupledChains < 1) {
        Msg::warning("chains must be >= 1; setting to 1.");
        numCoupledChains = 1;
    }

    if (deltaTemperatureProvided && numCoupledChains <= 1)
        Msg::warning("-delta_temperature has no effect without -coupled_chains > 1.");

    if (resampleEveryProvided && numCoupledChains <= 1)
        Msg::warning("-swap_interval has no effect without -coupled_chains > 1.");

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
    if (hessianFile.empty() == false && nstatesProvided && datatypeProvided) {
        int implied = (seqDataType == "aa" ? 20 : 4);
        if (implied != nStates)
            Msg::warning("-datatype implies " + std::to_string(implied) + " states but -n_states is " + std::to_string(nStates) + "; using " + std::to_string(nStates) + ".");
    }
    if (hessianFile.empty() == false && partitionFile.empty() == false)
        Msg::warning("-partition (NEXUS charset) applies to the sequence path only; ignoring it under -hessian.");

    if (ageOffset > 0.0) {
        double minCut = 0.0;
        bool haveCut = false;
        std::vector<std::vector<double>*> cutVecs = { &lambdaSkylineTimes, &muSkylineTimes, &psiSkylineTimes };
        for (std::map<std::string, std::vector<double>>::iterator it = psiTimesByName.begin(); it != psiTimesByName.end(); ++it)
            cutVecs.push_back(&it->second);
        for (std::vector<double>* v : cutVecs)
            for (double t : *v) { if (haveCut == false || t < minCut) { minCut = t; haveCut = true; } }
        if (haveCut && ageOffset >= minCut)
            Msg::error("-age_offset (" + std::to_string(ageOffset) + ") must be younger than the youngest rate change-time (" + std::to_string(minCut) + ").");
        for (std::vector<double>* v : cutVecs)
            for (double& t : *v) t -= ageOffset;
        if (conditionAgePriorSet) {
            if (conditionAgePrior == Probability::PriorFamily::FIXED) conditionAgePriorP1 -= ageOffset;
            else if (conditionAgePrior == Probability::PriorFamily::UNIFORM) { conditionAgePriorP1 -= ageOffset; conditionAgePriorP2 -= ageOffset; }
            else if (conditionAgePrior == Probability::PriorFamily::TRUNCATED_NORMAL) { conditionAgePriorP1 -= ageOffset; conditionAgePriorP3 -= ageOffset; }
            else conditionAgePriorP3 -= ageOffset;
        }
    }


    if (sbcMode == false) {
        if (treeOut.empty())
            Msg::warning("no tree output file specified (-tree_output).");
        if (parametersOut.empty())
            Msg::warning("no parameter output file specified (-log_output).");
    }
}

void UserSettings::parsePriorInto(const std::string& spec, Probability::PriorFamily& family, double& p1, double& p2, double& p3) {
    p3 = 0.0;
    const std::string help = "prior must be a fixed number or a distribution written improper[:offset], exp:rate[,offset], gamma:shape,rate[,offset], lognormal:mu,sigma[,offset], unif:a,b, or truncnormal:mean,sd[,offset]; got \"" + spec + "\".";
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
        if ((ps.size() != 1 && ps.size() != 2) || ps[0] <= 0.0) Msg::error("exp prior needs one positive rate and an optional offset: exp:rate[,offset].");
        family = Probability::PriorFamily::EXPONENTIAL; p1 = ps[0]; if (ps.size() == 2) p3 = ps[1];
    } else if (fam == "gamma") {
        if ((ps.size() != 2 && ps.size() != 3) || ps[0] <= 0.0 || ps[1] <= 0.0) Msg::error("gamma prior needs positive shape,rate and an optional offset: gamma:shape,rate[,offset].");
        family = Probability::PriorFamily::GAMMA; p1 = ps[0]; p2 = ps[1]; if (ps.size() == 3) p3 = ps[2];
    } else if (fam == "lognormal") {
        if ((ps.size() != 2 && ps.size() != 3) || ps[1] <= 0.0) Msg::error("lognormal prior needs mu,sigma with sigma>0 and an optional offset: lognormal:mu,sigma[,offset].");
        family = Probability::PriorFamily::LOGNORMAL; p1 = ps[0]; p2 = ps[1]; if (ps.size() == 3) p3 = ps[2];
    } else if (fam == "unif" || fam == "uniform") {
        if (ps.size() != 2 || ps[0] >= ps[1]) Msg::error("unif prior needs a<b: unif:a,b.");
        family = Probability::PriorFamily::UNIFORM; p1 = ps[0]; p2 = ps[1];
    } else if (fam == "truncnormal" || fam == "truncnorm" || fam == "normal") {
        if ((ps.size() != 2 && ps.size() != 3) || ps[1] <= 0.0) Msg::error("truncnormal prior needs mean,sd with sd>0 and an optional lower/offset: truncnormal:mean,sd[,offset].");
        family = Probability::PriorFamily::TRUNCATED_NORMAL; p1 = ps[0]; p2 = ps[1]; if (ps.size() == 3) p3 = ps[2];
    } else if (fam == "improper") {
        if (ps.size() > 1) Msg::error("improper prior takes at most an offset: improper[:offset].");
        family = Probability::PriorFamily::IMPROPER; if (ps.size() == 1) p3 = ps[0];
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
    if (numCoupledChains > 1)
        std::cout << "MC3 initial delta temperature: " << deltaTemperature << std::endl;
    if (numCoupledChains > 1)
        std::cout << "MC3 swap interval (gens):     " << resampleEvery << std::endl;
    std::cout << "Number of parallel chains:    " << numRuns << std::endl;
    std::cout << "Thinning:                     " << thinning << std::endl;
    if (sigma2Param == Sigma2Param::PNCP && clockPresent())
        std::cout << "PNCP tuning iterations:       " << pncpTuningGens << std::endl;
    std::cout << "Number of cores:              " << numCores << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;

}

void UserSettings::printHelp(void){
    std::cout << R"HELP(USAGE

  ufbd [flags]
  ufbd -config <file>     read settings from a config file (see README.md)

Every setting can be given as a flag (below) or in a config file. When a
setting appears in both, the config file is prioritized.

A prior is written one of these ways:
  improper | <number> (fixed value) | exp:rate | gamma:shape,rate |
  lognormal:mu,sigma | unif:a,b | truncnormal:mean,sd

REQUIRED
  -fossils <file>         tab-separated: taxon, min age, max age, clade,
                          CROWN|TOTAL|STEM, and an optional preservation type
  -conditioning <event> [prior]
                          event = crown | origin | anysample | extinct;
                          the optional prior is the root-age prior
                          (default improper)

INPUT
  -backbone_tree <file>   rooted binary NEWICK or NEXUS tree; if not given,
                          every taxon is treated as unresolved
  -clade_def <file>       name<TAB>taxon,taxon,... ; a clade is all taxa under
                          the MRCA of the listed ones

OUTPUT  (prefixes: the engine appends .log/.trees/.tree, _chainN per chain)
  -log_output <prefix>    parameter log: merged across parallel chains as
                          <prefix>.log, per chain as <prefix>_chainN.log
  -tree_output <prefix>   per-sample trees (same naming) plus the posterior
                          mean tree <prefix>.tree
                          The full console is also saved to
                          <prefix>.console.txt
  -no_latent_log          suppress the <prefix>_latent.log per-sample latent
                          trace log

MCMC
  -chain_length <N|auto>  number of generations, or 'auto' to stop on the
                          convergence targets below
  -thinning <N>           sample every N generations (default 1000)
  -burn_in <frac>         burn-in fraction (default 0.25)
  -parallel_chains <N>    independent replicate runs (needed for split R-hat) (default 4)
  -coupled_chains <N>     Metropolis-coupled chains per run (>1 enables MC3) (default 1)
  -delta_temperature <x>  MC3 initial temperature spacing (default 0.1); the range
                          self-tunes toward ~0.40 adjacent DEO swap acceptance
  -swap_interval <N>      MC3 generations between chain-swap attempts (default 1000)
  -cores <N>              number of threads (default 1)
  -seed <N>               random number seed (default random)

CONVERGENCE  (used only with -chain_length auto)
  -max_gen <N>            hard generation cap (default 1000000000)
  -min_ess <N>            pooled bulk-ESS and tail-ESS floor (default 100 x parallel_chains)
  -rhat <x>               target maximum rank-normalized split R-hat (default 1.01)

FBD RATES
  -rho <frac>             extant sampling fraction (default 1)
  -lambda_prior <prior>   speciation rate prior (default exp:5)
  -mu_prior <prior>       extinction rate prior (default exp:5)
  -psi_prior <prior>      fossil-sampling rate prior (default exp:5)

SKYLINE  (piecewise-constant rates; write lambda, mu, or psi for <r>)
  -<r>_skyline_times <t,...>
                          interior time boundaries (k boundaries -> k+1 bins)
  -<r>_groups <g,...>     group id per bin; bins sharing an id share a rate
  -<r>_prior_mode <indep|ou>
                          independent rates, or log-OU smoothing across bins (default indep)
  -<r>_prior <gid>:<prior>
                          per-bin prior for one bin group (mode indep only)
  -<r>_ou_theta <med,sd>  OU level: theta ~ Normal(log(med), sd) (mode ou) (default 0.2,0.5)
  -<r>_ou_sd <shape,rate> gamma prior on the OU log-SD (mode ou) (default 6,5)
  -<r>_ou_nu <shape,rate> gamma prior on the OU reversion rate /Myr (mode ou)
                          (default 4, rate ~11.2*median_bin_width)
  -age_offset <t_min>     shift the skyline support to [t_min, inf); using this implies
                          assuming that the clade went extinct by t_min if t_min>0 (default 0)

MULTIPLE PRESERVATION TYPES  (split psi by type)
  -psi_types <name,...>   declare types; names must match the fossil table's
                          6th column
  Give per-type flags a 'name:' prefix, e.g.
    -psi_prior type1:unif:0.01,0.04
    -psi_skyline_times type1:20,40
    -psi_prior_mode type1:indep

CLOCK & SUBSTITUTION MODEL
  -sequence <file>        alignment (FASTA / PHYLIP / NEXUS)
  -hessian <file>         Hessian for approximate dating (PAML in.BV format)
  -partition <file>       NEXUS partition file (-sequence path only)
  -datatype <nt|aa>       sequence data type (default nt)
  -n_states <N>           substitution state count (-hessian path only) (default 4)
  -ctmc_model <gtr|file>  'gtr', or a PAML matrix file for amino acids (default gtr)
  -ctmc_gamma_cat <N>     number of discrete gamma rate categories (default 4)
  -ctmc_inv <on|off>      add a proportion of invariant sites (default off)
  -ctmc_freq <model|empirical|estimated>
                          equilibrium frequencies (default model)
  -clock_model <ucln|gbm> uncorrelated lognormal, or geometric Brownian motion (default ucln)
  -clock_partitions <g,...>
                          clock group id per partition (omit = single clock)
  -rgene_gamma <a,b,c>    gamma prior on the mean clock rate (default 2,2000,1)
  -sigma2_gamma <a,b,c>   gamma prior on the clock rate variance (default 1,10,1)
  -sigma2_param <pncp|c|nc>
                          parameterization of the sigma2. pncp (default) adapts per branch to
                          how much the data constrain that branch; c is fully centered and nc fully
                          non-centered.

OTHER
  -config <file>          read a config file (its format is in README.md)
  -help, -h               show this help message

RESUME
  -resume                 Continue an interrupted run by re-running the SAME
                          command with -resume added.
)HELP";
}
