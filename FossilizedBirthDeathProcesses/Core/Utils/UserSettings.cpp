#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>
#include <cctype>

#include "Msg.hpp"
#include "Probability.hpp"
#include "UserSettings.hpp"


void UserSettings::checkSettings(void) {

    if (settingsInitialized == false)
        Msg::error("Settings are not initialized");
}

void UserSettings::initializeSettings(int argc, const char* argv[]) {
    if (settingsInitialized == true) {
        Msg::warning("Settings have already been initialized");
        return;
    }
        
    // Defaults
    treeOut         = "";
    parametersOut   = "";
    treeFile        = "";
    cladesFile      = "";
    fossilFile      = "";
    conditioningSet = false;
    conditionAgePriorSet = false;
    model           = Model::UFBD;
    rho             = 1.0;
    seed            = 0;
    seedSet         = false;
    chainLength     = 100;
    numChains       = 4;
    numThreads      = 4;
    printFrequency  = 1;
    sampleFrequency = 1;

    std::vector<std::string> arguments;
    for (int i = 0; i < argc; i++)
        arguments.push_back(std::string(argv[i]));

    executablePath = arguments[0];
    dumpLnp = false;

    // Known flags and whether they take a value
    std::set<std::string> knownFlags = {
        "-to", "-po", "-t", "-c", "-f", "-cond", "-model", "-rho", "-seed", "-n", "-p", "-s", "-nc", "-nt", "-help", "-h",
        "-lambda-prior", "-mu-prior", "-psi-prior", "-skyline-times", "-lnp"
    };
    std::set<std::string> valueFlags = {
        "-to", "-po", "-t", "-c", "-f", "-cond", "-model", "-rho", "-seed", "-n", "-p", "-s", "-nc", "-nt",
        "-lambda-prior", "-mu-prior", "-psi-prior", "-skyline-times"
    };

    for (int i = 1; i < (int)arguments.size(); i++) {
        std::string arg = arguments[i];

        // Check it looks like a flag
        if (arg.empty())
            Msg::error("Empty argument at position " + std::to_string(i));

        if (knownFlags.find(arg) == knownFlags.end())
            Msg::error("Unknown flag \"" + arg + "\". Use -help to see valid options.");

        // Help flag (no value)
        if (arg == "-help" || arg == "-h") {
            printHelp();
            return;
        }

        if (arg == "-lnp") {
            dumpLnp = true;
            continue;
        }


        // All remaining flags require a value — check it exists
        if (valueFlags.count(arg)) {
            if (i + 1 >= (int)arguments.size())
                Msg::error("Flag \"" + arg + "\" requires a value but none was provided.");

            std::string val = arguments[++i];

            // Catch accidentally passing another flag as a value
            if (knownFlags.count(val))
                Msg::error("Flag \"" + arg + "\" expects a value, but got another flag \"" + val + "\".");

            if (arg == "-to") {
                treeOut = val;
            } else if (arg == "-po") {
                parametersOut = val;
            } else if (arg == "-t") {
                treeFile = val;
            } else if (arg == "-c") {
                cladesFile = val;
            } else if (arg == "-f") {
                fossilFile = val;
            } else if (arg == "-cond") {
                std::string v = val;
                for (char& ch : v) ch = std::toupper((unsigned char)ch);
                if (v == "CROWN")       conditioning = Conditioning::CROWN;
                else if (v == "ORIGIN") conditioning = Conditioning::ORIGIN;
                else Msg::error("Flag \"-cond\" expects crown or origin, but got \"" + val + "\".");
                conditioningSet = true;
                if (i + 1 < (int)arguments.size() && knownFlags.find(arguments[i + 1]) == knownFlags.end()) {
                    parsePriorInto(arguments[++i], conditionAgePrior, conditionAgePriorP1, conditionAgePriorP2);
                    conditionAgePriorSet = true;
                }
            } else if (arg == "-model") {
                std::string v = val;
                for (char& ch : v) ch = std::toupper((unsigned char)ch);
                if (v == "FBD")         model = Model::FBD;
                else if (v == "HEA14")  model = Model::HEA14;
                else if (v == "UFBD")   model = Model::UFBD;
                else Msg::error("Flag \"-model\" expects FBD, HEA14, or UFBD, but got \"" + val + "\".");
            } else if (arg == "-rho") {
                try {
                    rho = std::stod(val);
                } catch (...) {
                    Msg::error("Flag \"-rho\" expects a number, but got \"" + val + "\".");
                }
                if (rho <= 0.0 || rho > 1.0)
                    Msg::error("Flag \"-rho\" must be in (0, 1].");
            } else if (arg == "-seed") {
                try {
                    seed = (unsigned int)std::stoul(val);
                } catch (...) {
                    Msg::error("Flag \"-seed\" expects a non-negative integer, but got \"" + val + "\".");
                }
                seedSet = true;
            } else if (arg == "-lambda-prior") {
                parsePriorInto(val, lambdaPrior.family, lambdaPrior.p1, lambdaPrior.p2); lambdaPrior.set = true;
            } else if (arg == "-mu-prior") {
                parsePriorInto(val, muPrior.family, muPrior.p1, muPrior.p2); muPrior.set = true;
            } else if (arg == "-psi-prior") {
                parsePriorInto(val, psiPrior.family, psiPrior.p1, psiPrior.p2); psiPrior.set = true;
            } else if (arg == "-skyline-times") {
                std::stringstream ss(val);
                std::string tok;
                while (std::getline(ss, tok, ','))
                    if (tok.empty() == false) skylineTimes.push_back(std::stod(tok));
                std::sort(skylineTimes.begin(), skylineTimes.end());
            }else {
                // Integer-valued flags
                // Check all characters are digits (allowing leading minus for negative detection)
                bool isNegative = (val[0] == '-');
                std::string digits = isNegative ? val.substr(1) : val;
                bool isInt = !digits.empty() && std::all_of(digits.begin(), digits.end(), ::isdigit);

                if (!isInt)
                    Msg::error("Flag \"" + arg + "\" expects an integer, but got \"" + val + "\".");

                int intVal = std::stoi(val);

                if (arg == "-n")        chainLength     = intVal;
                else if (arg == "-p")   printFrequency  = intVal;
                else if (arg == "-s")   sampleFrequency = intVal;
                else if (arg == "-nc")  numChains       = intVal;
                else if (arg == "-nt")  numThreads      = intVal;
            }
        }
    }

    settingsInitialized = true;

    // ── Post-parse validation  ──────────────────────────────

    if (conditioningSet == false)
        Msg::error("Flag \"-cond\" is required (crown or origin).");

    int maxNumThreads = (int)std::thread::hardware_concurrency();
    if (maxNumThreads <= 0) {
        Msg::warning("Could not determine hardware thread count; defaulting to 1 thread.");
        maxNumThreads = 1;
    }

    if (numChains < 1) {
        Msg::warning("Chains must be >= 1; resetting to 1.");
        numChains = 1;
    }

    if (chainLength < 10) {
        Msg::warning("Chain length " + std::to_string(chainLength) + " is very short; resetting to 10. Expect non-convergence!");
        chainLength = 10;
    }

    if (printFrequency < 1) {
        Msg::warning("Print frequency must be >= 1; resetting to 1.");
        printFrequency = 1;
    }

    if (sampleFrequency < 1) {
        Msg::warning("Sample frequency must be >= 1; resetting to 1.");
        sampleFrequency = 1;
    }

    if (numThreads >= maxNumThreads) {
        Msg::warning("Requested " + std::to_string(numThreads) +
                     " threads, but only " + std::to_string(maxNumThreads) +
                     " available; capping at " + std::to_string(maxNumThreads - 1) + ".");
        numThreads = maxNumThreads - 1;
    }

    if (numThreads > numChains) {
        Msg::warning("Threads (" + std::to_string(numThreads) +
                     ") cannot exceed chains (" + std::to_string(numChains) +
                     "); reducing threads to match.");
        numThreads = numChains;
    }

    if (numThreads < 1) {
        Msg::warning("Threads must be >= 1; resetting to 1.");
        numThreads = 1;
    }

    if (treeOut.empty())
        Msg::warning("No tree output file specified (-to). Use -help for usage.");
    if (parametersOut.empty())
        Msg::warning("No parameter output file specified (-po). Use -help for usage.");
}

void UserSettings::parsePriorInto(const std::string& spec, Probability::PriorFamily& family, double& p1, double& p2) {
    size_t lp = spec.find('(');
    size_t rp = spec.find(')');
    if (lp == std::string::npos || rp == std::string::npos || rp <= lp)
        Msg::error("prior must look like exp(rate), gamma(shape,rate), lognormal(mu,sigma), unif(a,b), or truncnormal(mean,sd); got \"" + spec + "\".");
    std::string fam = spec.substr(0, lp);
    for (char& ch : fam) ch = std::tolower((unsigned char)ch);
    std::vector<double> ps;
    std::stringstream ss(spec.substr(lp + 1, rp - lp - 1));
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        try { ps.push_back(std::stod(tok)); }
        catch (...) { Msg::error("prior parameter \"" + tok + "\" is not a number in \"" + spec + "\"."); }
    }
    if (fam == "exp" || fam == "exponential") {
        if (ps.size() != 1 || ps[0] <= 0.0) Msg::error("exp prior needs one positive rate: exp(rate).");
        family = Probability::PriorFamily::EXPONENTIAL; p1 = ps[0];
    } else if (fam == "gamma") {
        if (ps.size() != 2 || ps[0] <= 0.0 || ps[1] <= 0.0) Msg::error("gamma prior needs positive shape,rate: gamma(shape,rate).");
        family = Probability::PriorFamily::GAMMA; p1 = ps[0]; p2 = ps[1];
    } else if (fam == "lognormal" || fam == "lnorm") {
        if (ps.size() != 2 || ps[1] <= 0.0) Msg::error("lognormal prior needs mu,sigma with sigma>0: lognormal(mu,sigma).");
        family = Probability::PriorFamily::LOGNORMAL; p1 = ps[0]; p2 = ps[1];
    } else if (fam == "unif" || fam == "uniform") {
        if (ps.size() != 2 || ps[0] >= ps[1]) Msg::error("unif prior needs a<b: unif(a,b).");
        family = Probability::PriorFamily::UNIFORM; p1 = ps[0]; p2 = ps[1];
    } else if (fam == "truncnormal" || fam == "tn" || fam == "normal") {
        if (ps.size() != 2 || ps[1] <= 0.0) Msg::error("truncnormal prior needs mean,sd with sd>0: truncnormal(mean,sd).");
        family = Probability::PriorFamily::TRUNCATED_NORMAL; p1 = ps[0]; p2 = ps[1];
    } else {
        Msg::error("Unknown prior family \"" + fam + "\". Use exp / gamma / lognormal / unif / truncnormal.");
    }
}

void UserSettings::print(void) {

    checkSettings();
    std::cout << "Tree input file name:                  " << treeFile << std::endl;
    std::cout << "Clades input file name:                " << cladesFile << std::endl;
    std::cout << "Fossils input file name:               " << fossilFile << std::endl;
    std::cout << "Conditioning scheme:                   " << (conditioning == Conditioning::CROWN ? "crown" : "origin") << std::endl;
    std::cout << "Conditioning age prior:                ";
    if (!conditionAgePriorSet) std::cout << "improper uniform";
    else switch (conditionAgePrior) {
        case Probability::PriorFamily::EXPONENTIAL:      std::cout << "exp(" << conditionAgePriorP1 << ")"; break;
        case Probability::PriorFamily::GAMMA:            std::cout << "gamma(" << conditionAgePriorP1 << "," << conditionAgePriorP2 << ")"; break;
        case Probability::PriorFamily::LOGNORMAL:        std::cout << "lognormal(" << conditionAgePriorP1 << "," << conditionAgePriorP2 << ")"; break;
        case Probability::PriorFamily::UNIFORM:          std::cout << "unif(" << conditionAgePriorP1 << "," << conditionAgePriorP2 << ")"; break;
        case Probability::PriorFamily::TRUNCATED_NORMAL: std::cout << "truncnormal(" << conditionAgePriorP1 << "," << conditionAgePriorP2 << ")"; break;
        case Probability::PriorFamily::IMPROPER:         std::cout << "improper"; break;
    }
    std::cout << std::endl;
    std::cout << "Model:                                 " << (model == Model::FBD ? "FBD" : (model == Model::HEA14 ? "HEA14" : "UFBD")) << std::endl;
    std::cout << "Tree output file name:                 " << treeOut << std::endl;
    std::cout << "Parameter output file name:            " << parametersOut << std::endl;
    std::cout << "Extant sampling fraction (rho):        " << rho << std::endl;
    if (seedSet)
        std::cout << "Random number seed:                    " << seed << std::endl;
    else
        std::cout << "Random number seed:                    (time-based)" << std::endl;
    std::cout << "Chain length:                          " << chainLength << std::endl;
    std::cout << "Number of chains:                      " << numChains << std::endl;
    std::cout << "Number of threads:                     " << numThreads << std::endl;
    std::cout << "Print-to-screen frequency:             " << printFrequency << std::endl;
    std::cout << "Chain sampling frequency:              " << sampleFrequency << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;

}

void UserSettings::printHelp(void){
    std::cout << "Usage: XXXX [options]\n";
    std::cout << "Options:\n";
    std::cout << "  TBD \n";
}
