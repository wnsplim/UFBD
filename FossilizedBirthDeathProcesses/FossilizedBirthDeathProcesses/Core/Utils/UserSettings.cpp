#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>

#include "Msg.hpp"
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
    chainLength     = 100;
    numChains       = 4;
    numThreads      = 4;
    printFrequency  = 1;
    sampleFrequency = 1;

    std::vector<std::string> arguments;
    for (int i = 0; i < argc; i++)
        arguments.push_back(std::string(argv[i]));

    executablePath = arguments[0];

    // Known flags and whether they take a value
    std::set<std::string> knownFlags = {
        "-to", "-po", "-t", "-c", "-f", "-n", "-p", "-s", "-nc", "-nt", "-help", "-h"
    };
    std::set<std::string> valueFlags = {
        "-to", "-po", "-t", "-c", "-f", "-n", "-p", "-s", "-nc", "-nt"
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

void UserSettings::print(void) {

    checkSettings();
    std::cout << "Tree input file name:                  " << treeFile << std::endl;
    std::cout << "Clades input file name:                " << cladesFile << std::endl;
    std::cout << "Fossils input file name:               " << fossilFile << std::endl;
    std::cout << "Tree output file name:                 " << treeOut << std::endl;
    std::cout << "Parameter output file name:            " << parametersOut << std::endl;
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
