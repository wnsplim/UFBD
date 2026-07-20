#include <sstream>

#include "ChainRunner.hpp"
#include "Msg.hpp"
#include "RandomVariable.hpp"
#include "UserSettings.hpp"

std::ifstream ChainRunner::openCheckpoint(const std::string& path) {
    std::ifstream is(path);
    if(is.is_open() == false)
        Msg::error("could not open checkpoint file for -resume: " + path);
    return is;
}

void ChainRunner::reconcileThinning(int storedSf, int& thinning) {
    if(storedSf != thinning){
        Msg::warning("-thinning " + std::to_string(thinning) + " differs from the pre-resume thinning " + std::to_string(storedSf) + "; forcing " + std::to_string(storedSf));
        thinning = storedSf;
    }
}

void ChainRunner::requireCheckpointIntact(std::istream& is, const std::string& path) {
    if(is.fail())
        Msg::error("checkpoint file is truncated or corrupt: " + path);
}

void ChainRunner::resumeOutputs(void) {
    RandomVariable::setActiveInstance(resumeRng());

    std::vector<std::string> headStr = resumeParameterNames();
    traceNms.clear();
    traceNms.push_back("posterior");
    traceNms.push_back("likelihood");
    traceNms.push_back("prior");
    for(const std::string& s : headStr)
        traceNms.push_back(s);
    traceCols.assign(traceNms.size(), std::vector<double>());

    std::ifstream pin(paramOut);
    std::string header;
    std::vector<std::string> keep;
    if(pin.is_open()){
        std::getline(pin, header);
        std::string line;
        while(std::getline(pin, line)){
            if(line.empty())
                continue;
            std::stringstream ss(line);
            double nval;
            ss >> nval;
            if((unsigned long)nval > gen)
                break;
            keep.push_back(line);
        }
        pin.close();
    }

    for(const std::string& line : keep){
        std::stringstream ss(line);
        double nval;
        ss >> nval;
        double v;
        for(size_t j = 0; j < traceCols.size() && (ss >> v); j++)
            traceCols[j].push_back(v);
    }

    std::ofstream pout(paramOut, std::ios::out | std::ios::trunc);
    pout << header << '\n';
    for(const std::string& line : keep)
        pout << line << '\n';
    pout.close();

    params.addFilepath(paramOut, false);
    if(writeTrees){
        std::ifstream tin(treeOut);
        std::vector<std::string> tkeep;
        std::string tline;
        size_t nTree = 0;
        while(std::getline(tin, tline)){
            bool isTreeLine = (tline.find("TREE STATE_") != std::string::npos);
            if(isTreeLine){
                if(nTree < keep.size()){ tkeep.push_back(tline); nTree++; }
            }else if(tline != "END;"){
                tkeep.push_back(tline);
            }
        }
        tin.close();
        std::ofstream tout(treeOut, std::ios::out | std::ios::trunc);
        for(const std::string& l : tkeep)
            tout << l << '\n';
        tout.close();
        trees.addFilepath(treeOut, false);
    }

    resumeLatentOutput();
}

void ChainRunner::resumeLatentOutput(void) {
    std::vector<std::string> latNames = resumeLatentNames();
    if(latNames.empty() || UserSettings::userSettings().getWriteLatentLog() == false)
        return;

    latentOut = paramOut;
    size_t dp = latentOut.rfind(".log");
    latentOut.insert(dp != std::string::npos ? dp : latentOut.size(), "_latent");
    latentNms = latNames;
    latentCols.assign(latNames.size(), std::vector<double>());

    std::ifstream lin(latentOut);
    std::string header;
    std::vector<std::string> keep;
    if(lin.is_open()){
        std::getline(lin, header);
        std::string line;
        while(std::getline(lin, line)){
            if(line.empty())
                continue;
            std::stringstream ss(line);
            double nval;
            ss >> nval;
            if((unsigned long)nval > gen)
                break;
            keep.push_back(line);
        }
        lin.close();
    }

    for(const std::string& line : keep){
        std::stringstream ss(line);
        double nval;
        ss >> nval;
        double v;
        for(size_t j = 0; j < latentCols.size() && (ss >> v); j++)
            latentCols[j].push_back(v);
    }

    std::ofstream lout(latentOut, std::ios::out | std::ios::trunc);
    lout << header << '\n';
    for(const std::string& line : keep)
        lout << line << '\n';
    lout.close();

    latent.addFilepath(latentOut, false);
    writeLatent = true;
}
