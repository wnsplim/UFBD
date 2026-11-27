#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "ChainRunner.hpp"
#include "Msg.hpp"
#include "RandomVariable.hpp"
#include "PhylogeneticModel.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

std::ifstream ChainRunner::openCheckpoint(const std::string& path) {
    std::ifstream is(path);
    if(is.is_open() == false)
        Msg::error("-resume: cannot open the checkpoint " + path + ", so there is nothing to resume from.");
    return is;
}

bool ChainRunner::thinningWarned = false;

void ChainRunner::reconcileThinning(int storedSf, int& thinning) {
    if(storedSf != thinning){
        if(thinningWarned == false){
            Msg::warning("-thinning " + std::to_string(thinning) + " differs from the -thinning " + std::to_string(storedSf) + " the checkpoint was written with; the resumed run uses -thinning " + std::to_string(storedSf) + ".");
            thinningWarned = true;
        }
        thinning = storedSf;
        UserSettings::userSettings().setThinning(storedSf);
    }
}

void ChainRunner::writeDataShape(std::ostream& os, PhylogeneticModel* m) {
    os << "shape " << m->getTree()->getNumNodes() << ' ' << m->getNumFossils() << ' ' << m->getNumClockGroups() << '\n';
}

void ChainRunner::checkDataShape(std::istream& is, PhylogeneticModel* m) {
    std::streampos here = is.tellg();
    std::string tag;
    is >> tag;
    if(tag != "shape"){
        is.clear();
        is.seekg(here);
        return;
    }
    int nNode, nFossil, nGroup;
    is >> nNode >> nFossil >> nGroup;
    if(nNode != m->getTree()->getNumNodes() || nFossil != m->getNumFossils() || nGroup != m->getNumClockGroups())
        Msg::error("the input files do not match the checkpoint");
}

bool ChainRunner::checkpointIsCoupled(std::istream& is) {
    std::streampos here = is.tellg();
    std::string tag;
    is >> tag;
    is.clear();
    is.seekg(here);
    return tag != "shape";
}

void ChainRunner::requireCheckpointIntact(std::istream& is, const std::string& path) {
    if(is.fail())
        Msg::error("the checkpoint file " + path + " is truncated or corrupt.");
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
    fixedMask.assign(3, false);
    for(bool fb : resumeFixedMask())
        fixedMask.push_back(fb);

    std::ifstream pin(paramOut);
    std::string header;
    std::vector<std::string> keep;
    if(pin.is_open()){
        std::getline(pin, header);
        std::string line;
        while(std::getline(pin, line) && pin.eof() == false){
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
        while(std::getline(tin, tline) && tin.eof() == false){
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

// std::string ChainRunner::tunePath(const std::string& paramOut) {
//     std::string p = paramOut;
//     size_t dp = p.rfind(".log");
//     p.insert(dp != std::string::npos ? dp : p.size(), "_tune");
//     return p;
// }
// 
// void ChainRunner::runPsiRjTuning(PhylogeneticModel* m, unsigned long every, int label) {
//     rjTuning = true;
//     tuneOut = tunePath(paramOut);
//     bool resumed = UserSettings::userSettings().getResume() && std::ifstream(checkpointPath()).good();
//     if(resumed)
//         loadCheckpoint();
//     openTuneLog(m, resumed);
//     while(m->psiRjTuning())
//         advance(every);
//     tune.closeTSV();
//     std::ostringstream os;
//     os << std::setprecision(10) << "chain " << label << ": psi_rj tuning finished after " << gen << " generations, log prior odds " << m->getPsiRjState()[0] << "\n";
//     logLine(os.str());
//     rjTuning = false;
//     gen = 0;
// }
// 
// void ChainRunner::openTuneLog(PhylogeneticModel* m, bool resumed) {
//     std::vector<std::string> names = m->getParameterNames();
//     tuneIdx.clear();
//     for(const char* c : {"model", "delta", "log_prior_odds"})
//         tuneIdx.push_back((int)(std::find(names.begin(), names.end(), c) - names.begin()));
//     if(resumed == false){
//         tune.addFilepath(tuneOut, true);
//         tune.addColumnNamesTSV({"n", "posterior", "model", "delta", "log_prior_odds", "modification_factor"});
//         return;
//     }
//     std::ifstream in(tuneOut);
//     std::string header, line;
//     std::vector<std::string> keep;
//     std::getline(in, header);
//     while(std::getline(in, line) && in.eof() == false){
//         if(line.empty())
//             continue;
//         if((unsigned long)std::stod(line.substr(0, line.find('\t'))) > gen)
//             break;
//         keep.push_back(line);
//     }
//     in.close();
//     std::ofstream out(tuneOut, std::ios::out | std::ios::trunc);
//     out << header << '\n';
//     for(const std::string& l : keep)
//         out << l << '\n';
//     out.close();
//     tune.addFilepath(tuneOut, false);
// }
// 
// void ChainRunner::writeTuneRow(PhylogeneticModel* m, unsigned long n, double lnPost) {
//     std::vector<double> v = m->getParameterString(), st = m->getPsiRjState();
//     tune.appendDataTSV(std::vector<double>{(double)n, lnPost, v[tuneIdx[0]], v[tuneIdx[1]], v[tuneIdx[2]], std::expm1(st[1])});
// }
// 
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
        while(std::getline(lin, line) && lin.eof() == false){
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
