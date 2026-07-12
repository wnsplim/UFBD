#include <cstdio>
#include "FBDInput.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "ReadTSV.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#include <cctype>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

FBDInput::FBDInput(std::string treePath, std::string cladesPath, std::string fossilPath){
    if(treePath.empty()){
        Msg::warning("no backbone tree supplied (-backbone_tree): running with an empty backbone (every taxa unresolved).");
        tree = new Tree("");
        clades.push_back(Clade("whole", std::vector<std::string>(), tree->getCrown(), tree->getCrown()->getAncestor()));
    }else{
        tree = readTree(treePath);
        if(cladesPath.empty() == false)
            readClades(cladesPath);
    }
    if(fossilPath.empty() == false){
        readFossils(fossilPath);
        assignFossilAwareAges();
    }
}

static void stripBracketComments(std::string& s){
    std::string out;
    int depth = 0;
    for(char c : s){
        if(c == '[') depth++;
        else if(c == ']'){ if(depth > 0) depth--; }
        else if(depth == 0) out += c;
    }
    s.swap(out);
}

static std::string applyTranslate(const std::string& nwk, const std::map<std::string,std::string>& tr){
    if(tr.empty()) return nwk;
    std::string out, tok;
    for(size_t i = 0; i <= nwk.size(); i++){
        char c = (i < nwk.size()) ? nwk[i] : '\0';
        if(c == '(' || c == ')' || c == ',' || c == ':' || c == ';' || c == '\0'){
            if(tok.empty() == false){
                std::map<std::string,std::string>::const_iterator it = tr.find(tok);
                out += (it != tr.end()) ? it->second : tok;
                tok.clear();
            }
            if(c != '\0') out += c;
        }else
            tok += c;
    }
    return out;
}

static std::string extractNewick(std::string content){
    stripBracketComments(content);
    std::string lower(content);
    for(char& ch : lower) ch = (char)std::tolower((unsigned char)ch);
    bool nexus = (lower.find("#nexus") != std::string::npos) || (lower.find("begin trees") != std::string::npos);

    std::map<std::string,std::string> tr;
    size_t searchStart = 0;
    if(nexus){
        searchStart = lower.find("begin trees");
        if(searchStart == std::string::npos) searchStart = 0;
        size_t tpos = lower.find("translate", searchStart);
        if(tpos != std::string::npos){
            size_t semi = content.find(';', tpos);
            std::string body = content.substr(tpos + 9, (semi == std::string::npos ? content.size() : semi) - (tpos + 9));
            std::stringstream bs(body);
            std::string entry;
            while(std::getline(bs, entry, ',')){
                std::stringstream es(entry);
                std::string key, name;
                es >> key >> name;
                if(key.empty() == false && name.empty() == false) tr[key] = name;
            }
            if(semi != std::string::npos) searchStart = semi + 1;
        }
    }
    size_t lp = content.find('(', searchStart);
    size_t semi = (lp == std::string::npos) ? std::string::npos : content.find(';', lp);
    if(lp == std::string::npos || semi == std::string::npos)
        return "";
    return applyTranslate(content.substr(lp, semi - lp + 1), tr);
}

Tree* FBDInput::readTree(std::string path){
    std::ifstream file(path);
    if(file.is_open() == false)
        Msg::error("could not open backbone tree file '" + path + "'");
    std::stringstream ss;
    ss << file.rdbuf();
    std::string newick = extractNewick(ss.str());
    if(Tree::isValidNewick(newick) == false)
        Msg::error("backbone tree file '" + path + "' has no valid tree (use NEWICK or NEXUS)");
    Tree* t = new Tree(newick);
    t->validateBackbone();
    return t;
}

static std::vector<std::string> splitOnComma(const std::string& s){
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string cell;
    while(std::getline(ss, cell, ','))
        out.push_back(cell);
    return out;
}

void FBDInput::readClades(std::string path){
    ReadTSV reader(path, false, false, true);
    std::vector<std::vector<std::string>> rows = reader.getReadStringData();
    for(std::vector<std::string>& row : rows){
        if(row.size() < 2)
            Msg::error("clade row needs a name and a comma-separated taxon list");
        std::string name = row[0];
        std::vector<std::string> taxa = splitOnComma(row[1]);
        if(taxa.empty()){
            Msg::warning("clade '" + name + "' has no taxa; ignoring it");
            continue;
        }
        Node* crown = tree->getMRCA(taxa);
        for(Clade& c : clades)
            if(c.getCrown() == crown)
                Msg::error("clade '" + name + "' has the same MRCA as clade '" + c.getName() + "'");
        clades.push_back(Clade(name, taxa, crown, crown->getAncestor()));
    }
}

void FBDInput::readFossils(std::string path){
    ReadTSV reader(path, false, false, true);
    std::vector<std::vector<std::string>> rows = reader.getReadStringData();
    for(std::vector<std::string>& row : rows){
        if(row.size() < 5)
            Msg::error("fossil row needs: taxon, min_age, max_age, clade, assignment");
        std::string taxon = row[0];
        double minAge = 0.0;
        double maxAge = 0.0;
        try{
            minAge = std::stod(row[1]);
            maxAge = std::stod(row[2]);
        }catch(...){
            Msg::error("fossil '" + taxon + "' has a non-numeric age");
        }
        std::string cladeName = row[3];
        std::string assignStr = row[4];
        for(char& ch : assignStr)
            ch = std::toupper((unsigned char)ch);

        if(minAge > maxAge)
            Msg::error("fossil '" + taxon + "' has min_age greater than max_age");
        if(minAge < 0.0)
            Msg::error("fossil '" + taxon + "' has a negative age");

        Clade* clade = nullptr;
        for(Clade& c : clades)
            if(c.getName() == cladeName){
                clade = &c;
                break;
            }
        if(clade == nullptr)
            Msg::error("fossil '" + taxon + "' references undefined clade '" + cladeName + "'");

        if(assignStr != "CROWN" && assignStr != "TOTAL" && assignStr != "STEM")
            Msg::error("fossil '" + taxon + "' has unknown assignment '" + assignStr + "'");
        Assignment assignment = (assignStr == "CROWN") ? Assignment::CROWN
                              : (assignStr == "STEM")  ? Assignment::STEM
                              :                          Assignment::TOTAL;

        if(assignment == Assignment::STEM && maxAge == 0.0)
            Msg::error("fossil '" + taxon + "' is an unresolved extant assigned STEM; UE can only be CROWN or TOTAL.");

        if(assignment == Assignment::CROWN && clade->getCrown()->getIsTip()){
            Msg::warning("fossil '" + taxon + "' assigned CROWN to singleton clade '" + cladeName + "'; treating as TOTAL");
            assignment = Assignment::TOTAL;
        }
        if(assignment == Assignment::STEM && clade->getCrown()->getIsTip()){
            Msg::warning("fossil '" + taxon + "' assigned STEM to singleton clade '" + cladeName + "'; treating as TOTAL");
            assignment = Assignment::TOTAL;
        }
        if(assignment == Assignment::STEM && clade->getCrown() == tree->getCrown()
           && UserSettings::userSettings().getConditioning() == Conditioning::CROWN)
            Msg::error("conditioning on crown node but fossil '" + taxon + "' is STEM on the whole-tree clade '" + cladeName + "'.");
        if(UserSettings::userSettings().getConditioning() == Conditioning::CROWN
           && assignment == Assignment::TOTAL && clade->getCrown() == tree->getCrown()){
            Msg::warning("conditioning on crown node but fossil '" + taxon + "' is TOTAL on the whole-tree clade '" + cladeName + "'; treating as CROWN");
            assignment = Assignment::CROWN;
        }

        std::string typeName = (row.size() >= 6) ? row[5] : "";
        fossils.push_back(Fossil(taxon, minAge, maxAge, cladeName, assignment, typeName));
    }
}

void FBDInput::assignFossilAwareAges(void){
    double maxBound = 0.0;
    for(Fossil& f : fossils)
        if(f.getMaxAge() > maxBound)
            maxBound = f.getMaxAge();
    int numInternal = tree->getNumNodes() - tree->getNumBackbone();
    double unit = (maxBound > 0.0) ? maxBound / (numInternal + 1) : 1.0;

    std::map<Node*,double> minAges;
    for(Fossil& f : fossils){
        Clade* clade = nullptr;
        for(Clade& c : clades)
            if(c.getName() == f.getClade()){
                clade = &c;
                break;
            }
        Node* anchor = (f.getAssignment() == Assignment::CROWN) ? clade->getCrown() : clade->getOrigin();
        double bound = f.getMaxAge() * 1.05;
        if(minAges.count(anchor) == 0 || bound > minAges[anchor])
            minAges[anchor] = bound;
    }
    tree->assignStartingAges(minAges, unit);

    if(getenv("FBD_CHK_INIT") != nullptr){
        fprintf(stderr, "[init] maxFossil=%g numInternal=%d unit=%g  minAges:%d anchors\n", maxBound, numInternal, unit, (int)minAges.size());
        for(std::map<Node*,double>::iterator it = minAges.begin(); it != minAges.end(); ++it)
            fprintf(stderr, "[init]   anchor off=%d  floor=%g  -> age after assign=%g\n", it->first->getOffset(), it->second, it->first->getTime());
        fprintf(stderr, "[init] crown off=%d age=%g | crown->anc off=%d age=%g | root off=%d age=%g\n",
                tree->getCrown()->getOffset(), tree->getCrown()->getTime(),
                (tree->getCrown()->getAncestor() != nullptr) ? tree->getCrown()->getAncestor()->getOffset() : -1,
                (tree->getCrown()->getAncestor() != nullptr) ? tree->getCrown()->getAncestor()->getTime() : -1.0,
                tree->getRoot()->getOffset(), tree->getRoot()->getTime());
    }

    UserSettings& us = UserSettings::userSettings();
    if(us.getConditionAgePriorSet()){
        bool origin  = (us.getConditioning() == Conditioning::ORIGIN);
        bool bounded = (us.getConditionAgePrior() == Probability::PriorFamily::UNIFORM);
        double hi = bounded ? us.getConditionAgePriorP2() : std::numeric_limits<double>::max();
        if(bounded && hi < maxBound)
            Msg::error("conditioning age prior lower bound (" + std::to_string(hi) + ") is younger than the oldest fossil (" + std::to_string(maxBound) + ")");
        double mean = Probability::priorMean(us.getConditionAgePrior(), us.getConditionAgePriorP1(), us.getConditionAgePriorP2());
        double target = origin ? 0.9 * mean : mean;
        if(target < maxBound)
            target = bounded ? hi : (maxBound * 1.05);
        double cur = tree->getCrown()->getTime();
        if(cur > 0.0)
            tree->scaleInternalAges(target / cur);
    }
}
