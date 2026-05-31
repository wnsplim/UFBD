#include "FBDInput.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ReadTSV.hpp"
#include "Tree.hpp"

#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

FBDInput::FBDInput(std::string treePath, std::string cladesPath, std::string fossilPath){
    tree = readTree(treePath);
    readClades(cladesPath);
    readFossils(fossilPath);
}

Tree* FBDInput::readTree(std::string path){
    std::ifstream file(path);
    if(file.is_open() == false)
        Msg::error("could not open tree file: " + path);
    std::stringstream ss;
    ss << file.rdbuf();
    std::string newick = ss.str();
    if(Tree::isValidNewick(newick) == false)
        Msg::error("tree file is not valid newick: " + path);
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
        if(taxa.size() < 2)
            Msg::error("clade '" + name + "' must be defined by at least two taxa");
        Node* crown = tree->getMRCA(taxa);
        for(Clade& c : clades)
            if(c.getCrown() == crown)
                Msg::error("clade '" + name + "' has the same MRCA as clade '" + c.getName() + "'");
        clades.push_back(Clade(name, crown, crown->getAncestor()));
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

        bool cladeFound = false;
        for(Clade& c : clades)
            if(c.getName() == cladeName)
                cladeFound = true;
        if(cladeFound == false)
            Msg::error("fossil '" + taxon + "' references undefined clade '" + cladeName + "'");

        if(assignStr != "CROWN" && assignStr != "TOTAL")
            Msg::error("fossil '" + taxon + "' has unknown assignment '" + assignStr + "'");
        Assignment assignment = (assignStr == "CROWN") ? Assignment::CROWN : Assignment::TOTAL;

        fossils.push_back(Fossil(taxon, minAge, maxAge, cladeName, assignment));
    }
}
