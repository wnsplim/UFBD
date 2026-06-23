#include "AlignmentReader.hpp"
#include "Msg.hpp"

#include <cctype>
#include <fstream>
#include <map>
#include <sstream>

static int nucleotideBitmask(char c){
    switch(std::toupper((unsigned char)c)){
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 4;
        case 'T': return 8;
        case 'R': return 5;
        case 'Y': return 10;
        case 'S': return 6;
        case 'W': return 9;
        case 'K': return 12;
        case 'M': return 3;
        case 'B': return 14;
        case 'D': return 13;
        case 'H': return 11;
        case 'V': return 7;
        default:  return 15;
    }
}

AlignmentReader::AlignmentReader(const std::string& alignmentFile, const std::string& partitionFile, int numStates){
    readAlignment(alignmentFile);
    int nsite = matrix.empty() ? 0 : (int)matrix[0].size();
    if(partitionFile.empty()){
        partitionSites.resize(1);
        for(int s = 0; s < nsite; s++)
            partitionSites[0].push_back(s);
    }else{
        readPartitions(partitionFile);
    }
    compress();
}

void AlignmentReader::readAlignment(const std::string& file){
    std::ifstream f(file);
    if(f.is_open() == false)
        Msg::error("could not open alignment file: " + file);
    int c = f.peek();
    while(c == ' ' || c == '\t' || c == '\n' || c == '\r'){ f.get(); c = f.peek(); }
    if(c == '>')
        readFasta(f);
    else if(c == '#')
        readNexus(f);
    else if(std::isdigit(c))
        readPhylip(f);
    else
        Msg::error("unrecognized alignment format (use PHYLIP, FASTA, or NEXUS): " + file);

    int nsite = matrix.empty() ? 0 : (int)matrix[0].size();
    for(const std::vector<int>& row : matrix)
        if((int)row.size() != nsite)
            Msg::error("alignment rows have unequal length (not aligned): " + file);
}

void AlignmentReader::readPhylip(std::ifstream& f){
    int ntax = 0, nsite = 0;
    f >> ntax >> nsite;
    taxa.resize(ntax);
    matrix.assign(ntax, std::vector<int>(nsite));
    for(int i = 0; i < ntax; i++){
        std::string seq;
        f >> taxa[i] >> seq;
        if((int)seq.size() != nsite)
            Msg::error("alignment row for '" + taxa[i] + "' has wrong length");
        for(int s = 0; s < nsite; s++)
            matrix[i][s] = nucleotideBitmask(seq[s]);
    }
}

void AlignmentReader::readFasta(std::ifstream& f){
    std::string line, name, seq;
    auto flush = [&](){
        if(name.empty()) return;
        taxa.push_back(name);
        std::vector<int> row(seq.size());
        for(size_t s = 0; s < seq.size(); s++)
            row[s] = nucleotideBitmask(seq[s]);
        matrix.push_back(row);
        seq.clear();
    };
    while(std::getline(f, line)){
        if(line.empty()) continue;
        if(line[0] == '>'){
            flush();
            name = line.substr(1);
            size_t sp = name.find_first_of(" \t\r");
            if(sp != std::string::npos) name = name.substr(0, sp);
        }else{
            for(char c : line)
                if(std::isspace((unsigned char)c) == false) seq += c;
        }
    }
    flush();
}

void AlignmentReader::readNexus(std::ifstream& f){
    std::stringstream ss;
    ss << f.rdbuf();
    std::string text = ss.str();
    std::string clean;
    clean.reserve(text.size());
    int depth = 0;
    for(char c : text){
        if(c == '[') depth++;
        else if(c == ']'){ if(depth > 0) depth--; }
        else if(depth == 0) clean += c;
    }
    std::string low = clean;
    for(char& c : low) c = std::tolower((unsigned char)c);
    size_t mpos = low.find("matrix");
    if(mpos == std::string::npos)
        Msg::error("nexus alignment has no matrix block");
    size_t start = mpos + 6;
    size_t endpos = clean.find(';', start);
    if(endpos == std::string::npos)
        Msg::error("nexus matrix block is not terminated");
    std::string block = clean.substr(start, endpos - start);

    std::stringstream bs(block);
    std::string line;
    std::map<std::string, std::string> seqs;
    std::vector<std::string> order;
    while(std::getline(bs, line)){
        std::stringstream ls(line);
        std::string name;
        if(!(ls >> name)) continue;
        std::string chunk, part;
        while(ls >> part) chunk += part;
        if(chunk.empty()) continue;
        if(seqs.find(name) == seqs.end()) order.push_back(name);
        seqs[name] += chunk;
    }
    for(const std::string& nm : order){
        taxa.push_back(nm);
        const std::string& s = seqs[nm];
        std::vector<int> row(s.size());
        for(size_t i = 0; i < s.size(); i++)
            row[i] = nucleotideBitmask(s[i]);
        matrix.push_back(row);
    }
}

static void parseRange(const std::string& r, int nsite, std::vector<int>& sites){
    size_t dash = r.find('-');
    if(dash == std::string::npos){
        int v = std::atoi(r.c_str());
        if(v >= 1 && v <= nsite) sites.push_back(v - 1);
        return;
    }
    int start = std::atoi(r.substr(0, dash).c_str());
    std::string rest = r.substr(dash + 1);
    int stride = 1;
    size_t bs = rest.find('\\');
    int end;
    if(bs != std::string::npos){
        end = std::atoi(rest.substr(0, bs).c_str());
        stride = std::atoi(rest.substr(bs + 1).c_str());
    }else{
        end = std::atoi(rest.c_str());
    }
    if(stride < 1) stride = 1;
    for(int s = start; s <= end; s += stride)
        if(s >= 1 && s <= nsite) sites.push_back(s - 1);
}

void AlignmentReader::readPartitions(const std::string& file){
    std::ifstream f(file);
    if(f.is_open() == false)
        Msg::error("could not open partition file: " + file);
    std::stringstream ss;
    ss << f.rdbuf();
    std::string text = ss.str();
    int nsite = matrix.empty() ? 0 : (int)matrix[0].size();
    for(char& c : text)
        if(c == ';' || c == '=' || c == ',') c = ' ';
    std::stringstream ts(text);
    std::vector<std::string> tok;
    std::string t;
    while(ts >> t) tok.push_back(t);
    for(size_t i = 0; i < tok.size(); i++){
        std::string low = tok[i];
        for(char& c : low) c = std::tolower((unsigned char)c);
        if(low != "charset") continue;
        std::vector<int> sites;
        size_t j = i + 2;
        for(; j < tok.size(); j++){
            std::string lj = tok[j];
            for(char& c : lj) c = std::tolower((unsigned char)c);
            if(lj == "charset" || lj == "charpartition" || lj == "end" || lj == "begin")
                break;
            if(tok[j].find_first_of("0123456789") != std::string::npos)
                parseRange(tok[j], nsite, sites);
        }
        partitionSites.push_back(sites);
        i = j - 1;
    }
    if(partitionSites.empty())
        Msg::error("no charset found in partition file: " + file);
}

void AlignmentReader::compress(void){
    int ntax = (int)taxa.size();
    int nP = (int)partitionSites.size();
    partitionPatterns.assign(nP, std::vector<std::vector<int>>(ntax));
    partitionWeights.assign(nP, std::vector<int>());
    for(int p = 0; p < nP; p++){
        std::map<std::vector<int>, int> seen;
        for(int site : partitionSites[p]){
            std::vector<int> col(ntax);
            for(int t = 0; t < ntax; t++)
                col[t] = matrix[t][site];
            std::map<std::vector<int>, int>::iterator it = seen.find(col);
            if(it == seen.end()){
                seen[col] = (int)partitionWeights[p].size();
                for(int t = 0; t < ntax; t++)
                    partitionPatterns[p][t].push_back(col[t]);
                partitionWeights[p].push_back(1);
            }else{
                partitionWeights[p][it->second]++;
            }
        }
    }
}
