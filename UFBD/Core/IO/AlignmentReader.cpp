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
        case 'U': return 8;
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

static int aminoAcidBitmask(char c){
    static const char order[21] = "ARNDCQEGHILKMFPSTWYV";
    char u = (char)std::toupper((unsigned char)c);
    if(u == 'B') return aminoAcidBitmask('N') | aminoAcidBitmask('D');
    if(u == 'Z') return aminoAcidBitmask('Q') | aminoAcidBitmask('E');
    if(u == 'J') return aminoAcidBitmask('I') | aminoAcidBitmask('L');
    for(int i = 0; i < 20; i++)
        if(order[i] == u) return 1 << i;
    return (1 << 20) - 1;
}

int AlignmentReader::charBitmask(char c) const {
    return numStates == 20 ? aminoAcidBitmask(c) : nucleotideBitmask(c);
}

static std::string takeLabel(const std::string& s, size_t& pos){
    while(pos < s.size() && std::isspace((unsigned char)s[pos])) pos++;
    std::string label;
    if(pos < s.size() && s[pos] == '\''){
        pos++;
        while(pos < s.size()){
            if(s[pos] == '\''){
                if(pos + 1 < s.size() && s[pos + 1] == '\''){ label += '\''; pos += 2; }
                else { pos++; break; }
            }else
                label += s[pos++];
        }
    }else{
        while(pos < s.size() && std::isspace((unsigned char)s[pos]) == false)
            label += s[pos++];
    }
    return label;
}

static void appendSeqChars(std::string& dest, const std::string& s, size_t from){
    for(size_t i = from; i < s.size(); i++)
        if(std::isspace((unsigned char)s[i]) == false)
            dest += s[i];
}

AlignmentReader::AlignmentReader(const std::string& alignmentFile, const std::string& partitionFile, int numStates){
    this->numStates = numStates;
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
            Msg::error("alignment rows have unequal length: " + file);
}

void AlignmentReader::readPhylip(std::ifstream& f){
    int ntax = 0, nsite = 0;
    std::string line;
    while(std::getline(f, line)){
        std::stringstream hs(line);
        if(hs >> ntax >> nsite) break;
        for(char c : line) if(std::isspace((unsigned char)c) == false){ ntax = 0; break; }
        if(ntax == 0) Msg::error("invalid PHYLIP header: expected '<ntax> <nsite>'");
    }
    if(ntax <= 0 || nsite <= 0)
        Msg::error("invalid PHYLIP header: expected '<ntax> <nsite>'");

    std::vector<std::string> lines;
    while(std::getline(f, line)){
        size_t p = 0;
        while(p < line.size() && std::isspace((unsigned char)line[p])) p++;
        if(p < line.size()) lines.push_back(line);
    }

    bool sequential = ((int)lines.size() == ntax);
    if(sequential == false){
        bool block = ((int)lines.size() >= ntax) && ((int)lines.size() % ntax == 0);
        bool uniform = true;
        int frag0 = -1;
        if(block){
            for(int i = 0; i < ntax; i++){
                size_t pos = 0;
                takeLabel(lines[i], pos);
                int len = 0;
                for(size_t k = pos; k < lines[i].size(); k++)
                    if(std::isspace((unsigned char)lines[i][k]) == false) len++;
                if(i == 0) frag0 = len;
                else if(len != frag0){ uniform = false; break; }
            }
        }
        sequential = (block && uniform && frag0 < nsite) == false;
    }

    taxa.assign(ntax, "");
    std::vector<std::string> seq(ntax);
    if(sequential){
        int id = 0;
        for(const std::string& l : lines){
            if(id >= ntax)
                Msg::error("PHYLIP alignment has more sequence data than " + std::to_string(ntax) + " taxa");
            size_t pos = 0;
            if(seq[id].empty()) taxa[id] = takeLabel(l, pos);
            appendSeqChars(seq[id], l, pos);
            if((int)seq[id].size() >= nsite) id++;
        }
    }else{
        int nblock = (int)lines.size() / ntax;
        for(int b = 0; b < nblock; b++)
            for(int i = 0; i < ntax; i++){
                const std::string& l = lines[b * ntax + i];
                size_t pos = 0;
                if(b == 0) taxa[i] = takeLabel(l, pos);
                appendSeqChars(seq[i], l, pos);
            }
    }

    matrix.assign(ntax, std::vector<int>(nsite));
    for(int i = 0; i < ntax; i++){
        if((int)seq[i].size() != nsite)
            Msg::error("PHYLIP sequence for '" + taxa[i] + "' has length " + std::to_string(seq[i].size()) + " (expected " + std::to_string(nsite) + ")");
        for(int s = 0; s < nsite; s++)
            matrix[i][s] = charBitmask(seq[i][s]);
    }
}

void AlignmentReader::readFasta(std::ifstream& f){
    std::string line, name, seq;
    auto flush = [&](){
        if(name.empty()) return;
        taxa.push_back(name);
        std::vector<int> row(seq.size());
        for(size_t s = 0; s < seq.size(); s++)
            row[s] = charBitmask(seq[s]);
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

    int nchar = 0;
    size_t dpos = low.find("nchar");
    if(dpos != std::string::npos){
        size_t eq = clean.find('=', dpos);
        if(eq != std::string::npos) nchar = std::atoi(clean.c_str() + eq + 1);
    }

    char matchChar = 0;
    bool interleave = false;
    size_t fpos = low.find("format");
    if(fpos != std::string::npos){
        size_t fend = low.find(';', fpos);
        std::string fseg = clean.substr(fpos, (fend == std::string::npos ? clean.size() : fend) - fpos);
        std::string flow = fseg;
        for(char& c : flow) c = std::tolower((unsigned char)c);
        size_t mc = flow.find("matchchar");
        if(mc != std::string::npos){
            size_t eq = fseg.find('=', mc);
            if(eq != std::string::npos){
                size_t p = eq + 1;
                while(p < fseg.size() && std::isspace((unsigned char)fseg[p])) p++;
                if(p < fseg.size())
                    matchChar = (fseg[p] == '\'' && p + 1 < fseg.size()) ? fseg[p + 1] : fseg[p];
            }
        }
        size_t ip = flow.find("interleave");
        if(ip != std::string::npos){
            interleave = true;
            size_t eq = flow.find('=', ip);
            if(eq != std::string::npos){
                size_t p = eq + 1;
                while(p < flow.size() && std::isspace((unsigned char)flow[p])) p++;
                if(flow.compare(p, 2, "no") == 0) interleave = false;
            }
        }
    }

    size_t mpos = low.find("matrix");
    if(mpos == std::string::npos)
        Msg::error("nexus alignment has no matrix block");
    size_t start = mpos + 6;
    size_t endpos = clean.find(';', start);
    if(endpos == std::string::npos)
        Msg::error("nexus matrix block is not terminated");
    std::string block = clean.substr(start, endpos - start);

    std::map<std::string, std::string> seqs;
    std::vector<std::string> order;
    if(interleave == false && nchar > 0){
        size_t pos = 0;
        while(true){
            std::string name = takeLabel(block, pos);
            if(name.empty()) break;
            std::string chunk;
            while((int)chunk.size() < nchar && pos < block.size()){
                if(std::isspace((unsigned char)block[pos]) == false) chunk += block[pos];
                pos++;
            }
            order.push_back(name);
            seqs[name] = chunk;
        }
    }else{
        std::stringstream bs(block);
        std::string line;
        while(std::getline(bs, line)){
            size_t pos = 0;
            std::string name = takeLabel(line, pos);
            if(name.empty()) continue;
            std::string chunk;
            appendSeqChars(chunk, line, pos);
            if(chunk.empty()) continue;
            if(seqs.find(name) == seqs.end()) order.push_back(name);
            seqs[name] += chunk;
        }
    }

    if(matchChar != 0 && order.empty() == false){
        const std::string first = seqs[order[0]];
        for(size_t k = 1; k < order.size(); k++){
            std::string& s = seqs[order[k]];
            for(size_t j = 0; j < s.size() && j < first.size(); j++)
                if(s[j] == matchChar) s[j] = first[j];
        }
    }

    for(const std::string& nm : order){
        taxa.push_back(nm);
        const std::string& s = seqs[nm];
        std::vector<int> row(s.size());
        for(size_t i = 0; i < s.size(); i++)
            row[i] = charBitmask(s[i]);
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
        std::string name = (i + 1 < tok.size()) ? tok[i + 1] : "";
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
        partitionNames.push_back(name);
        partitionSites.push_back(sites);
        i = j - 1;
    }
    if(partitionSites.empty())
        Msg::error("no charset found in partition file: " + file);
    parseClockPartition(tok);
}

void AlignmentReader::parseClockPartition(const std::vector<std::string>& tok){
    int nPart = (int)partitionSites.size();
    for(size_t i = 0; i + 1 < tok.size(); i++){
        std::string low = tok[i];
        for(char& c : low) c = std::tolower((unsigned char)c);
        if(low != "charpartition") continue;
        std::string scheme = tok[i + 1];
        for(char& c : scheme) c = std::tolower((unsigned char)c);
        if(scheme != "clock" && scheme != "clocks") continue;
        partitionGroup.assign(nPart, -1);
        std::map<std::string, int> labelToGroup;
        int curGroup = -1;
        for(size_t j = i + 2; j < tok.size(); j++){
            std::string ll = tok[j];
            for(char& c : ll) c = std::tolower((unsigned char)c);
            if(ll == "charset" || ll == "charpartition" || ll == "end" || ll == "begin")
                break;
            std::string label, part;
            size_t colon = tok[j].find(':');
            if(colon != std::string::npos){
                label = tok[j].substr(0, colon);
                part  = tok[j].substr(colon + 1);
            }else{
                part = tok[j];
            }
            if(label.empty() == false){
                if(labelToGroup.count(label) == 0)
                    labelToGroup[label] = (int)labelToGroup.size();
                curGroup = labelToGroup[label];
            }
            if(part.empty() == false){
                int idx = -1;
                for(int k = 0; k < nPart; k++)
                    if(partitionNames[k] == part){ idx = k; break; }
                if(idx < 0)
                    Msg::error("charpartition clock references undefined charset '" + part + "'.");
                if(curGroup < 0)
                    Msg::error("charpartition clock: charset '" + part + "' precedes any 'label:' clock-group.");
                partitionGroup[idx] = curGroup;
            }
        }
        int nextGroup = (int)labelToGroup.size();
        for(int k = 0; k < nPart; k++)
            if(partitionGroup[k] < 0)
                partitionGroup[k] = nextGroup++;
        std::map<int, int> remap;
        for(int k = 0; k < nPart; k++){
            if(remap.count(partitionGroup[k]) == 0)
                remap[partitionGroup[k]] = (int)remap.size();
            partitionGroup[k] = remap[partitionGroup[k]];
        }
        break;
    }
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
