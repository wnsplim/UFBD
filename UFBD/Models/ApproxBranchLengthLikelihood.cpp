#include "ApproxBranchLengthLikelihood.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ThreadPool.hpp"
#include "Tree.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>

ApproxBranchLengthLikelihood::ApproxBranchLengthLikelihood(const std::string& hessianFile, const std::string& mlTreeFile, const std::vector<std::string>& rogue, int partitionIndex, int nStates) :
    nb(0),
    nPartitions(0),
    crownBranchIdx(-1),
    partitionIndex(partitionIndex),
    cJc((nStates - 1.0) / nStates),
    rogueTaxa(rogue.begin(), rogue.end()),
    cachedRoot(nullptr)
{
    if(nStates < 2)
        Msg::error("ApproxBranchLengthLikelihood: nStates must be >= 2");

    readHessianFile(hessianFile);

    std::ifstream tf(mlTreeFile);
    if(tf.is_open() == false)
        Msg::error("ApproxBranchLengthLikelihood: cannot open ML tree file '" + mlTreeFile + "'");
    std::string line, newick;
    while(std::getline(tf, line))
        if(line.find('(') != std::string::npos){ newick = line; break; }
    tf.close();
    if(newick.empty())
        Msg::error("ApproxBranchLengthLikelihood: no Newick in ML tree file '" + mlTreeFile + "'");
    size_t semi = newick.rfind(';');
    if(semi != std::string::npos)
        newick = newick.substr(0, semi + 1);

    buildBranchOrder(newick);
}

std::set<std::string> ApproxBranchLengthLikelihood::backboneTipsBelow(Node* n, Tree* tree){
    std::set<std::string> names;
    if(n->getIsTip()){
        if(rogueTaxa.find(n->getName()) == rogueTaxa.end())
            names.insert(n->getName());
        return names;
    }
    for(Node* d : tree->getAllDescendants(n))
        if(d->getIsTip() && rogueTaxa.find(d->getName()) == rogueTaxa.end())
            names.insert(d->getName());
    return names;
}

std::set<std::string> ApproxBranchLengthLikelihood::canonicalize(const std::set<std::string>& clade){
    if(clade.size() <= backboneTaxa.size() / 2)
        return clade;
    std::set<std::string> comp;
    for(const std::string& t : backboneTaxa)
        if(clade.find(t) == clade.end())
            comp.insert(t);
    return comp;
}

void ApproxBranchLengthLikelihood::readHessianFile(const std::string& fn){
    std::ifstream f(fn);
    if(f.is_open() == false)
        Msg::error("ApproxBranchLengthLikelihood: cannot open hessian file '" + fn + "'");
    std::vector<std::string> lines;
    std::string line;
    while(std::getline(f, line))
        lines.push_back(line);
    f.close();

    auto nextNonEmpty = [&](size_t& idx) -> bool {
        while(idx < lines.size()){
            if(lines[idx].find_first_not_of(" \t\r\n") != std::string::npos)
                return true;
            idx++;
        }
        return false;
    };
    auto isHeaderAt = [&](size_t idx) -> bool {
        size_t i = idx;
        if(nextNonEmpty(i) == false) return false;
        std::istringstream iss(lines[i]);
        std::string tok, extra;
        if(!(iss >> tok)) return false;
        if(iss >> extra) return false;
        for(char c : tok)
            if(std::isdigit((unsigned char)c) == false) return false;
        size_t j = i + 1;
        if(nextNonEmpty(j) == false) return false;
        return lines[j].find('(') != std::string::npos;
    };
    auto parseNumbers = [&](const std::string& s, std::vector<double>& out){
        std::istringstream iss(s);
        std::string tok;
        while(iss >> tok){
            char* end = nullptr;
            double v = std::strtod(tok.c_str(), &end);
            if(end != tok.c_str() && *end == '\0')
                out.push_back(v);
        }
    };
    auto parseHeader = [&](size_t& idx, int& nsOut, std::string& newickOut){
        if(nextNonEmpty(idx) == false)
            Msg::error("ApproxBranchLengthLikelihood: missing ns in hessian file");
        std::istringstream iss(lines[idx]);
        if(!(iss >> nsOut))
            Msg::error("ApproxBranchLengthLikelihood: invalid ns in hessian file");
        idx++;
        newickOut.clear();
        while(idx < lines.size()){
            if(lines[idx].find('(') != std::string::npos){ newickOut = lines[idx]; idx++; break; }
            idx++;
        }
        if(newickOut.empty())
            Msg::error("ApproxBranchLengthLikelihood: no Newick in hessian file");
    };

    size_t idx = 0;
    int ns = 0;
    std::string newick;
    parseHeader(idx, ns, newick);
    if(ns < 3)
        Msg::error("ApproxBranchLengthLikelihood: invalid ns in hessian file");
    nb = 2 * ns - 3;

    hessianNewick = newick;

    blMle.clear();
    gradient.clear();
    hessian.clear();
    nPartitions = 0;

    while(idx < lines.size()){
        std::vector<double> bl, grad, hess;
        while((int)bl.size() < nb && idx < lines.size()){
            if(isHeaderAt(idx)) break;
            parseNumbers(lines[idx++], bl);
        }
        if(bl.empty()) break;
        if((int)bl.size() != nb)
            Msg::error("ApproxBranchLengthLikelihood: wrong MLE branch length count");
        while((int)grad.size() < nb && idx < lines.size()){
            if(isHeaderAt(idx)) break;
            parseNumbers(lines[idx++], grad);
        }
        if((int)grad.size() != nb)
            Msg::error("ApproxBranchLengthLikelihood: wrong gradient count");
        while((int)hess.size() < nb * nb && idx < lines.size()){
            if(isHeaderAt(idx)) break;
            parseNumbers(lines[idx++], hess);
        }
        if((int)hess.size() != nb * nb)
            Msg::error("ApproxBranchLengthLikelihood: wrong Hessian count");

        blMle.push_back(bl);
        gradient.push_back(grad);
        hessian.push_back(hess);
        applyArcsinTransform(nPartitions);
        nPartitions++;

        if(isHeaderAt(idx)){
            int nsNext = 0;
            std::string newickNext;
            parseHeader(idx, nsNext, newickNext);
            if(nsNext != ns)
                Msg::error("ApproxBranchLengthLikelihood: ns mismatch between partitions");
        }
    }

    if(partitionIndex > 0){
        int p = partitionIndex - 1;
        if(p >= nPartitions)
            Msg::error("ApproxBranchLengthLikelihood: partition index out of range");
        blMle = { blMle[p] };
        gradient = { gradient[p] };
        hessian = { hessian[p] };
        nPartitions = 1;
    }
}

void ApproxBranchLengthLikelihood::applyArcsinTransform(int p){
    std::vector<double> dbu(nb), dbu2(nb);
    for(int i = 0; i < nb; i++){
        double b = blMle[p][i];
        double pDist = cJc - cJc * std::exp(-b / cJc);
        double denom = 1 - pDist / cJc;
        dbu[i] = std::sqrt(pDist * (1 - pDist)) / denom;
        dbu2[i] = (1 - 2 * pDist) / (2 * denom) + dbu[i] * dbu[i] / cJc;
    }
    for(int i = 0; i < nb; i++){
        for(int j = 0; j < i; j++){
            hessian[p][i*nb+j] *= dbu[i] * dbu[j];
            hessian[p][j*nb+i] = hessian[p][i*nb+j];
        }
        hessian[p][i*nb+i] = hessian[p][i*nb+i] * dbu[i] * dbu[i] + gradient[p][i] * dbu2[i];
    }
    for(int i = 0; i < nb; i++)
        gradient[p][i] *= dbu[i];
    for(int i = 0; i < nb; i++){
        double b = blMle[p][i];
        blMle[p][i] = 2.0 * std::asin(std::sqrt(cJc - cJc * std::exp(-b / cJc)));
    }
}

void ApproxBranchLengthLikelihood::newickCanonBiparts(const std::string& nwk, std::set<std::set<std::string>>& out){
    struct PNode { std::vector<int> children; int parent; std::string name; std::set<std::string> tips; };
    std::vector<PNode> nd;
    auto newNode = [&](int par) -> int { PNode p; p.parent = par; nd.push_back(p); return (int)nd.size() - 1; };
    int root = newNode(-1);
    int cur = root;
    size_t i = 0;
    while(i < nwk.size()){
        char c = nwk[i];
        if(c == '('){ int ch = newNode(cur); nd[cur].children.push_back(ch); cur = ch; i++; }
        else if(c == ','){ int sib = newNode(nd[cur].parent); nd[nd[cur].parent].children.push_back(sib); cur = sib; i++; }
        else if(c == ')'){ cur = nd[cur].parent; i++; }
        else if(c == ':'){ i++; while(i < nwk.size() && nwk[i] != '(' && nwk[i] != ')' && nwk[i] != ',' && nwk[i] != ';') i++; }
        else if(c == ';'){ break; }
        else { size_t j = i; while(j < nwk.size() && nwk[j] != '(' && nwk[j] != ')' && nwk[j] != ',' && nwk[j] != ':' && nwk[j] != ';') j++; nd[cur].name = nwk.substr(i, j - i); i = j; }
    }
    for(int k = (int)nd.size() - 1; k >= 0; k--){
        if(nd[k].children.empty()) nd[k].tips.insert(nd[k].name);
        else for(int ch : nd[k].children) nd[k].tips.insert(nd[ch].tips.begin(), nd[ch].tips.end());
    }
    for(int k = 0; k < (int)nd.size(); k++)
        if(nd[k].parent != -1)
            out.insert(canonicalize(nd[k].tips));
}

void ApproxBranchLengthLikelihood::buildBranchOrder(const std::string& backboneNewick){
    struct PNode { std::vector<int> children; int parent; std::string name; std::set<std::string> tips; int ibranch; bool dead; };
    std::vector<PNode> nd;
    auto newNode = [&](int par) -> int { PNode p; p.parent = par; p.ibranch = -1; p.dead = false; nd.push_back(p); return (int)nd.size() - 1; };
    int root = newNode(-1);
    int cur = root;
    size_t i = 0;
    while(i < backboneNewick.size()){
        char c = backboneNewick[i];
        if(c == '('){ int ch = newNode(cur); nd[cur].children.push_back(ch); cur = ch; i++; }
        else if(c == ','){ int sib = newNode(nd[cur].parent); nd[nd[cur].parent].children.push_back(sib); cur = sib; i++; }
        else if(c == ')'){ cur = nd[cur].parent; i++; }
        else if(c == ':'){ i++; while(i < backboneNewick.size() && backboneNewick[i] != '(' && backboneNewick[i] != ')' && backboneNewick[i] != ',' && backboneNewick[i] != ';') i++; }
        else if(c == ';'){ break; }
        else { size_t j = i; while(j < backboneNewick.size() && backboneNewick[j] != '(' && backboneNewick[j] != ')' && backboneNewick[j] != ',' && backboneNewick[j] != ':' && backboneNewick[j] != ';') j++; nd[cur].name = backboneNewick.substr(i, j - i); i = j; }
    }

    for(int k = 0; k < (int)nd.size(); k++)
        if(nd[k].children.empty() && rogueTaxa.find(nd[k].name) != rogueTaxa.end())
            nd[k].dead = true;
    bool changed = true;
    while(changed){
        changed = false;
        for(int k = 0; k < (int)nd.size(); k++){
            if(nd[k].dead) continue;
            std::vector<int> live;
            for(int ch : nd[k].children)
                if(nd[ch].dead == false) live.push_back(ch);
            if(live.size() != nd[k].children.size()){ nd[k].children = live; changed = true; }
            if(nd[k].children.empty() && nd[k].name.empty()){ nd[k].dead = true; changed = true; }
            else if(nd[k].children.size() == 1 && k != root){
                int only = nd[k].children[0];
                int par = nd[k].parent;
                nd[only].parent = par;
                for(int& s : nd[par].children) if(s == k) s = only;
                nd[k].children.clear();
                nd[k].dead = true;
                changed = true;
            }
        }
    }
    while(nd[root].children.size() == 1){ int only = nd[root].children[0]; nd[only].parent = -1; nd[root].dead = true; root = only; }

    backboneTaxa.clear();
    for(int k = (int)nd.size() - 1; k >= 0; k--){
        if(nd[k].dead) continue;
        if(nd[k].children.empty()){ nd[k].tips.insert(nd[k].name); backboneTaxa.insert(nd[k].name); }
        else for(int ch : nd[k].children) nd[k].tips.insert(nd[ch].tips.begin(), nd[ch].tips.end());
    }

    if((int)nd[root].children.size() != 2)
        Msg::error("ApproxBranchLengthLikelihood: backbone tree must have a binary root to order branches");

    int counter = 0;
    std::vector<int> st;
    for(int k = (int)nd[root].children.size() - 1; k >= 0; k--) st.push_back(nd[root].children[k]);
    while(st.empty() == false){
        int n = st.back(); st.pop_back();
        nd[n].ibranch = counter++;
        for(int k = (int)nd[n].children.size() - 1; k >= 0; k--) st.push_back(nd[n].children[k]);
    }

    int son1 = nd[root].children[0];
    int son2 = nd[root].children[1];
    if(nd[son1].children.empty()){
        nd[son1].ibranch = nb - 1;
        for(int k = 0; k < (int)nd.size(); k++)
            if(k != root && k != son1 && k != son2 && nd[k].ibranch >= 0) nd[k].ibranch -= 2;
    }else{
        nd[son1].ibranch = nd[son2].ibranch - 1;
        for(int k = 0; k < (int)nd.size(); k++)
            if(k != root && k != son1 && k != son2 && nd[k].ibranch >= 0) nd[k].ibranch -= 1;
    }
    nd[son2].ibranch = -1;
    nd[root].ibranch = -1;

    bipartitions.assign(nb, std::set<std::string>());
    for(int k = 0; k < (int)nd.size(); k++){
        int b = nd[k].ibranch;
        if(b < 0) continue;
        bipartitions[b] = canonicalize(nd[k].tips);
    }
    crownBranchIdx = nd[son1].ibranch;

    std::set<std::set<std::string>> hessBip;
    newickCanonBiparts(hessianNewick, hessBip);
    for(int b = 0; b < nb; b++)
        if(hessBip.find(bipartitions[b]) == hessBip.end())
            Msg::error("ApproxBranchLengthLikelihood: backbone topology does not match the Hessian tree");
}

Node* ApproxBranchLengthLikelihood::findNodeByBipartition(const std::set<std::string>& bp, Tree* tree){
    if(bp.size() == 1)
        return tree->getTaxonNode(*bp.begin());
    for(Node* n : tree->getDownPassSequence()){
        if(n == tree->getRoot() || n->getIsTip() || tree->isBackboneNode(n) == false) continue;
        if(canonicalize(backboneTipsBelow(n, tree)) == bp)
            return n;
    }
    return nullptr;
}

double ApproxBranchLengthLikelihood::computeLnL(Tree* tree, const std::vector<std::vector<double>>& branchRates){
    tree->ensureBackboneCache();
    Node* curRoot = tree->getRoot();
    if(curRoot != cachedRoot){
        branchNodeIdx.assign(nb, -1);
        for(int i = 0; i < nb; i++){
            if(i == crownBranchIdx) continue;
            Node* n = findNodeByBipartition(bipartitions[i], tree);
            if(n == nullptr){ cachedRoot = nullptr; return -INFINITY; }
            branchNodeIdx[i] = n->getOffset();
        }
        cachedRoot = curRoot;
    }

    Node* mrca = tree->getCrown();
    const std::vector<Node*>& mrcaChildren = tree->getBackboneChildren(mrca);

    if(crownBranchIdx >= 0 && mrcaChildren.size() != 2)
        return -INFINITY;

    const double ninf = -std::numeric_limits<double>::infinity();
    std::vector<double> partLnL(nPartitions, 0.0);
    ThreadPool::current().parallelFor(OP_CTMC, nPartitions, [&](int q0, int q1){
        std::vector<double> z(nb);
        for(int p = q0; p < q1; p++){
            const std::vector<double>& br = branchRates[p];
            bool bad = false;
            for(int i = 0; i < nb; i++){
                double predBl;
                if(i == crownBranchIdx){
                    double d0 = mrca->getTime() - mrcaChildren[0]->getTime();
                    double d1 = mrca->getTime() - mrcaChildren[1]->getTime();
                    predBl = br[mrcaChildren[0]->getOffset()] * d0 + br[mrcaChildren[1]->getOffset()] * d1;
                }else{
                    Node* n = tree->getNodeByOffset(branchNodeIdx[i]);
                    double d = tree->getBackboneParent(n)->getTime() - n->getTime();
                    predBl = br[branchNodeIdx[i]] * d;
                }
                if(predBl <= 0.0){ bad = true; break; }
                z[i] = 2.0 * std::asin(std::sqrt(cJc - cJc * std::exp(-predBl / cJc))) - blMle[p][i];
            }
            if(bad){ partLnL[p] = ninf; continue; }
            double pl = 0.0;
            for(int i = 0; i < nb; i++)
                pl += z[i] * gradient[p][i];
            for(int i = 0; i < nb; i++)
                for(int j = 0; j < nb; j++)
                    pl += 0.5 * z[i] * hessian[p][i*nb+j] * z[j];
            partLnL[p] = pl;
        }
    });
    double lnL = 0.0;
    for(double v : partLnL){
        if(v == ninf){ lnL = ninf; break; }
        lnL += v;
    }
    return lnL;
}