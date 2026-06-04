#include "ApproxBranchLengthLikelihood.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "Tree.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>

ApproxBranchLengthLikelihood::ApproxBranchLengthLikelihood(const std::string& hessianFile, const std::string& mlTreeFile, const std::vector<std::string>& rogue, int partitionIndex, int nStates) :
    nb(0),
    nPartitions(0),
    rootBranchIdx(-1),
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

    Tree* rootedMl = new Tree(newick);
    identifyRootBranch(rootedMl);
    delete rootedMl;
}

std::set<std::string> ApproxBranchLengthLikelihood::molecularDescendants(Node* n, Tree* tree){
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
            if(tok == "Hessian" || tok == "HESSIAN" || tok == "hessian") continue;
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

    Tree* molTree = new Tree(newick);
    backboneTaxa.clear();
    for(Node* n : molTree->getDownPassSequence())
        if(n->getIsTip())
            backboneTaxa.insert(n->getName());
    buildBipartitions(molTree);
    delete molTree;

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
        double u = 2.0 * std::asin(std::sqrt(cJc - cJc * std::exp(-b / cJc)));
        double sin2 = std::sin(u / 2.0);
        double cos2 = std::cos(u / 2.0);
        double denom = 1.0 - sin2 * sin2 / cJc;
        dbu[i] = sin2 * cos2 / denom;
        dbu2[i] = (cos2 * cos2 - sin2 * sin2) / 2.0 / denom + dbu[i] * dbu[i] / cJc;
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

void ApproxBranchLengthLikelihood::buildBipartitions(Tree* molTree){
    std::vector<std::set<std::string>> entries;
    for(Node* n : molTree->getDownPassSequence()){
        if(n == molTree->getRoot()) continue;
        std::set<std::string> desc;
        for(Node* d : molTree->getAllDescendants(n))
            if(d->getIsTip()) desc.insert(d->getName());
        if(n->getIsTip()) desc.insert(n->getName());
        entries.push_back(canonicalize(desc));
    }
    if((int)entries.size() != nb)
        Msg::error("ApproxBranchLengthLikelihood: hessian tree has " + std::to_string(entries.size()) + " branches, expected " + std::to_string(nb));
    bipartitions = entries;
}

void ApproxBranchLengthLikelihood::identifyRootBranch(Tree* rootedMlTree){
    Node* root = rootedMlTree->getRoot();
    Node* child0 = nullptr;
    for(Node* c : root->getNeighbors())
        if(c != root->getAncestor()){ child0 = c; break; }
    if(child0 == nullptr)
        Msg::error("ApproxBranchLengthLikelihood: ML root has no children");
    std::set<std::string> canon = canonicalize(molecularDescendants(child0, rootedMlTree));
    for(int i = 0; i < nb; i++)
        if(bipartitions[i] == canon){ rootBranchIdx = i; return; }
    Msg::error("ApproxBranchLengthLikelihood: cannot identify root branch; hessian topology must match inference tree");
}

Node* ApproxBranchLengthLikelihood::findNodeByBipartition(const std::set<std::string>& bp, Tree* tree){
    if(bp.size() == 1)
        return tree->getTaxonNode(*bp.begin());
    for(Node* n : tree->getDownPassSequence()){
        if(n == tree->getRoot() || n->getIsTip()) continue;
        if(canonicalize(molecularDescendants(n, tree)) == bp)
            return n;
    }
    return nullptr;
}

double ApproxBranchLengthLikelihood::computeLnL(Tree* tree, double rate, const std::vector<double>* branchRates){
    if(branchRates == nullptr && rate <= 0.0)
        return -INFINITY;

    Node* curRoot = tree->getRoot();
    if(curRoot != cachedRoot){
        branchNodeIdx.assign(nb, -1);
        for(int i = 0; i < nb; i++){
            if(i == rootBranchIdx) continue;
            Node* n = findNodeByBipartition(bipartitions[i], tree);
            if(n == nullptr){ cachedRoot = nullptr; return -INFINITY; }
            branchNodeIdx[i] = n->getOffset();
        }
        cachedRoot = curRoot;
    }

    Node* root = tree->getRoot();
    std::vector<Node*> rootChildren;
    for(Node* c : root->getNeighbors())
        if(c != root->getAncestor())
            rootChildren.push_back(c);

    std::vector<double> z(nb);
    double lnL = 0.0;
    for(int p = 0; p < nPartitions; p++){
        for(int i = 0; i < nb; i++){
            double predBl;
            if(i == rootBranchIdx){
                if(rootChildren.size() != 2) return -INFINITY;
                double d0 = rootChildren[0]->getAncestor()->getTime() - rootChildren[0]->getTime();
                double d1 = rootChildren[1]->getAncestor()->getTime() - rootChildren[1]->getTime();
                if(branchRates)
                    predBl = rate * ((*branchRates)[rootChildren[0]->getOffset()] * d0 + (*branchRates)[rootChildren[1]->getOffset()] * d1);
                else
                    predBl = rate * (d0 + d1);
            }else{
                Node* n = tree->getNodeByOffset(branchNodeIdx[i]);
                double d = n->getAncestor()->getTime() - n->getTime();
                if(branchRates)
                    predBl = rate * (*branchRates)[branchNodeIdx[i]] * d;
                else
                    predBl = rate * d;
            }
            if(predBl <= 0.0) return -INFINITY;
            z[i] = 2.0 * std::asin(std::sqrt(cJc - cJc * std::exp(-predBl / cJc))) - blMle[p][i];
        }
        for(int i = 0; i < nb; i++)
            lnL += z[i] * gradient[p][i];
        for(int i = 0; i < nb; i++)
            for(int j = 0; j < nb; j++)
                lnL += 0.5 * z[i] * hessian[p][i*nb+j] * z[j];
    }
    return lnL;
}
