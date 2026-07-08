#include "Eigen/Dense"
#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <limits>
#include <set>
#include <vector>

Tree::Tree(const Tree& t) {

    clone(t);
}


Tree::Tree(std::string newick){
    numTaxa = 0;
    std::vector<std::string> newickTokens = parseNewickString(newick);
    if(newickTokens.empty()){
        Node* r = addNode();
        r->setIsTip(false);
        r->setTime(0.0);
        crown = r;
        initializeDownPassSequence();
        return;
    }
    Node* p = nullptr;
    bool readingBl = false;
    for(int i = 0; i < newickTokens.size(); i++){
        std::string token = newickTokens[i];
        if(token == "("){
            Node* newNode = addNode();
            if(p == nullptr){
                crown = newNode;
            }else{
                p->addNeighbor(newNode);
                newNode->addNeighbor(p);
                newNode->setAncestor(p);
            }
            p = newNode;
        }else if (token == ")" || token == ","){
            if(p->getAncestor() == nullptr)
                Msg::error("no anc found for p");
            p = p->getAncestor();
        }else if (token == ";"){
            if(p != crown)
                Msg::error("expecting to be at crown");
        }else if (token == ":"){
            readingBl = true;
        }else{
            if(readingBl == false){
                Node* newNode = addNode();
                if(p == nullptr){
                    crown = newNode;
                }else{
                    newNode->addNeighbor(p);
                    newNode->setAncestor(p);
                    p->addNeighbor(newNode);
                }
                newNode->setName(token);
                newNode->setIsTip(true);
                numTaxa++;
                p = newNode;
            }else{
                readingBl = false;
            }
        }
    }
    initializeDownPassSequence();
    int idx = numTaxa;
    int tIdx = 0;
    for (Node* p : downPassSequence)
        {
            if(p->getIsTip() == true)
                p->setIndex(tIdx++);
            if (p->getIsTip() == false)
                p->setIndex(idx++);
        }
    initializeTimes();
}

Tree::~Tree(void) {

    deleteNodes();
}

Tree& Tree::operator=(const Tree& t) {

    if (this != &t)
        clone(t);
    return *this;
}

Node* Tree::addNode(void) {

    Node* newNode = new Node;
    newNode->setOffset((int)nodes.size());
    nodes.push_back(newNode);
    return newNode;
}

Node* Tree::addOriginNode(double x0) {

    Node* originNode = addNode();
    originNode->setIsTip(false);
    originNode->setIsFossil(false);
    originNode->addNeighbor(crown);
    crown->addNeighbor(originNode);
    origin = originNode;
    initializeDownPassSequence();
    originNode->setTime(x0);
    reindexNodes();
    return originNode;
}

Node* Tree::insertFossilTip(Node* hostChild, std::string name, double y, double z){
    Node* hostParent = hostChild->getAncestor();

    Node* split = addNode();
    split->setIsTip(false);
    split->setIsFossil(false);
    split->setTime(z);

    Node* fossil = addNode();
    fossil->setIsTip(true);
    fossil->setIsFossil(true);
    fossil->setName(name);
    fossil->setTime(y);

    hostParent->removeNeighbor(hostChild);
    hostChild->removeNeighbor(hostParent);

    hostParent->addNeighbor(split);
    split->addNeighbor(hostParent);
    split->setAncestor(hostParent);

    split->addNeighbor(hostChild);
    hostChild->addNeighbor(split);
    hostChild->setAncestor(split);

    split->addNeighbor(fossil);
    fossil->addNeighbor(split);
    fossil->setAncestor(split);

    initializeDownPassSequence();
    reindexNodes();
    return fossil;
}
 

void Tree::clone(const Tree& t) {
    
    if (this->nodes.size() != t.nodes.size())
        {
        deleteNodes();
        for (int i=0; i<t.nodes.size(); i++)
            addNode();
        }
        
    this->numTaxa = t.numTaxa;
    this->crown = this->nodes[t.crown->getOffset()];
    this->origin = (t.origin != nullptr) ? this->nodes[t.origin->getOffset()] : nullptr;
    
    for (int i=0; i<t.nodes.size(); i++)
        {
        Node* q = t.nodes[i];
        Node* p = this->nodes[i];
        p->setIndex(q->getIndex());
        p->setIsTip(q->getIsTip());
        p->setName(q->getName());
        p->setFlag(q->getFlag());
        p->setTime(q->getTime());
        p->setIsFossil(q->getIsFossil());
        if (q->getAncestor() != nullptr)
            p->setAncestor( this->nodes[q->getAncestor()->getOffset()] );
        else
            p->setAncestor(nullptr);
        std::set<Node*>& qNeighbors = q->getNeighbors();
        p->removeAllNeighbors();
        for (Node* qn : qNeighbors)
            p->addNeighbor( this->nodes[qn->getOffset()] );
        }
        
    initializeDownPassSequence();
}


void Tree::deleteNodes(void) {

    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
    nodes.clear();
}
    
std::vector<Node*> Tree::getAllDescendants(Node* n){
    std::vector<Node*> descendants;
    for (Node* child : n->getDescendants()) {
        descendants.push_back(child);
        if (!child->getIsTip()) {
            std::vector<Node*> childDescendants = getAllDescendants(child);
            descendants.insert(descendants.end(), childDescendants.begin(), childDescendants.end());
        }
    }
    return descendants;
}


double Tree::getBranchLength(Node* e1, Node* e2) {

    return std::fabs(e1->getTime() - e2->getTime());
}


// slow O(N) reference gamma-count, kept only for ref
int Tree::getNumLineagesAtTime(double t){
    std::set<Node*> branchesContainingFossils;
    branchesContainingFossils.clear();
    for(Node* n :downPassSequence){
        if(n != getRoot() && n->getTime() < t && n->getAncestor()->getTime() > t){
            branchesContainingFossils.insert(n);
        }
    }
    return((int)branchesContainingFossils.size());
}

std::string Tree::getNewickString(void) {

    std::stringstream strm;
    strm << std::setprecision(17);
    writeTree(getRoot(), strm);
    strm << ";";
    return strm.str();
}

static bool markBackbone(Node* n, std::set<Node*>& keep, std::map<Node*,std::string>& minName){
    if(n->getIsTip()){
        if(n->getIsFossil() == false){ keep.insert(n); minName[n] = n->getName(); return true; }
        return false;
    }
    bool any = false;
    std::string mn;
    for(Node* c : n->getNeighbors())
        if(c != n->getAncestor() && markBackbone(c, keep, minName)){
            any = true;
            const std::string& cm = minName[c];
            if(mn.empty() || cm < mn) mn = cm;
        }
    if(any){ keep.insert(n); minName[n] = mn; }
    return any;
}

static int keptChildren(Node* n, std::set<Node*>& keep, std::map<Node*,std::string>& minName, std::vector<Node*>& out){
    out.clear();
    for(Node* c : n->getNeighbors())
        if(c != n->getAncestor() && keep.count(c) > 0)
            out.push_back(c);
    std::sort(out.begin(), out.end(), [&](Node* a, Node* b){ return minName[a] < minName[b]; });
    return (int)out.size();
}

static Node* descendToBackbone(Node* n, std::set<Node*>& keep, std::map<Node*,std::string>& minName){
    std::vector<Node*> kc;
    while(n->getIsTip() == false && keptChildren(n, keep, minName, kc) == 1)
        n = kc[0];
    return n;
}

static void collectBbNodes(Node* p, std::set<Node*>& keep, std::map<Node*,std::string>& minName, std::vector<Node*>& out){
    if(p->getIsTip())
        return;
    out.push_back(p);
    std::vector<Node*> kc;
    keptChildren(p, keep, minName, kc);
    for(Node* c : kc)
        collectBbNodes(descendToBackbone(c, keep, minName), keep, minName, out);
}

static void writeBackbone(Node* p, std::set<Node*>& keep, std::map<Node*,std::string>& minName, std::map<Node*,int>& xidx, std::stringstream& strm){
    if(p->getIsTip()){
        std::string nm = p->getName();
        std::replace(nm.begin(), nm.end(), ' ', '_');
        strm << nm;
        return;
    }
    std::vector<Node*> kc;
    keptChildren(p, keep, minName, kc);
    strm << "(";
    for(size_t i = 0; i < kc.size(); i++){
        Node* t = descendToBackbone(kc[i], keep, minName);
        if(i > 0) strm << ",";
        writeBackbone(t, keep, minName, xidx, strm);
        strm << ":" << std::fabs(p->getTime() - t->getTime());
    }
    strm << ")";
    std::map<Node*,int>::iterator it = xidx.find(p);
    if(it != xidx.end())
        strm << "x" << it->second;
}

std::string Tree::getBackboneNewickString(void) {

    std::set<Node*> keep;
    std::map<Node*,std::string> minName;
    markBackbone(getRoot(), keep, minName);
    if(keep.empty())
        return "";
    Node* bbRoot = descendToBackbone(getRoot(), keep, minName);
    std::vector<Node*> bb;
    collectBbNodes(bbRoot, keep, minName, bb);
    std::map<Node*,int> xidx;
    for(size_t i = 0; i < bb.size(); i++)
        xidx[bb[i]] = (int)(i + 1);
    std::stringstream strm;
    strm << std::setprecision(17);
    writeBackbone(bbRoot, keep, minName, xidx, strm);
    strm << ";";
    return strm.str();
}

std::vector<Node*> Tree::getBackboneAgeNodes(void) {

    std::set<Node*> keep;
    std::map<Node*,std::string> minName;
    markBackbone(getRoot(), keep, minName);
    std::vector<Node*> out;
    if(keep.empty())
        return out;
    collectBbNodes(descendToBackbone(getRoot(), keep, minName), keep, minName, out);
    return out;
}

static double ageOf(Node* n, const std::map<Node*,double>& age){
    std::map<Node*,double>::const_iterator it = age.find(n);
    return it != age.end() ? it->second : n->getTime();
}

static std::pair<double,double> hpd95(std::vector<double> s){
    std::sort(s.begin(), s.end());
    int n = (int)s.size();
    int m = (int)std::ceil(0.95 * n);
    if(m < 1) m = 1;
    if(m >= n) return std::make_pair(s.front(), s.back());
    double best = std::numeric_limits<double>::max();
    int bi = 0;
    for(int i = 0; i + m - 1 < n; i++){
        double w = s[i + m - 1] - s[i];
        if(w < best){ best = w; bi = i; }
    }
    return std::make_pair(s[bi], s[bi + m - 1]);
}

static void writeBackboneSummary(Node* p, std::set<Node*>& keep, std::map<Node*,std::string>& minName, std::map<Node*,int>& xidx, const std::map<Node*,double>& age, const std::map<Node*,std::pair<double,double>>& hpd, std::stringstream& strm){
    if(p->getIsTip()){
        std::string nm = p->getName();
        std::replace(nm.begin(), nm.end(), ' ', '_');
        strm << nm;
        return;
    }
    std::vector<Node*> kc;
    keptChildren(p, keep, minName, kc);
    strm << "(";
    for(size_t i = 0; i < kc.size(); i++){
        Node* t = descendToBackbone(kc[i], keep, minName);
        if(i > 0) strm << ",";
        writeBackboneSummary(t, keep, minName, xidx, age, hpd, strm);
        strm << ":" << std::fabs(ageOf(p, age) - ageOf(t, age));
    }
    strm << ")";
    std::map<Node*,int>::iterator it = xidx.find(p);
    if(it != xidx.end())
        strm << "x" << it->second;
    std::map<Node*,std::pair<double,double>>::const_iterator h = hpd.find(p);
    if(h != hpd.end()){
        char buf[96];
        std::snprintf(buf, sizeof(buf), "[&95%%HPD={%.6g, %.6g}]", h->second.first, h->second.second);
        strm << buf;
    }
}

std::string Tree::getSummaryNewickString(const std::map<Node*,double>& age, const std::map<Node*,std::pair<double,double>>& hpd){
    std::set<Node*> keep;
    std::map<Node*,std::string> minName;
    markBackbone(getRoot(), keep, minName);
    if(keep.empty())
        return "";
    Node* bbRoot = descendToBackbone(getRoot(), keep, minName);
    std::vector<Node*> bb;
    collectBbNodes(bbRoot, keep, minName, bb);
    std::map<Node*,int> xidx;
    for(size_t i = 0; i < bb.size(); i++)
        xidx[bb[i]] = (int)(i + 1);
    std::stringstream strm;
    strm << std::setprecision(17);
    writeBackboneSummary(bbRoot, keep, minName, xidx, age, hpd, strm);
    strm << ";";
    return strm.str();
}

bool writeSummaryTree(Tree* tree, const std::vector<std::string>& names, const std::vector<std::vector<double>>& cols, double burninFrac, const std::string& path){
    std::vector<Node*> bb = tree->getBackboneAgeNodes();
    if(bb.empty())
        return false;
    std::map<std::string,int> colOf;
    for(size_t j = 0; j < names.size(); j++)
        colOf[names[j]] = (int)j;
    std::map<Node*,double> age;
    std::map<Node*,std::pair<double,double>> hpd;
    for(size_t i = 0; i < bb.size(); i++){
        std::map<std::string,int>::iterator it = colOf.find("x" + std::to_string(i + 1));
        if(it == colOf.end() || (size_t)it->second >= cols.size())
            return false;
        const std::vector<double>& col = cols[it->second];
        size_t b = (size_t)(burninFrac * col.size());
        std::vector<double> s(col.begin() + b, col.end());
        if(s.empty())
            return false;
        double sum = 0.0;
        for(double v : s) sum += v;
        age[bb[i]] = sum / (double)s.size();
        hpd[bb[i]] = hpd95(s);
    }
    std::string nwk = tree->getSummaryNewickString(age, hpd);
    if(nwk.empty())
        return false;
    std::ofstream os(path);
    os << nwk << "\n";
    return true;
}

void Tree::buildBackboneCache(void){
    int N = (int)nodes.size();
    bbParentByOffset.assign(N, nullptr);
    bbChildrenByOffset.assign(N, std::vector<Node*>());
    bbRateNodes.clear();
    bbRootNode = getRoot();
    bbCacheValid = true;

    std::set<Node*> keep;
    std::map<Node*,std::string> minName;
    markBackbone(bbRootNode, keep, minName);
    if(keep.empty())
        return;

    std::vector<Node*> stack;
    stack.push_back(bbRootNode);
    while(stack.empty() == false){
        Node* p = stack.back();
        stack.pop_back();
        std::vector<Node*> kc;
        keptChildren(p, keep, minName, kc);
        std::vector<Node*> bbch;
        for(Node* c : kc){
            Node* child = descendToBackbone(c, keep, minName);
            bbch.push_back(child);
            bbParentByOffset[child->getOffset()] = p;
            stack.push_back(child);
        }
        std::sort(bbch.begin(), bbch.end(), [](Node* a, Node* b){ return a->getOffset() < b->getOffset(); });
        bbChildrenByOffset[p->getOffset()] = bbch;
    }
    for(Node* n : downPassSequence)
        if(n != bbRootNode && bbParentByOffset[n->getOffset()] != nullptr)
            bbRateNodes.push_back(n);
}

void Tree::ensureBackboneCache(void){
    if(bbCacheValid == false)
        buildBackboneCache();
}

const std::vector<Node*>& Tree::getBackboneRateNodes(void){
    if(bbCacheValid == false)
        buildBackboneCache();
    return bbRateNodes;
}

Node* Tree::getBackboneRoot(void){
    if(bbCacheValid == false)
        buildBackboneCache();
    return bbRootNode;
}

Node* Tree::getBackboneParent(Node* n){
    if(bbCacheValid == false)
        buildBackboneCache();
    return bbParentByOffset[n->getOffset()];
}

const std::vector<Node*>& Tree::getBackboneChildren(Node* n){
    if(bbCacheValid == false)
        buildBackboneCache();
    return bbChildrenByOffset[n->getOffset()];
}

bool Tree::isBackboneNode(Node* n){
    if(bbCacheValid == false)
        buildBackboneCache();
    return n == bbRootNode || bbParentByOffset[n->getOffset()] != nullptr;
}

int Tree::getNumBackbone(void){
    int idx = 0;
    for(Node* n : downPassSequence)
        if(n->getIsTip() == true)
            idx++;
    return idx;
}

bool Tree::isBinary(void){
    for(Node* n : downPassSequence){
        if(n->getIsTip())
            continue;
        int numChildren = (int)n->getNeighbors().size();
        if(n != getRoot())
            numChildren -= 1;
        if(numChildren != 2)
            return false;
    }
    return true;
}

void Tree::validateBackbone(void){
    if(isBinary() == false)
        Msg::error("input tree must be fully resolved and rooted");
}

bool Tree::isValidNewick(const std::string& s){
    size_t begin = s.find_first_not_of(" \t\n\r");
    size_t end = s.find_last_not_of(" \t\n\r");
    if(begin == std::string::npos)
        return false;
    if(s[begin] != '(' || s[end] != ';')
        return false;
    int depth = 0;
    for(char c : s){
        if(c == '(')
            depth++;
        else if(c == ')'){
            depth--;
            if(depth < 0)
                return false;
        }
    }
    return depth == 0;
}

Node* Tree::getTaxonNode(std::string name) {

    for (Node* p : downPassSequence)
        {
        if (p->getName() == name)
            return p;
        }
    return nullptr;
}

Node* Tree::getMRCA(const std::vector<std::string>& taxonNames){
    Node* mrca = nullptr;
    for(const std::string& name : taxonNames){
        Node* n = getTaxonNode(name);
        if(n == nullptr)
            Msg::error("taxon not found in tree: " + name);

        if(mrca == nullptr){
            mrca = n;
            continue;
        }

        std::set<Node*> ancestors;
        for(Node* p = mrca; ; p = p->getAncestor()){
            ancestors.insert(p);
            if(p->getAncestor() == p) // root's ancestor is itself in this tree
                break;
        }

        Node* common = nullptr;
        for(Node* p = n; common == nullptr; p = p->getAncestor()){
            if(ancestors.count(p) > 0)
                common = p;
            else if(p->getAncestor() == p)
                break;
        }
        mrca = common;
    }
    return mrca;
}

void Tree::initializeDownPassSequence(void) {
    downPassSequence.clear();
    passDown(getRoot(), getRoot());
    bbCacheValid = false;
}

void Tree::initializeTimes(void){
    initializeDownPassSequence();
    std::map<Node*,double> noMinAges;
    assignStartingAges(noMinAges, 1.0);
}

void Tree::assignStartingAges(const std::map<Node*,double>& minAges, double unit){
    for(Node* n : downPassSequence){
        if(n->getIsTip()){
            n->setTime(0.0);
            n->setIsFossil(false);
            continue;
        }
        double maxChild = 0.0;
        for(Node* c : n->getNeighbors())
            if(c != n->getAncestor() && c->getTime() > maxChild)
                maxChild = c->getTime();
        double age = maxChild + unit;
        std::map<Node*,double>::const_iterator it = minAges.find(n);
        if(it != minAges.end() && it->second > age)
            age = it->second;
        n->setTime(age);
    }
}


std::vector<std::string>  Tree::parseNewickString(std::string newickStr){
    
    std::vector<std::string> tokens;
    
    std::string token = "";
    for(int i = 0; i < newickStr.length(); i++){
        char c = newickStr[i];
        if(c == ' ' || c == '\t' || c == '\n' || c == '\r')
            continue;
        if(c == '(' || c == ')' || c==',' || c==':' ||c == ';'){
            if(token != ""){
                tokens.push_back(token);
                token = "";
            }
            tokens.push_back(std::string(1,c));
        }else{
            token += std::string(1,c);
        }
    }
    
    if(token != "")
        tokens.push_back(token);

    return tokens;
}

void Tree::passDown(Node* p, Node* from) {

    if (p != nullptr)
        {
        std::vector<Node*> children;
        for (Node* d : p->getNeighbors())
            if (d != from)
                children.push_back(d);
        std::sort(children.begin(), children.end(), [](Node* a, Node* b){ return a->getOffset() < b->getOffset(); });
        for (Node* d : children)
            passDown(d, p);
        p->setAncestor(from);
        downPassSequence.push_back(p);
        }
}

void Tree::print(void) {

    showNode(getRoot(), 0);
    std::cout << "Postorder sequence: ";
    for (int i=0; i<downPassSequence.size(); i++)
        std::cout << downPassSequence[i]->getIndex() << " ";
    std::cout << std::endl;
}

void Tree::print(std::string header) {

    std::cout << header << std::endl;
    print();
}


void Tree::reindexNodes(void){
    int idx = 0;
    for (Node* p : downPassSequence)
        {
            if(p->getIsTip() == true)
                p->setIndex(idx++);
        }
    for (Node* p : downPassSequence)
        {
            if(p->getIsTip() == false)
                p->setIndex(idx++);
        }
}


void Tree::showNode(Node* p, int indent) {

    if (p != nullptr)
        {
        std::set<Node*>& pNeighbors = p->getNeighbors();
        for (int i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex() << " [" << p << "] ( ";
        for (Node* n : pNeighbors)
            {
            if (n == p->getAncestor())
                std::cout << "a_";
            std::cout << n->getIndex() << " ";
            }
        std::cout << ") ";
        std::cout << p->getName() << " " << p->getIsTip() << " ";
        std::cout << std::fixed << std::setprecision(6);
        if (p != getRoot())
            std::cout << getBranchLength(p, p->getAncestor()) << " ";
        else
            std::cout << "--- ";
        if (p == getRoot())
            std::cout << "<- Root ";
        std::cout << std::endl;

        for (Node* n : pNeighbors)
            {
            if (n != p->getAncestor())
                showNode(n, indent + 3);
            }
            
        }
}

double Tree::update(double scaleLambda){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    int numSlideable = 0;
    for(Node* n : downPassSequence)
        if(n != getRoot() && n->getIsTip() == false)
            numSlideable++;
    double crownWeight = 3.0;
    if(rng.uniformRv() * (numSlideable + crownWeight) >= numSlideable)
        return updateCrownAge(scaleLambda);
    lastUpdateWasScale = false;
    return updateNodeAge();
}

bool Tree::isSATip(Node* n){
    return n->getIsTip() && n->getIsFossil() && n->getAncestor() != n && n->getAncestor()->getTime() == n->getTime();
}

bool Tree::isFakeSplit(Node* n){
    if(n->getIsTip())
        return false;
    for(Node* c : n->getNeighbors())
        if(c != n->getAncestor() && isSATip(c))
            return true;
    return false;
}

int Tree::scaleInternalAges(double m){
    int numScaled = 0;
    for(Node* n : downPassSequence)
        if(n->getIsTip() == false && isFakeSplit(n) == false){
            n->setTime(n->getTime() * m);
            numScaled++;
        }
    lastUpdateWasScale = true;
    return numScaled;
}

int Tree::scaleSubtreeAges(Node* subtreeCrown, double m){
    int numScaled = 0;
    subtreeCrown->setTime(subtreeCrown->getTime() * m);
    numScaled++;
    for(Node* d : getAllDescendants(subtreeCrown))
        if(d->getIsTip() == false){
            d->setTime(d->getTime() * m);
            numScaled++;
        }
    return numScaled;
}

std::vector<Node*> Tree::getInternalAgeNodes(void){
    std::vector<Node*> candidates;
    for(Node* n : downPassSequence)
        if(n != getRoot() && n->getIsTip() == false && isFakeSplit(n) == false)
            candidates.push_back(n);
    return candidates;
}

double Tree::updateNodeAgeOnNode(Node* n){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double parentAge = n->getAncestor()->getTime();
    double maxChildAge = 0.0;
    for(Node* nb : n->getNeighbors())
        if(nb != n->getAncestor() && nb->getTime() > maxChildAge)
            maxChildAge = nb->getTime();

    std::map<Node*,double>::iterator it = ageFloors.find(n);
    if(it != ageFloors.end() && it->second > maxChildAge)
        maxChildAge = it->second;

    n->setTime(maxChildAge + rng.uniformRv() * (parentAge - maxChildAge));

    return 0.0;
}

double Tree::updateNodeAge(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    std::vector<Node*> candidates = getInternalAgeNodes();
    if(candidates.size() == 0)
        return 0.0;
    Node* n = candidates[(int)(rng.uniformRv() * candidates.size())];
    return updateNodeAgeOnNode(n);
}

double Tree::updateCrownAge(double scaleLambda){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double maxChildAge = 0.0;
    for(Node* c : getRoot()->getNeighbors())
        if(c->getTime() > maxChildAge)
            maxChildAge = c->getTime();

    std::map<Node*,double>::iterator it = ageFloors.find(getRoot());
    if(it != ageFloors.end() && it->second > maxChildAge)
        maxChildAge = it->second;

    double m = std::exp( scaleLambda * (rng.uniformRv() - 0.5) );
    double newCrownAge = getRoot()->getTime() * m;
    lastUpdateWasScale = true;
    if(newCrownAge <= maxChildAge)
        return -INFINITY;

    getRoot()->setTime(newCrownAge);

    return std::log(m);
}

void Tree::writeTree(Node* p, std::stringstream& strm) {

    if (p != nullptr)
        {
        std::set<Node*>& pNeighbors = p->getNeighbors();
        if (p->getIsTip() == false)
            strm << "(";
        else
            {
            std::string tName = p->getName();
            std::replace( tName.begin(), tName.end(), ' ', '_');
            strm << tName;
            }
        bool firstNode = true;
        for (Node* n : pNeighbors)
            {
            if (n != p->getAncestor())
                {
                if (firstNode == false)
                    strm << ",";
                writeTree(n, strm);
                strm << ":" << getBranchLength(p, n);
                firstNode = false;
                }
            }
        if (p->getIsTip() == false)
            strm << ")";

        }
}

void Tree::writeState(std::ostream& os) {

    os << nodes.size() << ' ' << crown->getOffset() << ' ' << (origin != nullptr ? origin->getOffset() : -1) << ' ' << numTaxa << '\n';
    for (Node* n : nodes)
        {
        os << (n->getIsTip() ? 1 : 0) << ' ' << (n->getIsFossil() ? 1 : 0) << ' ' << (n->getFlag() ? 1 : 0)
           << ' ' << n->getIndex() << ' ' << n->getTime();
        std::set<Node*>& nb = n->getNeighbors();
        os << ' ' << nb.size();
        std::vector<int> nbOff;
        for (Node* m : nb)
            nbOff.push_back(m->getOffset());
        std::sort(nbOff.begin(), nbOff.end());
        for (int off : nbOff)
            os << ' ' << off;
        const std::string& nm = n->getName();
        os << ' ' << (nm.empty() ? "*" : nm) << '\n';
        }
}

void Tree::readState(std::istream& is) {

    size_t nn;
    int crownOff, originOff;
    is >> nn >> crownOff >> originOff >> numTaxa;
    if (nodes.size() != nn)
        {
        deleteNodes();
        for (size_t i = 0; i < nn; i++)
            addNode();
        }
    std::vector<std::vector<int>> nbrs(nn);
    for (size_t i = 0; i < nn; i++)
        {
        int tip, foss, fl, idx;
        double tm;
        size_t nc;
        is >> tip >> foss >> fl >> idx >> tm >> nc;
        Node* p = nodes[i];
        p->setIsTip(tip != 0);
        p->setIsFossil(foss != 0);
        p->setFlag(fl != 0);
        p->setIndex(idx);
        p->setTime(tm);
        nbrs[i].resize(nc);
        for (size_t k = 0; k < nc; k++)
            is >> nbrs[i][k];
        std::string nm;
        is >> nm;
        p->setName(nm == "*" ? "" : nm);
        }
    for (size_t i = 0; i < nn; i++)
        {
        Node* p = nodes[i];
        p->removeAllNeighbors();
        for (int o : nbrs[i])
            p->addNeighbor(nodes[o]);
        }
    crown = nodes[crownOff];
    origin = (originOff >= 0) ? nodes[originOff] : nullptr;
    initializeDownPassSequence();
}
