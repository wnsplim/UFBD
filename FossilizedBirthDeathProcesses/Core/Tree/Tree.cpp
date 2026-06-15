#include "Eigen/Dense"
#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

Tree::Tree(const Tree& t) {

    clone(t);
}


Tree::Tree(std::string newick){
    numTaxa = 0;
    std::vector<std::string> newickTokens = parseNewickString(newick);
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
                tips.push_back(newNode);
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


int Tree::getNumLineagesAtTime(double t){
    std::set<Node*> branchesContainingFossils;
    branchesContainingFossils.clear();
    for(Node* n :downPassSequence){
        if(n != crown && n->getTime() < t && n->getAncestor()->getTime() > t){
            branchesContainingFossils.insert(n);
        }
    }
    return((int)branchesContainingFossils.size());
}

std::string Tree::getNewickString(void) {

    std::stringstream strm;
    strm << std::setprecision(17);
    writeTree(crown, strm);
    strm << ";";
    return strm.str();
}

int Tree::getNumTaxa(void){
    initializeDownPassSequence();
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
        if(n != crown)
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
            if(p->getAncestor() == p) // crown's ancestor is itself in this tree
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
    if(crown == nullptr)
        Msg::error("crown is nullptr");
    downPassSequence.clear();
    passDown(crown, crown);
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
    treeHeight = crown->getTime();
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
    
//    for(int i = 0; i < tokens.size(); i++){
//        std::cout << i << " " << tokens[i] << std::endl;
//    }
//        
    
    return tokens;
}

void Tree::passDown(Node* p, Node* from) {

    if (p != nullptr)
        {
        std::set<Node*>& pNeighbors = p->getNeighbors();
        for (Node* d : pNeighbors)
            {
            if (d != from)
                passDown(d, p);
            }
        p->setAncestor(from);
        downPassSequence.push_back(p);
        }
}

void Tree::print(void) {
    
    showNode(crown, 0);
    std::cout << "Postorder sequence: ";
    for (int i=0; i<downPassSequence.size(); i++)
        std::cout << downPassSequence[i]->getIndex() << " ";
    std::cout << std::endl;
    //printVCV();
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
        if (p != crown)
            std::cout << getBranchLength(p, p->getAncestor()) << " ";
        else
            std::cout << "--- ";
        if (p == crown)
            std::cout << "<- Crown ";
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
        if(n != crown && n->getIsTip() == false)
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

double Tree::updateNodeAge(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    std::vector<Node*> candidates;
    for(Node* n : downPassSequence)
        if(n != crown && n->getIsTip() == false && isFakeSplit(n) == false)
            candidates.push_back(n);
    if(candidates.size() == 0)
        return 0.0;

    Node* n = candidates[(int)(rng.uniformRv() * candidates.size())];
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

double Tree::updateCrownAge(double scaleLambda){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double maxChildAge = 0.0;
    for(Node* c : crown->getNeighbors())
        if(c->getTime() > maxChildAge)
            maxChildAge = c->getTime();

    std::map<Node*,double>::iterator it = ageFloors.find(crown);
    if(it != ageFloors.end() && it->second > maxChildAge)
        maxChildAge = it->second;

    double m = std::exp( scaleLambda * (rng.uniformRv() - 0.5) );
    double newCrownAge = crown->getTime() * m;
    lastUpdateWasScale = true;
    if(newCrownAge <= maxChildAge)
        return -INFINITY;

    crown->setTime(newCrownAge);

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
