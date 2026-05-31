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

Tree::Tree(std::vector<std::string> taxonNames, double lambda) {

    numTaxa = (int)taxonNames.size();

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    // start with three-species tree
    root = addNode();
    for (int i=0; i<3; i++)
        {
        Node* p = addNode();
        p->setIsTip(true);
        p->setIndex(i);
        p->setName(taxonNames[i]);
        p->addNeighbor(root);
        p->setAncestor(root);
        root->addNeighbor(p);
        }
        
    // randomly add the remaining taxa to the branches of the tree
    for (int i=3; i<numTaxa; i++)
        {
        Node* p = nullptr;
        do {
            double u = rng.uniformRv();
            int whichNode = (int)(u*nodes.size());
            p = nodes[whichNode];
            } while (p == root);
        Node* pAnc = p->getAncestor();
        p->removeNeighbor(pAnc);
        pAnc->removeNeighbor(p);

        Node* newTip = addNode();
        newTip->setIndex(i);
        newTip->setIsTip(true);
        newTip->setName(taxonNames[i]);
        Node* newInt = addNode();

        p->addNeighbor(newInt);
        p->setAncestor(newInt);
        pAnc->addNeighbor(newInt);
        newInt->addNeighbor(pAnc);
        newInt->addNeighbor(p);
        newInt->addNeighbor(newTip);
        newInt->setAncestor(pAnc);
        newTip->addNeighbor(newInt);
        newTip->setAncestor(newInt);
        }
        
    // initialize the post-order traversal sequence
    initializeDownPassSequence();
    
    // index interior nodes
    int idx = numTaxa;
    for (Node* p : downPassSequence)
        {
        if (p->getIsTip() == false)
            p->setIndex(idx++);
        }
    
    // add the branch lengths to the tree
    for (Node* p : downPassSequence)
        {
        if (p != root)
            setBranch(p, p->getAncestor(), Probability::Exponential::rv(&rng, lambda));
        }
    initializeTimes();
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
                root = newNode;
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
            if(p != root)
                Msg::error("expecting to be at root");
        }else if (token == ":"){
            readingBl = true;
        }else{
            if(readingBl == false){
                Node* newNode = addNode();
                newNode->addNeighbor(p);
                newNode->setAncestor(p);
                p->addNeighbor(newNode);
                newNode->setName(token);
                newNode->setIsTip(true);
                tips.push_back(newNode);
                numTaxa++;
                p = newNode;
            }else{
                double x = stod(token);
                if(p->getAncestor() != nullptr)
                    setBranch(p, p->getAncestor(), x);
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
 
void Tree::calculateTreeHeight(void){
    treeHeight = 0.0;
    double maxHeight = 0.0;
    for(Node* n : downPassSequence){
        if(n->getIsTip() == true){
            Node* p = n;
            Node* pAnc = p->getAncestor();
            double height = 0.0;
            while(p != root){
                pAnc = p->getAncestor();
                height += branchLengthFromMap(p, pAnc);
                p = pAnc;
            };
            if(height > maxHeight)
                maxHeight = height;
        }
    }
    treeHeight = maxHeight;
}

void Tree::checkBranchLengthsNeg(void){
    for(auto const& x : branchLengths)
        if(x.second < 0)
            Msg::error("negative branch lengths");
}

void Tree::clone(const Tree& t) {
    
    if (this->nodes.size() != t.nodes.size())
        {
        deleteNodes();
        for (int i=0; i<t.nodes.size(); i++)
            addNode();
        }
        
    this->numTaxa = t.numTaxa;
    this->root = this->nodes[t.root->getOffset()];
    
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
    this->branchLengths.clear();
    for (BranchLengths::const_iterator it = t.branchLengths.begin(); it != t.branchLengths.end(); it++)
        {
        if(it->first.first == nullptr)
            Msg::error("stop");
        Node* e1 = this->nodes[it->first.first->getOffset()];
        Node* e2 = this->nodes[it->first.second->getOffset()];
        setBranch(e1, e2, it->second);
        }
}

void Tree::collapseNode(Node* n){
    if(n->getNeighbors().size() != 2)
        Msg::error("cannot collapse this node because it does not have two neighbors");
        
//    std::cout << "deleting " << n->getIndex() << std::endl;
    
    std::vector<Node*> neighbors(n->getNeighbors().begin(), n->getNeighbors().end());
    double bl = getBranchLength(neighbors[0], n) + getBranchLength(neighbors[1], n);
    if(n == root){
        //need to select new root
        Node* foundNewRoot = nullptr;
        for(int i = 0; i < neighbors.size(); i++){
            if(neighbors[i]->getNeighbors().size() == 3){
                foundNewRoot = neighbors[i];
                neighbors.erase(neighbors.begin() + i);
            } //beacuse of how we set up the rot class, there should only be one node wiht this condition
        }
        
        if(neighbors.size() != 1)
            Msg::error("wonky neighbor sizing");
        if(foundNewRoot == nullptr)
            Msg::error("could not find suitable new root");
        
//        std::cout << "New root = " << foundNewRoot->getIndex() << std::endl;

        root = foundNewRoot;
        root->setAncestor(nullptr);
        Node* b = neighbors[0];
        b->addNeighbor(root);
        setBranch(root, b, bl);
        
        b->removeNeighbor(n);
        root->removeNeighbor(n);
        removeBranch(n, b);
        removeBranch(n, root);
        
    }else{
        neighbors[0]->addNeighbor(neighbors[1]);
        setBranch(neighbors[0], neighbors[1], bl);
        
        neighbors[0]->removeNeighbor(n);
        neighbors[1]->removeNeighbor(n);
        removeBranch(neighbors[0], n);
        removeBranch(neighbors[1], n);
    }
    
    n->removeAllNeighbors();
    n->setAncestor(nullptr);
    freeNode = n;
    freeNode->setIndex(-1);
//    for(int i = 0; i < nodes.size(); i++){
//        if(nodes[i] == n){
//            nodes.erase(nodes.begin() + i);
//            delete n;
//        }
//    }
    
    initializeDownPassSequence(); // resolves issues with DownPassSeq
}

void Tree::dropTip(std::string tip){
    Node* p = nullptr;
    for(Node* n : downPassSequence)
        if(n->getName() == tip)
            if(n->getIsTip() == true)
                p = n;
    if(p == nullptr)
        Msg::error("could not find tip in tree");
        
    // starting: (((p, q), a),b) p is tip to drop, a is int node to delete
    // ending:   ((q, some other tip that is already b descendent), b)
    Node* a = p->getAncestor();
    Node* b = a->getAncestor();
    Node* q = nullptr;
    std::vector<Node*> aDesc = a->getDescendants();
    if(aDesc[0] == p){
        q = aDesc[1];
    }else if (aDesc[1] == p){
        q = aDesc[0];
    }else{
        Msg::error("could not find tip in descendants");
    }
       
    if(a == root){
        //std::cout << "pongo pyg case" << std::endl;
        //q becomes new root
        removeBranch(q, a);
        q->removeNeighbor(a);
        root = q;
        initializeDownPassSequence();
    }else{
        double newBl = getBranchLength(q, a) + getBranchLength(a, b);
        //remove branches from branch map
        removeBranch(p, a);
        removeBranch(a, b);
        p->removeAllNeighbors();
        a->removeAllNeighbors();
        q->removeNeighbor(a);
        b->removeNeighbor(a);
        q->addNeighbor(b);
        q->setAncestor(b);
        b->addNeighbor(q);
        setBranch(q, b, newBl);
        initializeDownPassSequence();
    }
    
    
}

void Tree::dropTip(int idx){
    Node* p = nullptr;
    for(Node* n : downPassSequence)
        if(n->getIndex() == idx)
            if(n->getIsTip() == true)
                p = n;
    if(p == nullptr)
        Msg::error("could not find tip in tree");
        
    // starting: (((p, q), a),b) p is tip to drop, a is int node to delete
    // ending:   ((q, some other tip that is already b descendent), b)
    Node* a = p->getAncestor();
    Node* b = a->getAncestor();
    Node* q = nullptr;
    std::vector<Node*> aDesc = a->getDescendants();
    if(aDesc[0] == p){
        q = aDesc[1];
    }else if (aDesc[1] == p){
        q = aDesc[0];
    }else{
        Msg::error("could not find tip in descendants");
    }
       
    if(a == root){
        //std::cout << "pongo pyg case" << std::endl;
        //q becomes new root
        removeBranch(q, a);
        q->removeNeighbor(a);
        root = q;
        initializeDownPassSequence();
    }else{
        double newBl = getBranchLength(q, a) + getBranchLength(a, b);
        //remove branches from branch map
        removeBranch(p, a);
        removeBranch(a, b);
        p->removeAllNeighbors();
        a->removeAllNeighbors();
        q->removeNeighbor(a);
        b->removeNeighbor(a);
        q->addNeighbor(b);
        q->setAncestor(b);
        b->addNeighbor(q);
        setBranch(q, b, newBl);
        initializeDownPassSequence();
    }
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

void Tree::forceBinary(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    std::vector<Node*> rootDesc = root->getDescendants();
    if(rootDesc.size() != 2){
        Node* newOutgroup = rootDesc[(int)(rng.uniformRv() * rootDesc.size())];
        Node* newRoot = addNode();
        newRoot->setIndex((int)nodes.size() + 1);
        double blHalf = getBranchLength(newOutgroup, root) / 2;
        
        removeBranch(newOutgroup, root);
        newOutgroup->removeNeighbor(root);
        
        newOutgroup->addNeighbor(newRoot);
        root->addNeighbor(newRoot);
        setBranch(newRoot, root, blHalf);
        setBranch(newOutgroup, newRoot, blHalf);
        reroot(newRoot);
    }
    
    for(Node* n : downPassSequence)
        if(n->getDescendants().size() != 2 && n->getIsTip() == false){
            print();
            Msg::error("Tree is not binary");
        }
        
}

double Tree::getBranchLength(Node* e1, Node* e2) {

    return std::fabs(e1->getTime() - e2->getTime());
}

double Tree::branchLengthFromMap(Node* e1, Node* e2) {

    std::pair<Node*,Node*> key;
    initializeBranchLengthKey(key, e1, e2);
    BranchLengths::iterator it = branchLengths.find(key);
    if (it == branchLengths.end()){
        print("Tree with issue");
        std::cout << "for nodes: " << e1->getIndex() << " " << e2->getIndex() << std::endl;
        std::cout << "for nodes: " << e1 << " " << e2 << std::endl;
        std::cout << "e1 anc: " << e1->getAncestor()->getIndex() << std::endl;
        Msg::error("Could not find branch to return");
    }
    return it->second;
}

int Tree::getNumLineagesAtTime(double t){
    std::set<Node*> branchesContainingFossils;
    branchesContainingFossils.clear();
    for(Node* n :downPassSequence){
        if(n != root && n->getTime() < t && n->getAncestor()->getTime() > t){
            branchesContainingFossils.insert(n);
        }
    }
    return((int)branchesContainingFossils.size());
}

std::string Tree::getNewickString(void) {

    std::stringstream strm;
    writeTree(root, strm);
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

Node* Tree::getTaxonNode(std::string name) {

    for (Node* p : downPassSequence)
        {
        if (p->getName() == name)
            return p;
        }
    return nullptr;
}

Node* Tree::getMRCA(const std::vector<std::string>& taxonNames){
    if(taxonNames.size() < 2)
        Msg::error("a clade must be defined by at least two taxa");

    Node* mrca = getTaxonNode(taxonNames[0]);
    if(mrca == nullptr)
        Msg::error("taxon not found in tree: " + taxonNames[0]);

    for(int i = 1; i < (int)taxonNames.size(); i++){
        Node* n = getTaxonNode(taxonNames[i]);
        if(n == nullptr)
            Msg::error("taxon not found in tree: " + taxonNames[i]);

        std::set<Node*> ancestors;
        for(Node* p = mrca; p != nullptr; p = p->getAncestor())
            ancestors.insert(p);

        Node* common = nullptr;
        for(Node* p = n; p != nullptr && common == nullptr; p = p->getAncestor())
            if(ancestors.count(p) > 0)
                common = p;
        mrca = common;
    }
    return mrca;
}

void Tree::initializeBranchLengthKey(std::pair<Node*,Node*>& key, Node* e1, Node* e2) {

    if (e1 < e2)
        key = std::make_pair(e1,e2);
    else
        key = std::make_pair(e2,e1);
}

void Tree::initializeDownPassSequence(void) {
    if(root == nullptr)
        Msg::error("root is nullptr");
    downPassSequence.clear();
    passDown(root, root);
    checkBranchLengthsNeg();
}

void Tree::initializeTimes(void){
    calculateTreeHeight();
    
    //instantiate node times
    for (std::vector<Node*>::reverse_iterator rit = downPassSequence.rbegin(); rit != downPassSequence.rend(); rit++){
        Node* n = *rit;
        if(n==root)
            n->setTime(treeHeight);
        else{
            Node* p = n;
            Node* pAnc = p->getAncestor();
            double heightFromRoot = 0.0;
            while(p != root){
                pAnc = p->getAncestor();
                heightFromRoot += branchLengthFromMap(p, pAnc);
                p = pAnc;
            };
            n->setTime(treeHeight - heightFromRoot);
        }
    }
    
    //instantiate fossils
    for(Node* n : downPassSequence){
        if(n->getIsTip() == true){
            Node* p = n;
            Node* pAnc = p->getAncestor();
            double height = 0.0;
            while(p != root){
                pAnc = p->getAncestor();
                height += branchLengthFromMap(p, pAnc);
                p = pAnc;
            };
            if(height < treeHeight)
                n->setIsFossil(true);
            else
                n->setIsFossil(false);
        }
    }
}

void Tree::keepTips(std::vector<std::string> t){
    //function drops all tips but those in the vector tips
    initializeDownPassSequence();
    
    for(Node* n : tips){
        //std::cout << n->getName() << std::endl;
        
        int cnt = (int)count(t.begin(), t.end(), n->getName());

        // Check if the target value was found
        if (cnt == 0)
            dropTip(n->getName());
        
    }
    
}

std::vector<std::string>  Tree::parseNewickString(std::string newickStr){
    
    std::vector<std::string> tokens;
    
    std::string token = "";
    for(int i = 0; i < newickStr.length(); i++){
        char c = newickStr[i];
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
    
    showNode(root, 0);
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

std::pair<Node*,Node*> Tree::randomlyChooseBranch(void) {
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    std::pair<Node*,Node*> nodePair;
    
    Node* p = nullptr;
    do {
        p = nodes[(int)(rng.uniformRv()*nodes.size())];
        } while(p == root);
    Node* pAnc = p->getAncestor();
    nodePair = std::make_pair(p, pAnc);
    
    return nodePair;
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

void Tree::removeAllBranches(void) {

    branchLengths.clear();
}

void Tree::removeBranch(Node* e1, Node* e2) {

    std::pair<Node*,Node*> key;
    initializeBranchLengthKey(key, e1, e2);
    BranchLengths::iterator it = branchLengths.find(key);
    if (it != branchLengths.end())
        branchLengths.erase(it);
    else
        Msg::error("Could not find branch to remove");
}

void Tree::reroot(Node* r) {
    if(r->getIsTip() == true)
        Msg::error("why are you rerooting on a tip?");
    root = r;
    downPassSequence.clear();
    passDown(root, root);
}

void Tree::rSPR(void){
    int startingNodeSize = nodes.size();
    
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    std::vector<Node*> candidateMoves = {};
    for(Node* n : downPassSequence)
        if(n->getIsTip() == false)
            candidateMoves.push_back(n);
    if(candidateMoves.size() == 0){
        print();
        std::cout << getNewickString() << std::endl;
        Msg::error("number of candidate SPR moves = 0");
    }
    
    bool isBinary = false;
    if(root->getNeighbors().size() == 2)
        isBinary = true;
    
    //which interior node to move around
    Node* p = candidateMoves[(int)(rng.uniformRv() * candidateMoves.size())];
//    Node* p = root;
    
    std::set<Node*> pNeighbors = p->getNeighbors();
    if(pNeighbors.size() != 3 && p != root){
        std::cout << p->getIndex() << std::endl;
        print();
        Msg::error("pNeighbor size is not 3");
    }
    
    //sample neighbor at random
    int idx = (int)(rng.uniformRv() * pNeighbors.size());
    Node* q = *next(pNeighbors.begin(), idx);
    
    if(q == nullptr)
        Msg::error("could not sample q");
    
//    std::cout << "Int node: " << p->getIndex() << " node to move: " << q->getIndex() << std::endl;
    if(q == root){
        Node* transfer = p;
        p = q;
        q = transfer;
        pNeighbors = p->getNeighbors();
    }
    
    //p is what we're rotating around and collapsing
    //q is qhat actually moves
    if(p == root){
        if(isBinary == false){
            //root has three descendants
            //both neighbors of q tips with q int node?
            bool moveRootSubtree = false;
            pNeighbors.erase(q);
            if((*next(pNeighbors.begin(), 0))->getIsTip() == true && (*next(pNeighbors.begin(), 1))->getIsTip() == true)
                moveRootSubtree = true;
            if(moveRootSubtree == false){
    //            std::cout << "Move root subtree false " << std::endl;
                double bl = getBranchLength(p, q);
                q->removeNeighbor(p);
                collapseNode(p);
                initializeDownPassSequence();
                //finding a reattachment point (cannot be root)
                std::vector<Node*> candidateRegraftingLocations = {};
                for(Node* n : downPassSequence)
                    if(n != root && n != freeNode)
                        candidateRegraftingLocations.push_back(n);
                if(candidateRegraftingLocations.size() == 0){
                    print();
                    std::cout << getNewickString() << std::endl;
                    Msg::error("number of candidate SPR reattachment moves = 0");
                }
                
    //            Node* regraft = addNode();
                Node* regraft = freeNode;
                regraft->setIndex(nodes.size()+1);
                q->addNeighbor(regraft);
                q->setAncestor(regraft);
                setBranch(q, regraft, bl);

                
                //reattach q on branch subtending a
                Node* a = candidateRegraftingLocations[(int)(rng.uniformRv() * candidateRegraftingLocations.size())];
                Node* b = a->getAncestor();
                
                double regraftBL = getBranchLength(a, b);
                double regraftBottomHalf = rng.uniformRv() * regraftBL;
                double regraftTopHalf = regraftBL-regraftBottomHalf;
                
                a->removeNeighbor(b);
                removeBranch(a, b);
                
                a->addNeighbor(regraft);
                a->setAncestor(regraft);
                b->addNeighbor(regraft);
                regraft->setAncestor(b);
                setBranch(a, regraft, regraftTopHalf);
                setBranch(b, regraft, regraftBottomHalf);
            }else{
    //            std::cout << "Move root subtree true " << std::endl;
                //regrafting an internal node subtending the root-- really means we're moving the (tip1, tip2), root subtree
                double bl = getBranchLength(p, q);
                //first-- sample new root, reroot, and proceed as normal
                std::set<Node*> qNeighbors = q->getNeighbors();
                Node* newRoot = nullptr;
                for(Node* n : qNeighbors)
                    if(n->getNeighbors().size() == 3){
                        if(n != root){
                            newRoot = n;
                            break;
                        }
                    }
                if(newRoot == nullptr)
                    Msg::error("could not identify new root");
                
                reroot(newRoot);
                q->removeNeighbor(p);
                collapseNode(q);
                initializeDownPassSequence();
                
                //finding a reattachment point (cannot be root)
                std::vector<Node*> candidateRegraftingLocations = {};
                for(Node* n : downPassSequence)
                    if(n != root && n != freeNode)
                        candidateRegraftingLocations.push_back(n);
                if(candidateRegraftingLocations.size() == 0){
                    print();
                    std::cout << getNewickString() << std::endl;
                    Msg::error("number of candidate SPR reattachment moves = 0");
                }
                
                Node* regraft = freeNode;
    //            Node* regraft = addNode();
                regraft->setIndex(nodes.size()+1);
                p->addNeighbor(regraft);
                p->setAncestor(regraft);
                setBranch(p, regraft, bl);
                
                //reattach q on branch subtending a
                Node* a = candidateRegraftingLocations[(int)(rng.uniformRv() * candidateRegraftingLocations.size())];
                Node* b = a->getAncestor();
                
                double regraftBL = getBranchLength(a, b);
                double regraftBottomHalf = rng.uniformRv() * regraftBL;
                double regraftTopHalf = regraftBL-regraftBottomHalf;
                
                a->removeNeighbor(b);
                removeBranch(a, b);
                
                a->addNeighbor(regraft);
                a->setAncestor(regraft);
                b->addNeighbor(regraft);
                regraft->setAncestor(b);
                setBranch(a, regraft, regraftTopHalf);
                setBranch(b, regraft, regraftBottomHalf);
            }

        }else{
            //we need to remove p (which is binary), reroot on an internal node, and then force binary again
            
            Node* newRoot = nullptr;
            for(Node* n : pNeighbors){
                if(n->getNeighbors().size() == 3){
                    newRoot = n;
                    break;
                }
            }
            reroot(newRoot);
            
            Node* pAbove = p->getDescendants()[0];
            double blAbove = getBranchLength(pAbove, p);
            double blBelow = getBranchLength(p, newRoot);
            
            removeBranch(pAbove, p);
            removeBranch(p, newRoot);
            p->removeNeighbor(pAbove);
            p->removeNeighbor(newRoot);
            
            pAbove->addNeighbor(newRoot);
            pAbove->setAncestor(newRoot);
            setBranch(pAbove, newRoot, blAbove+blBelow);
            
            nodes.erase(find(nodes.begin(), nodes.end(), p));
            delete p;
            
            initializeDownPassSequence();
            pNeighbors.clear();
            p = root;
            pNeighbors = p->getNeighbors();
            
            //sample new q:
            idx = (int)(rng.uniformRv() * pNeighbors.size());
            q = *next(pNeighbors.begin(), idx);
            
            //root now has three descendants
            //both neighbors of q tips with q int node?
            bool moveRootSubtree = false;
            pNeighbors.erase(q);
            if((*next(pNeighbors.begin(), 0))->getIsTip() == true && (*next(pNeighbors.begin(), 1))->getIsTip() == true)
                moveRootSubtree = true;
            if(moveRootSubtree == false){
    //            std::cout << "Move root subtree false " << std::endl;
                double bl = getBranchLength(p, q);
                q->removeNeighbor(p);
                collapseNode(p);
                initializeDownPassSequence();
                //finding a reattachment point (cannot be root)
                std::vector<Node*> candidateRegraftingLocations = {};
                for(Node* n : downPassSequence)
                    if(n != root && n != freeNode)
                        candidateRegraftingLocations.push_back(n);
                if(candidateRegraftingLocations.size() == 0){
                    print();
                    std::cout << getNewickString() << std::endl;
                    Msg::error("number of candidate SPR reattachment moves = 0");
                }
                
    //            Node* regraft = addNode();
                Node* regraft = freeNode;
                regraft->setIndex(nodes.size()+1);
                q->addNeighbor(regraft);
                q->setAncestor(regraft);
                setBranch(q, regraft, bl);

                
                //reattach q on branch subtending a
                Node* a = candidateRegraftingLocations[(int)(rng.uniformRv() * candidateRegraftingLocations.size())];
                Node* b = a->getAncestor();
                
                double regraftBL = getBranchLength(a, b);
                double regraftBottomHalf = rng.uniformRv() * regraftBL;
                double regraftTopHalf = regraftBL-regraftBottomHalf;
                
                a->removeNeighbor(b);
                removeBranch(a, b);
                
                a->addNeighbor(regraft);
                a->setAncestor(regraft);
                b->addNeighbor(regraft);
                regraft->setAncestor(b);
                setBranch(a, regraft, regraftTopHalf);
                setBranch(b, regraft, regraftBottomHalf);
            }else{
    //            std::cout << "Move root subtree true " << std::endl;
                //regrafting an internal node subtending the root-- really means we're moving the (tip1, tip2), root subtree
                double bl = getBranchLength(p, q);
                //first-- sample new root, reroot, and proceed as normal
                std::set<Node*> qNeighbors = q->getNeighbors();
                Node* newRoot = nullptr;
                for(Node* n : qNeighbors)
                    if(n->getNeighbors().size() == 3){
                        if(n != root){
                            newRoot = n;
                            break;
                        }
                    }
                if(newRoot == nullptr)
                    Msg::error("could not identify new root");
                
                reroot(newRoot);
                q->removeNeighbor(p);
                collapseNode(q);
                initializeDownPassSequence();
                
                //finding a reattachment point (cannot be root)
                std::vector<Node*> candidateRegraftingLocations = {};
                for(Node* n : downPassSequence)
                    if(n != root && n != freeNode)
                        candidateRegraftingLocations.push_back(n);
                if(candidateRegraftingLocations.size() == 0){
                    print();
                    std::cout << getNewickString() << std::endl;
                    Msg::error("number of candidate SPR reattachment moves = 0");
                }
                
                Node* regraft = freeNode;
    //            Node* regraft = addNode();
                regraft->setIndex(nodes.size()+1);
                p->addNeighbor(regraft);
                p->setAncestor(regraft);
                setBranch(p, regraft, bl);
                
                //reattach q on branch subtending a
                Node* a = candidateRegraftingLocations[(int)(rng.uniformRv() * candidateRegraftingLocations.size())];
                Node* b = a->getAncestor();
                
                double regraftBL = getBranchLength(a, b);
                double regraftBottomHalf = rng.uniformRv() * regraftBL;
                double regraftTopHalf = regraftBL-regraftBottomHalf;
                
                a->removeNeighbor(b);
                removeBranch(a, b);
                
                a->addNeighbor(regraft);
                a->setAncestor(regraft);
                b->addNeighbor(regraft);
                regraft->setAncestor(b);
                setBranch(a, regraft, regraftTopHalf);
                setBranch(b, regraft, regraftBottomHalf);
            }
            forceBinary();
//            print();
            initializeDownPassSequence();
        }
    }else{
//        std::cout << "p is not root" << std::endl;
        if(p->getAncestor() == q){
//            std::cout << "here" << std::endl;
            Node* transfer = p;
            p = q;
            q = transfer;
        }
        
        double bl = getBranchLength(p, q);
        q->removeNeighbor(p);
        collapseNode(p);
        //finding a reattachment point (cannot be root)
        std::vector<Node*> candidateRegraftingLocations = {};
        for(Node* n : downPassSequence)
            if(n != root && n != freeNode)
                candidateRegraftingLocations.push_back(n);
        if(candidateRegraftingLocations.size() == 0){
            print();
            std::cout << getNewickString() << std::endl;
            Msg::error("number of candidate SPR reattachment moves = 0");
        }
        
//        Node* regraft = addNode();
        Node* regraft = freeNode;
        regraft->setIndex(nodes.size()+1);
        q->addNeighbor(regraft);
        q->setAncestor(regraft);
        setBranch(q, regraft, bl);
        
        //reattach q on branch subtending a
        Node* a = candidateRegraftingLocations[(int)(rng.uniformRv() * candidateRegraftingLocations.size())];
        Node* b = a->getAncestor();
        
//        std::cout << "Regrafting location: " << a->getIndex() << " to " << b->getIndex() << std::endl;
        
        double regraftBL = getBranchLength(a, b);
        double regraftBottomHalf = rng.uniformRv() * regraftBL;
        double regraftTopHalf = regraftBL-regraftBottomHalf;
        
        a->removeNeighbor(b);
        removeBranch(a, b);
        
        a->addNeighbor(regraft);
        a->setAncestor(regraft);
        b->addNeighbor(regraft);
        regraft->setAncestor(b);
        setBranch(a, regraft, regraftTopHalf);
        setBranch(b, regraft, regraftBottomHalf);
    }
    
    int endingNodesSize = nodes.size();
//    std::cout << "endingNodesSize " <<  endingNodesSize << std::endl;
    if(endingNodesSize != startingNodeSize)
        Msg::error("check node arithmatic");

    freeNode = nullptr;
    
    initializeDownPassSequence();
    reindexNodes();
    
//    print("After move");
}

void Tree::rSPR(std::string s){
//    print("before move");
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    Node* tip = nullptr;
    for(Node* n : downPassSequence)
        if(n->getName() == s){
            tip = n;
            break;
        }
    Node* p = tip->getAncestor();

    std::set<Node*> pNeighbors = p->getNeighbors();
    if(pNeighbors.size() != 3){
        std::cout << p->getIndex() << std::endl;
        print();
        Msg::error("pNeighbor size is not 3");
    }
    
    //sample neighbor at random
    Node* q = tip;
    
    double bl = getBranchLength(p, q);
    q->removeNeighbor(p);
    collapseNode(p);
    //finding a reattachment point (cannot be root)
    std::vector<Node*> candidateRegraftingLocations = {};
    for(Node* n : downPassSequence)
        if(n != root && n != freeNode)
            candidateRegraftingLocations.push_back(n);
    if(candidateRegraftingLocations.size() == 0){
        print();
        std::cout << getNewickString() << std::endl;
        Msg::error("number of candidate SPR reattachment moves = 0");
    }
    
//        Node* regraft = addNode();
    Node* regraft = freeNode;
    regraft->setIndex(nodes.size()+1);
    q->addNeighbor(regraft);
    q->setAncestor(regraft);
    setBranch(q, regraft, bl);
    
    //reattach q on branch subtending a
    Node* a = candidateRegraftingLocations[(int)(rng.uniformRv() * candidateRegraftingLocations.size())];
    Node* b = a->getAncestor();
    
//        std::cout << "Regrafting location: " << a->getIndex() << " to " << b->getIndex() << std::endl;
    
    double regraftBL = getBranchLength(a, b);
    double regraftBottomHalf = rng.uniformRv() * regraftBL;
    double regraftTopHalf = regraftBL-regraftBottomHalf;
    
    a->removeNeighbor(b);
    removeBranch(a, b);
    
    a->addNeighbor(regraft);
    a->setAncestor(regraft);
    b->addNeighbor(regraft);
    regraft->setAncestor(b);
    setBranch(a, regraft, regraftTopHalf);
    setBranch(b, regraft, regraftBottomHalf);
    
    int endingNodesSize = nodes.size();

    freeNode = nullptr;
    
    initializeDownPassSequence();
    reindexNodes();
//    print("after move");
}

void Tree::setBranch(Node* e1, Node* e2, double x) {
    std::pair<Node*,Node*> key;
    initializeBranchLengthKey(key, e1, e2);
    BranchLengths::iterator it = branchLengths.find(key);
    if (it != branchLengths.end())
        {
        it->second = x;
        }
    else
        {
        branchLengths.insert( std::make_pair(key, x) );
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
        if (p != root)
            std::cout << getBranchLength(p, p->getAncestor()) << " ";
        else
            std::cout << "--- ";
        if (p == root)
            std::cout << "<- Root ";
        std::cout << std::endl;

        for (Node* n : pNeighbors)
            {
            if (n != p->getAncestor())
                showNode(n, indent + 3);
            }
            
        }
}

double Tree::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    if(rng.uniformRv() < 0.5)
        return updateRootAge();
    return updateNodeAge();
}

double Tree::updateRootAge(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    double maxChildAge = 0.0;
    for(Node* c : root->getNeighbors())
        if(c->getTime() > maxChildAge)
            maxChildAge = c->getTime();

    double m = std::exp( 4.0 * (rng.uniformRv() - 0.5) );
    double newRootAge = root->getTime() * m;
    if(newRootAge <= maxChildAge)
        return -INFINITY;

    root->setTime(newRootAge);

    return std::log(m);
}

double Tree::updateNodeAge(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    std::vector<Node*> candidates;
    for(Node* n : downPassSequence)
        if(n != root && n->getIsTip() == false)
            candidates.push_back(n);
    if(candidates.size() == 0)
        return 0.0;

    Node* n = candidates[(int)(rng.uniformRv() * candidates.size())];
    double parentAge = n->getAncestor()->getTime();
    double maxChildAge = 0.0;
    for(Node* nb : n->getNeighbors())
        if(nb != n->getAncestor() && nb->getTime() > maxChildAge)
            maxChildAge = nb->getTime();

    n->setTime(maxChildAge + rng.uniformRv() * (parentAge - maxChildAge));

    return 0.0;
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
