#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

ParameterTree::ParameterTree(double prob, PhylogeneticModel* m) :
    Parameter(prob, m,"Tree"),
    numRejections(0),
    numAcceptances(0),
    useCachedLnP(false){
    //constructor used in derived classes
    trees[0] = nullptr;
    trees[1] = nullptr;
}

ParameterTree::ParameterTree(double prob, PhylogeneticModel* m, std::vector<std::string> taxonNames, double lam) :
    Parameter(prob, m,"Tree"),
    lambda(lam), //exponential prior on branch length
    numRejections(0),
    numAcceptances(0),
    useCachedLnP(false){
    //ParameterTree base class is responsible for topology and branch length moves only

    trees[0] = new Tree(taxonNames, lambda);
    trees[1] = new Tree(*trees[0]);
}

double ParameterTree::lnProbability(void) {
    if(useCachedLnP == false){
        Tree* t = trees[0];
        cachedLnP = 0.0;
        for (Node* p : t->getDownPassSequence())
            {
            if (p != t->getRoot())
                {
                double v = t->getBranchLength(p, p->getAncestor());
                cachedLnP += Probability::Exponential::lnPdf(lambda, v);
                }
            }
        useCachedLnP = true;
    }
    return cachedLnP;
}

void ParameterTree::print(void) {

    trees[0]->print();
}

void ParameterTree::setTree(Tree* t){
    delete trees[0];
    delete trees[1];
    trees[0] = new Tree(*t);
    trees[1] = new Tree(*trees[0]);
}

double ParameterTree::update(void) {
    //By convention: Tree class is responsible for ALL topology and branch length moves; responsible for nodes
    useCachedLnP = false;
    return trees[0]->update(); //right now, is just simple branch length update
}

void ParameterTree::updateForAcceptance(void) {
    numAcceptances++;
    useCachedLnP = false;
    *(trees[1]) = *(trees[0]);
}

void ParameterTree::updateForRejection(void) {
    numRejections++;
    useCachedLnP = false;
    *(trees[0]) = *(trees[1]);
}
