#include <iomanip>
#include <iostream>
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "PhylogeneticModel.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Serialize.hpp"
#include "Tree.hpp"

ParameterTree::ParameterTree(double prob, PhylogeneticModel* m) :
    Parameter(prob, m,"Tree"),
    numRejections(0),
    numAcceptances(0),
    useCachedLnP(false),
    scaleLambda(1.0),
    numScaleMoves(0){
    //constructor used in derived classes
    trees[0] = nullptr;
    trees[1] = nullptr;
}

double ParameterTree::lnProbability(void) {
    return 0.0;
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
    return trees[0]->update(scaleLambda);
}

void ParameterTree::updateForAcceptance(void) {
    numAcceptances++;
    if(trees[0]->getLastUpdateWasScale())
        tuneScale(true);
    useCachedLnP = false;
    *(trees[1]) = *(trees[0]);
}

void ParameterTree::updateForRejection(void) {
    numRejections++;
    if(trees[0]->getLastUpdateWasScale())
        tuneScale(false);
    useCachedLnP = false;
    *(trees[0]) = *(trees[1]);
}

void ParameterTree::writeState(std::ostream& os) {
    trees[1]->writeState(os);
    os << scaleLambda << ' ' << numAcceptances << ' ' << numRejections << ' ' << numScaleMoves << '\n';
    Serialize::writeBoolDeque(os, recentScaleAcceptRej);
}

void ParameterTree::readState(std::istream& is) {
    trees[1]->readState(is);
    *(trees[0]) = *(trees[1]);
    is >> scaleLambda >> numAcceptances >> numRejections >> numScaleMoves;
    Serialize::readBoolDeque(is, recentScaleAcceptRej);
    useCachedLnP = false;
}

void ParameterTree::tuneScale(bool accepted) {
    recentScaleAcceptRej.push_back(accepted);
    if(recentScaleAcceptRej.size() > 1000)
        recentScaleAcceptRej.pop_front();
    numScaleMoves++;
    if(numScaleMoves % 100 == 0 && numScaleMoves < 10000){
        double numAccepted = 0.0;
        for(bool b : recentScaleAcceptRej)
            if(b == true)
                numAccepted++;
        double rate = numAccepted / recentScaleAcceptRej.size();
        if(rate < 0.44 - 0.2)
            scaleLambda /= 1.1;
        else if(rate > 0.44 + 0.2)
            scaleLambda *= 1.1;
    }
}
