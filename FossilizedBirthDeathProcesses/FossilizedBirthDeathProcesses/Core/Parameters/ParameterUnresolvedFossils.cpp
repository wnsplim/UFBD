#include "FBDInput.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterUnresolvedFossils.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

ParameterUnresolvedFossils::ParameterUnresolvedFossils(double prob, PhylogeneticModel* m, Tree* bb, std::vector<Clade>& clades, std::vector<Fossil>& fossils) : Parameter(prob, m, "unresolvedFossils"){
    backbone = bb;
    numFossils = (int)fossils.size();
    numAcceptances = 0;
    numRejections = 0;
    lastFossil = -1;
    lastWasAttach = false;

    yMin.resize(numFossils);
    yMax.resize(numFossils);
    crownNode.resize(numFossils);
    originNode.resize(numFossils);
    isCrown.resize(numFossils);
    y[0].resize(numFossils);
    y[1].resize(numFossils);
    z[0].resize(numFossils);
    z[1].resize(numFossils);

    for(int i = 0; i < numFossils; i++){
        Fossil& f = fossils[i];
        Clade* clade = nullptr;
        for(Clade& c : clades)
            if(c.getName() == f.getClade()){
                clade = &c;
                break;
            }
        if(clade == nullptr)
            Msg::error("unresolved fossil '" + f.getTaxon() + "' references undefined clade '" + f.getClade() + "'");
        Node* cr = backbone->getMRCA(clade->getTaxa());
        crownNode[i]  = cr;
        originNode[i] = cr->getAncestor();
        isCrown[i]    = (f.getAssignment() == Assignment::CROWN);
        yMin[i]       = f.getMinAge();
        yMax[i]       = f.getMaxAge();
    }

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    for(int i = 0; i < numFossils; i++){
        double yi = yMin[i] + rng.uniformRv() * (yMax[i] - yMin[i]);
        y[0][i] = yi;
        y[1][i] = yi;
        double lo = getMinAttachAge(i);
        double hi = getMaxAttachAge(i);
        double zi = lo + rng.uniformRv() * (hi - lo);
        z[0][i] = zi;
        z[1][i] = zi;
    }
}

double ParameterUnresolvedFossils::getMinAttachAge(int i){
    return y[0][i];
}

double ParameterUnresolvedFossils::getMaxAttachAge(int i){
    return isCrown[i] ? crownNode[i]->getTime() : originNode[i]->getTime();
}

double ParameterUnresolvedFossils::update(void){
    if(numFossils == 0)
        return 0.0;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    lastFossil = (int)(rng.uniformRv() * numFossils);
    if(rng.uniformRv() < 0.5){
        lastWasAttach = false;
        return updateFossilAge(lastFossil);
    }
    lastWasAttach = true;
    return updateAttachAge(lastFossil);
}

double ParameterUnresolvedFossils::updateFossilAge(int i){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double hi = (yMax[i] < z[0][i]) ? yMax[i] : z[0][i]; // keep z older than y
    y[0][i] = yMin[i] + rng.uniformRv() * (hi - yMin[i]);
    return 0.0;
}

double ParameterUnresolvedFossils::updateAttachAge(int i){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double lo = getMinAttachAge(i);
    double hi = getMaxAttachAge(i);
    z[0][i] = lo + rng.uniformRv() * (hi - lo);
    return 0.0;
}

void ParameterUnresolvedFossils::updateForAcceptance(void){
    numAcceptances++;
    if(lastWasAttach)
        z[1][lastFossil] = z[0][lastFossil];
    else
        y[1][lastFossil] = y[0][lastFossil];
}

void ParameterUnresolvedFossils::updateForRejection(void){
    numRejections++;
    if(lastWasAttach)
        z[0][lastFossil] = z[1][lastFossil];
    else
        y[0][lastFossil] = y[1][lastFossil];
}

void ParameterUnresolvedFossils::print(void){
}
