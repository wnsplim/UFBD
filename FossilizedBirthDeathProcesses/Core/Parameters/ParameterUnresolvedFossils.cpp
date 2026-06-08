#include <cmath>

#include "FBDInput.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterDouble.hpp"
#include "ParameterUnresolvedFossils.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

ParameterUnresolvedFossils::ParameterUnresolvedFossils(double prob, PhylogeneticModel* m, Tree* bb, std::vector<Clade>& clades, std::vector<Fossil>& fossils, ParameterDouble* oa) : Parameter(prob, m, "unresolvedFossils"){
    backbone = bb;
    originAge = oa;
    numFossils = (int)fossils.size();
    numAcceptances = 0;
    numRejections = 0;
    lastFossil = -1;
    lastWasBulk = false;

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
    if(isCrown[i])
        return crownNode[i]->getTime();
    if(crownNode[i] == backbone->getCrown())
        return (originAge != nullptr) ? originAge->getValue() : backbone->getCrown()->getTime();
    return originNode[i]->getTime();
}

double ParameterUnresolvedFossils::update(void){
    if(numFossils == 0)
        return 0.0;
    lastWasBulk = false;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    lastFossil = (int)(rng.uniformRv() * numFossils);
    double u = rng.uniformRv();
    if(u < 1.0/3.0)
        return updateSampledAncestor(lastFossil);
    if(u < 2.0/3.0)
        return updateFossilAge(lastFossil);
    return updateAttachAge(lastFossil);
}

double ParameterUnresolvedFossils::updateFossilAge(int i){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    if(z[0][i] == y[0][i]){
        double ceiling = getMaxAttachAge(i);
        double saHi = (yMax[i] < ceiling) ? yMax[i] : ceiling;
        y[0][i] = yMin[i] + rng.uniformRv() * (saHi - yMin[i]);
        z[0][i] = y[0][i];
        return 0.0;
    }
    double hi = (yMax[i] < z[0][i]) ? yMax[i] : z[0][i]; // keep z older than y
    y[0][i] = yMin[i] + rng.uniformRv() * (hi - yMin[i]);
    return 0.0;
}

double ParameterUnresolvedFossils::updateAttachAge(int i){
    if(z[0][i] == y[0][i])
        return 0.0;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double lo = getMinAttachAge(i);
    double hi = getMaxAttachAge(i);
    z[0][i] = lo + rng.uniformRv() * (hi - lo);
    return 0.0;
}

double ParameterUnresolvedFossils::updateSampledAncestor(int i){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double lo = getMinAttachAge(i);
    double hi = getMaxAttachAge(i);
    double range = hi - lo;
    if(range <= 0.0)
        return -INFINITY;
    if(z[0][i] == y[0][i]){
        z[0][i] = lo + rng.uniformRv() * range;
        return std::log(range);
    }
    z[0][i] = y[0][i];
    return -std::log(range);
}

double ParameterUnresolvedFossils::scaleAllAttachAges(double m){
    lastWasBulk = true;
    for(int i = 0; i < numFossils; i++){
        if(z[0][i] == y[0][i])
            continue;
        if(z[0][i] * m < y[0][i])
            return -INFINITY;
    }
    int count = 0;
    for(int i = 0; i < numFossils; i++){
        if(z[0][i] == y[0][i])
            continue;
        z[0][i] *= m;
        count++;
    }
    return count * std::log(m);
}

double ParameterUnresolvedFossils::scaleAttachAges(const std::vector<int>& indices, double m){
    lastWasBulk = true;
    for(int i : indices)
        if(z[0][i] * m < y[0][i])
            return -INFINITY;
    for(int i : indices)
        z[0][i] *= m;
    return (int)indices.size() * std::log(m);
}

void ParameterUnresolvedFossils::updateForAcceptance(void){
    numAcceptances++;
    if(lastWasBulk){
        for(int i = 0; i < numFossils; i++)
            z[1][i] = z[0][i];
    }else{
        y[1][lastFossil] = y[0][lastFossil];
        z[1][lastFossil] = z[0][lastFossil];
    }
}

void ParameterUnresolvedFossils::updateForRejection(void){
    numRejections++;
    if(lastWasBulk){
        for(int i = 0; i < numFossils; i++)
            z[0][i] = z[1][i];
    }else{
        y[0][lastFossil] = y[1][lastFossil];
        z[0][lastFossil] = z[1][lastFossil];
    }
}

void ParameterUnresolvedFossils::print(void){
}
