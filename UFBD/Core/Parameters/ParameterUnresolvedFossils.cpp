#include <cmath>

#include "FBDInput.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterDouble.hpp"
#include "ParameterUnresolvedFossils.hpp"
#include "PhylogeneticModel.hpp"
#include "RandomVariable.hpp"
#include "Serialize.hpp"
#include "Tree.hpp"

ParameterUnresolvedFossils::ParameterUnresolvedFossils(double prob, PhylogeneticModel* m, Tree* bb, std::vector<Clade>& clades, std::vector<Fossil>& fossils, ParameterDouble* oa) : Parameter(prob, m, "unresolvedFossils"){
    backbone = bb;
    originAge = oa;
    numFossils = (int)fossils.size();
    numAcceptances = 0;
    numRejections = 0;
    lastFossil = -1;
    lastMove = SINGLE;
    spineIdx = -1;

    yMin.resize(numFossils);
    yMax.resize(numFossils);
    crownNode.resize(numFossils);
    originNode.resize(numFossils);
    isCrown.resize(numFossils);
    isStem.resize(numFossils);
    ue.resize(numFossils);
    y[0].resize(numFossils);
    y[1].resize(numFossils);
    z[0].resize(numFossils);
    z[1].resize(numFossils);
    sa[0].assign(numFossils, 0);
    sa[1].assign(numFossils, 0);
    az[0].assign(numFossils, 0);
    az[1].assign(numFossils, 0);
    azDomain.assign(numFossils, std::vector<int>(1, 0));

    for(int i = 0; i < numFossils; i++){
        Fossil& f = fossils[i];
        Clade* clade = nullptr;
        for(Clade& c : clades)
            if(c.getName() == f.getClade()){
                clade = &c;
                break;
            }
        Node* cr = clade->getTaxa().empty() ? backbone->getCrown() : backbone->getMRCA(clade->getTaxa());
        crownNode[i]  = cr;
        originNode[i] = cr->getAncestor();
        isCrown[i]    = (f.getAssignment() == Assignment::CROWN);
        isStem[i]     = (f.getAssignment() == Assignment::STEM);
        ue[i]         = (f.getMaxAge() == 0.0);
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

    if(backbone->getNumBackbone() == 0 && numFossils > 0 && originAge != nullptr){
        spineIdx = 0;
        for(int i = 1; i < numFossils; i++)
            if(y[0][i] < y[0][spineIdx])
                spineIdx = i;
    }
}

void ParameterUnresolvedFossils::setAttachmentZoneDomain(std::vector<std::vector<int> >& domain){
    azDomain = domain;
    for(int i = 0; i < numFossils; i++){
        az[0][i] = azDomain[i][0];
        az[1][i] = azDomain[i][0];
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
    lastMove = SINGLE;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    lastFossil = (int)(rng.uniformRv() * numFossils);
    int i = lastFossil;
    if(ue[i])
        return updateAttachAge(i);                          // a UE has no y and no SA state
    double u = rng.uniformRv();
    if(u < 1.0/3.0)
        return updateSampledAncestor(i);
    if(u < 2.0/3.0)
        return updateFossilAge(i);
    return updateAttachAge(i);
}

double ParameterUnresolvedFossils::updateFossilAge(int i){
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    if(spineIdx < 0){
        if(sa[0][i]){
            double ceiling = getMaxAttachAge(i);
            double saHi = (yMax[i] < ceiling) ? yMax[i] : ceiling;
            y[0][i] = yMin[i] + rng.uniformRv() * (saHi - yMin[i]);
            z[0][i] = y[0][i];
            return 0.0;
        }
        double hi = (yMax[i] < z[0][i]) ? yMax[i] : z[0][i];
        y[0][i] = yMin[i] + rng.uniformRv() * (hi - yMin[i]);
        return 0.0;
    }

    int Sold = spineIdx;
    bool iWasSA = (sa[0][i] != 0);
    double yiNew = yMin[i] + rng.uniformRv() * (yMax[i] - yMin[i]);
    int Snew = i;
    double bestY = yiNew;
    for(int j = 0; j < numFossils; j++){
        double yj = (j == i) ? yiNew : y[0][j];
        if(yj < bestY){ bestY = yj; Snew = j; }
    }
    if(Snew == Sold){
        if(i != Sold){
            if(iWasSA){
                y[0][i] = yiNew;
                z[0][i] = yiNew;
            }else{
                if(z[0][i] < yiNew)
                    return -INFINITY;
                y[0][i] = yiNew;
            }
        }else{
            y[0][i] = yiNew;
        }
        return 0.0;
    }
    bool SnewWasSA = (sa[0][Snew] != 0);
    double ynewSold = (i == Sold) ? yiNew : y[0][Sold];
    double zSold;
    if(SnewWasSA){
        zSold = ynewSold;
    }else{
        zSold = z[0][Snew];
        if(zSold < ynewSold)
            return -INFINITY;
    }
    lastMove = FLIP;
    flipS = Sold;
    flipT = Snew;
    flipSy = y[0][Sold];
    flipSz = z[0][Sold];
    flipTy = y[0][Snew];
    flipTz = z[0][Snew];
    flipSsa = sa[0][Sold];
    flipTsa = sa[0][Snew];
    y[0][i] = yiNew;
    z[0][Sold] = zSold;
    sa[0][Sold] = SnewWasSA ? 1 : 0;
    sa[0][Snew] = 0;
    spineIdx = Snew;
    return 0.0;
}

double ParameterUnresolvedFossils::updateAttachAge(int i){
    if(i == spineIdx)
        return 0.0;
    if(sa[0][i])
        return 0.0;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double lo = getMinAttachAge(i);
    double hi = getMaxAttachAge(i);
    z[0][i] = lo + rng.uniformRv() * (hi - lo);
    return 0.0;
}

double ParameterUnresolvedFossils::updateSampledAncestor(int i){
    if(i == spineIdx)
        return 0.0;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double lo = getMinAttachAge(i);
    double hi = getMaxAttachAge(i);
    double range = hi - lo;
    if(range <= 0.0)
        return -INFINITY;
    if(sa[0][i]){
        z[0][i] = lo + rng.uniformRv() * range;
        sa[0][i] = 0;
        return std::log(range);
    }
    z[0][i] = y[0][i];
    sa[0][i] = 1;
    return -std::log(range);
}

double ParameterUnresolvedFossils::scaleAllAttachAges(double m){
    lastMove = BULK;
    for(int i = 0; i < numFossils; i++){
        if(i == spineIdx || sa[0][i])
            continue;
        if(z[0][i] * m < y[0][i])
            return -INFINITY;
    }
    int count = 0;
    for(int i = 0; i < numFossils; i++){
        if(i == spineIdx || sa[0][i])
            continue;
        z[0][i] *= m;
        count++;
    }
    return count * std::log(m);
}

double ParameterUnresolvedFossils::scaleAttachAges(const std::vector<int>& indices, double m){
    lastMove = BULK;
    for(int i : indices)
        if(i != spineIdx && z[0][i] * m < y[0][i])
            return -INFINITY;
    int count = 0;
    for(int i : indices)
        if(i != spineIdx){ z[0][i] *= m; count++; }
    return count * std::log(m);
}

void ParameterUnresolvedFossils::updateForAcceptance(void){
    numAcceptances++;
    if(numFossils == 0)
        return;
    if(lastMove == FLIP){
        y[1][flipS] = y[0][flipS]; z[1][flipS] = z[0][flipS]; sa[1][flipS] = sa[0][flipS];
        y[1][flipT] = y[0][flipT]; z[1][flipT] = z[0][flipT]; sa[1][flipT] = sa[0][flipT];
    }else if(lastMove == BULK){
        for(int i = 0; i < numFossils; i++)
            z[1][i] = z[0][i];
    }else{
        y[1][lastFossil] = y[0][lastFossil];
        z[1][lastFossil] = z[0][lastFossil];
        sa[1][lastFossil] = sa[0][lastFossil];
        az[1][lastFossil] = az[0][lastFossil];
    }
}

void ParameterUnresolvedFossils::updateForRejection(void){
    numRejections++;
    if(numFossils == 0)
        return;
    if(lastMove == FLIP){
        y[0][flipS] = flipSy; z[0][flipS] = flipSz; y[1][flipS] = flipSy; z[1][flipS] = flipSz;
        y[0][flipT] = flipTy; z[0][flipT] = flipTz; y[1][flipT] = flipTy; z[1][flipT] = flipTz;
        sa[0][flipS] = sa[1][flipS] = flipSsa;
        sa[0][flipT] = sa[1][flipT] = flipTsa;
        spineIdx = flipS;
    }else if(lastMove == BULK){
        for(int i = 0; i < numFossils; i++)
            z[0][i] = z[1][i];
    }else{
        y[0][lastFossil] = y[1][lastFossil];
        z[0][lastFossil] = z[1][lastFossil];
        sa[0][lastFossil] = sa[1][lastFossil];
        az[0][lastFossil] = az[1][lastFossil];
    }
}

void ParameterUnresolvedFossils::writeState(std::ostream& os){
    Serialize::writeVec(os, y[1]);
    Serialize::writeVec(os, z[1]);
    Serialize::writeCVec(os, sa[1]);
    Serialize::writeIVec(os, az[1]);
    os << spineIdx << ' ' << numAcceptances << ' ' << numRejections << '\n';
}

void ParameterUnresolvedFossils::readState(std::istream& is){
    Serialize::readVec(is, y[1]);
    Serialize::readVec(is, z[1]);
    Serialize::readCVec(is, sa[1]);
    sa[0] = sa[1];
    Serialize::readIVec(is, az[1]);
    y[0] = y[1];
    z[0] = z[1];
    az[0] = az[1];
    is >> spineIdx >> numAcceptances >> numRejections;
}

void ParameterUnresolvedFossils::print(void){
}
