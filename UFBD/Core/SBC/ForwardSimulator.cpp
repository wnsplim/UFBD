#include "ForwardSimulator.hpp"
#include "Msg.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"

#include <sstream>

namespace {

struct SimNode {
    double    age;
    SimNode*  parent;
    SimNode*  left;
    SimNode*  right;
    bool      isFossil;
    bool      extantSampled;
    bool      inBackbone;
    bool      keep;
    int       label;
    int       type;
};

SimNode* newNode(double age, std::vector<SimNode*>& all){
    SimNode* n = new SimNode();
    n->age = age;
    n->parent = nullptr;
    n->left = nullptr;
    n->right = nullptr;
    n->isFossil = false;
    n->extantSampled = false;
    n->inBackbone = false;
    n->keep = false;
    n->label = 0;
    n->type = 0;
    all.push_back(n);
    return n;
}

void markKeep(SimNode* n, bool useBackbone){
    if(n->left == nullptr && n->right == nullptr){
        n->keep = useBackbone ? n->inBackbone : n->extantSampled;
        return;
    }
    bool k = false;
    if(n->left != nullptr){ markKeep(n->left, useBackbone); k = n->left->keep || k; }
    if(n->right != nullptr){ markKeep(n->right, useBackbone); k = n->right->keep || k; }
    n->keep = k;
}

SimNode* extantMRCA(SimNode* n){
    while(true){
        if(n->left == nullptr && n->right == nullptr)
            return n;
        SimNode* kl = (n->left != nullptr && n->left->keep) ? n->left : nullptr;
        SimNode* kr = (n->right != nullptr && n->right->keep) ? n->right : nullptr;
        if(kl != nullptr && kr != nullptr)
            return n;
        if(kl != nullptr){ n = kl; continue; }
        if(kr != nullptr){ n = kr; continue; }
        return n;
    }
}

void assignLabels(SimNode* n, int& next){
    if(n->keep == false)
        return;
    if(n->left == nullptr && n->right == nullptr){
        n->label = next++;
        return;
    }
    if(n->left != nullptr)
        assignLabels(n->left, next);
    if(n->right != nullptr)
        assignLabels(n->right, next);
}

void writeNewick(SimNode* n, std::stringstream& s){
    SimNode* kl = (n->left != nullptr && n->left->keep) ? n->left : nullptr;
    SimNode* kr = (n->right != nullptr && n->right->keep) ? n->right : nullptr;
    if(kl == nullptr && kr == nullptr){
        s << "T" << n->label;
        return;
    }
    if(kl != nullptr && kr == nullptr){
        writeNewick(kl, s);
        return;
    }
    if(kl == nullptr && kr != nullptr){
        writeNewick(kr, s);
        return;
    }
    s << "(";
    writeNewick(kl, s);
    s << ",";
    writeNewick(kr, s);
    s << ")";
}

}

SimResult ForwardSimulator::simulate(const SimParams& p){

    size_t nInt = p.intervalStart.size();

    SimResult res;
    std::vector<SimNode*> allNodes;

    const long maxAttempts = 1000000;
    long attempts = 0;
    bool ok = false;

    while(ok == false){
        if(++attempts > maxAttempts)
            Msg::error("exceeded max attempts");

        for(SimNode* n : allNodes)
            delete n;
        allNodes.clear();

        SimNode* root = newNode(p.startAge, allNodes);
        std::vector<SimNode*> active;
        if(p.originConditioning == false){
            SimNode* l = newNode(-1.0, allNodes); l->parent = root; root->left = l;
            SimNode* r = newNode(-1.0, allNodes); r->parent = root; root->right = r;
            active.push_back(l);
            active.push_back(r);
        }else{
            SimNode* stem = newNode(-1.0, allNodes); stem->parent = root; root->left = stem;
            active.push_back(stem);
        }

        int ii = (int)nInt - 1;
        while(ii > 0 && p.intervalStart[ii] > p.startAge)
            ii--;

        int numTypes = (int)p.psi.size();
        double curAge = p.startAge;
        while(curAge > 0.0 && active.empty() == false){
            double lam = p.lambda[p.lambdaIdx[ii]];
            double m   = p.mu[p.muIdx[ii]];
            double ps  = 0.0;
            for(int t = 0; t < numTypes; t++) ps += p.psi[t][p.psiIdx[t][ii]];
            double sum = lam + m + ps;
            int N = (int)active.size();
            double younger = p.intervalStart[ii];
            if(sum <= 0.0){
                if(ii == 0){ curAge = 0.0; break; }
                curAge = younger; ii--; continue;
            }
            double nextAge = curAge - Probability::Exponential::rv(rng, N * sum);
            if(nextAge < younger){
                if(ii == 0){ curAge = 0.0; break; }
                curAge = younger; ii--; continue;
            }
            int k = (int)(rng->uniformRv() * N);
            SimNode* n = active[k];
            n->age = nextAge;
            double u = rng->uniformRv();
            if(u < lam / sum){
                SimNode* l = newNode(-1.0, allNodes); l->parent = n; n->left = l;
                SimNode* r = newNode(-1.0, allNodes); r->parent = n; n->right = r;
                active[k] = l;
                active.push_back(r);
            }else if(u < (lam + m) / sum){
                active[k] = active.back();
                active.pop_back();
            }else{
                int ftype = 0;
                if(numTypes > 1){
                    double uu = rng->uniformRv() * ps;
                    double acc = 0.0;
                    for(int t = 0; t < numTypes; t++){ acc += p.psi[t][p.psiIdx[t][ii]]; if(uu < acc){ ftype = t; break; } }
                }
                SimNode* f = newNode(nextAge, allNodes); f->parent = n; f->isFossil = true; f->type = ftype; n->left = f;
                SimNode* c = newNode(-1.0, allNodes); c->parent = n; n->right = c;
                active[k] = c;
            }
            curAge = nextAge;
            if((int)active.size() > maxActiveEdges)
                Msg::error("tree exceeded " + std::to_string(maxActiveEdges) + " edges (lambda=" + std::to_string(lam) + ", mu=" + std::to_string(m) + ", x0=" + std::to_string(p.startAge) + "); net diversification likely too high for the origin age.");
        }

        for(SimNode* n : active){
            n->age = 0.0;
            n->extantSampled = (rng->uniformRv() < p.rho);
        }

        markKeep(root, false);
        bool hasFossil = false;
        for(SimNode* n : allNodes)
            if(n->isFossil){ hasFossil = true; break; }
        bool pass;
        if(p.originConditioning == false){
            pass = (root->left->keep && root->right->keep);
        }else{
            bool hasExtant = root->keep;
            if(p.condEvent == ConditioningEvent::ANYSAMPLE)
                pass = hasExtant || hasFossil;
            else if(p.condEvent == ConditioningEvent::EXTINCT)
                pass = (hasExtant == false) && hasFossil;
            else
                pass = hasExtant;
        }
        if(pass == false)
            continue;

        std::vector<SimNode*> extants;
        for(SimNode* n : allNodes)
            if(n->extantSampled)
                extants.push_back(n);
        int inBB = 0;
        for(SimNode* n : extants){
            n->inBackbone = (rng->uniformRv() < p.bb);
            if(n->inBackbone)
                inBB++;
        }
        ok = true;

        for(SimNode* n : allNodes){
            if(n->isFossil){
                res.fossilAges.push_back(n->age);
                res.fossilTypes.push_back(n->type);
                res.numFossils++;
            }
        }
        res.numExtantSampled = (int)extants.size();
        res.numUE = res.numExtantSampled - inBB;

        markKeep(root, true);
        SimNode* mrca = extantMRCA(root);
        int next = 1;
        assignLabels(mrca, next);
        res.numBackbone = next - 1;
        if(res.numBackbone > 0){
            std::stringstream s;
            writeNewick(mrca, s);
            s << ";";
            res.backboneNewick = s.str();
        }else{
            res.backboneNewick = "";
        }
    }

    for(SimNode* n : allNodes)
        delete n;

    return res;
}
