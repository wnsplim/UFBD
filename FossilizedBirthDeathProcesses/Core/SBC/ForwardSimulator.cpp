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

    const long maxAttempts = 2000000;
    long attempts = 0;
    bool ok = false;

    while(ok == false){
        if(++attempts > maxAttempts)
            Msg::error("ForwardSimulator: exceeded max conditioning attempts (params likely make the conditioning unattainable)");

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

        double curAge = p.startAge;
        bool failed = false;
        while(curAge > 0.0 && active.empty() == false){
            double lam = p.lambda[ii];
            double m   = p.mu[ii];
            double ps  = p.psi[ii];
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
                SimNode* f = newNode(nextAge, allNodes); f->parent = n; f->isFossil = true; n->left = f;
                SimNode* c = newNode(-1.0, allNodes); c->parent = n; n->right = c;
                active[k] = c;
            }
            curAge = nextAge;
            if((int)active.size() > maxLineages){ failed = true; break; }
        }
        if(failed)
            continue;

        for(SimNode* n : active){
            n->age = 0.0;
            n->extantSampled = (rng->uniformRv() < p.rho);
        }

        markKeep(root, false);
        bool pass;
        if(p.originConditioning == false)
            pass = (root->left->keep && root->right->keep);
        else
            pass = root->keep;
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
        if(inBB == 0 && extants.empty() == false){
            extants[0]->inBackbone = true;
            inBB = 1;
        }

        ok = true;

        for(SimNode* n : allNodes){
            if(n->isFossil){
                res.fossilAges.push_back(n->age);
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
        res.crownAge = mrca->age;
        res.originAge = p.originConditioning ? p.startAge : mrca->age;

        std::stringstream s;
        writeNewick(mrca, s);
        s << ";";
        res.backboneNewick = s.str();
    }

    for(SimNode* n : allNodes)
        delete n;

    return res;
}
