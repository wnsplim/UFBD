#include "SequenceLikelihood.hpp"
#include "Tree.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "Msg.hpp"
#include "ThreadPool.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

SequenceLikelihood::SequenceLikelihood(int numStates, int numCats) :
    numStates(numStates),
    numCats(numCats),
    numPartitions(0),
    mappedRoot(nullptr),
    cacheValid(false)
{
}

void SequenceLikelihood::addPartition(const std::vector<std::string>& taxa, const std::vector<std::vector<int>>& patterns, const std::vector<int>& weight){
    int npat = (int)weight.size();
    int fullMask = (1 << numStates) - 1;
    std::vector<int> mask(npat, fullMask);
    for(int h = 0; h < npat; h++)
        for(int t = 0; t < (int)taxa.size(); t++)
            mask[h] &= patterns[t][h];

    taxonNames.push_back(taxa);
    patternState.push_back(patterns);
    patternWeight.push_back(weight);
    constantState.push_back(mask);
    rateModel.push_back(GTRrateModel(numStates));
    numPartitions++;
}

void SequenceLikelihood::mapTaxaToNodes(Tree* tree){
    if(tree->getRoot() == mappedRoot)
        return;
    int numNodes = tree->getNumNodes();
    tipStateByOffset.assign(numPartitions, std::vector<std::vector<int>>(numNodes));
    for(int p = 0; p < numPartitions; p++)
        for(int t = 0; t < (int)taxonNames[p].size(); t++){
            Node* nd = tree->getTaxonNode(taxonNames[p][t]);
            if(nd == nullptr)
                Msg::error("taxon '" + taxonNames[p][t] + "' not found in the tree");
            tipStateByOffset[p][nd->getOffset()] = patternState[p][t];
        }
    mappedRoot = tree->getRoot();
}

double SequenceLikelihood::computeLnL(Tree* tree,
                                      const std::vector<std::vector<double>>& branchRates,
                                      const std::vector<std::vector<double>>& exchangeability,
                                      const std::vector<std::vector<double>>& frequency,
                                      const std::vector<double>& alpha,
                                      const std::vector<double>& proportionInvariant,
                                      const std::vector<std::vector<BranchMGF>>& branchMGF){
    mapTaxaToNodes(tree);
    tree->ensureBackboneCache();
    int numNodes = tree->getNumNodes();
    if(cacheValid == false){
        conP.assign(numPartitions, std::vector<std::vector<double>>(numNodes));
        cumScale.assign(numPartitions, std::vector<std::vector<double>>(numNodes));
        lastBl.assign(numPartitions, std::vector<double>(numNodes, -1.0));
        lastMGF.assign(numPartitions, std::vector<BranchMGF>(numNodes, BranchMGF{-1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
        lastExch.assign(numPartitions, std::vector<double>());
        lastFreq.assign(numPartitions, std::vector<double>());
        lastAlpha.assign(numPartitions, -1.0);
        lastPinv.assign(numPartitions, -1.0);
    }

    double lnL = 0.0;
    const double ninf = -std::numeric_limits<double>::infinity();
    if(numPartitions <= 1){
        lnL = computePartitionLnL(0, tree, branchRates, exchangeability, frequency, alpha, proportionInvariant, branchMGF, true);
    }else{
        std::vector<double> partLnL(numPartitions, 0.0);
        ThreadPool::current().parallelFor(OP_CTMC, numPartitions, [&](int q0, int q1){
            for(int q = q0; q < q1; q++)
                partLnL[q] = computePartitionLnL(q, tree, branchRates, exchangeability, frequency, alpha, proportionInvariant, branchMGF, false);
        });
        for(double v : partLnL){
            if(v == ninf){ lnL = ninf; break; }
            lnL += v;
        }
    }
    cacheValid = true;
    return lnL;
}

double SequenceLikelihood::computePartitionLnL(int p, Tree* tree,
                                      const std::vector<std::vector<double>>& branchRates,
                                      const std::vector<std::vector<double>>& exchangeability,
                                      const std::vector<std::vector<double>>& frequency,
                                      const std::vector<double>& alpha,
                                      const std::vector<double>& proportionInvariant,
                                      const std::vector<std::vector<BranchMGF>>& branchMGF,
                                      bool parallelPatterns){
    std::vector<Node*>& downPass = tree->getDownPassSequence();
    Node* root = tree->getRoot();
    int n = numStates;
    int numNodes = tree->getNumNodes();
    int full = (1 << n) - 1;

    std::vector<double> curBl(numNodes, 0.0);
    for(Node* node : downPass){
        if(node == root || tree->isBackboneNode(node) == false) continue;
        int off = node->getOffset();
        curBl[off] = branchRates[p][off] * (tree->getBackboneParent(node)->getTime() - node->getTime());
        if(curBl[off] < 0.0){
            return -std::numeric_limits<double>::infinity();
        }
    }
    {
        bool substDirty = (cacheValid == false) || (exchangeability[p] != lastExch[p]) || (frequency[p] != lastFreq[p]) || (alpha[p] != lastAlpha[p]) || (proportionInvariant[p] != lastPinv[p]);
        if((cacheValid == false) || (exchangeability[p] != lastExch[p]) || (frequency[p] != lastFreq[p]))
            rateModel[p].setParameters(exchangeability[p], frequency[p]);

        std::vector<double> cat;
        if(numCats > 1 && alpha[p] > 0.0){
            cat.resize(numCats);
            Probability::Gamma::discretization(cat, alpha[p], alpha[p], numCats, false);
        }else{
            cat.assign(1, 1.0);
        }
        double pinv = proportionInvariant[p];
        if(pinv > 0.0)
            for(double& r : cat) r /= (1.0 - pinv);
        int K = (int)cat.size();
        int npat = (int)patternWeight[p].size();

        std::vector<char> dirty(numNodes, 0);
        std::vector<Node*> dirtyNodes;
        std::vector<std::vector<Node*> > dirtyChildren;
        std::vector<std::vector<double> > Pcache(numNodes);
        for(Node* node : downPass){
            if(node != root && tree->isBackboneNode(node) == false) continue;
            int off = node->getOffset();
            if(node->getIsTip()){
                if(conP[p][off].empty()){
                    const std::vector<int>& st = tipStateByOffset[p][off];
                    conP[p][off].assign(K * npat * n, 0.0);
                    cumScale[p][off].assign(npat, 0.0);
                    for(int k = 0; k < K; k++)
                        for(int h = 0; h < npat; h++){
                            int m = st.empty() ? full : st[h];
                            for(int a = 0; a < n; a++)
                                conP[p][off][(k * npat + h) * n + a] = ((m >> a) & 1) ? 1.0 : 0.0;
                        }
                }
                continue;
            }
            const std::vector<Node*>& children = tree->getBackboneChildren(node);
            bool nd = substDirty;
            for(Node* c : children){
                int coff = c->getOffset();
                if(dirty[coff] || curBl[coff] != lastBl[p][coff] || branchMGF[p][coff] != lastMGF[p][coff]) nd = true;
            }
            if(nd == false) continue;
            dirty[off] = 1;
            dirtyNodes.push_back(node);
            dirtyChildren.push_back(children);
            conP[p][off].assign(K * npat * n, 1.0);
            cumScale[p][off].assign(npat, 0.0);
            for(Node* c : children){
                int coff = c->getOffset();
                Pcache[coff].assign(K * n * n, 0.0);
                for(int k = 0; k < K; k++)
                    rateModel[p].transitionProbabilities(curBl[coff], cat[k], branchMGF[p][coff], &Pcache[coff][k * n * n]);
            }
        }

        int croff = root->getOffset();
        const std::vector<double>& rootConP = conP[p][croff];
        std::vector<double> siteLn(npat);
        const double ninf = -std::numeric_limits<double>::infinity();
        auto patBody = [&](int h0, int h1){
            for(size_t di = 0; di < dirtyNodes.size(); di++){
                Node* node = dirtyNodes[di];
                int noff = node->getOffset();
                std::vector<double>& cp = conP[p][noff];
                std::vector<double>& cs = cumScale[p][noff];
                const std::vector<Node*>& children = dirtyChildren[di];
                for(Node* c : children){
                    int coff = c->getOffset();
                    const std::vector<double>& ccp = conP[p][coff];
                    const double* Pm = Pcache[coff].data();
                    for(int k = 0; k < K; k++){
                        const double* Pk = Pm + k * n * n;
                        for(int h = h0; h < h1; h++)
                            for(int a = 0; a < n; a++){
                                double sum = 0.0;
                                for(int b = 0; b < n; b++)
                                    sum += Pk[a * n + b] * ccp[(k * npat + h) * n + b];
                                cp[(k * npat + h) * n + a] *= sum;
                            }
                    }
                }
                static const bool noScale = (getenv("FBD_NOSCALE") != nullptr);
                for(int h = h0; h < h1; h++){
                    double mx = 0.0;
                    for(int k = 0; k < K; k++)
                        for(int a = 0; a < n; a++){
                            double v = cp[(k * npat + h) * n + a];
                            if(v > mx) mx = v;
                        }
                    double s = 0.0;
                    if(mx > 0.0 && noScale == false){
                        double inv = 1.0 / mx;
                        for(int k = 0; k < K; k++)
                            for(int a = 0; a < n; a++)
                                cp[(k * npat + h) * n + a] *= inv;
                        s = std::log(mx);
                    }
                    for(Node* c : children)
                        s += cumScale[p][c->getOffset()][h];
                    cs[h] = s;
                }
            }
            const std::vector<double>& csRoot = cumScale[p][croff];
            for(int h = h0; h < h1; h++){
                double gammaLk = 0.0;
                for(int k = 0; k < K; k++){
                    double lk = 0.0;
                    for(int a = 0; a < n; a++)
                        lk += frequency[p][a] * rootConP[(k * npat + h) * n + a];
                    gammaLk += lk / K;
                }
                double pinvLk = 0.0;
                if(pinv > 0.0){
                    int mask = constantState[p][h];
                    for(int a = 0; a < n; a++)
                        if(mask & (1 << a)) pinvLk += frequency[p][a];
                }
                double lnInv = (pinv > 0.0 && pinvLk > 0.0) ? std::log(pinv) + std::log(pinvLk) : ninf;
                double lnVar = (pinv < 1.0 && gammaLk > 0.0) ? std::log(1.0 - pinv) + std::log(gammaLk) + csRoot[h] : ninf;
                if(lnInv == ninf && lnVar == ninf){
                    siteLn[h] = ninf;
                    continue;
                }
                double mx = (lnInv > lnVar) ? lnInv : lnVar;
                siteLn[h] = patternWeight[p][h] * (mx + std::log(std::exp(lnInv - mx) + std::exp(lnVar - mx)));
            }
        };
        if(parallelPatterns)
            ThreadPool::current().parallelFor(OP_CTMC, npat, patBody);
        else
            patBody(0, npat);
        lastBl[p] = curBl;
        lastMGF[p] = branchMGF[p];
        lastExch[p] = exchangeability[p];
        lastFreq[p] = frequency[p];
        lastAlpha[p] = alpha[p];
        lastPinv[p] = proportionInvariant[p];

        double lnL = 0.0;
        for(int h = 0; h < npat; h++){
            if(siteLn[h] == -std::numeric_limits<double>::infinity()){
                return -std::numeric_limits<double>::infinity();
            }
            lnL += siteLn[h];
        }
        return lnL;
    }
}
