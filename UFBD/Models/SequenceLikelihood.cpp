#include "SequenceLikelihood.hpp"
#include "Tree.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "Msg.hpp"
#include "ThreadPool.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

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
    int full = (1 << numStates) - 1;
    tipStateByOffset.assign(numPartitions, std::vector<std::vector<int>>(numNodes));
    tipMissing.assign(numPartitions, std::vector<char>(numNodes, 0));
    for(int p = 0; p < numPartitions; p++)
        for(int t = 0; t < (int)taxonNames[p].size(); t++){
            Node* nd = tree->getTaxonNode(taxonNames[p][t]);
            if(nd == nullptr)
                Msg::error("taxon '" + taxonNames[p][t] + "' not found in the tree");
            int off = nd->getOffset();
            tipStateByOffset[p][off] = patternState[p][t];
            bool allGap = true;
            for(int m : patternState[p][t])
                if(m != full){ allGap = false; break; }
            tipMissing[p][off] = allGap ? 1 : 0;
        }
    buildRepeats(tree);
    mappedRoot = tree->getRoot();
}

void SequenceLikelihood::buildRepeats(Tree* tree){
    tree->ensureBackboneCache();
    std::vector<Node*>& downPass = tree->getDownPassSequence();
    Node* root = tree->getRoot();
    int numNodes = tree->getNumNodes();
    int full = (1 << numStates) - 1;
    clsId.assign(numPartitions, std::vector<std::vector<int>>(numNodes));
    clsRep.assign(numPartitions, std::vector<std::vector<int>>(numNodes));
    for(int p = 0; p < numPartitions; p++){
        int npat = (int)patternWeight[p].size();
        for(Node* node : downPass){
            if(node != root && tree->isBackboneNode(node) == false) continue;
            int off = node->getOffset();
            std::vector<int>& cid = clsId[p][off];
            std::vector<int>& rep = clsRep[p][off];
            cid.assign(npat, 0);
            rep.clear();
            if(node->getIsTip()){
                const std::vector<int>& st = tipStateByOffset[p][off];
                std::map<int, int> intern;
                for(int h = 0; h < npat; h++){
                    int m = st.empty() ? full : st[h];
                    auto it = intern.find(m);
                    if(it == intern.end()){ cid[h] = (int)intern.size(); intern.emplace(m, cid[h]); rep.push_back(h); }
                    else cid[h] = it->second;
                }
                continue;
            }
            std::vector<int> kids;
            for(Node* c : tree->getBackboneChildren(node)){
                int coff = c->getOffset();
                if(c->getIsTip() && tipMissing[p][coff]) continue;
                kids.push_back(coff);
            }
            std::map<std::vector<int>, int> intern;
            std::vector<int> key(kids.size());
            for(int h = 0; h < npat; h++){
                for(size_t j = 0; j < kids.size(); j++) key[j] = clsId[p][kids[j]][h];
                auto it = intern.find(key);
                if(it == intern.end()){ cid[h] = (int)intern.size(); intern.emplace(key, cid[h]); rep.push_back(h); }
                else cid[h] = it->second;
            }
        }
    }
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
                    const std::vector<int>& rp = clsRep[p][off];
                    int nc = (int)rp.size();
                    conP[p][off].assign(K * nc * n, 0.0);
                    cumScale[p][off].assign(nc, 0.0);
                    for(int k = 0; k < K; k++)
                        for(int c = 0; c < nc; c++){
                            int m = st.empty() ? full : st[rp[c]];
                            for(int a = 0; a < n; a++)
                                conP[p][off][(k * nc + c) * n + a] = ((m >> a) & 1) ? 1.0 : 0.0;
                        }
                }
                continue;
            }
            const std::vector<Node*>& children = tree->getBackboneChildren(node);
            bool nd = substDirty;
            for(Node* c : children){
                int coff = c->getOffset();
                if(c->getIsTip() && tipMissing[p][coff]) continue;
                if(dirty[coff] || curBl[coff] != lastBl[p][coff] || branchMGF[p][coff] != lastMGF[p][coff]) nd = true;
            }
            if(nd == false) continue;
            dirty[off] = 1;
            dirtyNodes.push_back(node);
            dirtyChildren.push_back(children);
            int nc = (int)clsRep[p][off].size();
            conP[p][off].assign(K * nc * n, 1.0);
            cumScale[p][off].assign(nc, 0.0);
            for(Node* c : children){
                int coff = c->getOffset();
                if(c->getIsTip() && tipMissing[p][coff]) continue;
                Pcache[coff].assign(K * n * n, 0.0);
                for(int k = 0; k < K; k++)
                    rateModel[p].transitionProbabilities(curBl[coff], cat[k], branchMGF[p][coff], &Pcache[coff][k * n * n]);
            }
        }

        int croff = root->getOffset();
        const std::vector<double>& rootConP = conP[p][croff];
        int ncRoot = (int)clsRep[p][croff].size();
        std::vector<double> siteLn(npat);
        const double ninf = -std::numeric_limits<double>::infinity();
        auto node4 = [&](size_t di, int c0, int c1){
            Node* node = dirtyNodes[di];
            int noff = node->getOffset();
            int nc = (int)clsRep[p][noff].size();
            const std::vector<int>& rp = clsRep[p][noff];
            double* cp = conP[p][noff].data();
            std::vector<double>& cs = cumScale[p][noff];
            const std::vector<Node*>& children = dirtyChildren[di];
            for(Node* c : children){
                int coff = c->getOffset();
                if(c->getIsTip() && tipMissing[p][coff]) continue;
                const double* Pm = Pcache[coff].data();
                if(c->getIsTip()){
                    const std::vector<int>& st = tipStateByOffset[p][coff];
                    for(int k = 0; k < K; k++){
                        const double* Pk = Pm + k * 16;
                        for(int cc = c0; cc < c1; cc++){
                            int m = st.empty() ? 15 : st[rp[cc]];
                            double* o = cp + (k * nc + cc) * 4;
                            if((m & (m - 1)) == 0){
                                int s = (m == 1) ? 0 : (m == 2) ? 1 : (m == 4) ? 2 : 3;
                                o[0] *= Pk[s]; o[1] *= Pk[4 + s]; o[2] *= Pk[8 + s]; o[3] *= Pk[12 + s];
                            }else{
                                double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
                                if(m & 1){ s0 += Pk[0];  s1 += Pk[4];  s2 += Pk[8];  s3 += Pk[12]; }
                                if(m & 2){ s0 += Pk[1];  s1 += Pk[5];  s2 += Pk[9];  s3 += Pk[13]; }
                                if(m & 4){ s0 += Pk[2];  s1 += Pk[6];  s2 += Pk[10]; s3 += Pk[14]; }
                                if(m & 8){ s0 += Pk[3];  s1 += Pk[7];  s2 += Pk[11]; s3 += Pk[15]; }
                                o[0] *= s0; o[1] *= s1; o[2] *= s2; o[3] *= s3;
                            }
                        }
                    }
                }else{
                    const double* ccp = conP[p][coff].data();
                    int ncc = (int)clsRep[p][coff].size();
                    const std::vector<int>& cidc = clsId[p][coff];
                    for(int k = 0; k < K; k++){
                        const double* Pk = Pm + k * 16;
                        for(int cc = c0; cc < c1; cc++){
                            const double* in = ccp + (k * ncc + cidc[rp[cc]]) * 4;
                            double b0 = in[0], b1 = in[1], b2 = in[2], b3 = in[3];
                            double* o = cp + (k * nc + cc) * 4;
                            o[0] *= Pk[0]  * b0 + Pk[1]  * b1 + Pk[2]  * b2 + Pk[3]  * b3;
                            o[1] *= Pk[4]  * b0 + Pk[5]  * b1 + Pk[6]  * b2 + Pk[7]  * b3;
                            o[2] *= Pk[8]  * b0 + Pk[9]  * b1 + Pk[10] * b2 + Pk[11] * b3;
                            o[3] *= Pk[12] * b0 + Pk[13] * b1 + Pk[14] * b2 + Pk[15] * b3;
                        }
                    }
                }
            }
            for(int cc = c0; cc < c1; cc++){
                double mx = 0.0;
                for(int k = 0; k < K; k++){
                    const double* o = cp + (k * nc + cc) * 4;
                    if(o[0] > mx) mx = o[0]; if(o[1] > mx) mx = o[1];
                    if(o[2] > mx) mx = o[2]; if(o[3] > mx) mx = o[3];
                }
                double s = 0.0;
                if(mx > 0.0){
                    double inv = 1.0 / mx;
                    for(int k = 0; k < K; k++){
                        double* o = cp + (k * nc + cc) * 4;
                        o[0] *= inv; o[1] *= inv; o[2] *= inv; o[3] *= inv;
                    }
                    s = std::log(mx);
                }
                for(Node* c : children){
                    int coff = c->getOffset();
                    s += cumScale[p][coff][clsId[p][coff][rp[cc]]];
                }
                cs[cc] = s;
            }
        };
        auto nodeGeneric = [&](size_t di, int c0, int c1){
            Node* node = dirtyNodes[di];
            int noff = node->getOffset();
            int nc = (int)clsRep[p][noff].size();
            const std::vector<int>& rp = clsRep[p][noff];
            std::vector<double>& cp = conP[p][noff];
            std::vector<double>& cs = cumScale[p][noff];
            const std::vector<Node*>& children = dirtyChildren[di];
            for(Node* c : children){
                int coff = c->getOffset();
                if(c->getIsTip() && tipMissing[p][coff]) continue;
                const std::vector<double>& ccp = conP[p][coff];
                int ncc = (int)clsRep[p][coff].size();
                const std::vector<int>& cidc = clsId[p][coff];
                const double* Pm = Pcache[coff].data();
                for(int k = 0; k < K; k++){
                    const double* Pk = Pm + k * n * n;
                    for(int cc = c0; cc < c1; cc++){
                        int chc = cidc[rp[cc]];
                        for(int a = 0; a < n; a++){
                            double sum = 0.0;
                            for(int b = 0; b < n; b++)
                                sum += Pk[a * n + b] * ccp[(k * ncc + chc) * n + b];
                            cp[(k * nc + cc) * n + a] *= sum;
                        }
                    }
                }
            }
            for(int cc = c0; cc < c1; cc++){
                double mx = 0.0;
                for(int k = 0; k < K; k++)
                    for(int a = 0; a < n; a++){
                        double v = cp[(k * nc + cc) * n + a];
                        if(v > mx) mx = v;
                    }
                double s = 0.0;
                if(mx > 0.0){
                    double inv = 1.0 / mx;
                    for(int k = 0; k < K; k++)
                        for(int a = 0; a < n; a++)
                            cp[(k * nc + cc) * n + a] *= inv;
                    s = std::log(mx);
                }
                for(Node* c : children){
                    int coff = c->getOffset();
                    s += cumScale[p][coff][clsId[p][coff][rp[cc]]];
                }
                cs[cc] = s;
            }
        };
        auto rootSum = [&](int h0, int h1){
            const std::vector<double>& csRoot = cumScale[p][croff];
            const std::vector<int>& cidRoot = clsId[p][croff];
            const double* rootP = rootConP.data();
            for(int h = h0; h < h1; h++){
                int rc = cidRoot[h];
                double gammaLk = 0.0;
                for(int k = 0; k < K; k++){
                    double lk = 0.0;
                    for(int a = 0; a < n; a++)
                        lk += frequency[p][a] * rootP[(k * ncRoot + rc) * n + a];
                    gammaLk += lk / K;
                }
                double pinvLk = 0.0;
                if(pinv > 0.0){
                    int mask = constantState[p][h];
                    for(int a = 0; a < n; a++)
                        if(mask & (1 << a)) pinvLk += frequency[p][a];
                }
                double lnInv = (pinv > 0.0 && pinvLk > 0.0) ? std::log(pinv) + std::log(pinvLk) : ninf;
                double lnVar = (pinv < 1.0 && gammaLk > 0.0) ? std::log(1.0 - pinv) + std::log(gammaLk) + csRoot[rc] : ninf;
                if(lnInv == ninf && lnVar == ninf){ siteLn[h] = ninf; continue; }
                double mx = (lnInv > lnVar) ? lnInv : lnVar;
                siteLn[h] = patternWeight[p][h] * (mx + std::log(std::exp(lnInv - mx) + std::exp(lnVar - mx)));
            }
        };
        for(size_t di = 0; di < dirtyNodes.size(); di++){
            int nc = (int)clsRep[p][dirtyNodes[di]->getOffset()].size();
            if(parallelPatterns){
                if(n == 4) ThreadPool::current().parallelFor(OP_CTMC, nc, [&](int a, int b){ node4(di, a, b); });
                else       ThreadPool::current().parallelFor(OP_CTMC, nc, [&](int a, int b){ nodeGeneric(di, a, b); });
            }else{
                if(n == 4) node4(di, 0, nc);
                else       nodeGeneric(di, 0, nc);
            }
        }
        if(parallelPatterns)
            ThreadPool::current().parallelFor(OP_CTMC, npat, rootSum);
        else
            rootSum(0, npat);
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
