#include "SequenceLikelihood.hpp"
#include "Tree.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "Msg.hpp"

#include <cmath>
#include <limits>

SequenceLikelihood::SequenceLikelihood(int numStates, int numCats) :
    numStates(numStates),
    numCats(numCats),
    numPartitions(0),
    mappedCrown(nullptr),
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
    if(tree->getCrown() == mappedCrown)
        return;
    int numNodes = tree->getNumNodes();
    tipStateByOffset.assign(numPartitions, std::vector<std::vector<int>>(numNodes));
    for(int p = 0; p < numPartitions; p++)
        for(int t = 0; t < (int)taxonNames[p].size(); t++){
            Node* nd = tree->getTaxonNode(taxonNames[p][t]);
            if(nd == nullptr)
                Msg::error("SequenceLikelihood: alignment taxon '" + taxonNames[p][t] + "' not found in the tree");
            tipStateByOffset[p][nd->getOffset()] = patternState[p][t];
        }
    mappedCrown = tree->getCrown();
}

double SequenceLikelihood::computeLnL(Tree* tree,
                                      const std::vector<std::vector<double>>& branchRates,
                                      const std::vector<std::vector<double>>& exchangeability,
                                      const std::vector<std::vector<double>>& frequency,
                                      const std::vector<double>& alpha,
                                      const std::vector<double>& proportionInvariant){
    mapTaxaToNodes(tree);
    std::vector<Node*>& downPass = tree->getDownPassSequence();
    Node* crown = tree->getCrown();
    int n = numStates;
    int numNodes = tree->getNumNodes();
    int full = (1 << n) - 1;
    std::vector<double> P(n * n);

    if(cacheValid == false){
        conP.assign(numPartitions, std::vector<std::vector<double>>(numNodes));
        lastBl.assign(numPartitions, std::vector<double>(numNodes, -1.0));
        lastExch.assign(numPartitions, std::vector<double>());
        lastFreq.assign(numPartitions, std::vector<double>());
        lastAlpha.assign(numPartitions, -1.0);
        lastPinv.assign(numPartitions, -1.0);
    }

    double lnL = 0.0;
    for(int p = 0; p < numPartitions; p++){
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

        std::vector<double> curBl(numNodes, 0.0);
        for(Node* node : downPass){
            if(node == crown) continue;
            int off = node->getOffset();
            curBl[off] = branchRates[p][off] * (node->getAncestor()->getTime() - node->getTime());
            if(curBl[off] < 0.0)
                return -std::numeric_limits<double>::infinity();
        }

        std::vector<char> dirty(numNodes, 0);
        for(Node* node : downPass){
            int off = node->getOffset();
            if(node->getIsTip()){
                if(conP[p][off].empty()){
                    const std::vector<int>& st = tipStateByOffset[p][off];
                    conP[p][off].assign(K * npat * n, 0.0);
                    for(int k = 0; k < K; k++)
                        for(int h = 0; h < npat; h++){
                            int m = st.empty() ? full : st[h];
                            for(int a = 0; a < n; a++)
                                conP[p][off][(k * npat + h) * n + a] = ((m >> a) & 1) ? 1.0 : 0.0;
                        }
                }
                continue;
            }
            std::vector<Node*>& kids = node->getDescendants();
            bool nd = substDirty;
            for(Node* c : kids){
                int coff = c->getOffset();
                if(dirty[coff] || curBl[coff] != lastBl[p][coff]) nd = true;
            }
            if(nd == false) continue;
            dirty[off] = 1;
            std::vector<double>& cp = conP[p][off];
            cp.assign(K * npat * n, 1.0);
            for(Node* c : kids){
                int coff = c->getOffset();
                const std::vector<double>& ccp = conP[p][coff];
                for(int k = 0; k < K; k++){
                    rateModel[p].transitionProbabilities(curBl[coff] * cat[k], &P[0]);
                    for(int h = 0; h < npat; h++)
                        for(int a = 0; a < n; a++){
                            double sum = 0.0;
                            for(int b = 0; b < n; b++)
                                sum += P[a * n + b] * ccp[(k * npat + h) * n + b];
                            cp[(k * npat + h) * n + a] *= sum;
                        }
                }
            }
        }

        int croff = crown->getOffset();
        const std::vector<double>& root = conP[p][croff];
        for(int h = 0; h < npat; h++){
            double gammaLk = 0.0;
            for(int k = 0; k < K; k++){
                double lk = 0.0;
                for(int a = 0; a < n; a++)
                    lk += frequency[p][a] * root[(k * npat + h) * n + a];
                gammaLk += lk / K;
            }
            double pinvLk = 0.0;
            if(pinv > 0.0){
                int mask = constantState[p][h];
                for(int a = 0; a < n; a++)
                    if(mask & (1 << a)) pinvLk += frequency[p][a];
            }
            double site = pinv * pinvLk + (1.0 - pinv) * gammaLk;
            if(site <= 0.0)
                return -std::numeric_limits<double>::infinity();
            lnL += patternWeight[p][h] * std::log(site);
        }

        lastBl[p] = curBl;
        lastExch[p] = exchangeability[p];
        lastFreq[p] = frequency[p];
        lastAlpha[p] = alpha[p];
        lastPinv[p] = proportionInvariant[p];
    }
    cacheValid = true;
    return lnL;
}
