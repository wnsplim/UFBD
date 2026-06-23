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
    mappedCrown(nullptr)
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
    std::vector<double> P(n * n);
    std::vector<std::vector<double>> conP(numNodes);
    double lnL = 0.0;

    for(int p = 0; p < numPartitions; p++){
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
        std::vector<double> siteG(npat, 0.0);

        for(int k = 0; k < K; k++){
            for(Node* node : downPass){
                int off = node->getOffset();
                conP[off].assign(npat * n, 0.0);
                if(node->getIsTip()){
                    const std::vector<int>& st = tipStateByOffset[p][off];
                    int full = (1 << n) - 1;
                    for(int h = 0; h < npat; h++){
                        int m = st.empty() ? full : st[h];
                        for(int a = 0; a < n; a++)
                            conP[off][h * n + a] = ((m >> a) & 1) ? 1.0 : 0.0;
                    }
                }else{
                    for(int h = 0; h < npat * n; h++)
                        conP[off][h] = 1.0;
                    for(Node* c : node->getDescendants()){
                        double bl = branchRates[p][c->getOffset()] * (node->getTime() - c->getTime());
                        if(bl < 0.0)
                            return -std::numeric_limits<double>::infinity();
                        rateModel[p].transitionProbabilities(bl * cat[k], &P[0]);
                        int coff = c->getOffset();
                        for(int h = 0; h < npat; h++)
                            for(int a = 0; a < n; a++){
                                double sum = 0.0;
                                for(int b = 0; b < n; b++)
                                    sum += P[a * n + b] * conP[coff][h * n + b];
                                conP[off][h * n + a] *= sum;
                            }
                    }
                }
            }
            int croff = crown->getOffset();
            for(int h = 0; h < npat; h++){
                double lk = 0.0;
                for(int a = 0; a < n; a++)
                    lk += frequency[p][a] * conP[croff][h * n + a];
                siteG[h] += lk / K;
            }
        }

        for(int h = 0; h < npat; h++){
            double pinvLk = 0.0;
            if(pinv > 0.0){
                int mask = constantState[p][h];
                for(int a = 0; a < n; a++)
                    if(mask & (1 << a)) pinvLk += frequency[p][a];
            }
            double site = pinv * pinvLk + (1.0 - pinv) * siteG[h];
            if(site <= 0.0)
                return -std::numeric_limits<double>::infinity();
            lnL += patternWeight[p][h] * std::log(site);
        }
    }
    return lnL;
}
