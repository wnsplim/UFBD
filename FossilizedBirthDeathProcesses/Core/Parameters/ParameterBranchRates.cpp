#include <cmath>

#include "Node.hpp"
#include "ParameterBranchRates.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"

ParameterBranchRates::ParameterBranchRates(double prob, PhylogeneticModel* m, Tree* t, int L, ClockModel cm, const double* rg, const double* s2) : Parameter(prob, m, "branchRates"){
    tree = t;
    numLoci = L;
    clockModel = cm;
    numNodes = t->getNumNodes();
    lastMove = -1;
    lastLocus = -1;
    lastNode = -1;
    thetaPara[0] = 2.0;
    thetaPara[1] = 2.0;
    thetaPara[2] = 1.0;
    for(int i = 0; i < 3; i++){
        rgenePara[i] = rg[i];
        sigma2Para[i] = s2[i];
    }
    for(int i = 0; i < 4; i++){
        step[i] = 1.0;
        acc[i] = 0;
        rej[i] = 0;
    }

    Node* root = t->getRoot();
    for(Node* n : t->getDownPassSequence())
        if(n != root && n->getIsFossil() == false)
            branchNodes.push_back(n->getOffset());

    double muInit = 1.0;
    double s2Init = sigma2Para[0] / sigma2Para[1];
    double thInit = thetaPara[0] / thetaPara[1];
    for(int s = 0; s < 2; s++){
        mu[s].assign(numLoci, muInit);
        sigma2[s].assign(numLoci, s2Init);
        theta[s].assign(numLoci, thInit);
        rate[s].assign(numLoci, std::vector<double>(numNodes, 1.0));
    }
}

double ParameterBranchRates::getAcceptanceRatio(void){
    int a = 0, r = 0;
    for(int i = 0; i < 4; i++){
        a += acc[i];
        r += rej[i];
    }
    return ((double)a) / ((double)a + (double)r);
}

double ParameterBranchRates::gammaDirichletLnP(const std::vector<double>& v, const double* para){
    double a = para[0], b = para[1], conc = para[2];
    int L = (int)v.size();
    double sum = 0.0, lnprod = 0.0;
    for(double x : v){
        if(x <= 0.0)
            return -INFINITY;
        sum += x;
        lnprod += std::log(x);
    }
    double lnp = (a - conc * L) * std::log(sum) - (b / L) * sum + (conc - 1.0) * lnprod;
    lnp += a * std::log(b / L) - Probability::Helper::lnGamma(a) + Probability::Helper::lnGamma(conc * L) - L * Probability::Helper::lnGamma(conc);
    return lnp;
}

double ParameterBranchRates::gammaLnPdf(double a, double b, double x){
    if(a <= 0.0 || b <= 0.0 || x <= 0.0)
        return -INFINITY;
    return a * std::log(b) - Probability::Helper::lnGamma(a) + (a - 1.0) * std::log(x) - b * x;
}

double ParameterBranchRates::lognormalLnP(double r, double s2, double m){
    if(r <= 0.0 || s2 <= 0.0)
        return -INFINITY;
    double logr = std::log(r);
    return Probability::Normal::lnPdf(std::log(m) - 0.5 * s2, s2, logr) - logr;
}

double ParameterBranchRates::whiteNoiseLnP(double r, double s2, double t, double m){
    if(r <= 0.0 || s2 <= 0.0 || t <= 0.0)
        return -INFINITY;
    double alpha = m * m * t / s2;
    double beta = m * t / s2;
    return alpha * std::log(beta) - Probability::Helper::lnGamma(alpha) + (alpha - 1.0) * std::log(r) - beta * r;
}

double ParameterBranchRates::gbmLnP(void){
    double lnp = 0.0;
    Node* root = tree->getRoot();
    for(int p = 0; p < numLoci; p++){
        double s2 = sigma2[0][p];
        if(s2 <= 0.0)
            return -INFINITY;
        lnp -= 0.5 * std::log(2.0 * PI) * (double)branchNodes.size();
        for(Node* inode : tree->getDownPassSequence()){
            std::vector<Node*>& sons = inode->getDescendants();
            if(sons.size() != 2)
                continue;
            double t = inode->getTime();
            double tA = (inode == root) ? 0.0 : (inode->getAncestor()->getTime() - t) / 2.0;
            double t1 = (t - sons[0]->getTime()) / 2.0;
            double t2 = (t - sons[1]->getTime()) / 2.0;
            double detT = t1 * t2 + tA * (t1 + t2);
            double Tinv0 = (tA + t2) / detT;
            double Tinv1 = -tA / detT;
            double Tinv3 = (tA + t1) / detT;
            double rA = (inode == root) ? mu[0][p] : rate[0][p][inode->getOffset()];
            double r1 = rate[0][p][sons[0]->getOffset()];
            double r2 = rate[0][p][sons[1]->getOffset()];
            double y1 = std::log(r1 / rA) + (tA + t1) * s2 / 2.0;
            double y2 = std::log(r2 / rA) + (tA + t2) * s2 / 2.0;
            double quadForm = y1 * y1 * Tinv0 + 2.0 * y1 * y2 * Tinv1 + y2 * y2 * Tinv3;
            lnp -= quadForm / (2.0 * s2) + 0.5 * std::log(detT * s2 * s2) + std::log(r1 * r2);
        }
    }
    return lnp;
}

double ParameterBranchRates::cirLnP(void){
    double lnp = 0.0;
    for(int p = 0; p < numLoci; p++){
        double s2 = sigma2[0][p];
        double th = theta[0][p];
        if(s2 <= 0.0 || s2 >= 1.0 || th <= 0.0)
            return -INFINITY;
        for(int b : branchNodes){
            Node* n = tree->getNodeByOffset(b);
            double L = n->getAncestor()->getTime() - n->getTime();
            if(L <= 0.0)
                return -INFINITY;
            double rhoUp = rate[0][p][n->getAncestor()->getOffset()];
            double decay = std::exp(-th * L);
            double mean = 1.0 + (rhoUp - 1.0) * decay;
            double var = s2 * ((1.0 - decay) * (1.0 - decay) + 2.0 * rhoUp * (decay - decay * decay));
            if(var <= 0.0)
                return -INFINITY;
            double alpha = mean * mean / var;
            double beta = mean / var;
            lnp += gammaLnPdf(alpha, beta, rate[0][p][b]);
        }
    }
    return lnp;
}

double ParameterBranchRates::besselIRatio(double nu, double x){
    double epsilon = 1e-300;
    double ratio = epsilon;
    double cNumer = ratio;
    double dDenom = 0;
    for(int j = 1; j < 1000000; j++){
        double coef = 2 * (nu + j - 1) / x;
        dDenom = coef + dDenom;
        if(dDenom == 0)
            dDenom = epsilon;
        cNumer = coef + 1 / cNumer;
        if(cNumer == 0)
            cNumer = epsilon;
        dDenom = 1 / dDenom;
        double delta = cNumer * dDenom;
        ratio *= delta;
        if(std::fabs(delta - 1) < 1e-12)
            break;
    }
    return ratio;
}

double ParameterBranchRates::getMeanTau(double rho, double rhoUp, double t, double sigma, double theta){
    if(t * theta < 0.0001)
        return (rho + rhoUp) * t / 2;
    double decay = std::exp(-theta * t);
    double rootTerm = std::sqrt(rho * rhoUp * decay);
    double term1 = sigma / theta * (-1 / theta + t * decay / (1 - decay) + theta * t / sigma + (rho + rhoUp) / sigma);
    double term2 = 2 * (rho + rhoUp) / theta / ((1 - decay) * (1 - decay)) * (decay - decay * decay - t * theta * decay) - t + sigma * t / 2 / theta;
    double besselArg = 4 * theta / sigma / (1 - decay) * rootTerm;
    double besselOrder = 2 * theta / sigma;
    double besselRatio = besselIRatio(besselOrder, besselArg);
    double besselFactor = besselRatio + 0.25 * (2 * theta / sigma - 1) / theta * sigma * (1 - decay) / rootTerm;
    double derivFactor = -4 / theta / (1 - decay) * rootTerm + 4 / ((1 - decay) * (1 - decay)) * rootTerm * t * decay + 2 / (1 - decay) / rootTerm * rho * rhoUp * t * decay;
    double term3 = besselFactor * derivFactor;
    return term1 + term2 + term3;
}

double ParameterBranchRates::lnProbability(void){
    double lnp = gammaDirichletLnP(mu[0], rgenePara) + gammaDirichletLnP(sigma2[0], sigma2Para);
    if(clockModel == ClockModel::GBM)
        return lnp + gbmLnP();
    if(clockModel == ClockModel::CIR)
        return lnp + gammaDirichletLnP(theta[0], thetaPara) + cirLnP();
    for(int p = 0; p < numLoci; p++){
        for(int b : branchNodes){
            if(clockModel == ClockModel::WN){
                Node* n = tree->getNodeByOffset(b);
                lnp += whiteNoiseLnP(rate[0][p][b], sigma2[0][p], n->getAncestor()->getTime() - n->getTime(), mu[0][p]);
            }else{
                lnp += lognormalLnP(rate[0][p][b], sigma2[0][p], mu[0][p]);
            }
        }
    }
    return lnp;
}

std::vector<std::vector<double>> ParameterBranchRates::getAbsoluteRates(void){
    std::vector<std::vector<double>> a(numLoci, std::vector<double>(numNodes, 0.0));
    Node* root = tree->getRoot();
    for(int p = 0; p < numLoci; p++){
        for(int b = 0; b < numNodes; b++){
            Node* n = tree->getNodeByOffset(b);
            if(clockModel == ClockModel::CIR && n != root){
                double L = n->getAncestor()->getTime() - n->getTime();
                double sigmaPB = 2.0 * theta[0][p] * sigma2[0][p];
                a[p][b] = mu[0][p] * getMeanTau(rate[0][p][b], rate[0][p][n->getAncestor()->getOffset()], L, sigmaPB, theta[0][p]) / L;
            }else if(clockModel == ClockModel::UCLN || clockModel == ClockModel::WN || clockModel == ClockModel::GBM){
                a[p][b] = rate[0][p][b];
            }else{
                a[p][b] = mu[0][p] * rate[0][p][b];
            }
        }
    }
    return a;
}

double ParameterBranchRates::bactrianMultiplier(int mt){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double ar = 0.0;
    for(bool b : recentAR[mt])
        if(b)
            ar++;
    if(recentAR[mt].empty() == false)
        ar /= recentAR[mt].size();
    int total = acc[mt] + rej[mt];
    if(total > 0 && total % 100 == 0){
        double gain = 1.0 / std::sqrt((double)(total / 100));
        step[mt] *= std::exp(gain * (ar - 0.3));
    }
    double m = 0.95;
    double s = std::sqrt(1.0 - m * m);
    double delta = m + Probability::Normal::rv(&rng) * s;
    if(Probability::Uniform::rv(&rng, 0.0, 1.0) < 0.5)
        delta = -delta;
    return std::exp(step[mt] * delta);
}

double ParameterBranchRates::scaleLocusRate(int p){
    double c = bactrianMultiplier(0);
    mu[0][p] *= c;
    return std::log(c);
}

double ParameterBranchRates::scaleLocusSigma2(int p){
    double c = bactrianMultiplier(1);
    sigma2[0][p] *= c;
    return std::log(c);
}

double ParameterBranchRates::scaleLocusTheta(int p){
    double c = bactrianMultiplier(3);
    theta[0][p] *= c;
    return std::log(c);
}

double ParameterBranchRates::scaleBranchRate(int p, int b){
    double c = bactrianMultiplier(2);
    rate[0][p][b] *= c;
    return std::log(c);
}

double ParameterBranchRates::globalRateBranchRatesScale(int p){
    double sf = bactrianMultiplier(0);
    mu[0][p] *= sf;
    for(int b : branchNodes)
        rate[0][p][b] *= sf;
    return (1.0 + (double)branchNodes.size()) * std::log(sf);
}

double ParameterBranchRates::vectorScale(int p){
    double sf = bactrianMultiplier(0);
    for(int b : branchNodes)
        rate[0][p][b] *= sf;
    return (double)branchNodes.size() * std::log(sf);
}

double ParameterBranchRates::update(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    lastLocus = (int)(rng.uniformRv() * numLoci);
    double u = rng.uniformRv();
    double branchW = (clockModel == ClockModel::CIR) ? 0.7 : 0.8;
    if(u < branchW){
        lastMove = 2;
        lastNode = branchNodes[(int)(rng.uniformRv() * branchNodes.size())];
        return scaleBranchRate(lastLocus, lastNode);
    }
    if(u < branchW + 0.1){
        lastMove = 0;
        return scaleLocusRate(lastLocus);
    }
    if(u < branchW + 0.2){
        lastMove = 1;
        return scaleLocusSigma2(lastLocus);
    }
    lastMove = 3;
    return scaleLocusTheta(lastLocus);
}

void ParameterBranchRates::updateForAcceptance(void){
    if(lastMove == 4){
        for(int k = 0; k < (int)cdNodes.size(); k++)
            for(int p = 0; p < numLoci; p++)
                rate[1][p][cdNodes[k]] = rate[0][p][cdNodes[k]];
        return;
    }
    acc[lastMove]++;
    recentAR[lastMove].push_back(true);
    if(recentAR[lastMove].size() > 1000)
        recentAR[lastMove].pop_front();
    if(lastMove == 0)
        mu[1][lastLocus] = mu[0][lastLocus];
    else if(lastMove == 1)
        sigma2[1][lastLocus] = sigma2[0][lastLocus];
    else if(lastMove == 3)
        theta[1][lastLocus] = theta[0][lastLocus];
    else
        rate[1][lastLocus][lastNode] = rate[0][lastLocus][lastNode];
}

void ParameterBranchRates::updateForRejection(void){
    if(lastMove == 4){
        for(int k = 0; k < (int)cdNodes.size(); k++)
            for(int p = 0; p < numLoci; p++)
                rate[0][p][cdNodes[k]] = rate[1][p][cdNodes[k]];
        return;
    }
    rej[lastMove]++;
    recentAR[lastMove].push_back(false);
    if(recentAR[lastMove].size() > 1000)
        recentAR[lastMove].pop_front();
    if(lastMove == 0)
        mu[0][lastLocus] = mu[1][lastLocus];
    else if(lastMove == 1)
        sigma2[0][lastLocus] = sigma2[1][lastLocus];
    else if(lastMove == 3)
        theta[0][lastLocus] = theta[1][lastLocus];
    else
        rate[0][lastLocus][lastNode] = rate[1][lastLocus][lastNode];
}

double ParameterBranchRates::constantDistanceMove(void){
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    std::vector<Node*> internals;
    for(Node* n : tree->getDownPassSequence())
        if(n != tree->getRoot() && n->getIsTip() == false)
            internals.push_back(n);
    Node* node = internals[(int)(rng.uniformRv() * internals.size())];
    Node* parent = node->getAncestor();
    std::vector<Node*> children;
    for(Node* c : node->getNeighbors())
        if(c != parent)
            children.push_back(c);
    double parentAge = parent->getTime();
    double myAge = node->getTime();
    double maxChild = 0.0;
    for(Node* c : children)
        if(c->getTime() > maxChild)
            maxChild = c->getTime();
    double newAge = maxChild + rng.uniformRv() * (parentAge - maxChild);
    node->setTime(newAge);
    cdNodes.clear();
    cdNodes.push_back(node->getOffset());
    double lnNum = std::log(parentAge - myAge);
    double lnDen = std::log(parentAge - newAge);
    for(int p = 0; p < numLoci; p++)
        rate[0][p][node->getOffset()] *= (parentAge - myAge) / (parentAge - newAge);
    for(Node* c : children){
        double prevC = myAge - c->getTime();
        double newC = newAge - c->getTime();
        cdNodes.push_back(c->getOffset());
        lnNum += std::log(prevC);
        lnDen += std::log(newC);
        for(int p = 0; p < numLoci; p++)
            rate[0][p][c->getOffset()] *= prevC / newC;
    }
    lastMove = 4;
    tree->setLastUpdateWasScale(false);
    return numLoci * (lnNum - lnDen);
}

void ParameterBranchRates::print(void){
}
