#include "PathModel.hpp"

#include "FBDTreeModel.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterDouble.hpp"
#include "ParameterOUField.hpp"
#include "ParameterTree.hpp"
#include "ParameterUnresolvedFossils.hpp"
#include "RelaxedClockTreeModel.hpp"
#include "SequenceCTMCModel.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

static std::map<std::string, PseudoPrior> fitted;

static bool fromCache(PseudoPrior& d, const std::string& key){
    std::map<std::string, PseudoPrior>::iterator hit = fitted.find(key);
    if(hit == fitted.end())
        return false;
    d.mean = hit->second.mean;
    d.K = hit->second.K;
    d.wMean = hit->second.wMean;
    d.chol = hit->second.chol;
    d.cov = hit->second.cov;
    return true;
}

static void fitMvn(PseudoPrior& d, const std::vector<Eigen::VectorXd>& rows, const std::vector<Eigen::VectorXd>& ws, const std::string& what){
    int n = (int)rows.size(), k = (int)rows[0].size(), q = ws.empty() ? 0 : (int)ws[0].size();
    if(n <= k + q)
        Msg::error(what + ": " + std::to_string(n) + " rows for " + std::to_string(k + q) + " coordinates; it needs more than " + std::to_string(k + q) + ".");
    std::vector<Eigen::VectorXd> x(n, Eigen::VectorXd(k + q));
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(k + q);
    for(int i = 0; i < n; i++){
        x[i].head(k) = rows[i];
        if(q > 0) x[i].tail(q) = ws[i];
        mean += x[i];
    }
    mean /= (double)n;
    Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(k + q, k + q);
    for(const Eigen::VectorXd& y : x)
        cov += (y - mean) * (y - mean).transpose();
    cov /= (double)(n - 1);
    d.mean = mean.head(k);
    d.wMean = mean.tail(q);
    Eigen::MatrixXd czz = cov.topLeftCorner(k, k);
    d.K = Eigen::MatrixXd::Zero(k, q);
    if(q > 0){
        Eigen::LLT<Eigen::MatrixXd> cw(cov.bottomRightCorner(q, q));
        if(cw.info() != Eigen::Success)
            Msg::error(what + ": the covariance of the " + std::to_string(q) + " log node ages over " + std::to_string(n) + " rows is not positive definite.");
        d.K = cw.solve(cov.bottomLeftCorner(q, k)).transpose();
        czz -= d.K * cov.bottomLeftCorner(q, k);
    }
    d.cov = czz;
    d.chol.compute(czz);
    if(d.chol.info() != Eigen::Success)
        Msg::error(what + ": the covariance of the " + std::to_string(k) + " coordinates over " + std::to_string(n) + " rows is not positive definite.");
}

std::vector<std::vector<double>> PathModel::readColumns(const std::string& path, const std::vector<std::string>& names, const std::string& key){
    std::ifstream in(path);
    if(in.is_open() == false)
        Msg::error("could not open " + key + ".");
    std::string line, tok;
    std::getline(in, line);
    std::vector<std::string> header;
    std::istringstream hs(line);
    while(std::getline(hs, tok, '\t'))
        header.push_back(tok);
    std::vector<int> col;
    for(const std::string& n : names){
        std::vector<std::string>::iterator it = std::find(header.begin(), header.end(), n);
        if(it == header.end())
            Msg::error(key + " has no column '" + n + "'.");
        col.push_back((int)(it - header.begin()));
    }
    std::vector<std::vector<double>> out;
    while(std::getline(in, line)){
        if(line.empty()) continue;
        std::istringstream rs(line);
        std::vector<double> r;
        while(std::getline(rs, tok, '\t')){
            try { r.push_back(std::stod(tok)); }
            catch(...) { Msg::error(key + ": '" + tok + "' is not a number."); }
        }
        if(r.size() != header.size())
            Msg::error(key + ": a row has " + std::to_string(r.size()) + " values, the header has " + std::to_string(header.size()) + " names.");
        std::vector<double> v;
        for(int c : col)
            v.push_back(r[c]);
        out.push_back(v);
    }
    return out;
}

static double column(const std::map<std::string,double>& row, const std::string& name, const std::string& path){
    std::map<std::string,double>::const_iterator it = row.find(name);
    if(it == row.end())
        Msg::error("start log " + path + " has no column '" + name + "'.");
    return it->second;
}

double PathModel::posteriorBurn(const std::string& base){
    return std::ifstream(base + "_chain0.log").good() ? 0.0 : UserSettings::userSettings().getBurninFraction();
}

static double lnMvn(const PseudoPrior& d, const Eigen::VectorXd& z, const Eigen::VectorXd& w){
    Eigen::VectorXd y = d.chol.matrixL().solve(z - d.mean - d.K * (w - d.wMean));
    return -0.5 * y.squaredNorm() - d.chol.matrixLLT().diagonal().array().log().sum() - 0.5 * (double)z.size() * std::log(2.0 * M_PI);
}

PathModel::PathModel(PhylogeneticModel* a, PhylogeneticModel* b, std::vector<Fossil>& fossils, unsigned int seed) :
    PhylogeneticModel(){

    rng.setSeed(seed);
    m[0] = a; m[1] = b;
    for(int k = 0; k < 2; k++){
        RelaxedClockTreeModel* r = dynamic_cast<RelaxedClockTreeModel*>(m[k]);
        f[k] = (r != nullptr) ? r->fbd : static_cast<FBDTreeModel*>(m[k]);
    }
    FBDTreeModel* fa = f[0];
    FBDTreeModel* fb = f[1];
    UserSettings& us = UserSettings::userSettings();
    shared[0] = us.getSharedLambda();
    shared[1] = us.getSharedMu();
    shared[2] = us.getSharedPsi();
    sharedClock = us.getSharedClock();
    sharedSubst = us.getSharedSubst();
    pathT = us.getPathT();
    parameters.push_back(fa->getParameterTree());
    for(int k = 0; k < 2; k++){
        rc[k] = dynamic_cast<RelaxedClockTreeModel*>(m[k]);
        if(rc[k] != nullptr){
            rc[k]->pathTarget = [this](){ sync(mover); return lnLikelihood() + lnPriorProbability(); };
            f[k]->pathTarget = rc[k]->pathTarget;
        }
        f[k]->pathPartner = f[1 - k];
        f[k]->pathSelfW = (k == 0) ? 1.0 - pathT : pathT;
        f[k]->pathPartnerW = (k == 0) ? pathT : 1.0 - pathT;
    }

    SkylineRates& sa = fa->a;
    SkylineRates& sb = fb->a;
    std::vector<std::pair<std::string,std::pair<ParameterOUField*,ParameterOUField*>>> pairs;
    if(shared[0] == false) pairs.push_back({"lambda", {sa.lambdaField, sb.lambdaField}});
    if(shared[1] == false) pairs.push_back({"mu", {sa.muField, sb.muField}});
    if(shared[2] == false)
        for(int ta = 0; ta < sa.numPsiTypes; ta++)
            for(int tb = 0; tb < sb.numPsiTypes; tb++)
                if(sa.psiTypeNames.empty() == sb.psiTypeNames.empty() && (sa.psiTypeNames.empty() || sa.psiTypeNames[ta] == sb.psiTypeNames[tb]))
                    pairs.push_back({sa.psiTypeNames.empty() ? std::string("psi") : "psi_" + sa.psiTypeNames[ta], {sa.psiField[ta], sb.psiField[tb]}});
    for(auto& [rate, fields] : pairs){
        ParameterOUField* fa = fields.first;
        ParameterOUField* fb = fields.second;
        if(fa == nullptr || fb == nullptr) continue;
        if(fa->isThetaFixed() != fb->isThetaFixed() || fa->isSdFixed() != fb->isSdFixed() || fa->isNuFixed() != fb->isNuFixed()
           || (fa->isThetaFixed() && fa->getTheta() != fb->getTheta()) || (fa->isSdFixed() && fa->getSdEq() != fb->getSdEq()) || (fa->isNuFixed() && fa->getNu() != fb->getNu()))
                Msg::error("[model 0] and [model 1] both put an OU prior on " + rate + "; the two models then share its hyperparameters theta, sdEq and nu, so each must be free in both models or fixed to the same value in both.");
        fa->shareHyperparameters(fb);
    }

    for(int k = 0; k < 2; k++){
        const std::vector<std::string>& other = f[1 - k]->a.psiTypeNames;
        lab[k] = (f[k]->a.psiTypeNames.size() <= 1 && other.size() > 1) ? lnTypeLabels(other, fossils) : 0.0;
    }
    std::string base[2] = {us.getPosteriorBaseA(), us.getPosteriorBaseB()};
    double burn[2];
    for(int k = 0; k < 2; k++)
        burn[k] = posteriorBurn(base[k]);
    for(int k = 0; k < 2; k++)
        fitRatePseudoPrior(f[k], pp[k], specificMask(f[k]), base[k] + ".log", burn[k], "model " + std::to_string(k));
    if(sharedClock == false)
        for(int k = 0; k < 2; k++)
            readClockPseudoPrior(k, base[k], burn[k]);
    if(sharedSubst == false)
        for(int k = 0; k < 2; k++)
            if(rc[k] != nullptr && rc[k]->ctmc != nullptr)
                readSubstPseudoPrior(k, base[k] + ".log", burn[k]);
    sync(0);
}

std::vector<std::string> PathModel::clockNames(int k){
    BranchRateModel* c = rc[k]->clock;
    std::vector<std::string> n, ln = rc[k]->getLatentNames();
    for(int g = 0; g < c->getNumClockGroups(); g++){
        std::string suf = (c->getNumClockGroups() > 1) ? ("_" + std::to_string(g)) : "";
        n.push_back("clockMean" + suf);
        n.push_back("clockSigma2" + suf);
    }
    n.insert(n.end(), ln.begin() + (long)f[k]->getLatentNames().size(), ln.end());
    return n;
}

std::vector<double> PathModel::clockValues(int k){
    BranchRateModel* c = rc[k]->clock;
    std::vector<double> v;
    for(int g = 0; g < c->getNumClockGroups(); g++){
        v.push_back(c->getClockGroupRate(g));
        v.push_back(c->getClockGroupSigma2(g));
    }
    for(int g = 0; g < c->getNumClockGroups(); g++)
        for(int i = 0; i < c->getNumBranchNodes(); i++)
            v.push_back(c->getBranchRate(g, i));
    return v;
}

void PathModel::setClockValues(int k, const std::vector<double>& v){
    BranchRateModel* c = rc[k]->clock;
    int ng = c->getNumClockGroups(), nb = c->getNumBranchNodes();
    for(int g = 0; g < ng; g++)
        c->setState(g, v[2 * g], v[2 * g + 1], std::vector<double>(v.begin() + 2 * ng + g * nb, v.begin() + 2 * ng + (g + 1) * nb));
    m[k]->invalidateLikelihoodCache();
}

double PathModel::specificPrior(int k){
    double p = lnRatePrior(f[k]);
    if(rc[k] != nullptr && sharedClock == false)
        p += rc[k]->clock->lnProbability();
    if(rc[k] != nullptr && sharedSubst == false && rc[k]->ctmc != nullptr)
        p += rc[k]->ctmc->lnPrior();
    return p;
}

std::vector<Node*> PathModel::varyingAgeNodes(int k){
    std::vector<Node*> nodes = f[k]->getAgeLogNodes();
    for(const FBDTreeModel::BackboneFossil& bf : f[k]->backboneFossils)
        if(bf.yMin < bf.yMax)
            nodes.push_back(bf.tip);
    return nodes;
}

Eigen::VectorXd PathModel::logAges(int k){
    std::vector<Node*> nodes = varyingAgeNodes(k);
    Eigen::VectorXd w((int)nodes.size());
    for(size_t i = 0; i < nodes.size(); i++)
        w((int)i) = std::log(nodes[i]->getTime());
    return w;
}

double PathModel::clockCoords(int k, const std::vector<double>& muV, const std::vector<double>& s2V, const std::vector<std::vector<double>>& r, Eigen::VectorXd& z){
    int ng = (int)muV.size();
    std::vector<double> u;
    double jac = static_cast<ParameterBranchRates*>(rc[k]->clock)->standardizedDeviations(muV, s2V, r, u);
    z.resize(2 * ng + (int)u.size());
    for(int g = 0; g < ng; g++){
        z(2 * g) = std::log(muV[g]);
        z(2 * g + 1) = std::log(s2V[g]);
        jac -= z(2 * g) + z(2 * g + 1);
    }
    for(size_t i = 0; i < u.size(); i++)
        z(2 * ng + (int)i) = u[i];
    return jac;
}

double PathModel::lnClockPseudoPrior(int k){
    const PseudoPrior& d = ppc[k];
    BranchRateModel* c = rc[k]->clock;
    int ng = c->getNumClockGroups(), nb = c->getNumBranchNodes();
    std::vector<double> muV(ng), s2V(ng);
    std::vector<std::vector<double>> r(ng, std::vector<double>(nb));
    for(int g = 0; g < ng; g++){
        muV[g] = c->getClockGroupRate(g);
        s2V[g] = c->getClockGroupSigma2(g);
        for(int i = 0; i < nb; i++)
            r[g][i] = c->getBranchRate(g, i);
    }
    Eigen::VectorXd z;
    double jac = clockCoords(k, muV, s2V, r, z);
    return lnMvn(d, z, logAges(k)) + jac;
}

void PathModel::readClockPseudoPrior(int k, const std::string& base, double burn){
    PseudoPrior& d = ppc[k];
    d.names = clockNames(k);
    int n = (int)d.names.size();
    std::string model = std::to_string(k);
    if(fromCache(d, "clock\t" + model + "\t" + base))
        return;
    int ng = rc[k]->clock->getNumClockGroups(), nb = rc[k]->clock->getNumBranchNodes();
    FBDTreeModel* fm = f[k];
    double off = UserSettings::userSettings().getAgeOffset();
    std::vector<Node*> ageNodes = fm->getAgeLogNodes();
    std::vector<Node*>& dp = fm->parameterTree->getTree()->getDownPassSequence();
    std::vector<double> t0;
    for(Node* nd : dp)
        t0.push_back(nd->getTime());
    std::vector<std::string> mainCols(d.names.begin(), d.names.begin() + 2 * ng), latCols(d.names.begin() + 2 * ng, d.names.end());
    for(size_t i = 0; i < ageNodes.size(); i++)
        mainCols.push_back("x" + std::to_string(i + 1));
    std::vector<Node*> yNodes;
    for(const FBDTreeModel::BackboneFossil& bf : fm->backboneFossils)
        if(bf.yMin < bf.yMax){
            latCols.push_back("y_" + bf.tip->getName());
            yNodes.push_back(bf.tip);
        }
    std::vector<Eigen::VectorXd> rows, ws;
    auto addLogs = [&](const std::string& mainPath, const std::string& latPath, double drop){
        if(std::ifstream(latPath).good() == false)
            return false;
        std::vector<std::vector<double>> rm = readColumns(mainPath, mainCols, "model " + model + " posterior log " + mainPath);
        std::vector<std::vector<double>> rl = readColumns(latPath, latCols, "model " + model + " posterior log " + latPath);
        if(rm.size() != rl.size())
            Msg::error("model " + model + " posterior logs " + mainPath + " and " + latPath + " have " + std::to_string(rm.size()) + " and " + std::to_string(rl.size()) + " rows; each row needs its latent row.");
        for(size_t j = (size_t)(drop * (double)rm.size()); j < rm.size(); j++){
            for(size_t i = 0; i < ageNodes.size(); i++)
                ageNodes[i]->setTime(rm[j][2 * ng + i] - off);
            for(size_t i = 0; i < yNodes.size(); i++)
                yNodes[i]->setTime(rl[j][n - 2 * ng + i] - off);
            std::vector<double> muV(ng), s2V(ng);
            std::vector<std::vector<double>> r(ng, std::vector<double>(nb));
            for(int g = 0; g < ng; g++){
                muV[g] = rm[j][2 * g];
                s2V[g] = rm[j][2 * g + 1];
                for(int i = 0; i < nb; i++)
                    r[g][i] = rl[j][g * nb + i];
            }
            Eigen::VectorXd z;
            clockCoords(k, muV, s2V, r, z);
            rows.push_back(z);
            ws.push_back(logAges(k));
        }
        return true;
    };
    if(addLogs(base + ".log", base + "_latent.log", burn) == false)
        for(int r = 0; addLogs(base + "_chain" + std::to_string(r) + ".log", base + "_chain" + std::to_string(r) + "_latent.log", UserSettings::userSettings().getBurninFraction()); r++) {}
    for(size_t i = 0; i < dp.size(); i++)
        dp[i]->setTime(t0[i]);
    fitMvn(d, rows, ws, "model " + model + " clock pseudo-prior from posterior logs " + base);
    fitted["clock\t" + model + "\t" + base] = d;
    std::cout << "Pseudo-prior of model " << model << ": " << n << " clock coordinates, conditional on " << ws[0].size() << " log node ages, fitted to " << rows.size() << " rows of its posterior logs\n";
}

double PathModel::lnTypeLabels(const std::vector<std::string>& types, std::vector<Fossil>& fossils){
    std::vector<int> count(types.size(), 0);
    int n = 0;
    for(Fossil& fo : fossils)
        if(fo.getMaxAge() > 0.0){
            count[std::find(types.begin(), types.end(), fo.getType()) - types.begin()]++;
            n++;
        }
    double lp = std::lgamma((double)types.size()) - std::lgamma((double)(n + types.size()));
    for(int c : count)
        lp += std::lgamma(c + 1.0);
    return lp;
}

std::vector<bool> PathModel::specificMask(FBDTreeModel* fm){
    std::vector<size_t> blocks = fm->rateBlockSizes(fm->a);
    std::vector<bool> spec;
    for(size_t j = 0; j < blocks.size(); j++)
        spec.insert(spec.end(), blocks[j], shared[std::min(j, (size_t)2)] == false);
    return spec;
}

void PathModel::setRateString(SkylineRates& s, const std::vector<double>& v){
    size_t k = 0;
    auto setField = [&](ParameterOUField* fl){
        std::vector<double> r(v.begin() + k, v.begin() + k + fl->getNumBins());
        k += fl->getNumBins();
        double th = v[k], sd = v[k + 1], nuv = v[k + 2];
        k += 3;
        fl->setState(r, th, sd, nuv);
    };
    if(s.lambdaField != nullptr) setField(s.lambdaField);
    else for(ParameterDouble* p : s.lambda) p->setValue(v[k++]);
    if(s.muField != nullptr) setField(s.muField);
    else for(ParameterDouble* p : s.mu) p->setValue(v[k++]);
    for(int tp = 0; tp < s.numPsiTypes; tp++){
        if(s.psiField[tp] != nullptr) setField(s.psiField[tp]);
        else for(ParameterDouble* p : s.psi[tp]) p->setValue(v[k++]);
    }
}

void PathModel::sync(int from){
    FBDTreeModel* src = f[from];
    FBDTreeModel* dst = f[1 - from];
    dst->parameterTree->copyState(*src->parameterTree);
    if(src->originAge != nullptr)
        dst->originAge->copyState(*src->originAge);
    if(src->unresolvedFossils != nullptr)
        dst->unresolvedFossils->copyState(*src->unresolvedFossils);
    if(rc[from] != nullptr){
        if(sharedClock)
            rc[1 - from]->clock->copyState(*rc[from]->clock);
        if(sharedSubst && rc[from]->ctmc != nullptr){
            std::vector<double> cv;
            rc[from]->ctmc->appendParameterValues(cv);
            rc[1 - from]->ctmc->setParameterValues(cv);
        }
    }
    SkylineRates& s = src->a;
    SkylineRates& d = dst->a;
    if(shared[0]){
        if(s.lambdaField != nullptr) d.lambdaField->copyState(*s.lambdaField);
        else for(size_t i = 0; i < s.lambda.size(); i++) d.lambda[i]->copyState(*s.lambda[i]);
    }
    if(shared[1]){
        if(s.muField != nullptr) d.muField->copyState(*s.muField);
        else for(size_t i = 0; i < s.mu.size(); i++) d.mu[i]->copyState(*s.mu[i]);
    }
    if(shared[2])
        for(int tp = 0; tp < s.numPsiTypes; tp++){
            if(s.psiField[tp] != nullptr) d.psiField[tp]->copyState(*s.psiField[tp]);
            else for(size_t i = 0; i < s.psi[tp].size(); i++) d.psi[tp][i]->copyState(*s.psi[tp][i]);
        }
    dst->zoneInit = true;
}

double PathModel::lnRatePrior(FBDTreeModel* fm){
    SkylineRates& s = fm->a;
    double lp = 0.0;
    if(shared[0] == false){
        for(ParameterDouble* p : s.lambda) lp += p->lnProbability();
        if(s.lambdaField != nullptr) lp += s.lambdaField->lnProbability();
    }
    if(shared[1] == false){
        for(ParameterDouble* p : s.mu) lp += p->lnProbability();
        if(s.muField != nullptr) lp += s.muField->lnProbability();
    }
    if(shared[2] == false){
        for(const std::vector<ParameterDouble*>& pv : s.psi) for(ParameterDouble* p : pv) lp += p->lnProbability();
        for(ParameterOUField* fl : s.psiField) if(fl != nullptr) lp += fl->lnProbability();
    }
    return lp;
}

double PathModel::lnPseudoPriorTerm(int k){
    double x = pp[k].names.empty() ? 0.0 : lnPseudoPrior(f[k], pp[k]) - lnRatePrior(f[k]);
    if(rc[k] != nullptr && sharedClock == false)
        x += lnClockPseudoPrior(k) - rc[k]->clock->lnProbability();
    if(rc[k] != nullptr && sharedSubst == false && rc[k]->ctmc != nullptr)
        x += lnSubstPseudoPrior(k) - rc[k]->ctmc->lnPrior();
    return x;
}

double PathModel::lnPseudoPrior(FBDTreeModel* fm, const PseudoPrior& d){
    std::vector<double> v = fm->getRateString(fm->a), u;
    double jac = 0.0;
    for(const std::pair<ParameterOUField*,int>& fl : d.fields){
        jac += fl.first->standardizedDeviations(u);
        std::copy(u.begin(), u.end(), v.begin() + fl.second);
    }
    int k = (int)d.names.size();
    Eigen::VectorXd z(k);
    for(int i = 0; i < k; i++){
        double x = v[d.index[i]];
        if(d.coord[i] == PseudoPriorCoord::LOG){
            x = std::log(x);
            jac -= x;
        }
        z(i) = x;
    }
    return lnMvn(d, z, Eigen::VectorXd()) + jac;
}

void PathModel::writePseudoPriorState(std::ostream& os, const PseudoPrior& d){
    int k = (int)d.mean.size();
    os << "pseudo_prior " << k;
    for(int i = 0; i < k; i++)
        os << ' ' << d.mean(i);
    for(int i = 0; i < k; i++)
        for(int j = 0; j < k; j++)
            os << ' ' << d.cov(i, j);
    os << '\n';
}

int PathModel::readPseudoPriorState(std::istream& is, PseudoPrior& d){
    std::string tag;
    int k = -1;
    is >> tag >> k;
    if(tag != "pseudo_prior" || k != (int)d.mean.size())
        return -1;
    Eigen::VectorXd mean(k);
    Eigen::MatrixXd cov(k, k);
    for(int i = 0; i < k; i++)
        is >> mean(i);
    for(int i = 0; i < k; i++)
        for(int j = 0; j < k; j++)
            is >> cov(i, j);
    if(mean == d.mean && cov == d.cov)
        return 0;
    d.mean = mean;
    d.cov = cov;
    d.chol.compute(cov);
    return 1;
}

void PathModel::fitRatePseudoPrior(FBDTreeModel* fm, PseudoPrior& d, const std::vector<bool>& spec, const std::string& path, double burn, const std::string& model){
    SkylineRates& s = fm->a;
    std::vector<std::string> own = fm->getRateNames(s);
    std::vector<bool> fixed, sharedHyper(own.size(), false);
    fm->appendRateMask(fixed, s);
    auto addField = [&](ParameterOUField* fld, const std::vector<std::string>& binNames){
        if(fld == nullptr || fld->sharesHyperparameters() == false) return;
        int j = (int)(std::find(own.begin(), own.end(), binNames[0]) - own.begin());
        d.fields.push_back({fld, j});
        std::fill(sharedHyper.begin() + j + fld->getNumBins(), sharedHyper.begin() + j + fld->getNumBins() + 3, true);
    };
    addField(s.lambdaField, s.lambdaName);
    addField(s.muField, s.muName);
    for(int tp = 0; tp < s.numPsiTypes; tp++) addField(s.psiField[tp], s.psiName[tp]);
    for(size_t j = 0; j < own.size(); j++){
        if(fixed[j] || sharedHyper[j] || spec[j] == false) continue;
        const std::string& n = own[j];
        PseudoPriorCoord c = (n.rfind("ou_", 0) == 0 && n.size() >= 6 && n.compare(n.size() - 6, 6, "_theta") == 0) ? PseudoPriorCoord::IDENTITY : PseudoPriorCoord::LOG;
        for(const std::pair<ParameterOUField*,int>& fp : d.fields)
            if((int)j >= fp.second && (int)j < fp.second + fp.first->getNumBins()) c = PseudoPriorCoord::STANDARDIZED_DEVIATION;
        d.names.push_back(n);
        d.index.push_back((int)j);
        d.coord.push_back(c);
    }
    int k = (int)d.names.size();
    if(k == 0)
        return;
    std::vector<bool> needed(own.size(), false);
    for(int j : d.index) needed[j] = true;
    for(const std::pair<ParameterOUField*,int>& fp : d.fields)
        std::fill(needed.begin() + fp.second, needed.begin() + fp.second + fp.first->getNumBins() + 3, true);
    if(fromCache(d, "rates\t" + model + "\t" + path))
        return;
    std::vector<std::string> cols;
    for(size_t j = 0; j < own.size(); j++)
        if(needed[j]) cols.push_back(own[j]);
    std::vector<std::vector<double>> table = readColumns(path, cols, model + " posterior log " + path);
    table.erase(table.begin(), table.begin() + (size_t)(burn * (double)table.size()));
    std::vector<Eigen::VectorXd> rows;
    for(const std::vector<double>& r : table){
        std::vector<double> v(own.size(), 0.0), u;
        size_t c = 0;
        for(size_t j = 0; j < own.size(); j++)
            if(needed[j]) v[j] = r[c++];
        for(const std::pair<ParameterOUField*,int>& fp : d.fields){
            int nb = fp.first->getNumBins(), j0 = fp.second;
            fp.first->standardizedDeviations(std::vector<double>(v.begin() + j0, v.begin() + j0 + nb), v[j0 + nb], v[j0 + nb + 1], v[j0 + nb + 2], u);
            std::copy(u.begin(), u.end(), v.begin() + j0);
        }
        Eigen::VectorXd z(k);
        for(int i = 0; i < k; i++)
            z(i) = d.coord[i] == PseudoPriorCoord::LOG ? std::log(v[d.index[i]]) : v[d.index[i]];
        rows.push_back(z);
    }
    fitMvn(d, rows, {}, model + " rate pseudo-prior from posterior log " + path);
    fitted["rates\t" + model + "\t" + path] = d;
    std::cout << "Pseudo-prior of " << model << ": " << k << " parameters fitted to the " << rows.size() << " rows of " << path << "\n";
}

double PathModel::lnSubstPseudoPrior(int k){
    std::vector<double> v, z;
    rc[k]->ctmc->appendParameterValues(v);
    double jac = rc[k]->ctmc->pseudoPriorCoords(v, z);
    return lnMvn(pps[k], Eigen::Map<const Eigen::VectorXd>(z.data(), (int)z.size()), Eigen::VectorXd()) + jac;
}

void PathModel::readSubstPseudoPrior(int k, const std::string& path, double burn){
    std::string model = std::to_string(k);
    if(fromCache(pps[k], "substitution\t" + model + "\t" + path))
        return;
    std::vector<std::string> names;
    rc[k]->ctmc->appendParameterNames(names);
    std::vector<std::vector<double>> table = readColumns(path, names, "model " + model + " posterior log " + path);
    table.erase(table.begin(), table.begin() + (size_t)(burn * (double)table.size()));
    std::vector<Eigen::VectorXd> rows;
    for(const std::vector<double>& v : table){
        std::vector<double> z;
        rc[k]->ctmc->pseudoPriorCoords(v, z);
        rows.push_back(Eigen::Map<const Eigen::VectorXd>(z.data(), (int)z.size()));
    }
    fitMvn(pps[k], rows, {}, "model " + model + " substitution pseudo-prior from posterior log " + path);
    fitted["substitution\t" + model + "\t" + path] = pps[k];
    std::cout << "Pseudo-prior of model " << model << ": " << rows[0].size() << " substitution coordinates fitted to the " << rows.size() << " rows of " << path << "\n";
}

std::map<std::string,double> PathModel::readLastRow(const std::string& path){
    std::ifstream in(path);
    if(in.is_open() == false)
        Msg::error("could not open start log: " + path);
    std::string line, last, tok;
    std::getline(in, line);
    std::vector<std::string> header;
    std::istringstream hs(line);
    while(std::getline(hs, tok, '\t'))
        header.push_back(tok);
    while(std::getline(in, line))
        if(line.empty() == false)
            last = line;
    if(last.empty())
        Msg::error("start log " + path + " has no data rows.");
    std::vector<double> vals;
    std::istringstream rs(last);
    while(std::getline(rs, tok, '\t')){
        try { vals.push_back(std::stod(tok)); }
        catch(...) { Msg::error("start log " + path + ": '" + tok + "' in the last row is not a number."); }
    }
    if(vals.size() != header.size())
        Msg::error("start log " + path + ": the last row has " + std::to_string(vals.size()) + " values, the header has " + std::to_string(header.size()) + " names.");
    std::map<std::string,double> row;
    for(size_t i = 0; i < header.size(); i++)
        row[header[i]] = vals[i];
    return row;
}

void PathModel::setState(int km, const std::map<std::string,double>& mainRow, const std::map<std::string,double>& latRow, const std::string& mainLog, const std::string& latentLog){
    FBDTreeModel* fm = f[km];
    double off = UserSettings::userSettings().getAgeOffset();
    std::vector<std::string> names = fm->getRateNames(fm->a);
    std::vector<double> v(names.size());
    for(size_t i = 0; i < names.size(); i++)
        v[i] = column(mainRow, names[i], mainLog);
    setRateString(fm->a, v);

    Tree* tr = fm->parameterTree->getTree();
    if(fm->originAge != nullptr){
        double x0 = column(mainRow, rc[km] != nullptr ? "x0" : "originAge", mainLog) - off;
        fm->originAge->setValue(x0);
        Node* rt = tr->getRoot();
        if(rt->getDescendants().size() == 1)
            rt->setTime(x0);
    }
    std::vector<Node*> ageNodes = fm->getAgeLogNodes();
    for(size_t i = 0; i < ageNodes.size(); i++)
        ageNodes[i]->setTime(column(mainRow, "x" + std::to_string(i + 1), mainLog) - off);
    for(const FBDTreeModel::BackboneFossil& bf : fm->backboneFossils)
        if(bf.yMin < bf.yMax)
            bf.tip->setTime(column(latRow, "y_" + bf.tip->getName(), latentLog) - off);
    fm->parameterTree->commitTree();

    RandomVariable* prevRng = RandomVariable::getActiveInstance();
    RandomVariable::setActiveInstance(&fm->rng);
    ParameterUnresolvedFossils* uf = fm->unresolvedFossils;
    if(uf != nullptr){
        int nf = uf->getNumFossils();
        std::vector<double> yv(nf);
        for(int i = 0; i < nf; i++)
            yv[i] = (uf->isUE(i) == false && uf->getYMin(i) < uf->getYMax(i)) ? column(latRow, "y_" + fm->unrFossilName[i], latentLog) - off : uf->getFossilAge(i);
        uf->setFossilAges(yv);
        for(int i = 0; i < nf; i++){
            const std::string& nm = fm->unrFossilName[i];
            if(uf->isUE(i) == false && column(latRow, "sa_" + nm, latentLog) != 0.0)
                uf->flipSA(i, true);
            else
                uf->setAttachAge(i, column(latRow, "z_" + nm, latentLog) - off);
            std::vector<int>& dom = uf->getAttachmentZoneDomain(i);
            int zone = dom[0];
            if(dom.size() > 1){
                zone = (int)column(latRow, "az_" + nm, latentLog);
                if(std::find(dom.begin(), dom.end(), zone) == dom.end())
                    Msg::error("start log " + latentLog + ": az_" + nm + " = " + std::to_string(zone) + " is not an attachment zone of this fossil.");
            }
            uf->initAttachmentZone(i, zone);
            uf->acceptFossil(i);
        }
    }
    RandomVariable::setActiveInstance(prevRng);
    fm->invalidateGammaCache();
    fm->zoneIndexBuilt = false;
    fm->eulerBuilt = false;
    fm->zoneInit = true;
}

void PathModel::setStart(int run){
    UserSettings& us = UserSettings::userSettings();
    int k = (pathT > 0.5) ? 1 : 0;
    std::string base = (k == 1) ? us.getPosteriorBaseB() : us.getPosteriorBaseA();
    int chains = 0;
    while(std::ifstream(base + "_chain" + std::to_string(chains) + ".log").good())
        chains++;
    std::string log = (chains > 0) ? base + "_chain" + std::to_string(run % chains) + ".log" : base + ".log";
    std::string latent = log.substr(0, log.size() - 4) + "_latent.log";
    std::map<std::string,double> mainRow = readLastRow(log), latRow = readLastRow(latent);
    if(f[k]->unresolvedFossils != nullptr && f[k]->zoneIndexBuilt == false)
        f[k]->buildZoneIndex();
    setState(k, mainRow, latRow, log, latent);
    if(rc[k] != nullptr){
        std::vector<std::string> cn = clockNames(k);
        std::vector<double> cv;
        for(const std::string& s : cn)
            cv.push_back(mainRow.count(s) > 0 ? mainRow[s] : column(latRow, s, latent));
        setClockValues(k, cv);
        if(rc[k]->ctmc != nullptr){
            std::vector<std::string> sn;
            rc[k]->ctmc->appendParameterNames(sn);
            std::vector<double> sv;
            for(const std::string& s : sn)
                sv.push_back(column(mainRow, s, log));
            rc[k]->ctmc->setParameterValues(sv);
        }
        m[k]->invalidatePriorCache();
    }
    sync(k);
    if(rc[1 - k] != nullptr && sharedClock == false){
        PseudoPrior& dc = ppc[1 - k];
        Eigen::MatrixXd Lc = dc.chol.matrixL();
        Eigen::VectorXd mc = dc.mean + dc.K * (logAges(1 - k) - dc.wMean);
        int ng = rc[1 - k]->clock->getNumClockGroups();
        std::vector<double> muV(ng), s2V(ng), u(mc.data() + 2 * ng, mc.data() + mc.size());
        for(int g = 0; g < ng; g++){
            muV[g] = std::exp(mc(2 * g) + 0.5 * Lc.row(2 * g).squaredNorm());
            s2V[g] = std::exp(mc(2 * g + 1) + 0.5 * Lc.row(2 * g + 1).squaredNorm());
        }
        static_cast<ParameterBranchRates*>(rc[1 - k]->clock)->setStandardizedDeviations(muV, s2V, u);
        m[1 - k]->invalidateLikelihoodCache();
    }
    if(rc[1 - k] != nullptr && sharedSubst == false && rc[1 - k]->ctmc != nullptr){
        const Eigen::VectorXd& ms = pps[1 - k].mean;
        rc[1 - k]->ctmc->setFromPseudoPriorCoords(std::vector<double>(ms.data(), ms.data() + ms.size()));
        m[1 - k]->invalidateLikelihoodCache();
    }
    PseudoPrior& d = pp[1 - k];
    if(d.names.empty() == false){
        FBDTreeModel* o = f[1 - k];
        std::vector<double> v = o->getRateString(o->a);
        Eigen::MatrixXd L = d.chol.matrixL();
        for(size_t i = 0; i < d.names.size(); i++){
            double mi = d.mean((int)i);
            v[d.index[i]] = d.coord[i] == PseudoPriorCoord::LOG ? std::exp(mi + 0.5 * L.row((int)i).squaredNorm()) : mi;
        }
        setRateString(o->a, v);
        for(const std::pair<ParameterOUField*,int>& fl : d.fields)
            fl.first->setStandardizedDeviations(std::vector<double>(v.begin() + fl.second, v.begin() + fl.second + fl.first->getNumBins()));
        o->invalidateGammaCache();
    }
}

double PathModel::dataPart(int k){
    double v = m[(rc[k] != nullptr && sharedClock && sharedSubst) ? mover : k]->lnLikelihood();
    if(m[k] != f[k])
        v += f[k]->lnLikelihood();
    return v + lab[k];
}

double PathModel::lnLikelihood(void){
    double v = (1.0 - pathT) * dataPart(0) + pathT * dataPart(1);
    if(std::isfinite(v) == false)
        return -INFINITY;
    return v;
}

double PathModel::lnPriorProbability(void){
    double p = f[0]->lnPriorProbability();
    if(rc[0] != nullptr){
        p += rc[0]->clock->lnProbability();
        if(rc[0]->ctmc != nullptr)
            p += rc[0]->ctmc->lnPrior();
    }
    return p + specificPrior(1) + pathT * lnPseudoPriorTerm(0) + (1.0 - pathT) * lnPseudoPriorTerm(1);
}

size_t PathModel::rateEnd(void){
    std::vector<std::string> names = m[0]->getParameterNames(), ra = f[0]->getRateNames(f[0]->a);
    return (size_t)(std::find(names.begin(), names.end(), ra[0]) - names.begin()) + ra.size();
}

std::vector<std::string> PathModel::getParameterNames(void){
    std::vector<std::string> names = m[0]->getParameterNames();
    size_t n = rateEnd();
    std::vector<std::string> rb = f[1]->getRateNames(f[1]->a), ins;
    std::vector<bool> spec = specificMask(f[1]);
    for(size_t i = 0; i < rb.size(); i++)
        if(spec[i]) ins.push_back("B_" + rb[i]);
    if(rc[1] != nullptr && sharedClock == false)
        for(int g = 0; g < 2 * rc[1]->clock->getNumClockGroups(); g++)
            ins.push_back("B_" + clockNames(1)[g]);
    ins.push_back("delta");
    names.insert(names.begin() + n, ins.begin(), ins.end());
    return names;
}

std::vector<double> PathModel::getParameterString(void){
    std::vector<double> vals = m[0]->getParameterString();
    size_t n = rateEnd();
    std::vector<double> rb = f[1]->getRateString(f[1]->a), ins;
    std::vector<bool> spec = specificMask(f[1]);
    for(size_t i = 0; i < rb.size(); i++)
        if(spec[i]) ins.push_back(rb[i]);
    if(rc[1] != nullptr && sharedClock == false){
        std::vector<double> cv = clockValues(1);
        ins.insert(ins.end(), cv.begin(), cv.begin() + 2 * rc[1]->clock->getNumClockGroups());
    }
    ins.push_back(dataPart(1) - dataPart(0) + lnPseudoPriorTerm(0) - lnPseudoPriorTerm(1));
    vals.insert(vals.begin() + n, ins.begin(), ins.end());
    return vals;
}

std::vector<bool> PathModel::getParameterFixedMask(void){
    std::vector<bool> mask = m[0]->getParameterFixedMask();
    size_t n = rateEnd();
    std::vector<bool> mb, ins;
    f[1]->appendRateMask(mb, f[1]->a);
    std::vector<bool> spec = specificMask(f[1]);
    for(size_t i = 0; i < mb.size(); i++)
        if(spec[i]) ins.push_back(mb[i]);
    if(rc[1] != nullptr && sharedClock == false)
        ins.insert(ins.end(), 2 * rc[1]->clock->getNumClockGroups(), false);
    ins.push_back(false);
    mask.insert(mask.begin() + n, ins.begin(), ins.end());
    return mask;
}

std::vector<std::string> PathModel::getLatentNames(void){ return m[0]->getLatentNames(); }

std::vector<double> PathModel::getLatentString(void){ return m[0]->getLatentString(); }

std::vector<std::string> PathModel::getZoneLegend(void){ return m[0]->getZoneLegend(); }

bool PathModel::treeIncludesFossils(void){ return m[0]->treeIncludesFossils(); }

int PathModel::getNumFossils(void){ return m[0]->getNumFossils(); }

void PathModel::invalidateLikelihoodCache(void){
    m[0]->invalidateLikelihoodCache();
    m[1]->invalidateLikelihoodCache();
}

void PathModel::invalidatePriorCache(void){
    m[0]->invalidatePriorCache();
    m[1]->invalidatePriorCache();
}

int PathModel::getLastMoveType(void){ return m[mover]->getLastMoveType(); }

void PathModel::setHeat(double h){
    heat = h;
    m[0]->setHeat(h);
    m[1]->setHeat(h);
}

void PathModel::print(void){ m[0]->print(); }

double PathModel::update(void){
    mover = (rng.uniformRv() < 0.5) ? 0 : 1;
    double r = m[mover]->update();
    sync(mover);
    return r;
}

void PathModel::updateForAcceptance(void){
    m[mover]->updateForAcceptance();
    sync(mover);
}

void PathModel::updateForRejection(void){
    m[mover]->updateForRejection();
    sync(mover);
}

void PathModel::writeState(std::ostream& os){
    for(int k = 0; k < 2; k++){
        m[k]->getRng()->writeState(os);
        m[k]->writeState(os);
    }
}

void PathModel::readState(std::istream& is){
    for(int k = 0; k < 2; k++){
        m[k]->getRng()->readState(is);
        m[k]->readState(is);
    }
    sync(0);
}
