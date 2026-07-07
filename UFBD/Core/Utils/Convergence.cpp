#include "Convergence.hpp"
#include "Probability.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <utility>

namespace {

typedef std::vector<std::vector<double>> Chains;

Chains splitChains(const Chains& chains){
    size_t minLen = 0;
    bool first = true;
    for(const std::vector<double>& c : chains){
        if(first || c.size() < minLen){ minLen = c.size(); first = false; }
    }
    size_t n = minLen / 2;
    Chains out;
    if(n < 2)
        return out;
    for(const std::vector<double>& c : chains){
        out.push_back(std::vector<double>(c.begin(), c.begin() + n));
        out.push_back(std::vector<double>(c.begin() + n, c.begin() + 2 * n));
    }
    return out;
}

void rankNormalize(Chains& chains){
    size_t S = 0;
    for(const std::vector<double>& c : chains)
        S += c.size();
    std::vector<std::pair<double, std::pair<int,int>>> all;
    all.reserve(S);
    for(int j = 0; j < (int)chains.size(); j++)
        for(int i = 0; i < (int)chains[j].size(); i++)
            all.push_back(std::make_pair(chains[j][i], std::make_pair(j, i)));
    std::sort(all.begin(), all.end(), [](const std::pair<double, std::pair<int,int>>& a, const std::pair<double, std::pair<int,int>>& b){ return a.first < b.first; });
    std::vector<double> ranks(S);
    size_t i = 0;
    while(i < S){
        size_t k = i;
        while(k + 1 < S && all[k+1].first == all[i].first) k++;
        double avg = (double)(i + k) / 2.0 + 1.0;
        for(size_t t = i; t <= k; t++) ranks[t] = avg;
        i = k + 1;
    }
    for(size_t t = 0; t < S; t++){
        double p = (ranks[t] - 0.375) / ((double)S - 0.25);
        double z = Probability::Normal::invCdfStandardNormal(p);
        chains[all[t].second.first][all[t].second.second] = z;
    }
}

double rhatOf(const Chains& chains){
    int m = (int)chains.size();
    int n = (int)chains[0].size();
    if(m < 2 || n < 2)
        return std::nan("");
    std::vector<double> means(m);
    double grand = 0.0;
    for(int j = 0; j < m; j++){
        double s = 0.0;
        for(double v : chains[j]) s += v;
        means[j] = s / n;
        grand += means[j];
    }
    grand /= m;
    double B = 0.0;
    for(int j = 0; j < m; j++) B += (means[j] - grand) * (means[j] - grand);
    B *= (double)n / (m - 1);
    double W = 0.0;
    for(int j = 0; j < m; j++){
        double s2 = 0.0;
        for(double v : chains[j]) s2 += (v - means[j]) * (v - means[j]);
        W += s2 / (n - 1);
    }
    W /= m;
    if(W <= 0.0)
        return std::nan("");
    double varPlus = ((double)(n - 1) * W + B) / n;
    return std::sqrt(varPlus / W);
}

void fft(std::vector<std::complex<double>>& a, bool invert){
    int n = (int)a.size();
    for(int i = 1, j = 0; i < n; i++){
        int bit = n >> 1;
        for(; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if(i < j) std::swap(a[i], a[j]);
    }
    for(int len = 2; len <= n; len <<= 1){
        double ang = 2.0 * std::acos(-1.0) / len * (invert ? 1.0 : -1.0);
        std::complex<double> wlen(std::cos(ang), std::sin(ang));
        for(int i = 0; i < n; i += len){
            std::complex<double> w(1.0, 0.0);
            for(int k = 0; k < len / 2; k++){
                std::complex<double> u = a[i + k];
                std::complex<double> v = a[i + k + len / 2] * w;
                a[i + k] = u + v;
                a[i + k + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if(invert)
        for(std::complex<double>& x : a) x /= (double)n;
}

std::vector<double> autocovSeq(const std::vector<double>& x, double mean){
    int n = (int)x.size();
    int L = 1;
    while(L < 2 * n) L <<= 1;
    std::vector<std::complex<double>> f(L, std::complex<double>(0.0, 0.0));
    for(int i = 0; i < n; i++) f[i] = std::complex<double>(x[i] - mean, 0.0);
    fft(f, false);
    for(int i = 0; i < L; i++) f[i] = std::complex<double>(std::norm(f[i]), 0.0);
    fft(f, true);
    std::vector<double> c(n);
    for(int t = 0; t < n; t++) c[t] = f[t].real() / (double)n;
    return c;
}

double multiChainEss(const Chains& chains){
    int m = (int)chains.size();
    if(m < 1) return 0.0;
    int n = (int)chains[0].size();
    if(n < 4) return (double)(m * n);
    std::vector<double> means(m), vars(m);
    std::vector<std::vector<double>> acov(m);
    for(int j = 0; j < m; j++){
        double s = 0.0;
        for(double v : chains[j]) s += v;
        means[j] = s / n;
        acov[j] = autocovSeq(chains[j], means[j]);
        vars[j] = acov[j][0] * (double)n / (n - 1);
    }
    double grand = 0.0;
    for(double mu : means) grand += mu;
    grand /= m;
    double B = 0.0;
    for(int j = 0; j < m; j++) B += (means[j] - grand) * (means[j] - grand);
    B *= (double)n / (m > 1 ? m - 1 : 1);
    double W = 0.0;
    for(double v : vars) W += v;
    W /= m;
    if(W <= 0.0) return (double)(m * n);
    double varPlus = (m > 1) ? (((double)(n - 1) * W + B) / n) : W;

    auto rhoHat = [&](int t) -> double {
        double meanAcov = 0.0;
        for(int j = 0; j < m; j++) meanAcov += acov[j][t];
        meanAcov /= m;
        return 1.0 - (W - meanAcov) / varPlus;
    };

    double sum = 1.0 + rhoHat(1);
    double prev = sum;
    for(int k = 1; 2 * k + 1 < n; k++){
        double g = rhoHat(2 * k) + rhoHat(2 * k + 1);
        if(g <= 0.0) break;
        if(g > prev) g = prev;
        prev = g;
        sum += g;
    }
    double tau = -1.0 + 2.0 * sum;
    if(tau < 1.0) tau = 1.0;
    return (double)(m * n) / tau;
}

double quantileSorted(const std::vector<double>& sorted, double p){
    double idx = p * (double)(sorted.size() - 1);
    size_t lo = (size_t)idx;
    if(lo + 1 >= sorted.size()) return sorted.back();
    double frac = idx - (double)lo;
    return sorted[lo] * (1.0 - frac) + sorted[lo+1] * frac;
}

}

namespace Convergence {

Diagnostic diagnose(const std::vector<std::vector<double>>& raw){
    Diagnostic d;
    d.rhat = std::nan("");
    d.bulkEss = std::nan("");
    d.tailEss = std::nan("");

    Chains split = splitChains(raw);
    if(split.size() < 2 || split[0].size() < 2)
        return d;

    std::vector<double> sorted;
    for(const std::vector<double>& c : split)
        for(double v : c) sorted.push_back(v);
    std::sort(sorted.begin(), sorted.end());

    if(sorted.back() - sorted.front() <= 1e-12 * (std::fabs(sorted.front()) + 1.0)){
        d.bulkEss = (double)sorted.size();
        d.tailEss = (double)sorted.size();
        return d;
    }

    double median = quantileSorted(sorted, 0.5);

    Chains bulk = split;
    rankNormalize(bulk);
    double rB = rhatOf(bulk);

    Chains fold = split;
    for(std::vector<double>& c : fold)
        for(double& v : c) v = std::fabs(v - median);
    rankNormalize(fold);
    double rF = rhatOf(fold);

    if(std::isnan(rB) && std::isnan(rF))
        d.rhat = std::nan("");
    else
        d.rhat = std::fmax(std::isnan(rB) ? 0.0 : rB, std::isnan(rF) ? 0.0 : rF);

    d.bulkEss = multiChainEss(bulk);

    double q05 = quantileSorted(sorted, 0.05);
    double q95 = quantileSorted(sorted, 0.95);
    Chains lo = split, hi = split;
    for(std::vector<double>& c : lo)
        for(double& v : c) v = (v <= q05) ? 1.0 : 0.0;
    for(std::vector<double>& c : hi)
        for(double& v : c) v = (v >= q95) ? 1.0 : 0.0;
    d.tailEss = std::fmin(multiChainEss(lo), multiChainEss(hi));

    return d;
}

bool meetsThresholds(const Diagnostic& d, double rhatMax, double essMin){
    if(std::isnan(d.rhat))
        return true;
    if(d.rhat > rhatMax)
        return false;
    if(d.bulkEss < essMin || d.tailEss < essMin)
        return false;
    return true;
}

}
