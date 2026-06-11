#include <cmath>

#include "AdaptiveOperatorSampler.hpp"

AdaptiveOperatorSampler::AdaptiveOperatorSampler(int nMoves){
    n = nMoves;
    adapting = false;
    adaptN = 0.0;
    picks = 0;
    weight.assign(n, 1.0 / n);
    rewSum.assign(n, 0.0);
    rewCnt.assign(n, 0.0);
    active.assign(n, 1);
    pinned.assign(n, 0);
}

int AdaptiveOperatorSampler::pick(double u){
    double cum = 0.0;
    for(int m = 0; m < n; m++){
        cum += weight[m];
        if(u < cum)
            return m;
    }
    for(int m = n - 1; m >= 0; m--)
        if(active[m])
            return m;
    return 0;
}

void AdaptiveOperatorSampler::reward(int move, double r){
    rewSum[move] += r;
    rewCnt[move] += 1.0;
    picks++;
    if(adapting && picks % 200 == 0)
        adapt();
}

void AdaptiveOperatorSampler::adapt(void){
    adaptN += 1.0;
    double g = 1.0 / std::sqrt(adaptN);
    double pinnedSum = 0.0;
    for(int m = 0; m < n; m++)
        if(pinned[m])
            pinnedSum += weight[m];
    double avail = 1.0 - pinnedSum;
    if(avail < 0.0)
        avail = 0.0;
    int na = 0;
    for(int m = 0; m < n; m++)
        if(active[m] && pinned[m] == 0)
            na++;
    if(na == 0)
        return;
    double fl = 0.05;
    if(fl * na > 0.8)
        fl = 0.8 / na;
    std::vector<double> r(n, 0.0);
    double sumR = 0.0;
    for(int m = 0; m < n; m++){
        if(active[m] && pinned[m] == 0 && rewCnt[m] > 0.0){
            r[m] = rewSum[m] / rewCnt[m];
            sumR += r[m];
        }
    }
    for(int m = 0; m < n; m++){
        if(pinned[m])
            continue;
        double target;
        if(active[m] == 0)
            target = 0.0;
        else if(sumR > 0.0)
            target = avail * (fl + (1.0 - fl * na) * (r[m] / sumR));
        else
            target = avail / na;
        weight[m] += g * (target - weight[m]);
    }
    for(int m = 0; m < n; m++){
        rewSum[m] *= 0.5;
        rewCnt[m] *= 0.5;
    }
}
