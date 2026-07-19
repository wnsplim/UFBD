#include "ConfigReader.hpp"
#include "Msg.hpp"

#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::string trim(const std::string& s){
    size_t a = 0, b = s.size();
    while(a < b && std::isspace((unsigned char)s[a])) a++;
    while(b > a && std::isspace((unsigned char)s[b - 1])) b--;
    return s.substr(a, b - a);
}

std::vector<std::string> splitChar(const std::string& s, char sep){
    std::vector<std::string> out;
    std::string cur;
    for(char c : s){ if(c == sep){ out.push_back(cur); cur.clear(); } else cur += c; }
    out.push_back(cur);
    return out;
}

bool isFamily(const std::string& k){
    return k == "exp" || k == "exponential" || k == "gamma" || k == "lognormal"
        || k == "unif" || k == "uniform" || k == "truncnormal" || k == "truncnorm"
        || k == "normal" || k == "improper";
}

double parseBound(const std::string& s){
    std::string l = s;
    for(char& c : l) c = (char)std::tolower((unsigned char)c);
    if(l == "inf" || l == "infinity") return std::numeric_limits<double>::infinity();
    try { return std::stod(s); }
    catch(...) { Msg::error("config time_bins: '" + s + "' is not a number or 'inf'"); }
    return 0.0;
}

std::string fmtNum(double x){
    if(std::isinf(x)) return "inf";
    std::ostringstream os; os << x; return os.str();
}

void emit(std::vector<std::string>& out, const std::string& flag, const std::string& value){
    out.push_back(flag);
    std::istringstream is(value);
    std::string tok;
    while(is >> tok) out.push_back(tok);
}

std::string getVal(const std::map<std::string,std::string>& kv, const std::string& k){
    std::map<std::string,std::string>::const_iterator it = kv.find(k);
    return it == kv.end() ? std::string() : it->second;
}

// grammar: name:(lo,hi)[+(lo,hi)...] | name:(lo,hi) | ...
void translateTimeBins(const std::string& spec, std::map<std::string,int>& nameToGid,
                       std::vector<double>& cuts, std::vector<int>& groups, double& outTmin){
    struct Iv { double lo, hi; int gid; };
    std::vector<Iv> ivs;
    for(std::string bt : splitChar(spec, '|')){
        bt = trim(bt);
        if(bt.empty()) continue;
        size_t c = bt.find(':');
        if(c == std::string::npos)
            Msg::error("config time_bins: expected 'name:(lo,hi)[+(lo,hi)]', got '" + bt + "'");
        std::string name = trim(bt.substr(0, c));
        if(nameToGid.count(name) == 0){ int g = (int)nameToGid.size(); nameToGid[name] = g; }
        int gid = nameToGid[name];
        for(std::string p : splitChar(bt.substr(c + 1), '+')){
            p = trim(p);
            size_t o = p.find('('), cl = p.find(')');
            if(o == std::string::npos || cl == std::string::npos || cl <= o)
                Msg::error("config time_bins: expected '(lo,hi)', got '" + p + "'");
            std::string inner = p.substr(o + 1, cl - o - 1);
            size_t comma = inner.find(',');
            if(comma == std::string::npos)
                Msg::error("config time_bins: expected '(lo,hi)', got '" + p + "'");
            double lo = parseBound(trim(inner.substr(0, comma)));
            double hi = parseBound(trim(inner.substr(comma + 1)));
            if(hi <= lo)
                Msg::error("config time_bins: need hi > lo in '" + p + "'");
            ivs.push_back({lo, hi, gid});
        }
    }
    double tmin = std::numeric_limits<double>::infinity();
    for(const Iv& iv : ivs)
        if(iv.lo < tmin) tmin = iv.lo;
    if(std::isinf(tmin)) tmin = 0.0;
    outTmin = tmin;
    std::set<double> cutSet;
    for(const Iv& iv : ivs){
        if(iv.lo > tmin) cutSet.insert(iv.lo);
        if(std::isinf(iv.hi) == false) cutSet.insert(iv.hi);
    }
    cuts.assign(cutSet.begin(), cutSet.end());
    std::vector<double> los, his;
    double prev = tmin;
    for(double cc : cuts){ los.push_back(prev); his.push_back(cc); prev = cc; }
    los.push_back(prev); his.push_back(std::numeric_limits<double>::infinity());
    groups.clear();
    for(size_t b = 0; b < los.size(); b++){
        double mid = std::isinf(his[b]) ? los[b] + 1.0 : 0.5 * (los[b] + his[b]);
        int found = -1, count = 0;
        for(const Iv& iv : ivs)
            if(iv.lo <= mid && mid < iv.hi){ found = iv.gid; count++; }
        if(count == 0)
            Msg::error("config time_bins leave a gap on [" + fmtNum(los[b]) + ", " + fmtNum(his[b]) + "); bins must tile [" + fmtNum(tmin) + ", inf).");
        if(count > 1)
            Msg::error("config time_bins overlap on [" + fmtNum(los[b]) + ", " + fmtNum(his[b]) + "); bins must be disjoint.");
        groups.push_back(found);
    }
}

void translatePrior(const std::string& ratePrefix, const std::string& typeLabel,
                    const std::string& spec, const std::map<std::string,int>& nameToGid,
                    std::vector<std::string>& out){
    std::string flag = "-" + ratePrefix + "_prior";
    std::vector<std::string> bins = splitChar(spec, '|');
    // <prior> (all bins) OR binName:<prior> | binName:<prior> (per chunk)
    if(bins.size() == 1){
        std::string t = trim(bins[0]);
        size_t c = t.find(':');
        std::string head = (c == std::string::npos) ? t : t.substr(0, c);
        if(c == std::string::npos || isFamily(head)){
            emit(out, flag, typeLabel.empty() ? t : typeLabel + ":" + t);
            return;
        }
    }
    for(std::string bt : bins){
        bt = trim(bt);
        if(bt.empty()) continue;
        size_t c = bt.find(':');
        if(c == std::string::npos)
            Msg::error("config prior: expected 'binName:prior', got '" + bt + "'");
        std::string name = trim(bt.substr(0, c));
        std::string pr = bt.substr(c + 1);
        std::map<std::string,int>::const_iterator it = nameToGid.find(name);
        if(it == nameToGid.end())
            Msg::error("config prior references bin '" + name + "' not defined in this rate's time_bins");
        std::string val = std::to_string(it->second) + ":" + pr;
        emit(out, flag, typeLabel.empty() ? val : typeLabel + ":" + val);
    }
}

void translateRateSection(const std::string& ratePrefix, const std::string& typeLabel,
                          const std::map<std::string,std::string>& kv, std::vector<std::string>& out, double& outTmin){
    std::map<std::string,int> nameToGid;
    outTmin = -1.0;
    std::string tb = getVal(kv, "time_bins");
    if(tb.empty() == false){
        std::vector<double> cuts; std::vector<int> groups;
        translateTimeBins(tb, nameToGid, cuts, groups, outTmin);
        if(cuts.empty() == false){
            std::string cutStr, grpStr;
            for(size_t i = 0; i < cuts.size(); i++){ if(i) cutStr += ","; cutStr += fmtNum(cuts[i]); }
            for(size_t i = 0; i < groups.size(); i++){ if(i) grpStr += ","; grpStr += std::to_string(groups[i]); }
            emit(out, "-" + ratePrefix + "_skyline_times", typeLabel.empty() ? cutStr : typeLabel + ":" + cutStr);
            emit(out, "-" + ratePrefix + "_groups", typeLabel.empty() ? grpStr : typeLabel + ":" + grpStr);
        }
    }
    std::string pr = getVal(kv, "prior");
    if(pr.empty() == false)
        translatePrior(ratePrefix, typeLabel, pr, nameToGid, out);
    std::string md = getVal(kv, "mode");
    if(md.empty() == false)
        emit(out, "-" + ratePrefix + "_prior_mode", typeLabel.empty() ? md : typeLabel + ":" + md);
}

} // namespace

std::vector<std::string> ConfigReader::translate(const std::string& path){
    std::ifstream f(path);
    if(f.is_open() == false)
        Msg::error("could not open config file: " + path);

    std::map<std::string,std::string> global, subst, clock, lambda, mu;
    std::vector<std::pair<std::string, std::map<std::string,std::string>>> psiSections;
    std::map<std::string,std::string>* cur = &global;

    std::string line;
    while(std::getline(f, line)){
        size_t h = line.find('#');
        if(h != std::string::npos) line = line.substr(0, h);
        std::string t = trim(line);
        if(t.empty()) continue;
        if(t[0] == '['){
            size_t close = t.find(']');
            if(close == std::string::npos) Msg::error("config: malformed section header '" + t + "'");
            std::istringstream is(trim(t.substr(1, close - 1)));
            std::string type, name;
            is >> type; is >> name;
            if(type == "substitution") cur = &subst;
            else if(type == "clock") cur = &clock;
            else if(type == "lambda") cur = &lambda;
            else if(type == "mu") cur = &mu;
            else if(type == "psi"){ psiSections.push_back(std::make_pair(name, std::map<std::string,std::string>())); cur = &psiSections.back().second; }
            else Msg::error("config: unknown section '[" + type + "]'");
            continue;
        }
        size_t eq = t.find('=');
        if(eq == std::string::npos) Msg::error("config: expected 'key = value', got '" + t + "'");
        std::string key = trim(t.substr(0, eq));
        std::string val = trim(t.substr(eq + 1));
        if(val.empty()) continue;
        (*cur)[key] = val;
    }

    static const std::set<std::string> globalKeys = {
        "backbone_tree","fossils","clade_def","sequence","hessian","log_output","tree_output",
        "chain_length","parallel_chains","thinning","burn_in","coupled_chains","cores","seed",
        "max_gen","min_ess","rhat","conditioning","rho","delta_temperature","swap_interval"
    };
    static const std::set<std::string> substKeys = {
        "partition","datatype","n_states","ctmc_model","ctmc_gamma_cat","ctmc_inv","ctmc_freq"
    };
    static const std::set<std::string> clockKeys = {
        "clock_partitions","clock_model","rgene_gamma","sigma2_gamma","sigma2_param","pncp_tuning"
    };
    static const std::set<std::string> rateKeys = { "time_bins","prior","mode" };

    std::vector<std::string> out;
    for(std::map<std::string,std::string>::iterator it = global.begin(); it != global.end(); ++it){
        if(globalKeys.count(it->first) == 0) Msg::error("config: unknown key '" + it->first + "' in the global section");
        emit(out, "-" + it->first, it->second);
    }
    for(std::map<std::string,std::string>::iterator it = subst.begin(); it != subst.end(); ++it){
        if(substKeys.count(it->first) == 0) Msg::error("config: unknown key '" + it->first + "' in [substitution]");
        emit(out, "-" + it->first, it->second);
    }
    for(std::map<std::string,std::string>::iterator it = clock.begin(); it != clock.end(); ++it){
        if(clockKeys.count(it->first) == 0) Msg::error("config: unknown key '" + it->first + "' in [clock]");
        emit(out, "-" + it->first, it->second);
    }
    for(std::map<std::string,std::string>::iterator it = lambda.begin(); it != lambda.end(); ++it)
        if(rateKeys.count(it->first) == 0) Msg::error("config: unknown key '" + it->first + "' in [lambda]");
    for(std::map<std::string,std::string>::iterator it = mu.begin(); it != mu.end(); ++it)
        if(rateKeys.count(it->first) == 0) Msg::error("config: unknown key '" + it->first + "' in [mu]");
    std::vector<double> sectionTmins;
    double secTmin = -1.0;
    translateRateSection("lambda", "", lambda, out, secTmin); if(secTmin >= 0.0) sectionTmins.push_back(secTmin);
    translateRateSection("mu", "", mu, out, secTmin); if(secTmin >= 0.0) sectionTmins.push_back(secTmin);

    std::vector<std::string> psiTypeOrder;
    for(size_t s = 0; s < psiSections.size(); s++){
        std::map<std::string,std::string>& kv = psiSections[s].second;
        for(std::map<std::string,std::string>::iterator it = kv.begin(); it != kv.end(); ++it)
            if(rateKeys.count(it->first) == 0) Msg::error("config: unknown key '" + it->first + "' in [psi " + psiSections[s].first + "]");
        secTmin = -1.0;
        translateRateSection("psi", psiSections[s].first, kv, out, secTmin);
        if(secTmin >= 0.0) sectionTmins.push_back(secTmin);
        if(psiSections[s].first.empty() == false) psiTypeOrder.push_back(psiSections[s].first);
    }
    if(psiTypeOrder.empty() == false){
        std::string joined;
        for(size_t i = 0; i < psiTypeOrder.size(); i++){ if(i) joined += ","; joined += psiTypeOrder[i]; }
        emit(out, "-psi_types", joined);
    }
    if(sectionTmins.empty() == false){
        double g = sectionTmins[0];
        for(double t : sectionTmins)
            if(t != g)
                Msg::error("all rate sections with time_bins must share the same youngest edge; got " + fmtNum(g) + " and " + fmtNum(t) + ".");
        if(g > 0.0)
            emit(out, "-age_offset", fmtNum(g));
    }
    return out;
}
