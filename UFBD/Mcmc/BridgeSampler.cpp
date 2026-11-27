// #include "BridgeSampler.hpp"
// #include "Msg.hpp"
// #include "PathModel.hpp"
// #include "UserSettings.hpp"
// 
// #include <algorithm>
// #include <cmath>
// #include <fstream>
// #include <iomanip>
// #include <iostream>
// #include <sstream>
// 
// static std::string pointText(double t){
//     std::ostringstream os;
//     os << std::setprecision(6) << t;
//     return os.str();
// }
// 
// static double mean(const std::vector<double>& x){
//     double m = 0.0;
//     for(double y : x) m += y / (double)x.size();
//     return m;
// }
// 
// static double variance(const std::vector<double>& x){
//     double m = mean(x), v = 0.0;
//     for(double y : x) v += (y - m) * (y - m) / (double)(x.size() - 1);
//     return v;
// }
// 
// static void gaussLegendre(int n, std::vector<double>& t, std::vector<double>& w){
//     t.assign(n, 0.0);
//     w.assign(n, 0.0);
//     for(int i = 0; i < n; i++){
//         double x = std::cos(M_PI * ((double)i + 0.75) / ((double)n + 0.5)), dp = 0.0, dx = 0.0;
//         do {
//             double p0 = 1.0, p1 = x;
//             for(int k = 2; k <= n; k++){
//                 double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / (double)k;
//                 p0 = p1;
//                 p1 = p2;
//             }
//             dp = (double)n * (x * p1 - p0) / (x * x - 1.0);
//             dx = p1 / dp;
//             x -= dx;
//         } while(std::abs(dx) > 1e-14);
//         t[n - 1 - i] = 0.5 * (1.0 + x);
//         w[n - 1 - i] = 1.0 / ((1.0 - x * x) * dp * dp);
//     }
// }
// 
// static double logMeanExp(const std::vector<double>& v){
//     double c = *std::max_element(v.begin(), v.end()), s = 0.0;
//     for(double x : v) s += std::exp(x - c);
//     return c + std::log(s / (double)v.size());
// }
// 
// static double logAddExp(double x, double y){
//     double m = std::max(x, y);
//     return m + std::log(std::exp(x - m) + std::exp(y - m));
// }
// 
// static double bridgeLogR(const std::vector<double>& lx, const std::vector<double>& ly){
//     double n1 = (double)lx.size(), n2 = (double)ly.size();
//     double ls1 = std::log(n1 / (n1 + n2)), ls2 = std::log(n2 / (n1 + n2)), rho = 0.0;
//     std::vector<double> a(lx.size()), c(ly.size());
//     UserSettings& us = UserSettings::userSettings();
//     for(int it = 0; it < us.getBfBridgeIter(); it++){
//         for(size_t i = 0; i < lx.size(); i++) a[i] = lx[i] - logAddExp(ls1 + rho, ls2 + lx[i]);
//         for(size_t i = 0; i < ly.size(); i++) c[i] = -logAddExp(ls1 + rho, ls2 + ly[i]);
//         double next = logMeanExp(a) - logMeanExp(c);
//         if(std::isfinite(next) == false)
//             Msg::error("bridge sampling: the bridge iteration gave a non-finite value.");
//         if(std::abs(next - rho) < us.getBfBridgeTol())
//             return next;
//         rho = next;
//     }
//     Msg::error("bridge sampling: the bridge iteration did not reach -bf_bridge_tol " + std::to_string(us.getBfBridgeTol()) + " in -bf_bridge_iter " + std::to_string(us.getBfBridgeIter()) + " steps.");
//     return 0.0;
// }
// 
// BridgeSampler::BridgeSampler(const std::string& base) : base(base) {
//     UserSettings& s = UserSettings::userSettings();
//     resumeAll = s.getResume();
//     perChain = s.getNumRuns() > 1;
//     burn = perChain ? 0.0 : s.getBurninFraction();
//     const std::vector<std::string>& cli = s.getCliArgs();
//     childArgs = {cli[0], "-config_bf", s.getBfConfigPath()};
//     for(size_t i = 1; i < cli.size(); i++){
//         if(cli[i] == "-bf_t" || cli[i] == "-bf_posterior"){ i++; continue; }
//         if(cli[i] == "-bf_pilot" || cli[i] == "-resume") continue;
//         childArgs.push_back(cli[i]);
//     }
// }
// 
// unsigned long BridgeSampler::pilotGenerations(void){
//     UserSettings& s = UserSettings::userSettings();
//     double gens = 0.0;
//     for(const std::string& b : {s.getPosteriorBaseA(), s.getPosteriorBaseB()})
//         gens = std::max(gens, PathModel::readLastRow(std::ifstream(b + "_chain0.log").good() ? b + "_chain0.log" : b + ".log")["n"]);
//     return (unsigned long)std::ceil(s.getBfPilotFrac() * gens);
// }
// 
// BridgeSampler::Job BridgeSampler::makeJob(const std::vector<std::string>& extra, const std::string& logPath, const std::string& label){
//     Job j{childArgs, label, logPath.substr(0, logPath.size() - 4) + ".console.txt"};
//     j.args.insert(j.args.end(), extra.begin(), extra.end());
//     std::string ckp = perChain ? logPath.substr(0, logPath.size() - 4) + "_chain0.log.ckp" : logPath + ".ckp";
//     if(resumeAll && std::ifstream(ckp).good())
//         j.args.push_back("-resume");
//     return j;
// }
// 
// void BridgeSampler::runJobs(const std::vector<Job>& jobs){
//     runChildJobs(jobs, UserSettings::userSettings().getBfParallelRuns(), "bridge sampling");
// }
// 
// std::vector<double> BridgeSampler::readPoint(double t, const std::string& text, bool pilot){
//     UserSettings& s = UserSettings::userSettings();
//     s.setPathPoint(t, text, pilot);
//     std::vector<double> d = readDelta();
//     int nb = s.getBfBatches();
//     if(pilot && d.size() < 2)
//         Msg::error("bridge sampling: " + s.getParamOutput() + " has " + std::to_string(d.size()) + " samples of delta after burn-in; the pilot standard deviation needs at least 2.");
//     if(pilot == false && d.size() < (size_t)nb)
//         Msg::error("bridge sampling: " + s.getParamOutput() + " has " + std::to_string(d.size()) + " samples of delta after burn-in; -bf_batches " + std::to_string(nb) + " needs at least " + std::to_string(nb) + ".");
//     return d;
// }
// 
// std::vector<BridgeSampler::Job> BridgeSampler::pointJobs(const std::vector<double>& pts, const std::vector<std::string>& text, bool pilot){
//     UserSettings& s = UserSettings::userSettings();
//     std::vector<Job> jobs;
//     for(size_t i = 0; i < pts.size(); i++){
//         s.setPathPoint(pts[i], text[i], pilot);
//         std::vector<std::string> extra = {"-bf_t", text[i]};
//         if(pilot)
//             extra.push_back("-bf_pilot");
//         jobs.push_back(makeJob(extra, s.getParamOutput(), (pilot ? "pilot run at t = " : "path point t = ") + text[i]));
//     }
//     return jobs;
// }
// 
// void BridgeSampler::run(void){
//     UserSettings& s = UserSettings::userSettings();
//     std::vector<Job> posteriors;
//     for(int k = 0; k < 2; k++){
//         std::string b = (k == 0) ? s.getPosteriorBaseA() : s.getPosteriorBaseB();
//         if(s.getPriorLogGiven(k)){
//             if(std::ifstream(b + ".log").good() == false)
//                 Msg::error("-config_bf: prior_log of [model " + std::to_string(k) + "] is " + b + ", but " + b + ".log does not exist.");
//             continue;
//         }
//         posteriors.push_back(makeJob({"-bf_posterior", std::to_string(k)}, b + ".log", "posterior of [model " + std::to_string(k) + "]"));
//     }
//     std::cout << "=======================================================================\n";
//     runJobs(posteriors);
//     std::vector<double> pts = s.getBfLadder();
//     std::vector<std::string> text = s.getBfLadderText();
//     if(pts.empty()){
//         pilot();
//         pts = placePoints();
//         for(double t : pts)
//             text.push_back(pointText(t));
//     }
//     ts = pts;
//     std::cout << "=======================================================================\n";
//     runJobs(pointJobs(pts, text, false));
//     for(size_t i = 0; i < pts.size(); i++)
//         deltas.push_back(readPoint(pts[i], text[i], false));
//     estimate();
//     writeTable();
//     print();
// }
// 
// std::vector<double> BridgeSampler::readDelta(void){
//     std::string path = UserSettings::userSettings().getParamOutput();
//     std::vector<double> d;
//     for(const std::vector<double>& r : PathModel::readColumns(path, {"delta"}, "path point log " + path))
//         d.push_back(r[0]);
//     d.erase(d.begin(), d.begin() + (size_t)(burn * (double)d.size()));
//     return d;
// }
// 
// void BridgeSampler::pilot(void){
//     UserSettings& s = UserSettings::userSettings();
//     gaussLegendre(s.getBfPilotNodes(), nodes, weights);
//     std::vector<std::string> text;
//     for(double& t : nodes){
//         t = std::stod(pointText(t));
//         text.push_back(pointText(t));
//     }
//     std::cout << "=======================================================================\n";
//     runJobs(pointJobs(nodes, text, true));
//     for(size_t i = 0; i < nodes.size(); i++){
//         std::vector<double> d = readPoint(nodes[i], text[i], true);
//         nodeSd.push_back(std::sqrt(variance(d)));
//         nodeN.push_back(d.size());
//     }
//     for(size_t i = 0; i < nodes.size(); i++)
//         length += weights[i] * nodeSd[i];
//     std::cout << "=======================================================================\n";
//     std::cout << "Pilot over " << nodes.size() << " Gauss-Legendre nodes, " << pilotGenerations() << " generations each (" << s.getBfPilotFrac()
//               << " of the generations of the longer posterior run)\n";
//     std::cout << std::setw(10) << "t" << std::setw(12) << "weight" << std::setw(8) << "n" << std::setw(14) << "sd_delta" << "\n";
//     for(size_t i = 0; i < nodes.size(); i++)
//         std::cout << std::setw(10) << pointText(nodes[i]) << std::setw(12) << std::setprecision(6) << weights[i] << std::setw(8) << nodeN[i]
//                   << std::setw(14) << std::setprecision(6) << nodeSd[i] << "\n";
//     std::cout << "thermodynamic length of the path = " << length << "\n" << std::defaultfloat;
//     if(length <= 0.0)
//         Msg::error("bridge sampling: delta does not vary at any pilot node, so the path has no length to place points along.");
// }
// 
// std::vector<double> BridgeSampler::placePoints(void){
//     std::vector<double> tt = {0.0}, sd = {nodeSd.front()};
//     tt.insert(tt.end(), nodes.begin(), nodes.end());
//     sd.insert(sd.end(), nodeSd.begin(), nodeSd.end());
//     tt.push_back(1.0);
//     sd.push_back(nodeSd.back());
//     std::vector<double> cum = {0.0};
//     for(size_t i = 1; i < tt.size(); i++)
//         cum.push_back(cum.back() + 0.5 * (tt[i] - tt[i - 1]) * (sd[i] + sd[i - 1]));
//     double scale = length / cum.back();
//     for(double& c : cum) c *= scale;
//     int k = (int)std::ceil(UserSettings::userSettings().getBfIntervals() * length) + 1;
//     std::vector<double> pts = {0.0};
//     size_t i = 0;
//     for(int j = 1; j + 1 < k; j++){
//         double target = length * (double)j / (double)(k - 1);
//         while(cum[i + 1] < target) i++;
//         double t = tt[i] + (target - cum[i]) / (cum[i + 1] - cum[i]) * (tt[i + 1] - tt[i]);
//         pts.push_back(std::stod(pointText(t)));
//     }
//     pts.push_back(1.0);
//     std::cout << k << " path points equidistant in thermodynamic length (" << UserSettings::userSettings().getBfIntervals() << " per unit length)\n";
//     return pts;
// }
// 
// double BridgeSampler::chained(const std::vector<std::vector<double>>& d){
//     double r = 0.0;
//     for(size_t j = 0; j + 1 < ts.size(); j++){
//         double h = ts[j + 1] - ts[j];
//         std::vector<double> lx(d[j].size()), ly(d[j + 1].size());
//         for(size_t i = 0; i < lx.size(); i++) lx[i] = h * d[j][i];
//         for(size_t i = 0; i < ly.size(); i++) ly[i] = h * d[j + 1][i];
//         r += bridgeLogR(lx, ly);
//     }
//     return r;
// }
// 
// void BridgeSampler::estimate(void){
//     lnBf = chained(deltas);
//     int nb = UserSettings::userSettings().getBfBatches();
//     std::vector<double> parts;
//     for(int b = 0; b < nb; b++){
//         std::vector<std::vector<double>> batch;
//         for(const std::vector<double>& d : deltas){
//             size_t q = d.size() / (size_t)nb, r = d.size() % (size_t)nb;
//             size_t from = (size_t)b * q + std::min((size_t)b, r);
//             batch.push_back(std::vector<double>(d.begin() + (long)from, d.begin() + (long)(from + q + ((size_t)b < r ? 1 : 0))));
//         }
//         parts.push_back(chained(batch));
//     }
//     se = std::sqrt(variance(parts) / (double)nb);
//     maxJ = 0.0;
//     for(size_t j = 0; j + 1 < ts.size(); j++)
//         maxJ = std::max(maxJ, (ts[j + 1] - ts[j]) * (mean(deltas[j + 1]) - mean(deltas[j])));
// }
// 
// void BridgeSampler::writeTable(void){
//     std::ofstream f(base + "_bf.tsv");
//     f << "t\tn\tmean_delta\tvar_delta\n" << std::setprecision(17);
//     for(size_t k = 0; k < ts.size(); k++)
//         f << ts[k] << "\t" << deltas[k].size() << "\t" << mean(deltas[k]) << "\t" << variance(deltas[k]) << "\n";
// }
// 
// void BridgeSampler::print(void){
//     int nb = UserSettings::userSettings().getBfBatches();
//     std::cout << "=======================================================================\n";
//     std::cout << "Bridge sampling over " << ts.size() << " path points\n";
//     std::cout << std::setw(10) << "t" << std::setw(8) << "n" << std::setw(16) << "mean_delta" << std::setw(16) << "var_delta" << "\n";
//     for(size_t k = 0; k < ts.size(); k++)
//         std::cout << std::setw(10) << pointText(ts[k]) << std::setw(8) << deltas[k].size() << std::setw(16) << std::fixed << std::setprecision(6) << mean(deltas[k])
//                   << std::setw(16) << std::setprecision(4) << variance(deltas[k]) << "\n" << std::defaultfloat;
//     std::cout << std::fixed << std::setprecision(6)
//               << "ln BF ([model 1] over [model 0]) = " << lnBf << " (SE " << se << ", batch means over " << nb << " batches per run)\n"
//               << "largest Jeffreys divergence between neighboring path points = " << maxJ << "\n" << std::defaultfloat;
// }
