// #ifndef BridgeSampler_hpp
// #define BridgeSampler_hpp
// 
// #include "ChildJobs.hpp"
// 
// #include <string>
// #include <vector>
// 
// class BridgeSampler {
// 
//     public:
//                                     BridgeSampler(void) = delete;
//                                     BridgeSampler(const std::string& base);
//         void                        run(void);
//         static unsigned long        pilotGenerations(void);
// 
//     private:
//         typedef ChildJob Job;
//         Job                        makeJob(const std::vector<std::string>& extra, const std::string& logPath, const std::string& label);
//         void                        runJobs(const std::vector<Job>& jobs);
//         std::vector<Job>            pointJobs(const std::vector<double>& pts, const std::vector<std::string>& text, bool pilot);
//         std::vector<double>         readPoint(double t, const std::string& text, bool pilot);
//         std::vector<double>         readDelta(void);
//         void                        pilot(void);
//         std::vector<double>         placePoints(void);
//         double                      chained(const std::vector<std::vector<double>>& d);
//         void                        estimate(void);
//         void                        writeTable(void);
//         void                        print(void);
//         std::string                 base;
//         std::vector<std::string>    childArgs;
//         bool                        resumeAll;
//         bool                        perChain;
//         double                      burn;
//         std::vector<double>         nodes, weights, nodeSd;
//         std::vector<size_t>         nodeN;
//         double                      length = 0.0;
//         std::vector<double>         ts;
//         std::vector<std::vector<double>> deltas;
//         double                      lnBf = 0.0, se = 0.0, maxJ = 0.0;
// };
// 
// #endif
