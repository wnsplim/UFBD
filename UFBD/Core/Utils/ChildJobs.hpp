#ifndef ChildJobs_hpp
#define ChildJobs_hpp

#include <string>
#include <vector>

struct ChildJob {
    std::vector<std::string> args;
    std::string             label;
    std::string             console;
    bool                    stdoutToConsole = false;
};

void runChildJobs(const std::vector<ChildJob>& jobs, int width, const std::string& what);

#endif
