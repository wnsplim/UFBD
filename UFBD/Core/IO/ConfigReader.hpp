#ifndef ConfigReader_hpp
#define ConfigReader_hpp

#include <string>
#include <vector>

// struct BfConfig {
//     std::vector<std::string>        shared;
//     std::vector<std::string>        model[2];
//     std::string                     priorLog[2];
// };
// 
class ConfigReader {

    public:
        static std::vector<std::string> translate(const std::string& path);
//         static BfConfig                 translateBf(const std::string& path);
//         static std::string              sectionOf(const std::string& flag);
};

#endif
