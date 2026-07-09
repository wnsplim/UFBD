#ifndef ConfigReader_hpp
#define ConfigReader_hpp

#include <string>
#include <vector>

class ConfigReader {

    public:
        static std::vector<std::string> translate(const std::string& path);
};

#endif
