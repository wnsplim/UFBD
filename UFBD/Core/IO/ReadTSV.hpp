#ifndef ReadTSV_hpp
#define ReadTSV_hpp

#include <string>
#include <vector>

class ReadTSV{
    public:
                                                ReadTSV(std::string filepath, bool rownames, bool colnames, bool string);
        std::vector<std::vector<std::string>>   getReadStringData(void){return stringData;}
    private:
        std::vector<std::vector<std::string>>   stringData;
};

#endif
