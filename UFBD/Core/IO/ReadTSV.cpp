#include "Msg.hpp"
#include "ReadTSV.hpp"

#include <fstream>
#include <sstream>

ReadTSV::ReadTSV(std::string filepath, bool, bool, bool){
    std::ifstream file(filepath);
    if(file.is_open() == false)
        Msg::error("cannot open input file '" + filepath + "'");

    std::string line;
    while(std::getline(file, line)){
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;
        while(std::getline(ss, cell, '\t'))
            row.push_back(cell);
        if(row.size() > 0)
            stringData.push_back(row);
    }
    file.close();
}
