#include "Msg.hpp"
#include "WriteTSV.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>

WriteTSV::WriteTSV(void){
}

WriteTSV::WriteTSV(std::string filepath, bool overwrite) : filepath(filepath), numRows(0), numCols(0){
    if(overwrite == true)
        std::remove(filepath.c_str());
    fout.open(filepath, std::ios::out | std::ios::app);
}

WriteTSV::~WriteTSV(void){
    fout.close();
}

void WriteTSV::addColumnNamesTSV(std::vector<std::string> cn){
    if(numRows == 0){
        for(int i = 0; i < cn.size(); i++){
            if(i < (cn.size() - 1))
                fout << cn[i] << "\t";
            else
                fout << cn[i];
        }
        fout << "\n";
        numRows++;
    }else{
        Msg::error("adding column names after data has alread been entered");
    }
    fout.flush();
}

void WriteTSV::addFilepath(std::string fp, bool overwrite){
    if (fout.is_open()) {
        fout.close();
    }
    filepath = fp;
    if(overwrite == true)
        std::remove(filepath.c_str());
    fout.open(filepath, std::ios::out | std::ios::app);
    numRows = 0;
    numCols = 0;
    if (!fout.is_open())
        Msg::error("File stream is not open: " + filepath);
}

void WriteTSV::appendDataTSV(std::vector<double> data){
    if (!fout.is_open())
        Msg::error("Attempting to write to TSV file before opening it.");

    for(int i = 0; i < data.size(); i++){
        if(i < (data.size()-1))
            fout << data[i] << "\t";
        else
            fout << data[i];
    }
    fout << "\n";
    numRows++;
    fout.flush();
}

void WriteTSV::appendDataTSV(std::string data){
    fout << data;
    fout << "\n";
    numRows++;
    fout.flush();
}

void WriteTSV::closeTSV(void){
    fout.close();
}
