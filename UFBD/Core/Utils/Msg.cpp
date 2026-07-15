#include <cstdlib>
#include <iostream>
#include <vector>
#include "Msg.hpp"

namespace {
    bool deferWarnings = false;
    std::vector<std::string> deferred;
}

void Msg::setDeferWarnings(bool b) {
    deferWarnings = b;
}

bool Msg::hasDeferredWarnings(void) {
    return deferred.empty() == false;
}

void Msg::flushWarnings(void) {
    for(const std::string& s : deferred)
        std::cout << "Warning: " << s << std::endl;
    deferred.clear();
}

void Msg::error(std::string s) {
    flushWarnings();
    std::cout << "Error: " << s << std::endl;
    std::cout << "Exiting program" << std::endl;
    std::exit(1);
}

void Msg::warning(std::string s) {
    if(deferWarnings){
        deferred.push_back(s);
        return;
    }
    std::cout << "Warning: " << s << std::endl;
}
