#include "Tree.hpp"

int main() {

    std::vector<std::string> taxonNames;
    for(int i = 0; i < 10; i++)
        taxonNames.push_back("t" + std::to_string(i));
        
    Tree pt = Tree(taxonNames, 10.0);
    pt.print();
    
    return 0;
}
