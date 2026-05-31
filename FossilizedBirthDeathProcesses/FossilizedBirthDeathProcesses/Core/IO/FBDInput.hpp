#ifndef FBDInput_hpp
#define FBDInput_hpp

#include <string>
#include <vector>

class Node;
class Tree;

enum class Assignment { CROWN, TOTAL };

class Clade {

    public:
                                            Clade(std::string n, std::vector<std::string> t, Node* c, Node* o) : name(n), taxa(t), crown(c), origin(o) {}
        std::string                         getName(void) { return name; }
        std::vector<std::string>&           getTaxa(void) { return taxa; }
        Node*                               getCrown(void) { return crown; }
        Node*                               getOrigin(void) { return origin; }

    private:
        std::string                         name;
        std::vector<std::string>            taxa;
        Node*                               crown;
        Node*                               origin;
};

class Fossil {

    public:
                                            Fossil(std::string t, double mn, double mx, std::string c, Assignment a) : taxon(t), minAge(mn), maxAge(mx), clade(c), assignment(a) {}
        std::string                         getTaxon(void) { return taxon; }
        double                              getMinAge(void) { return minAge; }
        double                              getMaxAge(void) { return maxAge; }
        std::string                         getClade(void) { return clade; }
        Assignment                          getAssignment(void) { return assignment; }

    private:
        std::string                         taxon;
        double                              minAge;
        double                              maxAge;
        std::string                         clade;
        Assignment                          assignment;
};

class FBDInput {

    public:
                                            FBDInput(std::string treePath, std::string cladesPath, std::string fossilPath);
        Tree*                               getTree(void) { return tree; }
        std::vector<Clade>&                 getClades(void) { return clades; }
        std::vector<Fossil>&                getFossils(void) { return fossils; }

    private:
        Tree*                               readTree(std::string path);
        void                                readClades(std::string path);
        void                                readFossils(std::string path);
        void                                assignFossilAwareAges(void);
        Tree*                               tree;
        std::vector<Clade>                  clades;
        std::vector<Fossil>                 fossils;
};

#endif
