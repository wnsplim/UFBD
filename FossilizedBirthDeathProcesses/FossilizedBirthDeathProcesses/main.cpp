#include "FBDTreeModel.hpp"
#include "Mcmc.hpp"
#include "MetropolisCoupledMcmc.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#include <iostream>

int main(int argc, const char* argv[]) {

    UserSettings& settings = UserSettings::userSettings();
    settings.initializeSettings(argc, argv);
    settings.print();

    std::vector<std::string> taxonNames;
    for(int i = 0; i < 10; i++)
        taxonNames.push_back("t" + std::to_string(i));
        
    Tree pt = Tree(taxonNames, 10.0);
    pt.print();
    
    int numChains = settings.getNumChains();
    if(numChains > 1){
        std::cout << "Running Metropolis-coupled MCMC with " << numChains << " chains parallelized across " << settings.getNumThreads() << " threads \n";
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::vector<PhylogeneticModel*> models;
        models.resize(numChains);
        for(int i = 0; i < numChains; i++)
            models[i] = new FBDTreeModel(&pt);
        MetropolisCoupledMcmc mcmc(settings.getChainLength(), settings.getPrintFrequency(), settings.getSampleFrequency(), models);
        mcmc.run();
    }else if (numChains == 1){
        std::cout << "Running standard MCMC \n";
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        PhylogeneticModel* model = new FBDTreeModel(&pt);
        Mcmc mcmc(settings.getChainLength(), settings.getPrintFrequency(), settings.getSampleFrequency(), model);
        mcmc.run();
    }

    return 0;
}
