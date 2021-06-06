#include "MyAnalyses/ZprimeSearch/interface/SelectEvent.h"
#include "FWCore/Framework/interface/MakerMacros.h"

// constructor -------------------------------------------------------------------------------------
SelectEvent::SelectEvent(const edm::ParameterSet &config) {
    runNrs = config.getParameter< std::vector<unsigned int> >("runNrs");
    lumiNrs = config.getParameter< std::vector<unsigned int> >("lumiNrs");
    eventNrs = config.getParameter< std::vector<unsigned int> >("eventNrs");
}

// destructor --------------------------------------------------------------------------------------
SelectEvent::~SelectEvent() {
    
}

// filter event ------------------------------------------------------------------------------------
bool SelectEvent::filter(edm::Event &event, const edm::EventSetup &setup) {
    for (unsigned int i = 0; i != runNrs.size(); ++i) {
        if (event.id().run() == runNrs.at(i) &&
            event.luminosityBlock() == lumiNrs.at(i) &&
            event.id().event() == eventNrs.at(i)) return true;
    }
    return false;
}

DEFINE_FWK_MODULE(SelectEvent);
