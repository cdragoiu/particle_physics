#ifndef SelectEvent_GUARD
#define SelectEvent_GUARD

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <vector>

class SelectEvent : public edm::EDFilter {
    
public:
    
    // constructor
	SelectEvent(const edm::ParameterSet &config);
    
    // destructor
	~SelectEvent();
    
private:
    
    // filter event
	bool filter(edm::Event &event, const edm::EventSetup &setup);
    
    // input variables
    std::vector<unsigned int> runNrs, lumiNrs, eventNrs;
    
};

#endif
