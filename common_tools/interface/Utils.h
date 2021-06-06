#ifndef Utils_GUARD
#define Utils_GUARD

#include "TH1.h"
#include "TH2.h"
#include <algorithm>
#include <sstream>
#include <string>

// book 1D histogram -------------------------------------------------------------------------------
TH1D* Book1D(std::string name, unsigned int nbins, double min, double max);

// book 2D histogram -------------------------------------------------------------------------------
TH2D* Book2D(std::string name, unsigned int nbinsX, double minX, double maxX,
             unsigned int nbinsY, double minY, double maxY);

// Convert mostly anything to string ---------------------------------------------------------------
template<typename T> std::string ToStr(T var) {
    std::ostringstream oss;
    oss << var;
    return oss.str();
}

// Order by pT -------------------------------------------------------------------------------------
template<typename T> void OrderByPt(std::vector<const T*> &vec) {
    if (vec.size() < 2) return;
    for (unsigned int i = 0; i != vec.size()-1; ++i) {
        for (unsigned int j = i+1; j != vec.size(); ++j) {
            if (vec.at(i)->p4.Pt() < vec.at(j)->p4.Pt()) {
                std::swap(vec.at(i), vec.at(j));
            }
        }
    }
}

// Percentage bar ----------------------------------------------------------------------------------
void ShowPercentage(unsigned int p);

#endif
