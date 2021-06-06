#include "MyAnalyses/CommonTools/interface/Utils.h"
#include <iostream>

// book 1D histogram -------------------------------------------------------------------------------
TH1D* Book1D(std::string name, unsigned int nbins, double min, double max) {
    TH1D *hist = new TH1D(name.c_str(), "", nbins, min, max);
    hist->Sumw2();
    return hist;
}

// book 2D histogram -------------------------------------------------------------------------------
TH2D* Book2D(std::string name, unsigned int nbinsX, double minX, double maxX,
             unsigned int nbinsY, double minY, double maxY) {
    TH2D *hist = new TH2D(name.c_str(), "", nbinsX, minX, maxX, nbinsY, minY, maxY);
    hist->Sumw2();
    return hist;
}

// Percentage bar ----------------------------------------------------------------------------------
void ShowPercentage(unsigned int p) {
    if (p == 0) {
        std::cout << std::endl;
        std::cout << "0% \033[0;47m                    \033[0m 100%";
        std::cout << std::flush;
    }
    else if (p == 10) {
        for (unsigned int i = 0; i != 28; ++i) std::cout << "\b";
        std::cout << "0% \033[0;42m  \033[0m\033[0;47m                  \033[0m 100%";
        std::cout << std::flush;
    }
    else if (p == 20) {
        for (unsigned int i = 0; i != 28; ++i) std::cout << "\b";
        std::cout << "0% \033[0;42m    \033[0m\033[0;47m                \033[0m 100%";
        std::cout << std::flush;
    }
    else if (p == 30) {
        for (unsigned int i = 0; i != 28; ++i) std::cout << "\b";
        std::cout << "0% \033[0;42m      \033[0m\033[0;47m              \033[0m 100%";
        std::cout << std::flush;
    }
    else if (p == 40) {
        for (unsigned int i = 0; i != 28; ++i) std::cout << "\b";
        std::cout << "0% \033[0;42m        \033[0m\033[0;47m            \033[0m 100%";
        std::cout << std::flush;
    }
    else if (p == 50) {
        for (unsigned int i = 0; i != 28; ++i) std::cout << "\b";
        std::cout << "0% \033[0;42m          \033[0m\033[0;47m          \033[0m 100%";
        std::cout << std::flush;
    }
    else if (p == 60) {
        for (unsigned int i = 0; i != 28; ++i) std::cout << "\b";
        std::cout << "0% \033[0;42m            \033[0m\033[0;47m        \033[0m 100%";
        std::cout << std::flush;
    }
    else if (p == 70) {
        for (unsigned int i = 0; i != 28; ++i) std::cout << "\b";
        std::cout << "0% \033[0;42m              \033[0m\033[0;47m      \033[0m 100%";
        std::cout << std::flush;
    }
    else if (p == 80) {
        for (unsigned int i = 0; i != 28; ++i) std::cout << "\b";
        std::cout << "0% \033[0;42m                \033[0m\033[0;47m    \033[0m 100%";
        std::cout << std::flush;
    }
    else if (p == 90) {
        for (unsigned int i = 0; i != 28; ++i) std::cout << "\b";
        std::cout << "0% \033[0;42m                  \033[0m\033[0;47m  \033[0m 100%";
        std::cout << std::flush;
    }
    else if (p == 100) {
        for (unsigned int i = 0; i != 28; ++i) std::cout << "\b";
        std::cout << "0% \033[0;42m                    \033[0m 100%";
        std::cout << std::endl;
    }
    return;
}
