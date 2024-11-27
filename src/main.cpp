#include <iostream>
#include "Kadet.hpp"

namespace Ct = Cantera;
namespace Kd = Kadet;


int demo_sh()
{
    std::shared_ptr<Ct::Solution> gas1 = Ct::newSolution("airNASA9.yaml");
    std::shared_ptr<Ct::Solution> gas2 = Ct::newSolution("airNASA9.yaml");

    std::string comp = "N2:78,O2:21,AR:1";
    double W = 2000;

    gas1->thermo()->setState_TP(250,40e+3);
    gas1->thermo()->setMassFractionsByName(comp);

    Kd::shock::shockEq(W,gas1,gas2,1e-4,1e-4);

    return 0;
}

int main()
{
    demo_sh();
    return 0;
}