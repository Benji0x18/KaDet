#include <iostream>
#include "Kadet.hpp"

namespace Kd = Kadet;
namespace Ct = Cantera;

int main()
{
    std::shared_ptr<Ct::Solution> gas1 = Ct::newSolution("airNASA9.yaml");
    std::shared_ptr<Ct::Solution> gas2 = Ct::newSolution("airNASA9.yaml");
    
    Kd::shock::calcEq(2000,gas1,gas2,1e-3,1e-3);

    return 0;
}