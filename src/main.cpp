#include <iostream>
#include "Kadet.hpp"

int main()
{
    std::shared_ptr<Cantera::Solution> gas1 = Cantera::newSolution("airNASA9.yaml");
    std::shared_ptr<Cantera::Solution> gas2 = Cantera::newSolution("airNASA9.yaml");
    
    KaDet::shock::eq_calc(gas1,gas2,1e-3,1e-3);

    return 0;
}