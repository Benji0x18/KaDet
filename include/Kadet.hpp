/*
    This is a header of the Shock and Detonation Library KaDet. Whenever using KaDet make sure to include this header in the file.
    
*/
#ifndef KADET
#define KADET
#include <iostream>
#include <cmath>
#include "cantera/core.h"

namespace Ct = Cantera;

namespace KaDet
{
    
    namespace shock
    {
        int eqCalc(std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, double errT, double errV);
    }

}

#endif