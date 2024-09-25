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
        std::array<double,6> LSQ_speedCJ(const std::vector<double> &x, const std::vector<double> &y);

        double hugFr(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Ct::Solution> gas);

        double hugEq(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Ct::Solution> gas);

        double calcCJ(std::shared_ptr<Ct::Solution> gas2, std::shared_ptr<Ct::Solution> gas1, const double &errV, const double &errT, const double &x);

        double speedCJ(const double &P1, const double &T1, const double &q, const std::string &mech);

        void FHFP(const double& vel, std::shared_ptr<Ct::Solution> gas2, std::shared_ptr<Ct::Solution> gas1 /*gas1*/,double &FH,double &FP);

        int eqCalc(const double &vel, std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, double errT, double errV);
    }

}

#endif