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

        double frHug(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Ct::Solution> gas);

        double eqHug(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Ct::Solution> gas);

        double calcCJ(std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, const double &errV, const double &errT, const double &x);

        double speedCJ(const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        std::array<std::shared_ptr<Ct::Solution>, 2> frPost(const double &U1, const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        std::array<std::shared_ptr<Ct::Solution>, 2> eqPost(const double &U1, const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        void FHFP(const double &U, std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, double &FH, double &FP);

        int frCalc(const double &U, std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, double errT, double errV);

        int eqCalc(const double &U, std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, double errT, double errV);

        double frRefl(std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, std::shared_ptr<Ct::Solution> gas3,const double &UI);

        double eqRefl(std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, std::shared_ptr<Ct::Solution> gas3,const double &UI);

        std::shared_ptr<Ct::Solution> frPostRefl(const double &U2, const std::string &mech);

        std::shared_ptr<Ct::Solution> eqPostRefl(const double &U2, const std::string &mech);

        void FHFP_Refl(const double &U2, std::shared_ptr<Ct::Solution> gas2, std::shared_ptr<Ct::Solution> gas3, double &FH, double &FP);

    }

    namespace utils
    {
        double getThermicity(std::shared_ptr<Ct::Solution> gas);
    }



}

#endif