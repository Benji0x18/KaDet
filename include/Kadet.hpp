/*
    This is a header of the Shock and Detonation Library KaDet. Whenever using KaDet make sure to include this header in the file.
    
*/
#ifndef KADET
#define KADET
#include <iostream>
#include <cmath>
#include "cantera/core.h"

namespace Ct = Cantera;

namespace Kadet
{

    
    class station
    {
        public:
            // Variables
            double                  T;
            double                  T_0;
            double                  P;
            double                  P_0;
            double                  D;      // Density
            double                  H;
            double                  H_0;
            double                  S;
            double                  cp;
            double                  ga;
            double                  la;
            double                  nu;
            int                     N;      // Species Number
            std::vector<double>     X;      // Mole Fractions
            std::vector<double>     Y;      // Mass Fractions
            double                  diff;   // Diffusion
            double                  U;      // Absolute velocity
            double                  Ma;     // Mach number
            double                  Pr;     // Prandtl number
            double                  Le;     // Lewis number

            // Methods
            station();
            ~station();
            /*!
            * Saves the gas properties and computes the stagnation propaerties based on velocity
            * @param gas Cantera solution object properties to be saved
            * @param Mach Speed of the gas at a given station
            */
            void storePropsU(std::shared_ptr<Ct::Solution> gas,  double velocity);
            void storePropsMa(std::shared_ptr<Ct::Solution> gas,  double Mach);

            void storePlasMa(std::shared_ptr<Ct::Solution> gas,  double Mach);
            void storePlasU(std::shared_ptr<Ct::Solution> gas,  double velocity);

            void setGas(std::shared_ptr<Ct::Solution> gas);
    };
    
    namespace shock
    {
        std::array<double,6> LSQ_speedCJ(const std::vector<double> &x, const std::vector<double> &y);

        double frHug(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Ct::Solution> gas);

        double eqHug(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Ct::Solution> gas);

        double calcCJ(std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, const double &errV, const double &errT, const double &x);

        double speedCJ(const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        std::array<std::shared_ptr<Ct::Solution>, 2> postShockFr(const double &U1, const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        std::array<std::shared_ptr<Ct::Solution>, 2> postShockEq(const double &U1, const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        void FHFP(const double &U, std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, double &FH, double &FP);

        int calcFr(const double &U, std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, double errT, double errV);

        int calcEq(const double &U, std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, double errT, double errV);

        double reflFr(std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, std::shared_ptr<Ct::Solution> gas3, const double &UI);

        double reflEq(std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Ct::Solution> gas2, std::shared_ptr<Ct::Solution> gas3, const double &UI);

        std::shared_ptr<Ct::Solution> postReflFr(const double &U2, const std::string &mech);

        std::shared_ptr<Ct::Solution> postReflEq(const double &U2, const std::string &mech);

        void FHFP_Refl(const double &U2, std::shared_ptr<Ct::Solution> gas2, std::shared_ptr<Ct::Solution> gas3, double &FH, double &FP);

    }

    namespace utils
    {
        double getThermicity(std::shared_ptr<Ct::Solution> gas);
        std::array<std::vector<double>,2> RK45();
    }

    namespace stagnation
    {
        class system
        {
            private:

            public:

        };
    }

    namespace detonation
    {
        void constP(std::shared_ptr<Ct::Solution> gas, double t, std::vector<double> y);

        void constV(std::shared_ptr<Ct::Solution> gas, double t, std::vector<double> y);

        void ZND(std::shared_ptr<Ct::Solution> gas, double t, std::vector<double> y);
    }



}

#endif