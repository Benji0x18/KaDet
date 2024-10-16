/*
    This is a header of the Shock and Detonation Library KaDet. Whenever using KaDet make sure to include this header in the file.
    
*/
#ifndef KADET
#define KADET
#include <iostream>
#include <cmath>
#include "cantera/core.h"

namespace Kadet
{
    
    namespace shock
    {
        std::array<double,6> LSQ_speedCJ(const std::vector<double> &x, const std::vector<double> &y);

        double frHug(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Cantera::Solution> gas);

        double eqHug(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Cantera::Solution> gas);

        double calcCJ(std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, const double &errV, const double &errT, const double &x);

        double speedCJ(const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        std::array<std::shared_ptr<Cantera::Solution>, 2> postShockFr(const double &U1, const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        std::array<std::shared_ptr<Cantera::Solution>, 2> postShockEq(const double &U1, const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        void FHFP(const double &U, std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, double &FH, double &FP);

        int calcFr(const double &U, std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, double errT, double errV);

        int calcEq(const double &U, std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, double errT, double errV);

        double reflFr(std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, std::shared_ptr<Cantera::Solution> gas3, const double &UI);

        double reflEq(std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, std::shared_ptr<Cantera::Solution> gas3, const double &UI);

        std::shared_ptr<Cantera::Solution> postReflFr(const double &U2, const std::string &mech);

        std::shared_ptr<Cantera::Solution> postReflEq(const double &U2, const std::string &mech);

        void FHFP_Refl(const double &U2, std::shared_ptr<Cantera::Solution> gas2, std::shared_ptr<Cantera::Solution> gas3, double &FH, double &FP);

    }

    namespace utils
    {
        double getThermicity(std::shared_ptr<Cantera::Solution> gas);
        std::vector<double> RK45();
        std::vector<double> LSODA();
        std::vector<double> Radau();
        std::shared_ptr<Cantera::Solution> initializeGas(const std::string &what,const std::string &mech, double propA, double propB,const std::vector<double> &comp);
    }

    namespace stagnation
    {
        class system
        {
            private:
            double t_0;
            double delta;
            std::shared_ptr<Cantera::Solution> gas;
            std::shared_ptr<Cantera::Solution> gas1;


            public:
            system();
            ~system();
            void solve();

        };

        std::shared_ptr<system> newSystem();
    }

    namespace detonation
    {
        class system
        {
            private:

            double t_0;
            double delta;
            std::vector<double> t;
            std::vector<double> y;
            std::shared_ptr<Cantera::Solution> gas;
            std::shared_ptr<Cantera::Solution> gas1;
            

            public:
            system();
            ~system();
            void solveCP(std::shared_ptr<Cantera::Solution> gas, double t, std::vector<double> y);

            void solceCV(std::shared_ptr<Cantera::Solution> gas, double t, std::vector<double> y);

            void solveZND(std::shared_ptr<Cantera::Solution> gas, double t, std::vector<double> y);
        };

        std::shared_ptr<system> newSystem();
    }



}

#endif