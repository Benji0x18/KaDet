/*
    This is a header of the Shock and Detonation Library KaDet. Whenever using KaDet make sure to include this header in the file.
    
*/
#ifndef KADET
#define KADET
#include <iostream>
#include <cmath>
#include <sstream>
#include "cantera/core.h"

namespace Kadet
{
    
    namespace shock
    {
        std::array<double,6> LSQ_speedCJ(const std::vector<double> &x, const std::vector<double> &y);

        double hugonoitFr(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Cantera::Solution> gas);

        double hugonoitEq(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Cantera::Solution> gas);

        double calculateCJ(std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, const double &errV, const double &errT, const double &x);

        double speedCJ(const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        std::array<std::shared_ptr<Cantera::Solution>, 2> postShockFr(const double &W1, const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        std::array<std::shared_ptr<Cantera::Solution>, 2> postShockEq(const double &W1, const double &P1, const double &T1, const std::string &comp, const std::string &mech);

        void FHFP(const double &W, std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, double &FH, double &FP);

        int shockFr(const double &W, std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, double errT, double errV);

        int shockEq(const double &W, std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, double errT, double errV);

        double reflectionFr(std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, std::shared_ptr<Cantera::Solution> gas3, const double &WI);

        double reflectionEq(std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, std::shared_ptr<Cantera::Solution> gas3, const double &WI);

        std::shared_ptr<Cantera::Solution> postReflectionFr(const double &W2, const std::string &mech);

        std::shared_ptr<Cantera::Solution> postReflectionEq(const double &W2, const std::string &mech);

        void FHFP_Refl(const double &W2, std::shared_ptr<Cantera::Solution> gas2, std::shared_ptr<Cantera::Solution> gas3, double &FH, double &FP);

    }

    namespace utils
    {
        double getThermicity(std::shared_ptr<Cantera::Solution> gas);
        std::vector< std::pair< double, std::vector<double> > > RK45(const double &y_0, const double &t_0, std::vector<double> (*)());
        std::vector< std::pair< double, std::vector<double> > > RK45(system sys);
        std::vector<std::array<double,2>> LSODA();
        std::vector<std::array<double,2>> Radau();
        std::shared_ptr<Cantera::Solution> initializeGas(const std::string &what_props, const std::string &mech, const double &property1, const double &property2, const std::vector<double> &comp);
        
        
        

        class system
        {
            private:

                double t_0;
                std::vector<double> y_0;
                double delta;
                std::vector<double> t;
                std::vector<std::vector<double>> y;
            
            public:
                system(const double &t_0, const std::vector<double> y_0, const double &delta);
                ~system();
        };
            
    }

    namespace stagnation
    {
        class system : private utils::system
        {
            private:
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
        class system : private utils::system
        {
            private:
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