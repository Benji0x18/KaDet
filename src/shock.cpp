#include "Kadet.hpp"

namespace KaDet
{
    namespace shock
    {
        void FHFP(const double& U, std::shared_ptr<Ct::Solution> gas2, std::shared_ptr<Ct::Solution> gas1, double &FH, double &FP )
        {
            double D = gas2->thermo()->density();
            double H = gas2->thermo()->enthalpy_mass();
            double P = gas2->thermo()->pressure();

            double D1 = gas1->thermo()->density();
            double H1 = gas1->thermo()->enthalpy_mass();
            double P1 = gas1->thermo()->pressure();

            double w1 = U *  U;
            double w2 = w1 * (D1/D) * (D1/D);

            FH = H + 0.5 * w2 - (H1 + 0.5*w1);
            FP = P + D * w2 - (P1 + D1 * w1);
        }


        int eqCalc(const double &U, std::shared_ptr<Ct::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, double errT, double errV)
        {
            // This is the function is the equivalent of shk_eq_calc.m from Shock Detonation Toolbox. https://shepherd.caltech.edu/EDL/PublicResources/sdt/ (retrieved 24.01.2024)
            
            int volumeBoundRatio = 5;
            double R1 = gas1->thermo()->density();
            double V1 = 1/R1;
            double P1 = gas1->thermo()->pressure();
            double T1 = gas1->thermo()->temperature();
            double deltaT = 1000;
            double deltaVo = 1000;

            //Preliminary Guess

            double V = V1/volumeBoundRatio;
            double R = 1/V;
            double P = P1 + R1*(U*U)*(1-R1/R);
            double T = T1*P*R1/(P1*R);
            gas2->thermo()->setState_TD(T, R);
            gas2->thermo()->equilibrate("TV");

            int j = 0;
            while (std::abs(deltaT) > ((errT)*T) || std::abs(deltaVo) > ((errV)*V))
            {
                j++;
                if(j==500)
                {
                    std::cout << "Shock Equilibrium Calculation did not converge for U = " << U << '\n';
                    break;
                    return 1;
                }
                // Calculate FH and FP for guess 1;
                double FH, FP;
                FHFP(U, gas2, gas1, FH, FP);

                // Calculate perturbations

                double DT, DV, Tper, Rper, Vper, FHX, FPX, DFHDT, DFHDV, DFPDT, DFPDV; 
                
                // Temperature Perturbation
                DT   = T*0.02;
                Tper = T+DT;
                Rper = R;

                gas2->thermo()->setState_TD(Tper, Rper);
                gas2->thermo()->equilibrate("TV");

                FHFP(U, gas2, gas1, FHX, FPX);

                DFHDT = (FHX - FH)/DT;
                DFPDT = (FPX - FP)/DT;

                // Volume perturbation
                DV = V*0.02;
                Tper = T;
                Vper = V+DV;
                Rper = 1/Vper;

                gas2->thermo()->setState_TD(Tper, Rper);
                gas2->thermo()->equilibrate("TV");

                FHFP(U, gas2, gas1, FHX, FPX);

                DFHDV = (FHX - FH)/DV;
                DFPDV = (FPX - FP)/DV;

                double J[] = {DFPDV, -DFHDV, -DFPDT, DFHDV};
                double detJ = DFHDT*DFPDV - DFHDV*DFPDT;
                double b[] = {-FH, -FP};
                deltaT = (J[0]*b[0] + J[1]*b[1])/detJ;
                deltaVo = (J[2]*b[0] + J[3]*b[1])/detJ;

                // Check and limit change in values

                // Temperature
                double DTM = 0.2*T;
                if (std::abs(deltaT) > DTM)
                {
                    deltaT = DTM*deltaT/std::abs(deltaT);
                }

                // Volume
                double v2x = V + deltaVo;
                double DVM;
                double deltaV = deltaVo;
                if (v2x > V1)
                {
                    DVM = 0.2*(V1-V);
                }else
                {
                    DVM = 0.2*V;
                }

                if(std::abs(deltaV) < 10)
                {
                    deltaV = 0.01*deltaV;
                }
                else if(std::abs(deltaV) > DVM)
                {
                    deltaV = 0.5*DVM*deltaV/std::abs(deltaV);
                }
                
                T = T + deltaT;
                V = V + deltaV;
                R = 1/V;

                gas2->thermo()->setState_TD(T,R);
                gas2->thermo()->equilibrate("TV");
                std::cout << gas2->thermo()->report() << '\n';
            }   
            return 0;
        }
    }
    
}