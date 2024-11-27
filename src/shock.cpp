#include "Kadet.hpp"

namespace Kadet::shock
{
    

    std::array<double,6> LSQ_speedCJ(const std::vector<double> &x, const std::vector<double> &y)
    {
        //Intialize values
        uint k = 0;
        double X = 0, X2 = 0, X3 = 0, X4 = 0;
        double Y = 0,Y1 = 0,Y2 = 0;
        double a = 0, b = 0, c = 0;
        uint n = x.size();

        while(k < n)
        {
            double x_k = x.at(k);
            double y_k = y.at(k);
            X = X + x_k;
            X2 = X2 + x_k*x_k;
            X3 = X3 + x_k*x_k*x_k;
            X4 = X4 + x_k*x_k*x_k*x_k;
            
            Y = Y + y_k;
            Y1 = Y1 + y_k*x_k;
            Y2 = Y2 + y_k*x_k*x_k;

            k++;
        }
        double m = Y/double(n);

        double den = (X3*(double(n)) - X2*X);

        double temp = (den*(X*X2 - X3*double(n)) + X2*X2*(X*X - double(n)*X2) - X4*double(n)*(X*X - X2*double(n)));
        double temp2 = (den*(Y*X2 - Y2*double(n)) + (Y1*double(n) - Y*X)*(X4*double(n) - X2*X2));

        std::array<double,3> coef;

        coef.at(0) = 1.0/den*(double(n)*Y1 - Y*X - b*(X2*double(n)-X*X));
        coef.at(1) = temp2/temp;
        coef.at(2) = 1/double(n)*(Y - a*X2 - b*X);

        k = 0;
        
        double SSE = 0; double SST = 0;

        while (k < n)
        {
            double x_k = x.at(k);
            double y_k = y.at(k);
            double f_k = a*x_k*x_k + b*x_k + c;
            SSE = SSE + (y_k - f_k)*(y_k - f_k);
            SST = SST + (y_k - m)*(y_k - m);
            k++;
        }

        std::array<double, 6> res;

        for(int n = 0; n < 3; n++)
        {
            res.at(n) = coef.at(n);
        }

        res.at(3) = 1-SSE/SST; // R2
        res.at(4) = SSE;
        res.at(5) = SST;

        return res;
    }

    double hugonoitFr(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Cantera::Solution> gas)
    {
        gas->thermo()->setState_TD(Ta,1/Vb);

        double Hb1 = gas->thermo()->enthalpy_mass();
        double Pb = Cantera::GasConstant * Ta / (gas->thermo()->meanMolecularWeight()*Vb);

        double Hb2 = H1 + 0.5*(Pb-P1)*(Vb + V1);

        return Hb2-Hb1;
    }

    double hugonoitEq(const double &Ta, const double &Vb, const double &H1, const double &P1, const double &V1, std::shared_ptr<Cantera::Solution> gas)
    {
        gas->thermo()->setState_TD(Ta,1/Vb);

        gas->thermo()->equilibrate("TV");

        double Hb1 = gas->thermo()->enthalpy_mass();
        double Pb = Cantera::GasConstant * Ta / (gas->thermo()->meanMolecularWeight()*Vb);

        double Hb2 = H1 + 0.5*(Pb-P1)*(Vb + V1);

        return Hb2-Hb1;
    }

    double calculateCJ(std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, const double &errV, const double &errT, const double &x)
    {
        double W1 = 2000;

        double R1 = gas1->thermo()->density();
        double V1 = 1/R1;

        double T = 2000; //Initial Guess
        double V = V1/x;

        double DT = 1000;
        double DW = 1000;

        uint i = 0;

        double Tper, Wper, Rper, Vper, FH, FP, FHX, FPX, DFHDT, DFHDW, DFPDT, DFPDW, DTM; 

        while(std::abs(DT) < errT*T && std::abs(DW) < errV*V)
        {
            if(i > 500)
            {
                std::cout << "CJ speed calculation didn't converge for specific volume ratio" << x << std::endl;
                return 0;
            }


            FHFP(W1,gas2,gas1,FH,FP);

            // Local temperature perturbation impact calculation
            DT = T*0.02;
            Tper = T + DT;

            Vper = V;
            Rper = 1/Vper;

            Wper = W1;

            gas2->thermo()->setState_TD(Tper,Rper);
            gas2->thermo()->equilibrate("TV");

            FHFP(Wper,gas2,gas1,FHX,FPX);

            DFHDT = FHX - FH;
            DFPDT = FPX - FP;

            // Local velocity perturbation inpact calculation
            Tper = T;

            Vper = V;
            Rper = 1/Vper;
            
            DW = W1*0.02;
            Wper = W1 + DW;

            gas2->thermo()->setState_TD(Tper,Rper);
            gas2->thermo()->equilibrate("TV");

            FHFP(Wper,gas2,gas1,FHX,FPX);

            DFHDW = FHX - FH;
            DFPDW = FPX - FP;

            double J[] = {DFPDW, -DFHDW, -DFPDT, DFHDT};
            double detJ = DFHDT*DFPDW - DFHDW*DFPDT;
            double b[] = {-FH, -FP};
            DT = (J[0]*b[0] + J[1]*b[1])/detJ;
            DW = (J[2]*b[0] + J[3]*b[1])/detJ;


            DTM = 0.02*T;
            if(std::abs(DT) > DTM)
            {
                DT = DTM*DT/std::abs(DT);
            }

            T = T + DT;
            W1 = W1 + DW;
        }

        return W1;
    }

    double speedCJ(const double &P1, const double &T1, const std::string &comp, const std::string &mech)
    {
        uint8_t steps = 20; double maxv = 2; double minv = 1.5;

        std::vector<double> W1;
        std::vector<double> rr;

        W1.reserve(21);
        rr.reserve(21);

        std::shared_ptr<Cantera::Solution> gas1 = Cantera::newSolution(mech);
        std::shared_ptr<Cantera::Solution> gas2 = Cantera::newSolution(mech);

        gas1->thermo()->setState_TPX(T1,P1,comp);

        double a,b,c,Rnew,cj_speed,R2 = 0; 
        std::array<double,6> results;
        for(int n = 0; n < 4 || R2 < 0.99999; n++)
        {
            double x_step = (maxv-minv)/steps ;
            double x = minv;

            for(int o = 0; x <= maxv; o++)
            {
                gas2->thermo()->setState_TPX(T1,P1,comp);
                W1.push_back(calculateCJ(gas1,gas2,1e-4,1e-4,x));
                rr.push_back(gas2->thermo()->density()/gas1->thermo()->density());
                x = x+x_step;
            }
            results = LSQ_speedCJ(rr,W1);
            a = results.at(0); b = results.at(1); c = results.at(2); R2 = results.at(3);
            Rnew = -results.at(1)/(2*results.at(0));
            minv = Rnew - Rnew*1e-3;
            maxv = Rnew + Rnew*1e-3;
        }

        return cj_speed;
    }

    void FHFP(const double& W, std::shared_ptr<Cantera::Solution> gas2, std::shared_ptr<Cantera::Solution> gas1, double &FH, double &FP )
    {
        double D = gas2->thermo()->density();
        double H = gas2->thermo()->enthalpy_mass();
        double P = gas2->thermo()->pressure();

        double D1 = gas1->thermo()->density();
        double H1 = gas1->thermo()->enthalpy_mass();
        double P1 = gas1->thermo()->pressure();

        double w1 = W *  W;
        double w2 = w1 * (D1/D) * (D1/D);

        FH = H + 0.5 * w2 - (H1 + 0.5*w1);
        FP = P + D * w2 - (P1 + D1 * w1);
    }


    int shockEq(const double &W, std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, double errT, double errV)
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
        double P = P1 + R1*(W*W)*(1-R1/R);
        double T = T1*P*R1/(P1*R);
        gas2->thermo()->setState_TD(T, R);
        gas2->thermo()->equilibrate("TV");
        double DT, DV, Tper, Rper, Vper, FH, FP, FHX, FPX, DFHDT, DFHDV, DFPDT, DFPDV; 
        uint j = 0;
        while (std::abs(deltaT) > ((errT)*T) || std::abs(deltaVo) > ((errV)*V))
        {
            j++;
            if(j==500)
            {
                std::cout << "Shock Equilibrium Calculation did not converge for W = " << W << '\n';
                break;
                return 1;
            }
            // Calculate FH and FP for guess 1;
            FHFP(W, gas2, gas1, FH, FP);

            // Calculate perturbations
            // Temperature Perturbation
            DT   = T*0.02;
            Tper = T+DT;
            Rper = R;

            gas2->thermo()->setState_TD(Tper, Rper);
            gas2->thermo()->equilibrate("TV");

            FHFP(W, gas2, gas1, FHX, FPX);

            DFHDT = (FHX - FH)/DT;
            DFPDT = (FPX - FP)/DT;

            // Volume perturbation
            DV = V*0.02;
            Tper = T;
            Vper = V+DV;
            Rper = 1/Vper;

            gas2->thermo()->setState_TD(Tper, Rper);
            gas2->thermo()->equilibrate("TV");

            FHFP(W, gas2, gas1, FHX, FPX);

            DFHDV = (FHX - FH)/DV;
            DFPDV = (FPX - FP)/DV;

            double J[] = {DFPDV, -DFHDV, -DFPDT, DFHDT};
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
        }   
        return 0;
    }

    int shockFr(const double &W, std::shared_ptr<Cantera::Solution> gas1, std::shared_ptr<Cantera::Solution> gas2, double errT, double errV)
    {
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
        double P = P1 + R1*(W*W)*(1-R1/R);
        double T = T1*P*R1/(P1*R);
        gas2->thermo()->setState_TD(T, R);
        gas2->thermo()->equilibrate("TV");

        int j = 0;
        while (std::abs(deltaT) > ((errT)*T) || std::abs(deltaVo) > ((errV)*V))
        {
            j++;
            if(j==500)
            {
                std::cout << "Shock Equilibrium Calculation did not converge for W = " << W << '\n';
                break;
                return 1;
            }
            // Calculate FH and FP for guess 1;
            double FH, FP;
            FHFP(W, gas2, gas1, FH, FP);

            // Calculate perturbations

            double DT, DV, Tper, Rper, Vper, FHX, FPX, DFHDT, DFHDV, DFPDT, DFPDV; 
            
            // Temperature Perturbation
            DT   = T*0.02;
            Tper = T+DT;
            Rper = R;

            gas2->thermo()->setState_TD(Tper, Rper);

            FHFP(W, gas2, gas1, FHX, FPX);

            DFHDT = (FHX - FH)/DT;
            DFPDT = (FPX - FP)/DT;

            // Volume perturbation
            DV = V*0.02;
            Tper = T;
            Vper = V+DV;
            Rper = 1/Vper;

            gas2->thermo()->setState_TD(Tper, Rper);

            FHFP(W, gas2, gas1, FHX, FPX);

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
        }   
        return 0;
    }
    
    
}