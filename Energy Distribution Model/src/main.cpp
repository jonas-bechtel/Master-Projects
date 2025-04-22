#include "pch.h"

#include "Application.h"

//#include "xForce.h"

int main(int argc, char** argv) 
{
    //xFrParam params;
    //params.mfield = 0.1;
    //params.Z = 1;
    //params.n_e = 1e11;
    //params.Smoos = 2;
    //params.tau = 1e-6;
    //params.V_tr_e = 10000;
    //params.V_tr_x = params.V_tr_e / sqrt(2);
    //params.V_tr_y = params.V_tr_e / sqrt(2);
    //params.V_long_e = 5000;
    //params.V_eff_e = 5000;
    //iForce.v[0] = 0;
    //iForce.v[1] = 0;
    //iForce.v[2] = -5000;
    //iForce.Vtr = 0;
    //iForce.Magnetized = 1;
    //iForce.Fast = 1;
    //iForce.Adiabatic = 1;
    //iForce.dt = 50;
    //iForce.dl = 50;
    //iForce.nfi = 50;
    //iForce.N_M = 50;
    //iForce.D3dl = 50;
    //iForce.D3dx = 50;
    //iForce.D3dy = 50;
    //// for toepper
    //iForce.TFast = 1;
    //iForce.Tight = 1;
    //iForce.Stretched = 1;
    //iForce.Tdl = 50;
    //iForce.Tdt = 50;
    //iForce.Tnfi = 50;


    ////iForce.DerSkr(params);
    ////iForce.NonMag(params);
    ////iForce.Parhom(params);
    ////iForce.D4(params);
    //iForce.Toepffer(params);
    //std::cout << iForce.Ftr.v << ", " << iForce.f[2].v << std::endl;

    return 0;
    Application app;
    app.Run();
    //try
    //{
    //    app.Run();
    //}
    //catch (const std::runtime_error& e)
    //{
    //    // Code to handle the exception
    //    std::cerr << "Caught an error: " << e.what() << std::endl;
    //}
    //catch (const std::exception& e)
    //{
    //    // Catches any other std::exception derived exceptions
    //    std::cerr << "Caught a generic error: " << e.what() << std::endl;
    //}

    return 0;
}
