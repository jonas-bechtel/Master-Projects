#include "pch.h"

#include "Application.h"

#include "force.h"
#include "beam.h"

int main(int argc, char** argv) 
{
    //UniformCylinder beam(1,1);
    //beam.set_v_rms(20000, 5000);
   // ForcePark model;
   // model.set_v_eff(0);
   // model.set_mag_field(0.01);
   // model.set_time_cooler(1e-6);
   // std::vector<double> v_tr = { 0 , 0 };
   // std::vector<double> v_l = { 0, 20000 };
   // std::vector<double> n_e = { 1e12, 1e12 };
   // std::vector<double> f_tr;
   // std::vector<double> f_l;
   // model.friction_force(1, 2, v_tr, v_l, n_e, beam, f_tr, f_l);
   //
   // std::cout << f_l.at(0) / k_e << " , " << f_l.at(1) / k_e << std::endl;

    //ForceNonMagNumeric3D model2;
    //model2.set_time_cooler(1e-6);
    ////model2.set_mag_field(0.01);
    //std::vector<double> v_tr = { 1000 , 0 };
    //std::vector<double> v_l = { 20000, 20000 };
    //std::vector<double> n_e = { 1e12, 1e12 };
    //std::vector<double> f_tr;
    //std::vector<double> f_l;
    //model2.friction_force(1, 2, v_tr, v_l, n_e, beam, f_tr, f_l);
    //
    //std::cout << f_l.at(0) / k_e << " , " << f_l.at(1) / k_e << std::endl;
    //return 0;

    //TApplication ROOTapp = TApplication("app1", nullptr, nullptr);
    //TApplication ROOTapp2 = TApplication("app1", nullptr, nullptr);
   
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
