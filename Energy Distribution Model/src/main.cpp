#include "pch.h"

#include "Application.h"


int main(int argc, char** argv) 
{
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
