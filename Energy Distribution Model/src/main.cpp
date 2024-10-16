#include "Application.h"

#include <vector>

void test();

int main(int argc, char** argv) {
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

    //test();

    return 0;
}
