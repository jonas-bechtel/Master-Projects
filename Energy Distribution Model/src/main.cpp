#include "Application.h"

#include "CL/cl2.hpp"
#include <vector>

void test();

int main(int argc, char** argv) {
    Application app;
    try
    {
        app.Run();
    }
    catch (const std::runtime_error& e)
    {
        // Code to handle the exception
        std::cerr << "Caught an error: " << e.what() << std::endl;
    }
    catch (const std::exception& e)
    {
        // Catches any other std::exception derived exceptions
        std::cerr << "Caught a generic error: " << e.what() << std::endl;
    }

    //test();

    return 0;
}

void test()
{
    const char* kernelSource = R"(
__kernel void vectorAdd(__global const float* A,
                        __global const float* B,
                        __global float* C,
                        const unsigned int n) {
    int i = get_global_id(0);
    if (i < n) {
        C[i] = A[i] + B[i];
    }
}
)";

    // Initialize vectors
    const int N = 1024;
    std::vector<float> A(N, 1.0f); // Vector A initialized to 1.0
    std::vector<float> B(N, 2.0f); // Vector B initialized to 2.0
    std::vector<float> C(N, 0.0f); // Result vector C

    // Get platforms and devices
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    cl::Platform platform = platforms[0]; // Use the first platform

    std::vector<cl::Device> devices;
    platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
    cl::Device device = devices[0];

    std::cout << device.getInfo<CL_DEVICE_NAME>() << device.getInfo<CL_DEVICE_VERSION>() << "\n";
    std::cout << "available: " << device.getInfo<CL_DEVICE_AVAILABLE>() << "\n";

    // Create context and command queue
    cl_int err;
    cl::Context context(device, NULL, NULL, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error creating context: " << err << std::endl;
        return ;
    }

    cl::CommandQueue queue(context, device, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error creating context: " << err << std::endl;
        return;
    }

    // Create buffers
    cl::Buffer bufferA(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * N, A.data());
    cl::Buffer bufferB(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * N, B.data());
    cl::Buffer bufferC(context, CL_MEM_READ_WRITE, sizeof(float) * N);

    // Create and build the program
    cl::Program program(context, kernelSource);
    if (program.build({ device }) != CL_SUCCESS) {
        std::cout << " Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << "\n";
        exit(1);
    }

    // Create the kernel
    cl::Kernel kernel(program, "vectorAdd");
    std::cout << kernel.getInfo<CL_KERNEL_FUNCTION_NAME>() << "\n";

    // Set kernel arguments
    kernel.setArg(0, bufferA);
    kernel.setArg(1, bufferB);
    kernel.setArg(2, bufferC);
    kernel.setArg(3, N);

    // Execute the kernel
    cl::NDRange global(N);
    auto result = queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, cl::NullRange);
    if (result != CL_SUCCESS)
    {
        std::cerr << "Error queuing up the kernel: " << result << "\n";
        return;
    }
    queue.finish(); // Wait for the kernel to finish
    
    
    // Read the result back to host
    queue.enqueueReadBuffer(bufferC, CL_TRUE, 0, sizeof(float) * N, C.data());

    // Verify the result
    for (int i = 0; i < N; ++i) {
        if (C[i] != 3.0f) {
            std::cerr << "Error at index " << i << ": expected 3.0, got " << C[i] << std::endl;
            return;
        }
    }

    std::cout << "Vector addition completed successfully!" << std::endl;
}
