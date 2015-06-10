#include "libStraticounter.h"
#include <iostream>

int main(int ac, const char *av[])
{

    // Initialize the MATLAB Compiler Runtime global state
    if (!mclInitializeApplication(NULL,0))
    {
        std::cerr << "Could not initialize the application properly."
                  << std::endl;
    	return -1;
    }

    // Initialize the library
    if( !libStraticounterInitialize() )
    {
        std::cerr << "Could not initialize the library properly."
                  << std::endl;
	return -1;
    }

    // Must declare all MATLAB data types after initializing the
    // application and the library, or their constructors will fail.
    mwArray settingsPath(av[1]);
    mwArray outputPath(av[2]);

    // Call the function
    straticounter_scibox(settingsPath, outputPath);

    // Shut down the library and the application global state.
    libStraticounterTerminate();
    mclTerminateApplication();
}
