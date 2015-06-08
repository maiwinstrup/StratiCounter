#include "IOStructlib.h"
#include <iostream>

#define NUM_IN 2
#define NUM_OUT 2

struct inStruct
{
  double var1;
  /* The following two fields will be entered as a substructure */
  const char *var2;
};

struct outStruct
{
  double var1;
  const char *var2;
};

int run_main(int argc, char **argv)
{
    // Call application and library initialization. Perform this 
    // initialization before calling any API functions or
    // Compiler-generated libraries.
    if (!mclInitializeApplication(NULL,0)) 
    {
        std::cerr << "could not initialize the application properly"
                   << std::endl;
    	return -1;
    }
    if( !IOStructlibInitialize() )
    {
        std::cerr << "could not initialize the library properly"
                   << std::endl;
	return -1;
    }
    else
    {
        try
        {
		
		const char *in_field_names[] = {"val1", "val2"};
		const char *out_field_names[] = {"val3", "val4"};

		mwArray inS(1,1,2,in_field_names);
		inS.Get("val1",1,1).Set(mwArray(2));
		inS.Get("val2",1,1).Set(mwArray(17));

		mwArray outS(1,1,2,out_field_names);

		inAndOut(1,outS,inS);
		

		outStruct myOut;
		for (int i=0; i<NUM_OUT; i++){
			myOut.var1 = int(outS.Get("val3",1,1));
			myOut.var2 = outS.Get("val4",1,1).ToString();
		}
		
		std::cout<< myOut.var1 << std::endl;
		std::cout<< myOut.var2 << std::endl;
		
        }
        catch (const mwException& e)
        {
          std::cerr << e.what() << std::endl;
          return -2;
        }
        catch (...)
        {
          std::cerr << "Unexpected error thrown" << std::endl;
          return -3;
        }     
		
        // Call the application and library termination routine
        IOStructlibTerminate();
    }
/* Note that you should call mclTerminate application at the end of
 * your application. 
 */
    mclTerminateApplication();
    return 0;
}

int main()
{
    mclmcrInitialize();
    return mclRunMain((mclMainFcnType)run_main,0,NULL);
}




