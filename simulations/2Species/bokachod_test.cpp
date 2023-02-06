#include "2species_stochastic.h"
#include "Debug.h"

void increase_stack_limit(long long stack_size){

	//Source: https://stackoverflow.com/questions/2279052/increase-stack-size-in-linux-with-setrlimit
	//Credit: https://stackoverflow.com/users/253056/paul-r
	//Used with modification

	// Necessary because of complexity of integrations.

	const rlim_t kStackSize = stack_size * 1024LL * 1024LL;   // min stack size = 128 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    }
}

bool maxis(int a, int b) 
{
	//Binary function necessary for determining max in a list.
 	return (a < b); 
}


void add_three(int a, int b, int c)
{
	int d = a + b +c; cout << d << endl;
}

void init_randframe(D2Vec_Double &array, int Sp, int size, double mean[], double sd[])
{
	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);
	//Returns array with all elements initialised to random values drawn from a Gaussian distribution with mean mu and SD sigma,
	// with all values 
	for(int s=0; s< Sp; s++) 
	{
		for(int i=0; i< size; i++)
		{	
		normal_distribution<double> norm;
		norm = normal_distribution<double>(mean[s], sd[s]);
		array[s][i] =norm(rng);
		if(array[s][i] < 0)
			array[s][i] = 0.0; //Values less than 0 are unrealistic and not permissible.
		}
	}
}