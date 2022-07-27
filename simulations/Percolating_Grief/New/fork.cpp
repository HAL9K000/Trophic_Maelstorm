#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>

#include <random>
#include <vector>
#include<bits/stdc++.h>
using namespace std;



int main()
{
    double a = 0.03; double dt = 0.1; double dx = 0.5; double D = 0.25; double sigma = 0.25; 
    double beta = a - (10.0/3.0)*D/(dx*dx);
    double lambda_exp = exp( (beta)*dt); double lambda = 2*(beta)/(sigma*sigma*(lambda_exp -1.0));

    int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
	//https://cplusplus.com/forum/beginner/220120/
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
    poisson_distribution<int> poisson; gamma_distribution <double> gamma;

    for(int i=0; i< 25; i++)
    {
        double mu = -1.0;
	    poisson = poisson_distribution<int>(lambda*lambda_exp*0.0);
        int sod= poisson(rng);	
	    gamma = gamma_distribution<double>(mu + 1.0 + sod, 1.0); 
        double fod = gamma(rng);

        cout << "Poisson Value:\t"<< sod << " & Gamma Value:\t"<< fod << " for Trial: " << i;
    }
}