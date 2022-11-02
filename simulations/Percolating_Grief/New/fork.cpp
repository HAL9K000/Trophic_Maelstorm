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

#include <sys/resource.h>
#include <omp.h>

void determine_neighbours_R2(int neighbours_R2[][2][5], int g);
void determine_neighbours_S8(int neighbours_R2[][3][3], int g);
void fuckwad(vector <vector<int>> &hic, int st, int thr);


int main()
{
    double a = 0.03; double dt = 0.1; double dx = 0.5; double D = 0.25; double sigma = 0.25; 
    double beta = a - (10.0/3.0)*D/(dx*dx);
    double lambda_exp = exp( (beta)*dt); double lambda = 2*(beta)/(sigma*sigma*(lambda_exp -1.0));
    int g =4;

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

        cout << "Poisson Value:\t"<< sod << "  Gamma Value:\t"<< fod << " for Trial: " << i;
    }
    cout << "Lambda:\t" << lambda << " with exponent:\t" << lambda_exp <<endl;

    int nS8[g*g][3][3] ={0}; int nR2[g*g][2][5] ={0};

    cout << "TESTING S8" <<endl;

    determine_neighbours_S8(nS8, g);

    for (int i =0; i<g*g; i++)
    {
        cout << "For index i= " << i << "\t the neighbours in a Norm 1 radius are:" <<endl;
        for(int y =0 ; y < 3; y++)
        {
            for(int x =0 ; x < 3; x++)
                cout << nS8[i][y][x] << " ";
            cout << endl;
        } 
    }

    cout << "TESTING R2" <<endl;
    determine_neighbours_R2(nR2, g);

    for (int i =0; i<g*g; i++)
    {
        cout << "For index i= " << i << "\t the neighbours in a Van Neumann Radius of 2 are:" <<endl;
        for(int y =0 ; y < 2; y++)
        {
            for(int x =0 ; x < 5; x++)
                cout << nR2[i][y][x] << " ";
            cout << endl;
        } 
    }
    std::vector<vector<int>> vec;
    #pragma omp parallel
    {
        std::vector<vector<int>> vec_private;

        #pragma omp for nowait
      for (int i=0; i < 200; i+=25)
      {
        int nom = omp_get_thread_num();
        std::vector<vector <int>> CExpRho_a;
        fuckwad(CExpRho_a, i , nom);
        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
        vector<vector <int>>().swap(CExpRho_a);
      }

      #pragma omp critical
      {
        vec.insert(vec.end(), vec_private.begin(), vec_private.end());
        // Inserting critical exponent data for each grid size in order.
        stringstream message3;
        message3 << "Is this happening?\n";
        cout << message3.str();

      
	    //vector<vector <int>>().swap(vec_private);
      }
    }

    ofstream output_dp;
    // Creating a file instance called output to store output data as CSV.
	output_dp.open("FORK_ULTIMATE.csv");

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
	output_dp << " Thread # , Num \n";
	cout << "The vector elements are: "<< endl;
  	cout << "  Thread # , Num \n";
    cout << "Vec size: " << vec.size() <<endl;
	for(int i=0; i< vec.size(); i++)
        {
            cout << vec[i][0] << "," << vec[i][1] <<endl;
            output_dp << vec[i][0] << "," << vec[i][1] <<endl;
        }
    output_dp.close();
}

void fuckwad(vector <vector<int>> &hic, int st, int thr)
{
    ofstream output_dp;
    // Creating a file instance called output to store output data as CSV.
    stringstream num; num << thr;
	output_dp.open("Fork_" + num.str() + ".csv");

    for(int i= st; i < st +25; i++)
    {    
        output_dp << thr << "," << i <<endl;
        hic.push_back({thr,i});
    }
    output_dp.close();
}

void determine_neighbours_R2(int neighbours_R2[][2][5], int g)
{
	//Stores neighbouring sites (Van-Neumann Radius of 2) of each cell site in a 
	// g^2*2*5 matrix.
	int a=0; int b= 0; //Useful for accounting purposes
	for (int i=0; i<g*g; i++)
	{
		int p = int(i/g); int q = i%g;
		//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.
		for(int j=0; j < 5; j++)
		{
			int k= j-2;
			a= ((p+k)%g + g)%g; // Eq to (p+k)mod g. -2 <= k <= 2
			b= ((q+k)%g + g)%g; // Eq to (q+k)mod g. -2 <= k <= 2

			//Recall, first row stores vertical neighbours, second row stores horizontal neighbours.
			neighbours_R2[i][0][j] = a*g + q; neighbours_R2[i][1][j] = p*g + b;
		}
	}

}

void determine_neighbours_S8(int neighbours_R2[][3][3], int g)
{
	//Stores neighbouring sites in a square (Norm 1) radius of 1
	int a=0; int b= 0; //Useful for accounting purposes
	for (int i=0; i<g*g; i++)
	{
		int p = int(i/g); int q = i%g;
		//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.
		for(int y=0; y < 3; y++)
		{	int ky= y-1; a= ((p+ky)%g + g)%g; // Eq to (p+ky)mod g. -1 <= ky <= 1
			for(int x=0; x < 3; x++)
			{
				int kx = x-1;
				b= ((q+kx)%g + g)%g; // Eq to (q+kx)mod g. -1 <= kx <= 1
				//Recall, first row stores vertical neighbours, second row stores horizontal neighbours.
				neighbours_R2[i][y][x] = a*g + b; 
			}	
		}
	}

}