
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <list>
#include <thread>
using namespace std::chrono;
//static std::mt19937_64 rng(time(NULL));
#include<bits/stdc++.h>
using namespace std;
#include <sys/resource.h>
#include <omp.h>

#include <map>

#ifndef tuSP_DYNAMICS_H

#define tuSP_DYNAMICS_H

inline const int Sp = 2; //Number of species in the system.
inline const int Sp4_1 = 4*Sp +1; // Used for generating statistics on surviving runs
inline const int Sp2 = 2*Sp; // Used for generating statistics on surviving runs

//------------------------- TYPEDEFS --------------------------------------------------------------//

typedef std::vector<std::vector <double>> D2Vec_Double;
typedef std::vector<std::vector <int>> D2Vec_Int;

typedef std::vector<std::vector <std::vector <double>>> D3Vec_Double;
typedef std::vector<std::vector <std::vector <int>>> D3Vec_Int;

typedef std::vector<std::vector <std::vector <std::vector <double>>>> D4Vec_Double;
typedef std::vector<std::vector <std::vector <std::vector <int>>>> D4Vec_Int;



//---------------------Some Primitives ------------------------------------------------------ //


//Credits: Mohmahkho, https://codeforces.com/blog/entry/76149

// Meta programming template used to create a vector of dimension D, with specified lengths of each dimension (given by Args)
// with default length of each dimension set to 0 if not specified. Last parameter

/**

template<int D, typename T> struct createVec : public vector<createVec<D - 1, T>> 
{
  static_assert(D >= 1, "Vector dimension must be > 0!");
  template<typename... Args> createVec(int n = 0, Args... args) : vector<createVec<D - 1, T>>(n, createVec<D - 1, T>(args...)) 
  {
  }
};
template<typename T> struct createVec<1, T> : public vector<T> 
{
  createVec(int n = 0, const T& val = T()) : vector<T>(n, val) 
  {
  }
};

*/



struct coordinates {
   int x;
   int y;
};

struct pde_val {
	int label;  //Stores grid location of PDE values
	std::vector<double> pde;   //PDE values at label index.
};

struct f_coordinates {
   float x;
   float y;
};

struct zd_coordinates {
   double x;
   double y;
   double z;
};



void increase_stack_limit(long long stack_size);
bool maxis(int a, int b);
void add_three(int a, int b, int c); //Test function.

void init_fullframe(D2Vec_Double &array, int Sp, int size);
void init_constframe(D2Vec_Double &array, int Sp, int size, double constant[]);
void init_randframe(D2Vec_Double &array, int Sp, int size, double mean[], double sd[]);
void init_quarterframe(D2Vec_Double &array, double c1, double c2, double c3, double c4, int L);
void init_solitarytear(D2Vec_Double &array, int length);

double mean_of_array(double array[],int size);
double standard_deviation_of_array(double array[],int size);
double mean_of_vector(vector<double> array,int size);
double standard_deviation_of_vector(vector<double> array,int size);
double occupied_sites_of_vector(vector<double> array,int size);
auto meansq_spread_of_vector(vector<double> array, int g, int c_x, int c_y);

void var_mean_incremental_surv_runs(double t_avg_var_rep_N[][Sp4_1], double X_curr[][Sp2], int size);

void var_mean_incremental_all_runs(double rep_avg_var[][3], double X_curr[][2], int size, int r);
void spreading_incremental_surv_runs(double t_N_R2_rep[][4], double X_curr[][3], int size);
void var_mean_incremental(double rep_avg_var[][2], vector<vector <double>> &X_curr, int size, int r);

// --------------------------- Partitioning functions to make life easier -----------------------------------//

template<typename T> std::vector<double> linspace(T start_in, T end_in, int num_in);
template<typename T>std::vector<double> lnspace(T start_in, T end_in, int log_points);
std::vector<double> logarithmic_time_bins(double t_max, double dt);
std::vector<double> logarithm10_time_bins(double t_max, double dt);

//----------------------------- Misc. Supporting Add-ons -------------------------------------------------

int theborderwall(D2Vec_Double &Rho_t, int g);
void determine_neighbours_R2( int g, int S, D3Vec_Int &neighbours_R2);
void determine_neighbours_Sq8(int g, int S, D3Vec_Int &neighbours_Sq8);

//----------------------------- Stochastic 2PAC Dornic Integration Machinery --------------------------------------------//

void RK4_Integrate(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, double a, double b, double D[], double A[Sp][Sp], 
	double H[Sp][Sp], double E[], double M[], double t, double dt, double dx, int g);
void tupac_percolationDornic_2D(D2Vec_Double &Rho, vector <double> &t_meas, D2Vec_Double &Rh0, double t_max, double a, double b, double c,
	  double D[], double sigma[], double A[Sp][Sp], double H[Sp][Sp], double E[], double M[], double dt, double dx, int r,  int g);

void first_order_critical_exp_delta(D2Vec_Double &Rh0, int div, double t_max, double a_start, double a_end, double b, double c, 
	double D[], double sigma[], double A[Sp][Sp], double H[Sp][Sp], double E[], double M[], double dt, double dx, int r,  int g);

#endif


// ALT INTEGRATION MACHINERY WHERE MIXING OF TIMESTEPS DUE TO LAPLACIAN MAY BE AN ISSUE.

/**

void f_2D(auto &f, auto &Rho_M, double a, double b, double D[], double A[][Sp], 
	double H[][Sp], double E[], double t, double dt, double dx, int g)
{
	//Vector function that updates an array containing (dR/dt, dC1/dt)

	for(int i=0; i < g*g; i++)
	{
		f[0][i] = -b*Rho_M[0][i]*Rho_M[0][i] 
				- (E[0]*A[0][1]*Rho_M[0][i]*Rho_M[1][i])/(1+ A[0][1]*H[0][1]*Rho_M[0][i]);
		//Equivalent to dR/dt
		f[1][i] = (E[1]*A[0][1]*Rho_M[0][i]*Rho_M[1][i])/(1+ A[0][1]*H[0][1]*Rho_M[0][i]);
	}
	
	

}

void RK4_Integrate(auto &Rho_t, auto &Rho_tsar, double a, double b, double D[], double A[][Sp], 
	double H[][Sp], double E[], double M[], double t, double dt, double dx, int g)
{
	// Rho_tstar stores value of rho*, provided by sampling Gaussian.
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	double dt6 = dt/6.0; double dt2 = dt/2.0;

	double K1[Sp][g*g] ={0.0}; double K2[Sp][g*g] ={0.0}; double K3[Sp][g*g] ={0.0}; double K4[Sp][g*g] ={0.0};
	double Rho_M[Sp] ={0.0, 0.0};


	f_2D(K1, Rho_tsar, a, b, D, A, H, E, t, dt, dx, g); //K1 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K1[s][i];
	}

	f_2D(K2, Rho_M, a, b, D, A, H, E, t + dt2, dt, dx, g); //K2 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K2[s][i];
	}
	f_2D(K3, Rho_M, a, b, D, A, H, E, t + dt2, dt, dx, g); //K3 updated.
	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt)*K3[s][i];
	}
	f_2D(K4, Rho_M, a, b, D, A, H, E, t + dt, dt, dx, g); //K4 updated.

		
		
	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
		{

		
			Rho_t[s][i]+= (dt6)*( K1[s][i] + 2.0*K2[s][i] + 2.0*K3[s][i] + K4[s][i]);
	

			if( Rho_t[s][i] < 0 || isfinite(Rho_t[s][i]) == false || Rho_t[s][i] > 5)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "RK4 WAS KO'ED WITH :\t" << Rho_t[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " For Species:\t:" << s << " with K1[s][i]:   " << K1[s][i] << "\t, K2[s][i]:\t" << K2[s][i] 
				<< "\t, K3[s][i]:\t" << K3[s][i] << "\t AND K4[s][i]:\t" << K4[s][i] << endl; //cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			
			
		}
	}
}

**/