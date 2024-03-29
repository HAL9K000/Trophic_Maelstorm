
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

#include <limits>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <list>
#include <thread>
#include <mutex>

using namespace std::chrono;
//static std::mt19937_64 rng(time(NULL));
#include<bits/stdc++.h>
using namespace std;
#include <sys/resource.h>
#include <omp.h>

#include <map>

#ifndef BERSERK_DYNAMICS_H

#define BERSERK_DYNAMICS_H



inline const int Sp = 4; //Number of species in the system.
inline const int Sp4_1 = 4*(Sp -2) +3; // Used for generating statistics on surviving runs
inline const int Sp2 = 2*Sp; // Used for generating statistics on surviving runs


//------------------------- TYPEDEFS --------------------------------------------------------------//

typedef std::vector<std::vector <double>> D2Vec_Double;
typedef std::vector<std::vector <int>> D2Vec_Int;

typedef std::vector<std::vector <std::vector <double>>> D3Vec_Double;
typedef std::vector<std::vector <std::vector <int>>> D3Vec_Int;

typedef std::vector<std::vector <std::vector <std::vector <double>>>> D4Vec_Double;
typedef std::vector<std::vector <std::vector <std::vector <int>>>> D4Vec_Int;





//---------------------Some Primitives ------------------------------------------------------ //

//Stolen from: https://web.archive.org/web/20200629174939/http://timmurphy.org/2013/06/27/template-metaprogramming-in-c/
// Generic exponent compile-time calculations. Pow<2,3>::result == 2^3
template < long long B, long long E>
struct Pow
{
    static const long long result = B * Pow<B, E-1>::result;
};

template <long long B>
struct Pow<B, 0>
{
    static const long long result = 1;
};

// Bailey-Borwein-Plouffe formula for calculating pi
// http://en.wikipedia.org/wiki/Bailey-Borwein-Plouffe_formula
template <long long int N>
struct CalculatePi
{
    static const constexpr double pi =
    (
        1.0/Pow<16,N>::result *
        (
            4.0/(8*N + 1.0) - 2.0/(8*N + 4.0) -
            1.0/(8*N + 5.0) - 1.0/(8*N + 6.0)
        )
    ) + CalculatePi<N-1>::pi;
};

template <>
struct CalculatePi<-1>
{
    static const constexpr  double pi = 0;
};

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

inline const double PI = CalculatePi<14>::pi;

template <typename T> int sgn(T val);
void increase_stack_limit(long long stack_size);
bool maxis(int a, int b);
void add_three(int a, int b, int c); //Test function.
// Custom comparator for sorting pairs based on the squared distance
bool comparePairs(const pair<int, int>& a, const pair<int, int>& b); 

void init_fullframe(D2Vec_Double &array, int Sp, int size);
void init_csvconstframe(D2Vec_Double &array, D2Vec_Double &const_ind_val, const std::string& filename, const vector<int> &columns, int size);
void init_randconstframe(D2Vec_Double &array, int Sp, int size, double perc,   double c_high[], double c_low[]);
void init_constframe(D2Vec_Double &array,  int Sp, int size, double constant[]);
void init_randframe(D2Vec_Double &array, int Sp, int size, double mean[], double sd[]);
void init_quarterframe(D2Vec_Double &array, double c1, double c2, double c3, double c4, int L);
void init_solitarytear(D2Vec_Double &array, int length);

double mean_of_array(double array[],int size);
double standard_deviation_of_array(double array[],int size);
double mean_of_vector(vector<double> array,int size);
double standard_deviation_of_vector(vector<double> array,int size);
double occupied_sites_of_vector(vector<double> array,int size);
auto meansq_spread_of_vector(vector<double> array, int g, int c_x, int c_y);

void var_mean_incremental_surv_runs(double t_avg_var_rep_N[][Sp4_1], double X_curr[][Sp2], int size, int j);

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
void determine_neighbours_Sq4(int g, D3Vec_Int &neighbours_Sq4);
void determine_neighbours_Sq8(int g, int S, D3Vec_Int &neighbours_Sq8);
std::vector<std::pair<int, int>> computeNeighboringSitesCentral(int R);
std::vector<int> getNeighborsInBall(const vector<int>& sortedNeighboringSites, double f);
std::vector<int> generateNeighboringSitesFromCentral(const vector<pair<int, int>>& centralNeighboringSites, int i, int j, int L); 

//----------------------------- Regular Rietkerk Integration Machinery --------------------------------------------//

void f_2D(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, double alpha, double d, double rW, double W0,  
	double D[], double K[],  double t, double dt, double dx, int g);
void RK4_Integrate(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, double a, double b, double D[], double A[Sp][Sp], 
	double H[Sp][Sp], double E[], double M[], double t, double dt, double dx, int g);
void RK4_Wrapper_2D(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, 
	  double alpha, double d, double rW, double W0, double D[], double K[], double a_st, double a_end,  double dt, double dx, double dP, int r, int g);

void first_order_critical_exp_delta(int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double d, double rW, double W0, double D[], double K[], double dt, double dx, double dP, int r,  int g);


//----------------------------- Stochastic Rietkerk-Dornic Integration Machinery --------------------------------------------//

//------------------- Only Vegetation (+ Soil Water + Surface Water) -------------------//

void f_2Dor(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, 
	double alpha, double d, double rW, double W0, double D[], double K[], double t, double dt, double dx1_2, double g);
void RK4_Integrate_Stochastic(D2Vec_Double &Rho_dt, D2Vec_Double &Rho_tsar, D3Vec_Int &nR2,double a,double c,double gmax,double alpha,
		double d, double rW, double W0, double D[], double K[],double t,double dt,double dx, int g);
void rietkerk_Dornic_2D(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double d, 
	  double rW, double W0, double D[], double K[], double sigma[], double a_st, double a_end,  double dt, double dx, double dP, int r, int g);
void first_order_critical_exp_delta_stochastic(int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double d, double rW, double W0,  double D[], double K[], double sigma[], double dt, double dx, double dP, int r,  int g);

//------------------- Vegetation + Grazer (+ Soil Water + Surface Water) -------------------//

void f_2Dor_2Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, 
	double alpha, double rW, double W0, double D[], double K[], double A[Sp][Sp], double H[Sp][Sp], double E[], double t, double dt, double dx1_2, double g);
void RK4_Integrate_Stochastic_2Sp(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D3Vec_Int &nR2,double a,double c,double gmax,double alpha,
	double rW, double W0, double D[], double K[], double A[Sp][Sp], double H[Sp][Sp], double E[], double t,double dt,double dx, int g);
void rietkerk_Dornic_2D_2Sp(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double rW, double W0, 
	double D[], double v[], double K[], double sigma[], double a_st, double a_end, double A[Sp][Sp], double H[Sp][Sp], double E[], double M[], double pR[], 
	double dt, double dx, double dP, int r, int g);
void first_order_critical_exp_delta_stochastic_2Sp(int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double rW, double W0,  double D[], double v[], double K[], double sigma[], double A[Sp][Sp], double H[Sp][Sp], double E[], double M[], double pR[],
	double dt, double dx, double dP, int r,  int g);

//------------------- Vegetation + Grazer + Predator (+ Soil Water + Surface Water) -------------------//

void f_2Dor_3Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, 
	double alpha, double rW, double W0, double D[], double K[], double A[Sp][Sp], double H[Sp][Sp], double E[], double t, double dt, double dx1_2, double g);
void RK4_Integrate_Stochastic_3Sp(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D3Vec_Int &nR2,double a,double c,double gmax,double alpha,
	double rW, double W0, double D[], double K[], double A[Sp][Sp], double H[Sp][Sp], double E[], double t,double dt,double dx, int g);
void rietkerk_Dornic_2D_3Sp(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double rW, double W0, 
	double D[], double v[], double K[], double sigma[], double a_st, double a_end, double A[Sp][Sp], double H[Sp][Sp], double E[], double M[], double pR[], 
	double dt, double dx, double dP, int r, int g);
void first_order_critical_exp_delta_stochastic_3Sp(int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double rW, double W0,  double D[], double v[], double K[], double sigma[], double A[Sp][Sp], double H[Sp][Sp], double E[], double M[], double pR[],
	double dt, double dx, double dP, int r,  int g);


#endif