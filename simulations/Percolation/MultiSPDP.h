#include <string>
//#include <format>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>


#include <regex>
#include <limits>
#include <random>
#include <vector>
#include <algorithm>
#include <ranges>
#include <span>
#include <chrono>
#include <cmath>
#include <list>
#include <thread>
#include <mutex>
/** NOTE: The following include is for the exprtk library. Download and add to your project/include path.
// if not already present. The library is used for parsing mathematical expressions. */
#include <exprtk.hpp>

using namespace std::chrono;
//static std::mt19937_64 rng(time(NULL));
#include<bits/stdc++.h>
using namespace std;
#include <sys/resource.h>
#include <omp.h>

#include <map>

#if defined(__GNUC__) && (__GNUC__ >= 9)
#include <filesystem>
namespace fs = std::filesystem;
#endif

#ifndef BERSERK_DYNAMICS_H

#define BERSERK_DYNAMICS_H

#ifndef SPB
	#define SPB 3 //Number of biota species in the system (default value used if NOT defined).
#endif

#ifndef INIT
	#define INIT 0 //Initial condition for the system (0 = Homogeneous, 1 = Random Heterogenous Speckles, 2 = Burn-in)
#endif

#if SPB == 3
	#include "MultiSPDP_constants_3Sp.h"
#elif SPB == 2
	#include "MultiSPDP_constants_2Sp.h"
#elif SPB == 1
	#include "MultiSPDP_constants_1Sp.h"
#else
	#error "Number of species not supported. Pass valid compiler flag as -DSPB=2 or -DSPB=3"
#endif

inline const int SpB = Sp; //Number of biota species in the system.
inline const int Sp_NV = Sp-1; //Number of grazer and predator species in system.
inline const int Sp4_1 = 4*(Sp) +1; // Used for generating statistics on surviving runs
inline const int Sp2 = 2*Sp; // Used for generating statistics on surviving runs

// Global user-defined parameters for the Rietkerk model.
inline double dt2 /** = dt/2.0*/; inline double dt6 /** = dt/6.0*/;
inline double dx2 /** = dx*dx*/; inline double dx1_2 /** = 1.0/(dx*dx)*/; inline double dxb2 /** = dx/2.0/ */; inline double dx1 /** =  1/dx */;

// Global user-defined parameters for the predator-prey trophic chain.
inline double K_P; /** Carrying capacity of predator */
inline double K_P1; /** Inverse of carrying capacity of predator */

// Rietkerk Model Parameters, set in set_global_user_Rietkerk_params() function.
//inline double cgmax /** = c*gmax */; inline double K2W0 /**= K[2]*W0 */;inline double A01H01 /** = A[0][1]*H[0][1] */;inline double A12H12 /**= A[1][2]*H[1][2] **/;

// Strings that are used for I/O operations and file naming in the simulations.
inline string prefix = "DiC-NREF-NI";
inline string frame_folder = "../Data/DP/Frames/Stochastic/"+ std::to_string(SpB) +"Sp/" + prefix + "_"; //Folder to store frames.
inline string prelim_folder = "../Data/DP/Prelims/Stochastic/"+ std::to_string(SpB) +"Sp/"+ prefix +"_"; //Folder to store preliminary data.
inline const string frame_prefix = "/FRAME_P_c_DP_G_"; //Prefix for frame files.
inline const string gamma_prefix = "/GAMMA_G_"; //Prefix for gamma files.
inline const string prelim_prefix = "/PRELIM_AGGRAND_P_c_ID_"; //Prefix for preliminary data files.
inline const string replicate_prefix = "/PRELIM_TSERIES_P_c_DP_G_"; //Prefix for replicate time-series data files.

inline const string input_prefix = "/FRAME_T_"; //Prefix for input files.
inline string stat_prefix = "../Data/DP/Stochastic/"+ std::to_string(SpB) +"Sp/1stCC_"  + prefix + "_P_c_G_";

inline string input_frame_parenfolder = "../Input/DP/"+ std::to_string(SpB) +"Sp"; //Parendirectory for input data.
inline string input_frame_subfolder =""; //Subfolder for input data.
// NOTE: Complete path to input data is input_frame_parenfolder + input_frame_subfolder + input_prefix +"*.csv"
inline std::map<string, string> input_keys = { {"outPRE", ""}, {"a", "0"}, {"a_c", "0"}, {"ascaled", "0"}, {"T", "0"}, {"dP", "0"}, {"Geq", "0"}, {"R", ""}, 
{"MAMP", "0"}, {"SEP", "0"}, {"WID", "0"}, {"MSF", "0"}};

//inline const string prelimheader = " a , r, L, t , <<P(x; t)>_x>_r, Var[<P(x; t)>_x]_r, # Surviving Runs P(x; t),"
//	" # Active Sites P(x; t), <<G(x; t)>_x>_r, Var[<G(x; t)>_x]_r, # Surviving Runs G(x; t), # Active Sites G(x; t), \n";

// Strings that are used for the Expertk library to parse mathematical expressions.
// These represent MFT-versions of the species Coexistance equilibria wrt to the order parameter (Rainfall, given by a).

inline vector<string> MFT_Vec_CoexExpr(2*Sp, ""); //Vector of MFT Coexistance expressions for all species.

//inline exprtk::symbol_table<double> global_symbol_table; //Symbol table for the Expertk library.


//------------------------- TYPEDEFS --------------------------------------------------------------//

typedef std::vector<std::vector <double>> D2Vec_Double;
typedef std::vector<std::vector <int>> D2Vec_Int;

typedef std::vector<std::vector <std::vector <double>>> D3Vec_Double;
typedef std::vector<std::vector <std::vector <int>>> D3Vec_Int;

typedef std::vector<std::vector <std::vector <std::vector <double>>>> D4Vec_Double;
typedef std::vector<std::vector <std::vector <std::vector <int>>>> D4Vec_Int;

typedef std::vector<double> Vec_Double;
typedef std::vector<int> Vec_Int;





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

// Generic template to 

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

//void set_global_user_Rietkerk_params(double c,double gmax,double alpha, double rW, double W0, double (&D)[Sp], 
//	double (&K)[3], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB])

template <typename T> int sgn(T val);
void increase_stack_limit(long long stack_size);
bool maxis(int a, int b);
string format_str(const std::string& s, const std::map<string, string>& values);
void add_three(int a, int b, int c); //Test function.
void set_Prefix(string& user_prefix);
void set_input_Prefix(string& user_inputframepath, string& user_prefix, double user_a_c, double user_dP, double user_Geq = -1, double user_input_T = -1);
void set_global_system_params(double dt, double dx);
void set_global_predator_params(double Km);

void display_symbol_table(const exprtk::symbol_table<double>& symbol_table);
// Custom comparator for sorting pairs based on the squared distance
bool comparePairs(const pair<int, int>& a, const pair<int, int>& b); 

void init_fullframe(D2Vec_Double &array, int Sp, int size);

void init_randbistableframe(D2Vec_Double &array, int size, double R, double R_c, double perc,  double c_high[], double c_low[]);
// NOTE: The following function REQUIRES the exprtk library to be included in the project.
void init_exprtk_randbiMFTframe(D2Vec_Double &array, int size, double R, double R_c, double dP, double perc,  double c_spread[]);
void init_exprtk_readCSVcolumns_frame(D2Vec_Double &array, vector<int> &const_index, const string &parendir, const string &filenamePattern, 
					const std::vector<std::string> &read_cols, double a, double a_c,  int L, int r, double c_spread[]);

void init_csvconstframe(D2Vec_Double &array, D2Vec_Double &const_ind_val, const std::string& filename, const vector<int> &columns, int size);
void init_burnin_wrapper(D2Vec_Double &Rho_dt, double a, double a_c, int a_scalingfactor, double dP, double perc,  int L, int r, double c_spread[] );
void init_randconstframe(D2Vec_Double &array, int Sp, int size, double perc,   double c_high[], double c_low[]);
void init_constframe(D2Vec_Double &array,  int Sp, int size, double constant[]);
void init_randframe(D2Vec_Double &array, int Sp, int size, double mean[], double sd[]);
void init_quarterframe(D2Vec_Double &array, double c1, double c2, double c3, double c4, int L);
void init_gaussframe(D2Vec_Double &array, int size, vector<double> &sd, vector<double> &amp, int s_gauss = SpB);
void init_gaussiantear(D2Vec_Double &array, int Sp_lim, int g, double cx[], double cy[], const vector<double> &sd,const vector<double> &amp);
void init_solitarytear(D2Vec_Double &array, int length);

double mean_of_array(double array[],int size);
double standard_deviation_of_array(double array[],int size);
double mean_of_vector(vector<double> array,int size);
double standard_deviation_of_vector(vector<double> array,int size);
double occupied_sites_of_vector(vector<double> array,int size);
auto meansq_spread_of_vector(vector<double> array, int g, int c_x, int c_y);

void var_mean_incremental_surv_runs(D2Vec_Double &t_avg_var_rep_N, const D2Vec_Double &X_curr, int size, int j);

void var_mean_incremental_all_runs(double rep_avg_var[][3], double X_curr[][2], int size, int r);
void spreading_incremental_surv_runs(double t_N_R2_rep[][4], double X_curr[][3], int size);
void var_mean_incremental(double rep_avg_var[][2], vector<vector <double>> &X_curr, int size, int r);

// --------------------------- Partitioning functions to make life easier -----------------------------------//

template<typename T> std::vector<double> linspace(T start_in, T end_in, int num_in);
template<typename T>std::vector<double> lnspace(T start_in, T end_in, int log_points);
std::vector<double> logarithmic_time_bins(double t_max, double dt);
std::vector<double> logarithm10_time_bins(double t_max, double dt);

//----------------------------- Misc. Supporting Add-ons -------------------------------------------------

void recursive_dir_create(const string& path_to_dir);
int findMaxRepNo(const string& parendir, const string& filenamePattern);
int theborderwall(D2Vec_Double &Rho_t, int g);
void determine_neighbours_R2( int g, int S, D3Vec_Int &neighbours_R2);
void determine_neighbours_Sq4(int g, D3Vec_Int &neighbours_Sq4);
void determine_neighbours_Sq8(int g, int S, D3Vec_Int &neighbours_Sq8);
void set_density_dependentmortality(double m);
std::vector<std::pair<int, int>> computeNeighboringSitesCentral(int R);
std::vector<int> getNeighborsInBall(const vector<int>& sortedNeighboringSites, double f);
void generateNeighboringSitesFromCentral(const auto range_CentralNeighboringSites, 
										 vector<int>& neighboringSites, int i, int j, int L);
//void generateNeighboringSitesFromCentral(const vector<pair<int, int>>& centralNeighboringSites, 
//                                         vector<int>& neighboringSites, int i, int j, int L); 

//----------------------------- Regular DP Integration Machinery --------------------------------------------//

//----------------------------- Stochastic DP-Dornic Integration Machinery --------------------------------------------//

void save_prelimframe(D2Vec_Double &Rho_t, const string &parendir, const string &filenamePattern, double a, double a_st, double a_end, 
		double t, double dt, double dx, double dP, int r, int g, string header ="", bool overwrite = false, bool delete_previous = false);
void save_frame(D2Vec_Double &Rho_t, const string &parendir, const string &filenamePattern, double a, double a_st, double a_end, 
		double t, double dt, double dx, double dP, int r, int g, string header ="", bool overwrite = false);

//------------------- Only Vegetation -------------------//

void f_DP_Dor_1Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double b, double c, double t, double dt, double g);

//------------------- Vegetation + Grazer-------------------//

void calc_gamma_2Sp_NonRefugia(const vector<pair<int, int>>& centralNeighboringSites, D2Vec_Double &Rho_t, D2Vec_Double &gamma,  
	double (&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, double nVeg_frac, int r_max, int L);
void f_DP_Dor_2Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double b, double c, 
	double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t, double dt, double g);


//------------------- Vegetation + Grazer + Predator -------------------//


//Calculates gamma for 3 Sp Rietkerk model (details in PDF, 
// ASSUMING PREDATORS AREN'T DETERRED BY HIGH LOCAL VEGETATION DENSITY)
void calc_gamma_3Sp_NonRefugia(const vector<pair<int, int>>& centralNeighboringSites, D2Vec_Double &Rho_t, D2Vec_Double &gamma,  
			double (&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, double nVeg_frac, int r_max, int L);
//Calculates gamma for 3 Sp Rietkerk model (details in PDF)
void calc_gamma_3Sp(const vector<pair<int, int>>& centralNeighboringSites, D2Vec_Double &Rho_t, D2Vec_Double &gamma,  
			double (&Rho_avg)[Sp], vector <std::pair<double, int>>& rfrac, double nVeg_frac, int r_max, int L);
void f_DP_Dor_3Sp(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double b, double c, 
	double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t, double dt, double g);
void RK4_Integrate_Stochastic_MultiSp(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D2Vec_Double &K1, D2Vec_Double &K2, D2Vec_Double &K3, D2Vec_Double &K4, D3Vec_Int &nR2,
	double b, double c, double (&Dxd2)[Sp], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double t,double dt,double dx, int g);
void dP_Dornic_2D_MultiSp(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double b, double c, double (&D)[Sp], double v[],
    double sigma[], double a_st, double a_end, double a_c, double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double pR[], 
	double chigh[], double clow[], double dt, double dx, double dP, int r, int g, double Gstar = -1, double Vstar = -1);
void first_order_critical_exp_delta_stochastic_MultiSp(int div, double t_max, double a_start, double a_end, double a_c, double b,  double c, 
	double (&D)[Sp], double v[], double sigma[], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double pR[], double ch[], double clo[],
	double dt, double dx, double dP, int r,  int g, double Gstar = -1, double Vstar = -1);

#endif