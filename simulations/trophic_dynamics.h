#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>

#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
using namespace std::chrono;
static std::mt19937_64 rng(time(NULL));
#include<bits/stdc++.h>
using namespace std;
#include <sys/resource.h>
#include <omp.h>

#include <map>

#ifndef TROPHIC_DYNAMICS_H
#define TROPHIC_DYNAMICS_H

//----------------------- Primitives ---------------------------------//

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

void increase_stack_limit(long stack_size);
double mean_of_array(double array[],int size);
double standard_deviation_of_array(double array[],int size);
double mean_of_vector(vector<double> array,int size);
double standard_deviation_of_vector(vector<double> array,int size);
template<typename T> std::vector<double> linspace(T start_in, T end_in, int num_in);

//----------------------------- Theoretical Constructs --------------------------------------------//

vector<double> rho_eqpk(map<string, double> &c);
vector<double> rhono_eq(map<string, double> &c);
int check_stability(map<string, double> &c, string flag);

//----------------------------------- Basic Integration Schema ------------------------------------//

void f(double f[], double Pijm[], double t, map<string, double> &c);
void simulateRK4(vector<vector <double>> &Pijm, vector <double> &t_list, float t, double dt, map<string, double> &c);

//------------------------- Basic Trophic Chain Simulation --------------------------------------//

void topical_trophics(map<string, double> &c, vector<vector<double>> &Pijm, double dt, float t);

//------------------------- Spatial Trophic Chain Simulation --------------------------------------//

void find_neighbours_R2(double neighbours_R2[][2][5], vector<vector<double>> &Pijm, int n, int p, int q, int g);
void f_2D(double f[], double Pijm_f[], vector<vector<double>> &Pijm, int n, int p, int q, int g, double t,double dt, map<string, double> &c);
void simulate_RK42D(vector<vector<double>> &Pijm,  vector<vector<double>> &Pij, double t, double dt, int p, int q, int g, map<string, double> &c);
void topical_trophics_2D(map<string, double> &c, vector<vector<double>> &Pij, vector<double> &t_stop, int div, double t_max, double dt,  int g);


#endif
