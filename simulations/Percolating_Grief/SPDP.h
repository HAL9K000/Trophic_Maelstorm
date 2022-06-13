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
#include <list>
#include <thread>
using namespace std::chrono;
//static std::mt19937_64 rng(time(NULL));
#include<bits/stdc++.h>
using namespace std;
#include <sys/resource.h>
#include <omp.h>

#include <map>

#ifndef PERCOL_DYNAMICS_H
#define PERCOL_DYNAMICS_H

//---------------------Some Primitives ------------------------------------------------------ //

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

//----------------------------- Stochastic SPDP RK4 Integration Machinery --------------------------------------------//

void find_neighbours_R2(double neighbours_R2[2][5], vector<double> &Rho_t, int p, int q, int g);
void f_2D(double f[], vector<double> &Rho_t, double p, double D, double t, double dt, double dx, int g);
void percolation_2D(vector<vector<double>> &Rho, vector<double> &t_stop, double Rh0[],
  double t_max, double p, double D, double sigma, double dt, double dx, int r,  int g);
void stoc_dem_percolation(double Rho_0[], vector<double> &t_stop, int div, double t_max,
	double p_start, double p_end, double D, double sigma, double dt, double dx, int r,  int g);

#endif
