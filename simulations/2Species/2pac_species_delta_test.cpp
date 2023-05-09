#include "2species_stochastic.h"

int main()
{
  increase_stack_limit(2048L); //Increase stack limit to 512 MB.



  double a_start, a_end; double r;
  double dt, dx, t_max; int g, div;

  double p0i = 1.2; double p0j= 0.05; // In kg/m^2
  
  /**
  double mR = 1; double mC = 10; double mP=100; //Mass of creatures. [in kg]
  double a = (1)*1.71*pow(10.0,-6.0)*pow(mR,-0.25)*3600;
  double b = a*1*pow(mR, 0.25); double aij =pow(10,-3.08)*pow(mC,-0.37)*3600;
  double hij = 1*pow(mC,-0.02)/3600.0; double m= 4.15*pow(10,-8.0)*pow(mC,-0.25)*3600;
  double ajm = pow(10,-3.08)*pow(mP,-0.37)*3600; double hjm = 1*pow(mP,-0.02)/3600.0;
  double mm= 4.15*pow(10,-8.0)*pow(mP,-0.25)*3600; double ej =0.45; double em = 0.85;
  double d0 = pow(10,-5.0); double d1 = pow(10,-3.0); double d2 = pow(10,-3.0);
  */

<<<<<<< HEAD
 double a =1; double b= 0.5/1.349; double aij, hij, ej, m; double d0, d1; //Slightly 
 aij= 1.161; hij = 0.000389; m = 0.000007605; ej =0.1; 
 
  d0 = 0.0012; d1 = 48420;
=======
 double a =1; double b=1; double aij, hij, ej, m; double d0, d1; //Slightly 
 aij= 0.25; hij = 0.25; m = 0.1; ej =0.5; d0 = 0.1; d1 = 0.5;
>>>>>>> c1d9a8cc3e7948e831a134e960943d98df4f6152

  double H[Sp][Sp] ={0, hij, 
                    hij, 0.0};    // Handling time matrix [T]. No canabilism, symmetric effects.
  double A[Sp][Sp] ={0, aij, 
                     aij, 0.0};  // Attack rate matrix []. No canabilism, symmetric effects.
  double M[Sp] = {0.0, m};     // Mortality Rate of Species.
  double E[Sp] ={1.0, ej}; //Efficiency of consumption.
  double D[Sp] ={d0, d1}; //Diffusion coefficients for species.

<<<<<<< HEAD
  double sigma[Sp] ={sqrt(4.8*pow(10.0, -7.0)),sqrt(2.0)}; dx=50.0;
=======
  double sigma[Sp] ={sqrt(1.0),sqrt(1.0)}; dx=1.0;
>>>>>>> c1d9a8cc3e7948e831a134e960943d98df4f6152

  cout << "Enter desired time-step: ";
  cin >> dt;

  cout << "Enter maximum duration of simulation (in hr): ";
  cin >> t_max;

  cout << "Enter the desired grid-size (L): ";
  cin >> g;

  cout << "Enter the desired number of replicates (r): ";
  cin >> r;

  cout << "Note that a relevant value of 'a' may be construed as:\t" << setprecision(10) << a <<endl;

  cout << "Enter starting a-value: ";
  cin >> a_start;

  cout << "Enter ending a-value: ";
  cin >> a_end;

  cout << "Enter number of divisions of p (n+1): ";
  cin >> div;

  //cout << "Enter 1 to read from csv or 0 to input standardised initial condions"

  /* cout << "Enter initial density of Primary Producer (kg/m^2): ";
  cin >> p0i;

  cout << "Enter initial density of Primary Consumer (kg/m^2): ";
  cin >> p0j;

	cout << "Enter initial density of Secondary Producer (kg/m^2): ";
  cin >> p0m; */

  //add_two(3.14, 1.78);
  cout << maxis(3.14, 1.78) <<endl;

  int c =8;
  add_three(1,2,c);

  D2Vec_Double Rho_0(Sp, vector<double> (g*g, 1.0));
  // Rho_0 is 2D [Spx(L*L)] vector initialised to 1.00.

  //init_fullframe(Rho_0, Sp, g*g); //Returns Rho_0 with a full initial frame filled with ones.
  //double p0i = 1.0; double p0j= 0.05;`
<<<<<<< HEAD
  //double mean[Sp] = {p0i, p0j}; double sd[Sp] = {p0i/4.0, p0j/4.0};
	//init_randframe(Rho_0, Sp,  g*g, mean, sd); //Returns Rho_0 with a full initial frame filled with 0.2.
  double mean[Sp] = {p0i, p0j};
  init_constframe(Rho_0, Sp,  g*g, mean); //Returns Rho_0 with a full initial frame filled with 0.2.
=======
  double mean[Sp] = {p0i, p0j}; double sd[Sp] = {p0i/4.0, p0j/4.0};
	init_randframe(Rho_0, Sp,  g*g, mean, sd); //Returns Rho_0 with a full initial frame filled with 0.2.

>>>>>>> c1d9a8cc3e7948e831a134e960943d98df4f6152
  

  cout << " A subset of the initial frame:\t" << endl;
  for(int s =0; s < Sp; s++)
  {
      cout << " For SPECIES:\t" << s << endl;
      for(int i = 0; i < 10; i++)
      {
	    for(int j =0; j < 10; j++)
		    cout << setprecision(3) << Rho_0[s][i*g + j] << " ";
	    cout << endl;
      }
  }


  first_order_critical_exp_delta(Rho_0, div, t_max, a_start, a_end, b, c, D, sigma, A,H,E,M, dt, dx, r, g);

  //tupac_percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, auto &Rh0,
	 //double t_max, double a[], double b[], double c[], double D[], double sigma[], double dt, double dx, int r,  int g)

  return 0;
}