#include "rietkerk_bjork_basic.h"

int main()
{
  increase_stack_limit(1024L); //Increase stack limit to 1024 MB.


  double a =0.7;
  double c, gmax, alpha, d, rW, W0, t_max, dt, dx; //int g, div;
  double a_start, a_end; double r; double dP; // Kick for high initial state
  int g, div;

  //c= 10; d =0.25; gmax= 0.05; alpha =0.2; W0= 0.2; rW = 0.2; // From Bonachela et al 2015
  /**
  dx= 0.2 ; //From Bonachela et al 2015 (in m)
  d0 = 0.001; d1=0; d2 = 0.001; d3= 0.1; //From Bonachela et al 2015 (in m^2/day)
  k0= 0; k1 = 5; k2 =5;
  s0 = sqrt(0.025); // ~ D/(dx)^2 (in g^{0.5}/(m day))
  */

  // Using units of mass = kg, length = km, time = hr.

  c = 10000; gmax = 0.05*pow(10, -3.0)/24.0; d = 0.25/24.0; alpha =0.2/24.0; W0 = 0.2; rW = 0.2/24.0; // From Bonachela et al 2015

  //double p0i = 0.5; double p0j= p0i/200; double p0j= 2.25; double p0m= 8; // In g/m^2

  double k0, k1, k2; double d0, d1, d2, d3; double s0, s1; double v1;//Slightly 
  
  dx= 0.1 ; //From Bonachela et al 2015 (in km)
  d0 = 0.00025/24.0; d1=0.0298; d2 = 0.00025/24.0; d3= 0.025/24.0; //From Bonachela et al 2015 (in km^2/hr)
  k0= 0; k1 = 5; k2 =5000;
  s0 = sqrt(d0/(dx*dx)); // ~ D0/(dx)^2 (in kg^{0.5}/(km hr))
  s1 = 1; // ~ D/(dx)^2 (in kg^{0.5}/(km hr))
  

  v1 = 0.416582; //In km/hr

  double D[Sp] ={d0, d1, d2, d3}; //Diffusion coefficients for species.
  double K[Sp] ={k0, k1, k2}; //Diffusion coefficients for species.
  double sigma[Sp] ={s0, s1, 0, 0}; //Demographic stochasticity coefficients for species.
  double v[Sp] ={0, v1, 0, 0}; //Velocity of species.

  /**

  double beta = c*gmax - d;
  double p0istar, p0jstar, p0mstar; // Analytic steady state values.

  p0jstar = d*K[1]/beta; 
  p0istar = (c/d)*(R - rW*p0jstar);
  p0mstar = (R/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
  */
 double aij, hij, ej, m;
  aij = 3.6*pow(10.0, -6.08)*pow(20.0, -0.37); // in km^2/(hr kg)
  hij = 1; // Handling time in hrs
  ej =0.45; m = 0.061609*pow(20.0, -0.25)/8760.0; // Mortality rate in hr^{-1}
  double H[Sp][Sp] ={{0, hij, 0, 0}, 
                     {hij, 0.0, 0, 0}, 
                     {0, 0, 0, 0}, 
                     {0, 0, 0, 0}};    // Handling time matrix [T]. No canabilism, symmetric effects.
  double A[Sp][Sp] ={{0, aij, 0, 0}, 
                      {aij, 0.0, 0, 0}, 
                      {0, 0, 0, 0}, 
                      {0, 0, 0, 0}};  // Attack rate matrix []. No canabilism, symmetric effects. 
  double M[Sp] = {d, m};     // Mortality Rate of Species. (in hr^{-1})
  double E[Sp] ={1.0, ej}; //Efficiency of consumption.
  //double D[Sp] ={d0, d1}; //Diffusion coefficients for species. (in km^2/hr)
  double pR[Sp] ={0.0, 0.51}; //Perception rate of species (in km)

  cout << "Presets are as follows:\n";
  cout << "For the vegetation:\n";
  cout << "c = " << c << "\t gmax = " << gmax << "\t d = " << d << "\t alpha = " << alpha << "\t W0 = " << W0 << "\t rW = " << rW << endl;
  cout << "For the grazer:\n";
  cout << "aij = " << aij << "\t hij = " << hij << "\t ej = " << ej << "\t m = " << m << endl;
  cout << "\n System-level details:\n";
  cout << "pi = " << PI  << endl;

  cout << "Diffusion Constants For Species:\n";
  cout << "D0 = " << setprecision(16) << D[0] << "\t D1 = " << D[1] << "\t D2 = " << setprecision(16) << D[2] << "\t D3 = " << setprecision(16) << D[3] << endl;

  cout << "Enter desired time-step (~0.1 is a good choice): ";
  cin >> dt;

  cout << "Enter maximum duration of simulation (in hrs): ";
  cin >> t_max;

  cout << "Enter the desired grid-size (L) (dx = 100 m): ";
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

  cout << "Enter kick for high state: ";
  cin >> dP;

  cout << "\n System-level details:\n";
  cout << "dx = " << dx << "\t dt = " << dt << "\t t_max = " << t_max << "\t g = " << g << "\t r = " << r << endl;

  if( 2*v[1]*dt/dx > 1)
  {
    cout << "Warning: CFL condition violated. The velocity of the grazer is too high for the given dx and dt. Please reduce the velocity and dt or increase dx ." << endl;
  }
  else
  {
    cout << "CFL condition satisfied." << endl;
  }

  //cout << "Enter 1 to read from csv or 0 to input standardised initial condions"

  /* cout << "Enter initial density of Primary Producer (kg/m^2): ";
  cin >> p0i;

  cout << "Enter initial density of Primary Consumer (kg/m^2): ";
  cin >> p0j;

	cout << "Enter initial density of Secondary Producer (kg/m^2): ";
  cin >> p0m; */

  //add_two(3.14, 1.78);
  cout << maxis(3.14, 1.78) <<endl;

  int go =8;
  add_three(1,2,go);

  //D2Vec_Double Rho_0(Sp, vector<double> (g*g, 1.0));
  // Rho_0 is 2D [Spx(L*L)] vector initialised to 1.00.

  //init_fullframe(Rho_0, Sp, g*g); //Returns Rho_0 with a full initial frame filled with ones.
  //double p0i = 1.0; double p0j= 0.05;`
  //double mean[Sp] = {p0i, p0j}; double sd[Sp] = {p0i/4.0, p0j/4.0};
	//init_randframe(Rho_0, Sp,  g*g, mean, sd); //Returns Rho_0 with a full initial frame filled with 0.2.
  //double mean[Sp] = {p0i, p0j, p0m}; double sd[Sp] = {p0i/2.0, p0j/2.0, p0m/2.0};
  //init_randframe(Rho_0, Sp,  g*g, mean, sd); //Returns Rho_0 with a full random frame.

  //double perc = 0.015; double c_high[Sp] ={p0i, p0j, p0m}; double c_low[Sp] ={p0i, p0j, p0m};

  //init_randconstframe(Rho_0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.

  //init_constframe(Rho_0, Sp,  g*g, mean); //Returns Rho_0 with a full initial frame filled with 0.2.

  stringstream ast, est, dPo; ast << a_start; est  << a_end; dPo << dP;

  stringstream foldername;
	foldername << "../Data/Rietkerk/Frames/Stochastic/2Sp/" << ast.str() << "-" << est.str() << "_dP_" << dPo.str() << "/";
	// Creating a string stream instance to store the values of the parameters in the file name.
	// This is done to avoid overwriting of files.

	struct stat info2;
	if( stat( foldername.str().c_str(), &info2 ) != 0 )
	{
		cout << "Cannot access " << foldername.str() << ". Creating directory." << endl;
		const int dir_err = system(("mkdir " + foldername.str()).c_str());
		if (-1 == dir_err)
		{
			printf("Error creating directory!n");
			exit(1);
		}
	}

  stringstream foldername2;
	foldername2 << "../Data/Rietkerk/Prelims/Stochastic/2Sp/" << ast.str() << "-" << est.str() << "_dP_" << dPo.str() << "/";
	// Creating a string stream instance to store the values of the parameters in the file name.
	// This is done to avoid overwriting of files.

	struct stat info1;
	if( stat( foldername2.str().c_str(), &info1 ) != 0 )
	{
		cout << "Cannot access " << foldername2.str() << ". Creating directory." << endl;
		const int dir_err = system(("mkdir " + foldername2.str()).c_str());
		if (-1 == dir_err)
		{
			printf("Error creating directory!n");
			exit(1);
		}
	}

    
  stringstream foldername3;
	foldername3 << "../Data/Rietkerk/Stochastic/2Sp/";
	// Creating a string stream instance to store the values of the parameters in the file name.
	// This is done to avoid overwriting of files.

	struct stat info0;
	if( stat( foldername3.str().c_str(), &info0 ) != 0 )
	{
		cout << "Cannot access " << foldername3.str() << ". Creating directory." << endl;
		const int dir_err = system(("mkdir " + foldername3.str()).c_str());
		if (-1 == dir_err)
		{
			printf("Error creating directory!n");
			exit(1);
		}
	}
  //first_order_critical_exp_delta_stochastic(div, t_max, a_start, a_end, c, gmax, alpha, d, rW, W0, D, K, sigma, dt, dx, dP, r, g);
  first_order_critical_exp_delta_stochastic_2Sp(div, t_max, a_start, a_end, c, gmax, alpha, rW, W0, D, v, K, sigma, A, H, E, M, pR,dt, dx, dP, r, g);

  //tupac_percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, auto &Rh0,
	 //double t_max, double a[], double b[], double c[], double D[], double sigma[], double dt, double dx, int r,  int g)

  return 0;
}