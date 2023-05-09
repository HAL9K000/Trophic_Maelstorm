#include "rietkerk_berserk_basic.h"

int main()
{
  increase_stack_limit(2048L); //Increase stack limit to 512 MB.


  double a =0.575;
  double c, gmax, alpha, d, rW, W0, t_max, dt, dx; //int g, div;
  double a_start, a_end; double r; double dP; // Kick for high initial state
  int g, div;

  c= 10; d =0.25; gmax= 0.05; alpha =0.2; W0= 0.2; rW = 0.2; // From Bonachela et al 2015

  double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2

 double k0, k1, k2; double d0, d1, d2; //Slightly 
  
  dx= 0.2 ; //From Bonachela et al 2015 (in m)
  d0 = 0.001; d1 = 0.001; d2= 0.1; //From Bonachela et al 2015 (in m^2/day)
  k0= 0; k1 = 5; k2 =5;

  double D[Sp] ={d0, d1, d2}; //Diffusion coefficients for species.
  double K[Sp] ={k0, k1, k2}; //Diffusion coefficients for species.
  

  cout << "Enter desired time-step: ";
  cin >> dt;

  cout << "Enter maximum duration of simulation (in days): ";
  cin >> t_max;

  cout << "Enter the desired grid-size (L) (dx = 0.2 m): ";
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

  D2Vec_Double Rho_0(Sp, vector<double> (g*g, 1.0));
  // Rho_0 is 2D [Spx(L*L)] vector initialised to 1.00.

  //init_fullframe(Rho_0, Sp, g*g); //Returns Rho_0 with a full initial frame filled with ones.
  //double p0i = 1.0; double p0j= 0.05;`
  //double mean[Sp] = {p0i, p0j}; double sd[Sp] = {p0i/4.0, p0j/4.0};
	//init_randframe(Rho_0, Sp,  g*g, mean, sd); //Returns Rho_0 with a full initial frame filled with 0.2.
  //double mean[Sp] = {p0i, p0j, p0m}; double sd[Sp] = {p0i/2.0, p0j/2.0, p0m/2.0};
  //init_randframe(Rho_0, Sp,  g*g, mean, sd); //Returns Rho_0 with a full random frame.

  double perc = 0.015; double c_high[Sp] ={dP, p0j, p0m}; double c_low[Sp] ={p0i, p0j, p0m};

  init_randconstframe(Rho_0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.

  //init_constframe(Rho_0, Sp,  g*g, mean); //Returns Rho_0 with a full initial frame filled with 0.2.

  stringstream ast, est, dPo; ast << a_start; est  << a_end; dPo << dP;

  stringstream foldername;
	foldername << "../Data/Rietkerk/Frames/" << ast.str() << "-" << est.str() << "_dP_" << dPo.str() << "/";
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
	foldername2 << "../Data/Rietkerk/Prelims/" << ast.str() << "-" << est.str() << "_dP_" << dPo.str() << "/";
	// Creating a string stream instance to store the values of the parameters in the file name.
	// This is done to avoid overwriting of files.

	struct stat info1;
	if( stat( foldername2.str().c_str(), &info2 ) != 0 )
	{
		cout << "Cannot access " << foldername2.str() << ". Creating directory." << endl;
		const int dir_err = system(("mkdir " + foldername2.str()).c_str());
		if (-1 == dir_err)
		{
			printf("Error creating directory!n");
			exit(1);
		}
	}
  

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

  first_order_critical_exp_delta(Rho_0, div, t_max, a_start, a_end, c, gmax, alpha, d, rW, W0, D, K, dt, dx, dP, r, g);

  //tupac_percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, auto &Rh0,
	 //double t_max, double a[], double b[], double c[], double D[], double sigma[], double dt, double dx, int r,  int g)

  return 0;
}