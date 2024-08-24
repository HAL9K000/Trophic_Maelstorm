#include "MultiSPDP.h"

int main(int argc, char *argv[])
{
  increase_stack_limit(1024L); //Increase stack limit to 1024 MB.

  string preFIX; // Prefix for the output files
  //double a_c =1/24.0;
  double a_c, b, c, gmax, alpha, d, rW, W0, t_max, dt, dx; //int g, div;
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
  // ASSUMING MASS OF  GRAZER = 20 kg, MASS OF PREDATOR = 100 kg.
  b = pow(10.0, -6.00); // ~ aij in km^2/(hr kg)

  //double p0i = 0.5; double p0j= p0i/200; double p0j= 2.25; double p0m= 8; // In g/m^2

  double k0, k1, k2; double d0, d1, d2, d3, d4; double s0, s1, s2; double v1, v2;//Slightly 
  
  dx= 0.1 ; //From Bonachela et al 2015 (in km)
  d0 = 0.00025/24.0; d1=0.0298; d2= 0.05221; d3 = 0.00025/24.0; d4= 0.025/24.0; //From Bonachela et al 2015 (in km^2/hr)
  s0 = sqrt(d0/(dx*dx)); // ~ D0/(dx)^2 (in kg^{0.5}/(km hr))
  s1 = 1; // ~ D/(dx)^2 (in kg^{0.5}/(km hr))
  s2 = 3; // ~ D/(dx)^2 (in kg^{0.5}/(km hr))
  
  // Allometric scaling for velocity: 
  // v (m/hr) = 43.706*M^1.772*(1 - e^{-14.27*(M)^(-1.865)}), where M is mass in kg.
  v1 = 0.45427; //In km/hr
  v2 = 0.40587; //In km/hr

  double D[Sp] ={d0, d1}; //Diffusion coefficients for species.
  double sigma[Sp] ={s0, s1}; //Demographic stochasticity coefficients for species.
  double v[Sp] ={0, v1}; //Velocity of species.

  /**

  double beta = c*gmax - d;
  double p0istar, p0jstar, p0mstar; // Analytic steady state values.

  p0jstar = d*K[1]/beta; 
  p0istar = (c/d)*(R - rW*p0jstar);
  p0mstar = (R/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
  */
  double aij, hij, ej, mj; double ajm, hjm, em, mm; // Parameters for the grazer and predator.
  
  aij = 3.6*pow(10.0, -6.08)*pow(20.0, -0.37); // in km^2/(hr kg)
  hij = 1; // Handling time in hrs
  ej =0.45; mj = 0.061609*pow(20.0, -0.25)/8760.0; // Mortality rate in hr^{-1}
  //Predator of mass 100 kg.
  ajm = 3.6*pow(10.0, -6.08)*pow(100.0, -0.37); // in km^2/(hr kg)
  hjm = 1; // Handling time in hrs
  em =0.85; mm = 0.061609*pow(100.0, -0.25)/8760.0; // Mortality rate in hr^{-1}

  

  double Km = pow(10.0, 1.22)*pow(100.0, -0.31); // Carrying capacity of predator in kg/km^2s

  //double Gstar = mm/((em -mm*hjm)*ajm); //Steady state for grazer.
  double H[SpB][SpB] ={{0, hij}, 
                     {hij, 0.0}};    // Handling time matrix [T]. No canabilism, symmetric effects.
  double A[SpB][SpB] ={{0, aij}, 
                      {aij, 0.0}};  // Attack rate matrix []. No canabilism, symmetric effects. 
  double M[SpB] = {d, mj};     // Mortality Rate of Species. (in hr^{-1})
  double E[SpB] ={1.0, ej}; //Efficiency of consumption.
  //double D[Sp] ={d0, d1}; //Diffusion coefficients for species. (in km^2/hr)
  double pR[Sp] ={0.0, 1.04285/2.0}; //Perception radius of species (in km)

  double kappa = (A[0][1]*(E[1] - M[1]*H[0][1])); // Check notation in MFT derivation PDF.
  a_c = M[1]*b/(kappa); // Critical value of a for grazer co-existance.

  
  double Vstar = M[1]/(kappa); // MFT value of V for grazer co-existance.
  double init_frac_grazpred = 1.0; //Initial fraction of grazers and predators in the high state.
  
  cout << "Presets are as follows:\n";
  cout << "For the vegetation:\n";
  cout << "c = " << c << "\t gmax = " << gmax << "\t d = " << d << "\t alpha = " << alpha << "\t W0 = " << W0 << "\t rW = " << rW << endl;
  cout << "\n MFT Vegetation Co-existance Equilibrium Density = " << Vstar << " kg/km^2\n" << endl;
  cout << "For the grazer:\n";
  cout << "aij = " << aij << "\t hij = " << hij << "\t ej = " << ej << "\t mj = " << mj << endl;
  cout << "\nMFT Grazer Co-existance Critical Threshold = " << a_c << " hr^{-1}\n" << endl;
  cout << "For the predator:\n";
  cout << "ajm = " << ajm << " " << A[1][2] << "\t hjm = " << hjm << " " << H[1][2] 
       << "\t em = " << em << " " << E[2] <<  "\t mm = " << mm <<  " " << M[2] << endl;
  cout << "\n System-level details:\n";
  cout << "pi = " << PI  << endl;

  cout << "Diffusion Constants For Species:\n";
  cout << "D0 = " << setprecision(16) << D[0] << "\t D1 = " << D[1] << "\t D2 = " << setprecision(16) << D[2] << "\t D3 = " << setprecision(16) << D[3] << endl;

  cout << "This is a 2Sp (Two) Stochastic DP Model Script WITHOUT DDM\n";

  cout << "\n============================================================\n";

  // There are eight user-defined parameters, as follow below: dt, t_max, g, r, a_start, a_end, div, dP.

  //Check if the user has entered the correct number of arguments. If not prompt the user to enter all the arguments within the script.
  // First check if no arguments are entered. In this case, the user will be prompted to enter all the arguments using cin.

  if(argc == 1)
  {
    cout << "Please enter all the arguments within the script.\n";
    cout << "Enter desired time-step (~0.1 hr is a good choice): ";
    cin >> dt;

    cout << "Enter maximum duration of simulation (in hrs): ";
    cin >> t_max;

    cout << "Enter the desired grid-size (L) (dx = 100 m): ";
    cin >> g;

    cout << "Enter the desired number of replicates (r): ";
    cin >> r;

    cout << "Note that a relevant value of 'a'(Grazer Critical Threshold) may be construed as:\t" << setprecision(10) << a_c <<endl;

    cout << "Enter starting a-value: ";
    cin >> a_start;

    cout << "Enter ending a-value: ";
    cin >> a_end;

    cout << "Enter number of divisions of p (n+1): ";
    cin >> div;

    cout << "Enter kick for high state: ";
    cin >> dP;

    cout << "Enter fractional deviation from MFT for initial conditions of Grazer and Predator above R_c (0.5-1.2 is a good choice): ";
    cin >> init_frac_grazpred;

    cout << "Enter Prefix (common choices include 'DiC-REF', 'DDM-DiC-BURNIN', 'nDiC' etc.) If not provided, default value: " 
        + prefix + " will be used: ";
    cin >> preFIX;

  }
  // If the user has entered the correct number of arguments, then the arguments are read from the command line.
  else if(argc == 11)
  {
    dt = atof(argv[1]);
    t_max = atof(argv[2]);
    g = atoi(argv[3]);
    r = atoi(argv[4]);
    a_start = atof(argv[5]);
    a_end = atof(argv[6]);
    div = atoi(argv[7]);
    dP = atof(argv[8]);
    init_frac_grazpred = atof(argv[9]);
    preFIX = argv[10];
  }
  // If there are otherwise an arbitrary number of arguements, terminate the program.
  else
  {
    cout << "Please enter the correct number of arguments.\n";
    cout << "The correct number of arguments is 10.\n";
    cout << "The arguments are as follows: dt, t_max, g, r,a_start, a_end, div, dP, init_frac_grazpred, PREFIX.\n";
    exit(1);
  }

  cout << "\n System-level details:\n";
  cout << "dx = " << dx << "\t dt = " << dt << "\t t_max = " << t_max << "\t g = " << g << "\t r = " << r << endl;

  if( 2*v[1]*dt/dx > 1)
  {
    cout << "Warning: CFL condition violated. The velocity of the grazer is too high for the given dx and dt. Please reduce the velocity and dt or increase dx ." << endl;
    exit(2);
  }
  else
  {
    cout << "CFL condition satisfied." << endl;
  }

  cout << "\n=============================================================\n";


  if(preFIX == "N/A" || preFIX == "n/a" || preFIX == "NA" || preFIX == "na" || preFIX == "N" || preFIX == "n" || preFIX == "")
  {
    cout << "No prefix provided. Using default prefix: " << prefix << endl;
    preFIX = prefix;
  }
  //Defining the save folders.
  
  set_Prefix(preFIX); //Set the prefix (and thus save folders) to the user-defined prefix.
  cout << "Save directory for frames: " << frame_folder << "\n";
  cout << "Save directory for preliminary data: " << prelim_folder << "\n";
  cout << "Save directory for final data: " << stat_prefix << "\n";

  set_global_system_params(dt, dx); //Set the global parameters for the simulation.
  cout << "Global parameters set, with dt/2.0 = " << dt2 << " and dx*dx = " << dx2 <<  " and 1/(dx*dx) = " << dx1_2 << "\n";
  //INITIAL CONDITIONS:

  // Equations for MFT E Eqilibrium values  as functions of a (Rainfall).

  
  double mG = E[1]/(kappa); double cG = -b*M[1]/(kappa*kappa); 
  // Coexistance equilibrium density of grazer as a function of a. (Geq = mG*a + cG)
  double mV = 1/b; // Vegetation only equilibrium density as a function of a. (Veq = a/b =mV*a)
  double mW = -3.90344309; double cW = 5.18;
  double mO = 120.4744181; double cO = 0.58243061;
  double  A_Pr = 4.49723728e+05; double b_Pr = -1.73381257;

  string MFT_PreV = std::to_string(mV) + " * a";
  string MFT_PreG = std::to_string(0.0);
  string MFT_V = std::to_string(Vstar);
  string MFT_G = std::to_string(mG) + " * a + " + std::to_string(cG);
  

  MFT_Vec_CoexExpr.assign({MFT_PreV, MFT_PreG, MFT_V, MFT_G});
  // NOTE: The following clow is used for analytic MFT based initial conditions. 
  double scaling_factor[2*Sp] = {5, dP/dP,  5, init_frac_grazpred};

  // LE ORIGINAL (NORMAL CLOSE TO MFT)
  //double clow[2*Sp] = {0, dP/50000.0, dP/500000.0, 4, 20, 10000.0, Gstar, 10, 4, 10};
  // HIGH INIT, > MFT
  //double clow[2*Sp] = {0, dP/dP, dP/(10.0*dP), 4, 20, dP, Gstar*1.5, 15, 4, 10};
  // LOW INIT, < MFT
  //double clow[2*Sp] = {0, dP/dP, dP/(10.0*dP), 4, 20, 10000.0, 2, 0.6, 4, 10};
  // CLOSE TO MFT
  //double clow[2*Sp] = {0, dP/dP, dP/(10.0*dP), 4, 20, 10000.0, Gstar, 10, 4, 10};

  // NORMAL INIT CLOSE TO MFT
  double chigh[2*Sp] = {dP, dP/10000.0,  10000.0 + dP, 10};
  // HIGH INIT > MFT
  //double chigh[2*Sp] = {dP, dP/dP, dP/(10.0*dP), 4, 20, 10000.0 + dP, Gstar*1.5, 15, 4, 10};
  // LOW INIT < MFT
  //double chigh[2*Sp] = {dP, dP/dP, dP/(10.0*dP), 4, 20, 10000.0 + dP, 2, 0.6, 4, 10};


  //NOTE: The following clow is used for Gaussian initial conditions for testing and is NOT MEANT FOR PRODUCTION RUNS.

  //double clow[2*Sp] = {500, Gstar*10, 40, 4, 20, 500, Gstar*10, 40, 4, 10};
  
  //c_high and c_low are arrays of size 2*Sp, with the values of the constants for each species.
	// If R < R_c, then the first Sp elements of c_high and c_low are passed to init_randconstframe() to initialise the frame.
	// If R >= R_c, then the last Sp elements of c_high and c_low are passed to init_randconstframe() to initialise the frame.




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

  stringstream ast, est, dPo, geq; ast << a_start; est  << a_end; dPo << dP;

  
	// This is done to avoid overwriting of files.
  string path_to_folder = frame_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str();
  recursive_dir_create(path_to_folder);

  /**
  stringstream foldername;
	foldername << frame_folder 
  << ast.str() << "-" << est.str() << "_dP_" << dPo.str() << "_Geq_" << geq.str() << "/";
	// Creating a string stream instance to store the values of the parameters in the file name.
	struct stat info2;
	if( stat( foldername.str().c_str(), &info2 ) != 0 )
	{
		cout << "Cannot access " << foldername.str() << ". Creating directory." << endl;
		const int dir_err = system(("mkdir " + foldername.str()).c_str());
		if (-1 == dir_err)
		{
			printf("Error creating directory!n");
			exit(3);
		}
	}
  */

  path_to_folder = prelim_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "_Geq_" + geq.str() +"/TimeSeries";
  recursive_dir_create(path_to_folder);
  
  /**
  stringstream foldername2;
	foldername2 << prelim_folder 
  << ast.str() << "-" << est.str() << "_dP_" << dPo.str() << "_Geq_" << geq.str() << "/";
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
			exit(3);
		}
	}
  */
  
  recursive_dir_create("../Data/Rietkerk/Stochastic/3Sp");
  
  /**
	// Creating a string stream instance to store the values of the parameters in the file name.
	// This is done to avoid overwriting of files.
  stringstream foldername3;
	foldername3 << "../Data/Rietkerk/Stochastic/3Sp/";
	struct stat info0;
	if( stat( foldername3.str().c_str(), &info0 ) != 0 )
	{
		cout << "Cannot access " << foldername3.str() << ". Creating directory." << endl;
		const int dir_err = system(("mkdir " + foldername3.str()).c_str());
		if (-1 == dir_err)
		{
			printf("Error creating directory!n");
			exit(3);
		}
	}
  */
  first_order_critical_exp_delta_stochastic_MultiSp(div, t_max, a_start, a_end, a_c, b, c, D, v, sigma, A, H, E, M, pR, chigh, scaling_factor, dt, dx, dP, r, g);
  
  //Finally recursively delete all the contents of the Temp folder.
  #ifdef _WIN32
  string command = "rmdir /s /q " + prelim_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "_Geq_" + geq.str() +"/Temp";
  #else
  string command = "rm -r " + prelim_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "_Geq_" + geq.str() +"/Temp";
  #endif
  system(command.c_str());

  return 0;
}