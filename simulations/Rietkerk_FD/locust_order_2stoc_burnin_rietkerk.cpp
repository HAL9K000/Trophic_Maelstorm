#include "rietkerk_bjork_basic.h"

int main(int argc, char *argv[])
{
  increase_stack_limit(1024L); //Increase stack limit to 1024 MB.

  #ifndef LOCUST_WIND 
    string error_string = "ERROR: The LOCUST_WIND flag MUST be set globally, as g++ -DLOCUST_WIND locust_order_unity_burnin_rietkerk.cpp ... \n";
    cout << error_string; cerr << error_string;
    exit(1);
  #endif

  #if defined(INIT) && INIT != 2
    string error_string = "ERROR: The initial condition MUST be set to 2 (Burn-in).\n";
    cout << error_string; cerr << error_string;
    exit(1);
  #endif

  // Set wind parameters.
  //wind_speed.first = -10.0; wind_speed.second = 0; // Wind speed range in km/hr.
  double w1, w2; // Wind parameters.

  string preFIX; // Prefix for the output files
  string input_preFIX= ""; // Prefix for the input files
  string input_frame_subdir; // Subdirectory for the input frames.
  double a =0.7; double a_c = 1/24.0; // a_c is the critical value of a for which vegetation appears in classic Rietkerk model.
  double c, gmax, alpha, d, rW, W0, t_max, dt, dx; //int g, div;
  double a_start, a_end; double r; double dP; // Kick for high initial state
  int g, div;

  double mGrazer; mGrazer = 0.001; // Mass of grazer and predator in kg.
  //(FOR DESERT LOCUST CASE, WITH GRAZER = 1g, PREDATOR (COMMON KERNEL) = 0.15 kg [Mullie et al 2021 10.1371/journal.pone.0244733])

  // ASSUMING MASS OF  GRAZER = 20 kg, MASS OF PREDATOR = 100 kg.
  //mGrazer = 20; mPredator = 100;
  // ASSUMING MASS OF  GRAZER = 1 g, MASS OF PREDATOR = 0.15 kg. 
  c = 10000; gmax = 0.05*pow(10, -3.0)/24.0; d = 0.25/24.0; alpha =0.2/24.0; W0 = 0.2; rW = 0.2/24.0; // From Bonachela et al 2015


  //double p0i = 0.5; double p0j= p0i/200; double p0j= 2.25; double p0m= 8; // In g/m^2

  double k0, k1, k2; double d0, d1, d2, d3; double s0, s1; double v1; double dtv1;
  
  // Diffusion coefficient allometries for consumer species:
  //General allometry: ln(D (in km^2/hr)) = 0.3473*ln(M) -4.15517, where M is mass in kg of general consumer.
  //Grazer allometry: ln(D (in km^2/hr)) = 0.5942*ln(M) -6.3578, where M is mass in kg of grazer.
  
  dx= 2 ; //From Bonachela et al 2015 (in km)
  //d0 = 0.00025/24.0; d1=0.0298; d2 = 0.00025/24.0; d3= 0.025/24.0; //From Bonachela et al 2015 (in km^2/hr) and general D allometry.
  d0 = 0.00025/24.0; d1 = 0.0009579; d2 = 0.00025/24.0; d3= 0.025/24.0; //From Bonachela et al 2015 (in km^2/hr) and specific D allometries.
  k0= 0; k1 = 5; k2 =5000;
  s0 = sqrt(d0/(dx*dx)); s1 = sqrt(d1/(dx*dx)); // ~ D0/(dx)^2 (in kg^{0.5}/(km hr))
  
  // Allometric scaling for velocity: 
  // ALTERNATIVELY FOR LOW BODY MASSES: Use Hirt allometric parameters fit to Ricardo data= v (m/hr) = 250*M^0.256*(1 - e^{-20*(M)^(-0.6)}), where M is mass in kg.
  v1 = 0.25*pow(mGrazer, 0.26)*(1.0 - exp(-20.0*pow(mGrazer, -0.6))); //In km/hr [Equivalent to ~ 0.041 km/hr for 1g]

  //Time-scale of advection.
  dtv1 = 0.016958; //In hr, based on relation ln(dtv_grazer) = 0.2596*ln(M) + 1.8106, where M is mass in kg, and dtv_grazer is in min.

  double D[Sp] ={d0, d1, d2, d3}; //Diffusion coefficients for species.
  double K[3] ={k0, k1, k2}; //Diffusion coefficients for species.
  double sigma[Sp] ={s0, s1, 0, 0}; //Demographic stochasticity coefficients for species.
  double v[Sp] ={0, v1, 0, 0}; //Velocity of species.

  /**

  double beta = c*gmax - d;
  double p0istar, p0jstar, p0mstar; // Analytic steady state values.

  p0jstar = d*K[1]/beta; 
  p0istar = (c/d)*(R - rW*p0jstar);
  p0mstar = (R/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
  */
  double aij, hij, ej, m, m_scale, aij_scale;
  m_scale = 1.0; aij_scale = 1.0; 
  // Attack rate, Handling time, Efficiency of consumption, Mortality rate, scaling factor for mortality rate.
  aij = 3.6*pow(10.0, -6.08)*pow(mGrazer, -0.37); // in km^2/(hr kg)
  hij = 1; // Handling time in hrs
  ej =0.45; m = 0.061609*pow(mGrazer, -0.25)/8760.0; // Mortality rate in hr^{-1}
  double H[SpB][SpB] ={{0, hij}, 
                     {hij, 0.0}};    // Handling time matrix [T]. No canabilism, symmetric effects.
  double A[SpB][SpB] ={{0, aij}, 
                      {aij, 0.0}};  // Attack rate matrix []. No canabilism, symmetric effects. 
  double M[SpB] = {d, m};     // Mortality Rate of Species. (in hr^{-1})
  double E[SpB] ={1.0, ej}; //Efficiency of consumption.
  double pR[Sp] ={0.0, 1.04285/2.0}; //Perception rate of species (in km)

  double Vstar = m/((ej -m*hij)*aij); //Steady state for grazer.
  double Wstar_veg = d*K[1]/(c*gmax - d); //Steady state for vegetation.

  double init_frac_graz = 1.0; //Initial fraction of grazers in the high state.

  cout << "Presets are as follows:\n";
  cout << "For the vegetation:\n";
  cout << "c = " << c << "\t gmax = " << gmax << "\t d = " << d << "\t alpha = " << alpha << "\t W0 = " << W0 << "\t rW = " << rW << endl;
  cout << "\nBase MFT biomass density of Vegetation = " << Vstar << " kg/km^2\n" << endl;
  cout << "For the grazer:\n";
  cout << "aij = " << aij << "\t hij = " << hij << "\t ej = " << ej << "\t m = " << m << endl;
  cout << "\n System-level details:\n";
  cout << "pi = " << PI  << endl;

  cout << "Diffusion Constants For Species:\n";
  cout << "D0 = " << setprecision(16) << D[0] << "\t D1 = " << D[1] << "\t D2 = " << setprecision(16) << D[2] << "\t D3 = " << setprecision(16) << D[3] << endl;

  cout << "This is a 2Sp (Two) Stochastic Rietkerk Model Script FOR LOCUSTS WITH BURN-IN FRAMES AND UNIT INITIALISATION OF GRAZERS\n";
  cout << "Header Species Check: " << std::to_string(SpB) << "\n";


  if(argc == 1)
  {
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

    cout << "Enter wind speed in x-direction (in km/hr): ";
    cin >> w1;

    cout << "Enter wind speed in y-direction (in km/hr): ";
    cin >> w2;

    cout << "Enter fractional deviation from unity for initial conditions of Grazer above R_c (: ";
    cin >> init_frac_graz;

    //cout << "Enter scaling factor for grazer nat. mortality (mj) which has a base value " << m << " \t:";
    //cin >> m_scale;

    cout << "Enter scaling factor for grazer attacking rate (aij) which has a base value " << aij << " \t:";
    cin >> aij_scale;

    cout << "Enter Prefix (common choices include 'DiC-REF', 'DDM-DiC-BURNIN', 'nDiC' etc.) If not provided, default value: " 
        + prefix + " will be used: ";
    cin >> preFIX;

    cout << "NOTE: Valid keys for BURN-IN frames include : [";
    // Iterate over the keys in the map input_keys.
    for (auto const& x : input_keys)
      cout << x.first << '\t'; // string (key) 
    cout << "]\n";
    
    cout << "Enter input subdirectory for BURN-IN frames (for example '/HEXBLADE/L_{L}_a_0/SEP_10_WID_4/MSF_0.25/MAMP_20000_MINVAL_0.0'): ";
    cin >> input_frame_subdir;

  }
  // If the user has entered the correct number of arguments, then the arguments are read from the command line.
  else if(argc == 15)
  {
    dt = atof(argv[1]);
    t_max = atof(argv[2]);
    g = atoi(argv[3]);
    r = atoi(argv[4]);
    a_start = atof(argv[5]);
    a_end = atof(argv[6]);
    div = atoi(argv[7]);
    dP = atof(argv[8]);
    w1 = atof(argv[9]); w2 = atof(argv[10]);
    init_frac_graz = atof(argv[11]);
    aij_scale = atof(argv[12]);
    preFIX = argv[13];
    input_frame_subdir = argv[14];
  }
  // If there are otherwise an arbitrary number of arguements, terminate the program.
  else
  {
    cout << "Please enter the correct number of arguments.\n";
    cout << "The correct number of arguments is 14.\n";
    cout << "The arguments are as follows: dt, t_max, g, r, a_start, a_end, div, dP, w1, w2, init_GP,  PREFIX, BRIN_DIR.\n";
    exit(1);
  }

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

  cout << "\n=============================================================\n";


  if(preFIX == "N/A" || preFIX == "n/a" || preFIX == "NA" || preFIX == "na" || preFIX == "N" || preFIX == "n" || preFIX == "")
  {
    cout << "No prefix provided. Using default prefix: " << prefix << endl;
    preFIX = prefix;
  }
  //Defining the save folders.
  
  set_LocustPrefix(preFIX); //Set the prefix (and thus save folders) to the user-defined prefix.
  cout << "Save directory for frames: " << frame_folder << "\n";
  cout << "Save directory for preliminary data: " << prelim_folder << "\n";
  cout << "Save directory for final data: " << stat_prefix << "\n";

  // Define input file folders.
  set_input_Prefix(input_frame_subdir, preFIX, a_c, dP, Vstar); 
  //Set the prefix (and thus save folders) to the user-defined prefix.
  cout << "Heuristic Input directory for frames: " << input_frame_parenfolder  + input_frame_subfolder << "\n";

  set_global_system_params(dt, dx); //Set the global parameters for the simulation.
  wind_speed.first = w1; wind_speed.second = w2; // Wind speed range in km/hr.
  cout << "Wind speed (x,y) set to: [ " << wind_speed.first << ", " << wind_speed.second << " ] km/hr.\n";
  cout << "Global parameters set, with dt/2.0 = " << dt2 << " and dx*dx = " << dx2 <<  " and 1/(dx*dx) = " << dx1_2 << "\n";

  // Setting m to the scaled value, and assosciated parameters (M[] and MFT expressions) to the scaled value.
  //m = m_scale*m; M[1] = m; // Mortality Rate of Species. (in hr^{-1)
  aij = aij_scale*aij; A[0][1] = aij; A[1][0] = aij; // Attack rate matrix []. No canabilism, symmetric effects.
  Vstar = m/((ej -m*hij)*aij); //Steady state for grazer.
  cout << "Global parameters set, with dt/2.0 = " << dt2 << " and dx*dx = " << dx2 <<  " and 1/(dx*dx) = " << dx1_2 << "\n";
  //cout << "Scaled grazer nat. mortality (mj) = " << M[1] << "\n";
  cout << "Scaled grazer attacking rate (aij) = " << A[0][1] << " " << A[1][0] << "\n";
  cout << "Base MFT biomass density of Vegetation (V*) = " << Vstar << " kg/km^2\n" << endl;
  
  // Next, create the integer array dtV[SpB], which is initialised to {0, round(dtv1/dt), round(dtv2/dt)}.
  int dtV[SpB] = {0, max(1, int(round(dtv1/dt)))};
  // Note that the dtV vector is used to store the time-steps for the grazer and predator, respectively.
  // and that the time-steps are rounded to the nearest integer value, with a minimum value of 1 (to avoid division by zero).

  cout << "Values of dtV for the grazer is: " << dtV[1] << "\n";
  ///** CORRECT MFT INITIALISATIONS
  //double mG = (6.76550102e+08)/m_scale; double cG =  (2.81895876e+07)/m_scale; 
  //double mW_Prev = w; double cW_Prev = 0.0;
  double mW = 1/rW; double cW = 0.0;
  double mW_Prev = 1/rW; double cW_Prev  = 0.0;
  double A_G = 7332.541; double b_G = -15.5375; // DEFAULTS when aij_scale = m_scale = e_scale = 1

  if (m_scale == 1.0)
  {
    if(aij_scale == 1.0)
    { A_G = 7332.541; b_G = -15.5375; }
    else if(aij_scale == 10)
    { A_G = 645.346; b_G = -18.94669; }
    else if(aij_scale == 25)
    { A_G = 240; b_G = -20.7705;      }
    else
    { std::cerr << "Scaling factor for aij not supported. Please use 1, 10, or 25." << endl; exit(1); }
  }
  else
  { std::cerr << "Scaling factor for mortality rate not supported. Please use 1." << endl; exit(1); }
  double mO = (Vstar + K[2] )/((Vstar + K[2]*W0)*alpha); double cO = 0.0; // O* = R*(V* + K2)/[(V* + K2*W0)*alpha]
  double mO_Prev = 1/(alpha*W0); double cO_Prev = 0.0; 
  string MFT_V = std::to_string(Vstar); 
  string MFT_V_Prev = std::to_string(0);
  string MFT_G = std::to_string(0);
  string MFT_G_Prev = std::to_string(0);
  string MFT_W = std::to_string(mW) + " * a - " + std::to_string(cW);
  string MFT_W_Prev = std::to_string(Wstar_veg);
  string MFT_O = std::to_string(mO) + " * a";
  string MFT_O_Prev = std::to_string(mO_Prev) + " * a";
  
  /** MFT INITIALISATIONS (OLD)
  double mG = (6.76550102e+08)/m_scale; double cG =  (2.81895876e+07)/m_scale;
  double mW = 60; double cW = 2.50;
  double mW_Prev = 1/rW; double cW_Prev = 0.0;
  double mO = (Vstar + K[2] )/((Vstar + K[2]*W0)*alpha); double cO = 0.0; // O* = R*(V* + K2)/[(V* + K2*W0)*alpha]
  double mO_Prev = 1/(alpha*W0); double cO_Prev = 0.0; 
  //double mO = 596.37; double cO = 0.0;
  string MFT_V = std::to_string(Vstar); string MFT_V_Prev = std::to_string(0);
  string MFT_G = std::to_string(mG) + " * a + " + std::to_string(cG);
  string MFT_G_Prev = std::to_string(0);
  string MFT_W = std::to_string(mW) + " * a - " + std::to_string(cW);
  string MFT_W_Prev = std::to_string(mW_Prev) + " * a";
  string MFT_O = std::to_string(mO) + " * a";
  string MFT_O_Prev = std::to_string(mO_Prev) + " * a";
  //*/

  // NOTE: The following clow is used for analytic MFT based initial conditions.
  #if defined(INIT) && INIT == 0
    double scaling_factor[2*Sp] = {Vstar, dP/dP, 1, 1,  1, init_frac_graz, 1, 1};
    // USED FOR HOMOGENEOUS MFT BASED INITIAL CONDITIONS.
  #else
    double scaling_factor[2*Sp] = {5, dP/dP*init_frac_graz, 1, 1,  5, init_frac_graz, 1, 1};
    // USED FOR PERIODIC MFT BASED INITIAL CONDITIONS.
  #endif

  MFT_Vec_CoexExpr.assign({MFT_V_Prev, MFT_G_Prev, MFT_W_Prev, MFT_O_Prev, MFT_V, MFT_G, MFT_W, MFT_O});
  double chigh[2*Sp]; 

  //NOTE: The following scaling_factor is used for Gaussian initial conditions for testing and is NOT MEANT FOR PRODUCTION RUNS.
  //double scaling_factor[2*Sp] = {500, 100, 4, 20, 500, 100, 4, 10};

  stringstream ast, est, veq, dPo; ast << a_start; est  << a_end; dPo << dP; veq << setprecision(5) << Vstar;

  // This is done to avoid overwriting of files.
  string path_to_folder = frame_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "_Veq_" + veq.str();
  recursive_dir_create(path_to_folder);
  
  path_to_folder = prelim_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "_Veq_" + veq.str() +"/TimeSeries";
  recursive_dir_create(path_to_folder);

  recursive_dir_create("../Data/LocRiet/Stochastic/"+ std::to_string(SpB) +"Sp");

  //first_order_critical_exp_delta_stochastic(div, t_max, a_start, a_end, c, gmax, alpha, d, rW, W0, D, K, sigma, dt, dx, dP, r, g);
  //first_order_critical_exp_delta_stochastic_2Sp(div, t_max, a_start, a_end, c, gmax, alpha, rW, W0, D, v, K, sigma, A, H, E, M, pR,dt, dx, dP, r, g);
  //first_order_critical_exp_delta_stochastic_2Sp(div, t_max, a_start, a_end, a_c, c, gmax, alpha, rW, W0, D, v, K, sigma, A, H, E, M, pR, chigh, scaling_factor, dt, dx, dP, r, g, -1, Vstar);
  first_order_critical_exp_delta_stochastic_MultiSp(div, t_max, a_start, a_end, a_c, c, gmax, alpha, rW, W0, D, v, K, sigma, A, H, E, M, pR, dtV, scaling_factor, dt, dx, dP, r, g, -1, Vstar);
  //tupac_percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, auto &Rh0,
	 //double t_max, double a[], double b[], double c[], double D[], double sigma[], double dt, double dx, int r,  int g)

  return 0;
}