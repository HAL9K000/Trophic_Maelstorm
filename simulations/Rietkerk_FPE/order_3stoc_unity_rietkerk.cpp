#include "rietkerk_bjork_basic.h"

int main(int argc, char *argv[])
{
  increase_stack_limit(1024L); //Increase stack limit to 1024 MB.

  string preFIX; // Prefix for the output files
  double a_c =1/24.0;
  double c, gmax, alpha, d, rW, W0, t_max, dt, dx; 
  double a_start, a_end; double r; double dP; // Kick for high initial state
  int g, div;
  // Using units of mass = kg, length = km, time = hr.

  double dt_analytical; //Analytical time-step for the system.

  dx= 0.1 ; //From Bonachela et al 2015 (in km)
  
  double mGrazer; double mPredator; // Mass of grazer and predator in kg.


  // ASSUMING MASS OF  GRAZER = 20 kg, MASS OF PREDATOR = 100 kg.
  mGrazer = 20; mPredator = 100;
  // ASSUMING MASS OF  GRAZER = 1 g, MASS OF PREDATOR = 0.15 kg. 
  //(FOR DESERT LOCUST CASE, WITH GRAZER = 1g, PREDATOR (COMMON KERNEL) = 0.15 kg [Mullie et al 2021 10.1371/journal.pone.0244733])
  //mGrazer = 0.001; mPredator = 0.15;

  // ASSUMING MASS OF  GRAZER = 150 g, MASS OF PREDATOR = 2 kg.
  //mGrazer = 0.15; mPredator = 2;

  c = 10000; gmax = 0.05*pow(10, -3.0)/24.0; d = 0.25/24.0; alpha =0.2/24.0; W0 = 0.2; rW = 0.2/24.0; // From Bonachela et al 2015

  //double p0i = 0.5; double p0j= p0i/200; double p0j= 2.25; double p0m= 8; // In g/m^2

  double k0, k1, k2; double d0, d1, d2, d3, d4; double s0, s1, s2; double v1, v2;  double dtv1, dtv2;

  k0= 0; k1 = 5; k2 =5000;
  
  // Diffusion coefficient allometries for consumer species:
  //General allometry: ln(D (in km^2/hr)) = 0.3473*ln(M) -4.15517, where M is mass in kg of general consumer.
  //Grazer allometry: ln(D (in km^2/hr)) = 0.5942*ln(M) -6.3578, where M is mass in kg of grazer.
  //Predator allometry: ln(D (in km^2/hr)) = 0.6072*ln(M) -4.3289, where M is mass in kg of predator.
  
  

  /** // Allometric scaling for diffusion for : Generic consumers - D [km^2/hr] = e^{-4.5517}*M^{0.3473}, where M is mass in kg.
  d0 = 0.00025/24.0; d1 = exp(-4.5517)*pow(mGrazer, 0.3473); d2 = exp(-4.5517)*pow(mPredator, 0.3473); d3 = 0.00025/24.0; d4= 0.025/24.0; s2 = 3; //From Bonachela et al 2015 (in km^2/hr) and general D allometry.
  //From Bonachela et al 2015 (in km^2/hr) and general D allometry.
  //d0 = 0.00025/24.0; d1 = 0.0009579; d2 = 0.00545856; d3 = 0.00025/24.0; d4= 0.025/24.0; s2 = 10; //From Bonachela et al 2015 (in km^2/hr) and specific D allometries.
  s0 = sqrt(d0/(dx*dx)); s1 = sqrt(d1/(dx*dx)); s2 = sqrt(d2/(dx*dx)); // ~ D0/(dx)^2 (in kg^{0.5}/(km hr))
  // Allometric scaling for velocity: 
  // v (m/hr) = 43.706*M^1.772*(1 - e^{-14.27*(M)^(-1.865)}), where M is mass in kg.
  //v1 = 43.706*pow(mGrazer, 1.772)*(1.0 - exp(-14.27*pow(mGrazer, -1.865)))/1000.0;
  //v2 = 43.706*pow(mPredator, 1.772)*(1.0 - exp(-14.27*pow(mPredator, -1.865)))/1000.0; 
  //In km/hr [Equivalent to ~ 0.000211125*0.001 km/hr for 1g, 1.515567*0.001 km/hr for 150g]
  // ALTERNATIVELY FOR LOW BODY MASSES: Use Hirt allometric parameters fit to Ricardo data= v (m/hr) = 250*M^0.256*(1 - e^{-20*(M)^(-0.6)}), where M is mass in kg.
  v1 = 0.25*pow(mGrazer, 0.26)*(1.0 - exp(-20.0*pow(mGrazer, -0.6))); //In km/hr [Equivalent to ~ 0.152 km/hr for 150g grazer ]
  v2 = 0.25*pow(mPredator, 0.26)*(1.0 - exp(-20.0*pow(mPredator, -0.6))); //In km/hr [Equivalent to ~ 0.299369 km/hr for 2kg predator]

  //Time-scale of advection.
  dtv1 = (exp(1.8106)*pow(mGrazer, 0.2596))/60.0; // 0.062272 In hr, based on relation ln(dtv_grazer) = 0.2596*ln(M) + 1.8106, where M is mass in kg, and dtv_grazer is in min.
  dtv2 = (exp(1.8106)*pow(mPredator, 0.2596))/60.0; // 0.12199 In hr, based on relation ln(dtv_pred) = 0.2596*ln(M) + 1.8106, where M is mass in kg, and dtv_pred is in min.
  //*/

  // SPECIFIC MASS OF  GRAZER = 20 kg, MASS OF PREDATOR = 100 kg.  
  d0 = 0.00025/24.0; d1=0.0298; d2= 0.05221; d3 = 0.00025/24.0; d4= 0.025/24.0; s2 = 3; //From Bonachela et al 2015 (in km^2/hr) and general D allometry.
  //d0 = 0.00025/24.0; d1=0.01028; d2= 0.215965; d3 = 0.00025/24.0; d4= 0.025/24.0; s2 = 10; //From Bonachela et al 2015 (in km^2/hr) and specific D allometries.
  
  s0 = sqrt(d0/(dx*dx)); // ~ D0/(dx)^2 (in kg^{0.5}/(km hr))
  s1 = 1; // ~ D/(dx)^2 (in kg^{0.5}/(km hr))
  //s2 = 3; // ~ D/(dx)^2 (in kg^{0.5}/(km hr))
  //s2 = 10; // ~ D/(dx)^2 (in kg^{0.5}/(km hr))
  
  // Allometric scaling for velocity: 
  // v (m/hr) = 43.706*M^1.772*(1 - e^{-14.27*(M)^(-1.865)}), where M is mass in kg.
  v1 = 0.45427; //In km/hr
  v2 = 0.40587; //In km/hr

  //Time-scale of advection.
  dtv1 = 0.11817455; //In hr, based on relation ln(dtv_grazer) = 0.2596*ln(M) + 1.8106, where M is mass in kg, and dtv_grazer is in min.
  dtv2 = 0.765270868; //In hr, based on relation ln(dtv_predator) = 0.7*ln(M) + 0.6032, where M is mass in kg, and dtv_predator is in min.
  //*/

  // Round the time-scales to the two decimal places.
  dtv1 = round(dtv1*100.0)/100.0; dtv2 = round(dtv2*100.0)/100.0;

  // Set dt_analytical = dt
  dt = dt_analytical = min(dtv1, dtv2); //Set dt_analytical to the minimum of the two time-scales.

  double D[Sp] ={d0, d1, d2, d3, d4}; //Diffusion coefficients for species.
  double K[3] ={k0, k1, k2}; //Diffusion coefficients for species.
  double sigma[Sp] ={s0, s1, s2, 0, 0}; //Demographic stochasticity coefficients for species.
  double v[SpB] ={0, v1, v2}; //Velocity of species.
  double dt_adv[SpB] = {dt_analytical, dtv1, dtv2}; //Time-scales of advection for species.


  #if defined(MV_INVARIANCE) && MV_INVARIANCE == 0
    cout << "Displacement invariant advection.\n";
    for(int s = 1; s < SpB; s++)
      v[s] *= dt_adv[s]/dt_analytical; // Scale the velocity by the ratio of the time-scales to ensure displacement invariance.
  #elif defined(MV_INVARIANCE) && MV_INVARIANCE == 1
    cout << "Velocity invariant advection.\n";
  #endif

  double aij, hij, ej, mj; double ajm, hjm, em, mm; // Parameters for the grazer and predator.
  double mj_scale, aij_scale, mm_scale, ajm_scale;
  mj_scale = 1.0; aij_scale = 1.0; mm_scale = 1.0; ajm_scale = 1.0;
  
  aij = 3.6*pow(10.0, -6.08)*pow(mGrazer, -0.37); // in km^2/(hr kg)
  hij = 1; // Handling time in hrs
  ej =0.45; mj = 0.061609*pow(mGrazer, -0.25)/8760.0; // Mortality rate in hr^{-1}
  //Predator of mass 100 kg.
  ajm = 3.6*pow(10.0, -6.08)*pow(mPredator, -0.37); // in km^2/(hr kg)
  hjm = 1; // Handling time in hrs
  em =0.85; mm = 0.061609*pow(mPredator, -0.25)/8760.0; // Mortality rate in hr^{-1}

  double Km = pow(10.0, 1.22)*pow(mPredator, -0.31); // Carrying capacity of predator in kg/km^2s

  double Gstar = mm/((em -mm*hjm)*ajm); //Steady state for grazer.
  double H[SpB][SpB] ={{0, hij, 0}, 
                     {hij, 0.0, hjm}, 
                     {0, hjm, 0}};    // Handling time matrix [T]. No canabilism, symmetric effects.
  double A[SpB][SpB] ={{0, aij, 0}, 
                      {aij, 0.0, ajm}, 
                      {0, ajm, 0}};  // Attack rate matrix []. No canabilism, symmetric effects. 
  double M[SpB] = {d, mj, mm};     // Mortality Rate of Species. (in hr^{-1})
  double E[SpB] ={1.0, ej, em}; //Efficiency of consumption.
  //double D[Sp] ={d0, d1}; //Diffusion coefficients for species. (in km^2/hr)
  double pR[Sp] ={0.0, 1.04285/2.0, 1.25536/2.0}; //Perception radius of species (in km)

  double init_frac_pred = 1.0; //Initial fraction of grazers and predators in the high state.
  
  cout << "Presets are as follows:\n";
  cout << "For the vegetation:\n";
  cout << "c = " << c << "\t gmax = " << gmax << "\t d = " << d << "\t alpha = " << alpha << "\t W0 = " << W0 << "\t rW = " << rW << endl;
  cout << "For the grazer of mass " << mGrazer << " kg:\n";
  cout << "aij = " << aij << "\t hij = " << hij << "\t ej = " << ej << "\t mj = " << mj << endl;
  cout << "\nMFT biomass density of Grazer = " << Gstar << " kg/km^2\n" << endl;
  cout << "For the predator of mass " << mPredator << " kg:\n";
  cout << "ajm = " << ajm << " " << A[1][2] << "\t hjm = " << hjm << " " << H[1][2] 
       << "\t em = " << em << " " << E[2] <<  "\t mm = " << mm <<  " " << M[2] << endl;
  cout << "\n System-level details:\n";
  cout << "pi = " << PI  << endl;

  cout << "Diffusion Constants For Species:\n";
  cout << "D0 = " << setprecision(16) << D[0] << "\t D1 = " << D[1] << "\t D2 = " << setprecision(16) << D[2] << "\t D3 = " << setprecision(16) << D[3] << endl;

  cout << "This is a 3Sp (Three) Stochastic FPE MOVEMENT Rietkerk Model Script AND UNIT INITIALISATION OF CONSUMERS\n";

  cout << "\n============================================================\n";

  // There are eight user-defined parameters, as follow below: dt, t_max, g, r, a_start, a_end, div, dP.

  //Check if the user has entered the correct number of arguments. If not prompt the user to enter all the arguments within the script.
  // First check if no arguments are entered. In this case, the user will be prompted to enter all the arguments using cin.

  if(argc == 1)
  {
    cout << "Please enter all the arguments within the script.\n";
    //cout << "Enter desired time-step (~0.1 is a good choice): ";
    cout << "dt set to " << dt_analytical << " hrs (Analytical time-step for the system)\n";
    //cin >> dt;

    cout << "Enter maximum duration of simulation (in hrs): ";
    cin >> t_max;

    cout << "Enter the desired grid-size (L) (dx = 100 m): ";
    cin >> g;

    cout << "Enter the desired number of replicates (r): ";
    cin >> r;

    cout << "Note that a relevant value of 'a' may be construed as:\t" << setprecision(10) << a_c <<endl;

    cout << "Enter starting a-value: ";
    cin >> a_start;

    cout << "Enter ending a-value: ";
    cin >> a_end;

    cout << "Enter number of divisions of p (n+1): ";
    cin >> div;

    cout << "Enter kick for high state: ";
    cin >> dP;

    cout << "Enter fractional deviation from 0.1 for initial conditions of Predator above R_c (0.5-1.2 is a good choice): ";
    cin >> init_frac_pred;

    cout << "Enter scaling factor for grazer attacking rate (aij) which has a base value " << aij << " \t:";
    cin >> aij_scale;

    cout << "Enter scaling factor for predator attacking rate (ajm) which has a base value " << ajm << " \t:";
    cin >> ajm_scale;

    cout << "Enter Prefix (common choices include 'DiC-REF', 'DDM-DiC-BURNIN', 'nDiC' etc.) If not provided, default value: " 
        + prefix + " will be used: ";
    cin >> preFIX;

  }
  // If the user has entered the correct number of arguments, then the arguments are read from the command line.
  else if(argc == 13)
  {
    dt = dt_analytical;
    //dt = atof(argv[1]);
    t_max = atof(argv[2]);
    g = atoi(argv[3]);
    r = atoi(argv[4]);
    a_start = atof(argv[5]);
    a_end = atof(argv[6]);
    div = atoi(argv[7]);
    dP = atof(argv[8]);
    init_frac_pred = atof(argv[9]);
    aij_scale = atof(argv[10]);
    ajm_scale = atof(argv[11]);
    preFIX = argv[12];
  }
  // If there are otherwise an arbitrary number of arguements, terminate the program.
  else
  {
    cout << "Please enter the correct number of arguments.\n";
    cout << "The correct number of arguments is 13.\n";
    cout << "The arguments are as follows: dt, t_max, g, r, a_start, a_end, div, dP, init_GP, PREFIX.\n";
    exit(1);
  }

  cout << "\n System-level details:\n";
  cout << "dx = " << dx << "\t dt = " << dt << "\t t_max = " << t_max << "\t g = " << g << "\t r = " << r << endl;

  cout << "\n=============================================================\n";


  if(preFIX == "N/A" || preFIX == "n/a" || preFIX == "NA" || preFIX == "na" || preFIX == "N" || preFIX == "n" || preFIX == "")
  {
    cout << "No prefix provided. Using default prefix: " << prefix << endl;
    preFIX = prefix;
  }
  //Defining the save folders.
  
  //IMP: Don't pass mGrazer and mPredator as arguments to set_Prefix if you don't want their masses to be included in the folder names.
  set_Prefix(preFIX, mGrazer, mPredator); //Set the prefix (and thus save folders) to the user-defined prefix.
  cout << "Save directory for frames: " << frame_folder << "\n";
  cout << "Save directory for preliminary data: " << prelim_folder << "\n";
  cout << "Save directory for final data: " << stat_prefix << "\n";
  cout << "=============================================================\n";

  set_global_FKE_params(D, v, dt_adv, dt_analytical, dx, g, 4); //Set the global parameters for system-wide and FKE movement parameters.
  cout << "Global parameters set, with dt/2.0 = " << dt2 << " and dx*dx = " << dx2 <<  " and 1/(dx*dx) = " << dx1_2 << "\n";

  cout << "FKE Movement Parameters set as follows:\n";
  for(int s=0; s < SpB; s++)
  {
    cout  << "For species " << s << " with [D, v, dtv]: " << D[s] << ", " << v[s] << ", " << dt_adv[s] << ", we have:\n";
    cout << "Beta_V = " << beta_vel[s] << "\t sigma_D = " << sigma_D[s] << " with scaled bound = " << sig_D_ScaledBounds[s] 
    << "\t sigma_vD( X, Y) = { " << sigma_vD[s].first << ", " << sigma_vD[s].second << " } with scaled bounds = { " 
    << sig_vD_ScaledBounds[s].first << ", " << sig_vD_ScaledBounds[s].second << " } and AVG Shift mu : " << mu_vel_prefactor[s]*v[s] << "\n";
    if ( sigma_vD[s].first == sigma_vD[s].second)
      cout << "Additionally velocity OU Process for species " << s << " IS SYMMETRIC.\n";
    else
      cout << "Additionally velocity OU Process for species " << s << " IS NOT SYMMETRIC.\n";
    cout << "____________________________________________________\n";
  }

  #ifdef BARRACUDA
    cout << "FKE Movement Parameters on CUDA device memory set as follows:\n";
    reportConstantMemory();
  #endif

  
  
  aij = aij_scale*aij; A[0][1] = aij; A[1][0] = aij; // Attack rate matrix []. No canabilism, symmetric effects.
  ajm = ajm_scale*ajm; A[1][2] = ajm; A[2][1] = ajm; // Attack rate matrix []. No canabilism, symmetric effects.
  Gstar = mm/((em -mm*hjm)*ajm); //Steady state for grazer.
  double Vstar_veg = mj/((ej -mj*hij)*aij); //Steady state for vegetation (Veg) IFF there is NO predator.
  //Gstar = mm/((em -mm*hjm)*ajm); //Steady state for grazer.

  // Next, create the integer array dtV[SpB], which is initialised to {0, round(dtv1/dt), round(dtv2/dt)}.
  int dtV[SpB] = {0, max(1, int(round(dtv1/dt))), max(1, int(round(dtv2/dt)))};
  //int dtV[SpB] = {0, 1,1}; // For testing purposes.
  // Note that the dtV vector is used to store the time-steps for the grazer and predator, respectively.
  // and that the time-steps are rounded to the nearest integer value, with a minimum value of 1 (to avoid division by zero).

  cout << "Values of dtV for the grazer and predator are: " << dtV[1] << " and " << dtV[2] << " respectively.\n";


  

  ///** CORRECT MFT INITIALISATIONS
  double mV = 960000.0; double cV = -40000.0; // GOOD APPROXIMATION FOR GENERAL RANGE OF A AND M VALUES CONSIDERED HERE.
  double mV_Prev = 0.0; double cV_Prev = 0.0;
  double mW = -3.90344309; double cW = 5.18; 
  double A_W = 0; double b_W = 0; 
  // NOTE: If these values are ! = 0, then the MFT for W is not a simple linear function, but a decaying exponential.
  // of the form: mW = A_W + b_W*exp(mW*(a - cW))) [mW < 0]
  double mW_Prev = 1/rW; double cW_Prev  = 0.0; 
  double mO = 120.0; double cO = 0.0; // GOOD APPROXIMATION FOR GENERAL RANGE OF A AND M VALUES CONSIDERED HERE.
  double mO_Prev = 1/(alpha*W0); double cO_Prev = 0.0; 
  double  A_Pr = 4.49723728e+05; double b_Pr = -1.73381257; 
  // DEFAULTS when aij_scale = mj_scale = ej_scale = mm_scale = ajm_scale = = em_scale = 1
  double A_G = 7332.541; double b_G = -15.5375; // DEFAULTS when aij_scale = mj_scale = ej_scale = 1

  double mG_Prev = 0.0; double cG_Prev = 0.0;
  double mPr_Prev = 0.0; double cPr_Prev = 0.0;
  if( mj_scale == 1.0)
  {
    if(aij_scale == 1.0)
    { A_G = 7332.541; b_G = -15.5375; }
    else if(aij_scale == 10)
    { A_G = 645.346; b_G = -18.94669; }
    else if(aij_scale == 25)
    { A_G = 240; b_G = -20.7705;      }
    else
    { std::cerr << "Scaling factor for aij not supported. Please use 1, 10, or 25." << endl; exit(1); }

    if(mm_scale == 1.0)
    {
      if(aij_scale == 1.0)
      { 
        if (ajm_scale == 1.0)
        { A_Pr = 4.49723728e+05; b_Pr = -1.73381257; mW = -3.90344309; cW = 5.18;  }
        else if(ajm_scale == 10.0)
        { A_Pr = 4.55158720e+04; b_Pr = -1.70802424; mW = -0.374252;  cW = 5.01475; }
        else if(ajm_scale == 25.0)
        { A_Pr = 2.93592196e+05; b_Pr = -0.0911713457; mW = -0.149703;  cW = 5.0059; }
        else
        { std::cerr << "Scaling factor for ajm not supported. Please use 1, 10, or 25." << endl; exit(1); }
      }
      else if(aij_scale == 25.0)
      {
        if(ajm_scale == 1.0)
        { A_Pr = 6.8933839e+05; b_Pr = -20.0807587; A_W = 1; b_W = 8.59213; mW = -32.0623;  cW = 0.0243051;}
        else if(ajm_scale == 10.0)
        { A_Pr = 6.8798713e+04; b_Pr = -20.4693606; A_W = 4.65526; b_W = 14.9438; mW = -55.1092;  cW = -0.0286443; }
        else if(ajm_scale == 25.0)
        { A_Pr = 2.75159916e+04; b_Pr = -20.4953090; A_W = 4.86169; b_W = 12.9258; mW = -53.9559;  cW = -0.044096; }
        else if(ajm_scale == 50.0)
        { A_Pr = 1.37574163e+04; b_Pr = -20.5039489; A_W = 4.92608802; b_W = 3.68792022; mW = -33.60729439;  cW = -0.07048515; }
        else if(ajm_scale == 100.0)
        { A_Pr = 6.87856303e+03; b_Pr = -20.5082747; A_W = 4.96537423; b_W = 3.19028452; mW = -53.4178045;  cW = -0.0445609397; }
        else if(ajm_scale == 150.0)
        { A_Pr = 4.58567622e+03; b_Pr = -20.5097222; A_W = 4.97691262; b_W = 2.68932522; mW = -53.9559;  cW = -0.044096; }
        else
        { std::cerr << "Scaling factor for ajm not supported. Please use 1, 10, 25, 50, 100 or 150." << endl; exit(1); }
      }
      else if (aij_scale == 10.0)
      {
        if(ajm_scale == 1.0)
        { A_Pr = 6.11008374e+05; b_Pr = -10.8608332; mW = -16.3029;  cW = 5.64117; }
        else if(ajm_scale == 10.0)
        { A_Pr = 6.0987417e+04; b_Pr = -10.9556230; A_W = 4.6581; b_W = 13.2612; mW = -11.889; cW = 0.263293; }
        else if(ajm_scale == 25.0)
        { A_Pr = 2.4392e+04; b_Pr = -10.9619; A_W = 4.8; b_W = 2.68932522; mW = -53.3597073;  cW = -0.0490352778; }
        else
        { std::cerr << "Scaling factor for ajm not supported. Please use 1, 10, or 25." << endl; exit(1); }
      }
      else
      { std::cerr << "Scaling factor for aij not supported. Please use 1, 10, or 25." << endl; exit(1); }
    }
  else
  { std::cerr << "Scaling factor for mm not supported. Please use 1." << endl; exit(1); }
  }
  else
  { std::cerr << "Scaling factor for mortality rate not supported. Please use 1." << endl; exit(1); }
  
  string MFT_V = std::to_string(mV) + " * a + " + std::to_string(cV);
  string MFT_V_Prev = std::to_string(mV_Prev) + " * a + " + std::to_string(cV_Prev);
  string MFT_G = std::to_string(0);
  // string MFT_G = std::to_string(A_G) + " * ( 1 - exp( " + std::to_string(b_G) + " * ( a - " + std::to_string(a_c) + " )) )";
  string MFT_G_Prev = std::to_string(0);
  string MFT_Pr = std::to_string(0);
  string MFT_Pr_Prev = std::to_string(0);
  string MFT_W;
  if(A_W == 0)
   MFT_W = std::to_string(mW) + " * a + " + std::to_string(cW);
  else
   MFT_W = std::to_string(A_W) + " + " + std::to_string(b_W) + 
   " * exp( " + std::to_string(mW) + " * ( a - " + std::to_string(cW) + " ))";
  string MFT_W_Prev = std::to_string(mW_Prev) + " * a";
  string MFT_O = std::to_string(mO) + " * a + " + std::to_string(cO);
  string MFT_O_Prev = std::to_string(mO_Prev) + " * a";

  MFT_Vec_CoexExpr.assign({ MFT_V_Prev, MFT_G_Prev, MFT_Pr_Prev, MFT_W_Prev, MFT_O_Prev, MFT_V, MFT_G, MFT_Pr, MFT_W, MFT_O});
  // NOTE: The following clow is used for analytic MFT based initial conditions.
  //double clow[2*Sp] = {0, dP/dP, dP/(10.0*dP), 1, 1, 5, init_frac_pred, init_frac_pred, 1, 1};
  // NOTE: The following scaling_factor is used for analytic MFT based initial conditions.
  #if defined(INIT) && INIT == 0
    double scaling_factor[2*Sp] = {10, dP/dP, dP/(10.0*dP), 1,  1, 1, 1, 0.1*init_frac_pred, 1, 1};
    // USED FOR HOMOGENEOUS MFT BASED INITIAL CONDITIONS.
  #else
    double scaling_factor[2*Sp] = {0, dP/dP, dP/(10.0*dP)*init_frac_pred, 1, 1, 5, 1, 0.1*init_frac_pred, 1, 1};
    // USED FOR PERIODIC MFT BASED INITIAL CONDITIONS.
  #endif

  //*/


  // LE ORIGINAL (NORMAL CLOSE TO MFT)
  //double clow[2*Sp] = {0, dP/50000.0, dP/500000.0, 4, 20, 10000.0, Gstar, 10, 4, 10};
  // HIGH INIT, > MFT
  //double clow[2*Sp] = {0, dP/dP, dP/(10.0*dP), 4, 20, dP, Gstar*1.5, 15, 4, 10};
  // LOW INIT, < MFT
  //double clow[2*Sp] = {0, dP/dP, dP/(10.0*dP), 4, 20, 10000.0, 2, 0.6, 4, 10};
  // CLOSE TO MFT
  //double clow[2*Sp] = {0, dP/dP, dP/(10.0*dP), 4, 20, 10000.0, Gstar, 10, 4, 10};

  // NORMAL INIT CLOSE TO MFT
  double chigh[2*Sp] = {dP, dP/10000.0, dP/100000.0, 4, 20, 10000.0 + dP, Gstar, 10, 4, 10};
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

  stringstream ast, est, dPo, geq; ast << a_start; est  << a_end; dPo << dP, geq << setprecision(5) << Gstar;

	// This is done to avoid overwriting of files.
  string path_to_folder = frame_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "_Geq_" + geq.str();
  recursive_dir_create(path_to_folder);

  path_to_folder = prelim_folder + ast.str() + "-" + est.str() + "_dP_" + dPo.str() + "_Geq_" + geq.str() +"/TimeSeries";
  recursive_dir_create(path_to_folder);
  
  recursive_dir_create("../Data/Rietkerk/Stochastic/"+ std::to_string(SpB) +"Sp");
  
 
  first_order_critical_exp_delta_stochastic_MultiSp(div, t_max, a_start, a_end, a_c, c, gmax, alpha, rW, W0, D, v, K, sigma, A, H, E, M, pR, dtV, scaling_factor, dt, dx, dP, r, g, Gstar);

  return 0;
}