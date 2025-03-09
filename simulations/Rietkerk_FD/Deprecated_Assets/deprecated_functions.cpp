#include "rietkerk_bjork_basic.h"
#include "Debug.h"


/** Old version of var_mean_incremental_surv_runs() function. Retained for reference.
void var_mean_incremental_surv_runs_MultiSp(double t_avg_var_rep_N[][Sp4_1], double X_curr[][Sp2], int size, int j)
{
	/** X_curr contains time-series data (indexed as (N(t), <rho(t)>x)) from the current surviving replicate (given by r) with a length "size".
	 *  rep_avg_var stores the time(t), running density average and variance of X_curr, as well as N(t) and number of survivng replicates 
	 * for each time-step from previous surviving replicates ( <= r-1) respectively.
	 * In this function, X_curr is used to update the running avg and variance of rep_avg_var to account for current surviving replicate r(t).
	 * Note rep_avg_Var[t][3] stores number of surviving runs at t.
	*/
/**
	double mean_prev[size];
	for(int s=0; s < Sp; s++)
	{
		int s_eff = s-2;
		if(s ==0)
		{  s_eff = 0; }
		else if(s_eff <= 0)
		{	continue; }
		for(int i=0; i<size; i++)
		{
			int r= t_avg_var_rep_N[i][4*s_eff+ 5]+1; //Number of surviving replicates
			//Number of surviving run at time t.
			if(X_curr[i][2*s+ 0] != 0.0)
			{
				// Only surviving runs are considered.
				mean_prev[i] = t_avg_var_rep_N[i][4*s_eff+ 3]; //Stores X_n-1 mean data as temp variable.
				if(t_avg_var_rep_N[i][4*s_eff + 5] == 0.0)
				{	//First surviving run encountered at time t, encountered here.
					t_avg_var_rep_N[i][4*s_eff+ 3] = X_curr[i][2*s+ 1]; t_avg_var_rep_N[i][4*s+ 4] = 0.0; 
					t_avg_var_rep_N[i][4*s_eff+ 6] = X_curr[i][2*s+ 0];
				}
				else
				{	//Not the first surviving run, which means r = rep_avg_var[i][3]+1 >=2
					t_avg_var_rep_N[i][4*s_eff+ 3] = (X_curr[i][2*s+1] + (r-1)*t_avg_var_rep_N[i][4*s_eff+ 3])/r; //Note r>= 1.
					t_avg_var_rep_N[i][4*s_eff+ 6] = (X_curr[i][2*s+0] + (r-1)*t_avg_var_rep_N[i][4*s_eff+6])/r; //Note r>= 1.
					//Calculating and storing incremental variance (V^2_n) [For formula derivation, refer to notes.]
					t_avg_var_rep_N[i][4*s_eff+4] = (t_avg_var_rep_N[i][4*s_eff+3] - mean_prev[i])*(t_avg_var_rep_N[i][4*s_eff+3] - mean_prev[i]) +
					(1/(r-1))*( (r-2)*t_avg_var_rep_N[i][4*s_eff+4] + 
					(X_curr[i][2*s+1] - t_avg_var_rep_N[i][4*s_eff+3])*(X_curr[i][2*s+1] - t_avg_var_rep_N[i][4*s_eff+3]));
				}
				t_avg_var_rep_N[i][4*s_eff+5] = r; //At the end, update number of succesful runs.
			}
		}
	}
	for(int s=1; s < 2; s++)
	{
		for(int i=0; i<size; i++)
		{
		//For soil water and surface water.
			t_avg_var_rep_N[i][s] = (X_curr[i][2*s+1] + (j-1)*t_avg_var_rep_N[i][s])/j; //Note r>= 1.
		}
	}

}
**/

/** //NON RANGE VERSION OF THE FUNCTION generateNeighboringSitesFromCentral(). Retained for reference.
void generateNeighboringSitesFromCentral(const vector<pair<int, int>>& centralNeighboringSites, 
                                         vector<int>& neighboringSites, int i, int j, int L)
{
    // Resize neighboringSites if necessary
	if (neighboringSites.size() != centralNeighboringSites.size())
    	neighboringSites.resize(centralNeighboringSites.size());

    // Generate neighboring sites for the lattice site (i, j) based on the central neighboring sites for the lattice site (0, 0)
    for (size_t k = 0; k < centralNeighboringSites.size(); ++k) {
        const auto& offset = centralNeighboringSites[k];
        int dx = offset.first;
        int dy = offset.second;

        int nx = (i + dx + L) % L;  // Equivalent to (i + dx) mod L with reflective boundary conditions.
        int ny = (j + dy + L) % L;  // Equivalent to (j + dy) mod L with reflective boundary conditions.
        neighboringSites[k] = nx * L + ny;
    }
}
// */


/** Old version of the function generateNeighboringSitesFromCentral(). Retained for reference.
std::vector<int> generateNeighboringSitesFromCentral(const vector<pair<int, int>>& centralNeighboringSites, int i, int j, int L)
{
	vector<int> neighboringSites;
    // Generate neighboring sites for the lattice site (i, j) based on the central neighboring sites for the lattice site (0, 0)
    for (const auto& offset : centralNeighboringSites) {
        int dx = offset.first;
        int dy = offset.second;

        int nx = (i + dx + L) % L;  // Equivalent to (i + dx) mod L with reflective boundary conditions.
        int ny = (j + dy + L) % L;  // Equivalent to (j + dy) mod L with reflective boundary conditions.
        neighboringSites.push_back(nx * L + ny);
    }

    return neighboringSites;
} 
*/

// --------------------- OLD NON-GENERALISED FUNCTIONS FOR RIETKERK VEG MODEL--------------------- //


void f_2Dor_OLD(D2Vec_Double &f, D2Vec_Double &Rho_M, D3Vec_Int &nR2, double a, double c, double gmax, 
	double alpha, double d, double rW, double W0, double D[], double K[], double t, double dt, double dx1_2, double g)
{
	//Vector function that updates an array containing (dP/dt, dW/dt, DO/dt) for each site in the lattice.
    //Based on the Rietkerk model for plant vegetation dynamics with the Dornic twist where linear and stoch term for vegetation are already taken care of.
	// The following values are precomputed to reduce time complexity.
	double cgmax = c*gmax; double K2W0 = K[2]*W0;
	for(int i=0; i < g*g; i++)
	{
        //Equations for the density of plants, soil water and surface water at each site.
        //Note that the Laplacian is calculated using reflective boundary conditions.
        //The Laplacian is calculated using the 5-point stencil method.
        /**
		 * NOTE!!!!!:  Equivalent to dP/dt = c*g_max*P*W/(W+K1)
		 *  Linear and stochastic terms taken care of by Dornic integration routine previously.
		**/
		
        f[0][i] = cgmax*Rho_M[1][i]*Rho_M[0][i]/(Rho_M[1][i] +K[1]);
        //Equivalent to dW/dt = alpha*(P+K2*W)/(P+K2) - rW*W + D*(Laplacian of W)
        f[1][i] = alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[2][i] -gmax*Rho_M[1][i]*Rho_M[0][i]/(Rho_M[1][i] +K[1]) - rW*Rho_M[1][i] 
        + (D[1]*dx1_2)*(Rho_M[1][nR2[i][0][0]]  + Rho_M[1][nR2[i][0][1]]  + Rho_M[1][nR2[i][1][0]]  + Rho_M[1][nR2[i][1][1]]  - 4*Rho_M[1][i]);
        //Equivalent to dO/dt = a - alpha*(P+K*W)/(P+K)*O + D*(Laplacian of O)
		f[2][i] = a - alpha*(Rho_M[0][i]+ K[2]*W0)/(Rho_M[0][i] +K[2])*Rho_M[2][i] + (D[2]*dx1_2)*(Rho_M[2][nR2[i][0][0]]  + Rho_M[2][nR2[i][0][1]] 
        + Rho_M[2][nR2[i][1][0]]  + Rho_M[2][nR2[i][1][1]]  - 4*Rho_M[2][i]);
		
	}

}


void RK4_Integrate_Stochastic(D2Vec_Double &Rho_t, D2Vec_Double &Rho_tsar, D3Vec_Int &nR2,double a,double c,double gmax,double alpha,
		double d, double rW, double W0, double D[], double K[],double t,double dt,double dx, int g)
{
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	double dt6 = dt/6.0; double dt2 = dt/2.0; double dx1_2 = 1/(dx*dx);

	D2Vec_Double K1(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K2(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double K3(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K4(Sp, vector<double> (g*g, 0.0));
    D2Vec_Double Rho_M(Sp, vector<double> (g*g, 0.0));

	f_2Dor_OLD(K1, Rho_tsar, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t, dt, dx1_2, g); //K1 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K1[s][i];
	}

	f_2Dor_OLD(K2, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt2, dt, dx1_2, g); //K2 updated.

	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt2)*K2[s][i];
	}
	f_2Dor_OLD(K3, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt2, dt, dx1_2, g); //K3 updated.
	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
			Rho_M[s][i] = Rho_tsar[s][i] + (dt)*K3[s][i];
	}
	f_2Dor_OLD(K4, Rho_M, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t + dt, dt, dx1_2, g); //K4 updated.
    
	for(int i=0; i < g*g; i++)
	{
		for(int s= 0; s <Sp; s++)
		{
			Rho_t[s][i]+= (dt6)*( K1[s][i] + 2.0*K2[s][i] + 2.0*K3[s][i] + K4[s][i]);
	

			if( Rho_t[s][i] < 0 || isfinite(Rho_t[s][i]) == false || isnan(Rho_t[s][i]) == true)
			{
				stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m6 << "RK4 WAS KO'ED WITH:\t" << Rho_t[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
				<< " at time:\t:" << t << " For Species:\t:" << s << " with K1[s][i]:   " << K1[s][i] << "\t, K2[s][i]:\t" << K2[s][i] 
				<< "\t, K3[s][i]:\t" << K3[s][i] << "\t AND K4[s][i]:\t" << K4[s][i] << endl; //cout << m6.str();
				errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
			}
			
		}
	}
}

void rietkerk_Dornic_2D(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double d, 
	  double rW, double W0, double D[], double K[], double sigma[], double a_st, double a_end,  double dt, double dx, double dP, int r, int g)
{

	D3Vec_Int nR2 (g*g, D2Vec_Int (2, vector<int> (2, 0)));
    //nR2 is 3D [(L*L)x2x2] vector initialised to 0, which stores neighbours of all sites i (second dim) in a ball of NORM 1 (SQUARE) radius 1 (2ND & 3RD dim).
	// nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 

	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_Sq4(g, nR2); //Assigning neigbours.
	int tot_iter = 1 + t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events

	//t_meas[0] = 100; t_meas[1] = 200; t_meas[2] = 300; 

	//Defining Dornic variables for integration schema.
	//double bdt = b*dt; double cdt = c*dt; 
	double dx2 = dx*dx;
	double diff_coefficient[Sp]={0.0}; double beta[Sp]={0.0}; double lambda_exp[Sp]={0.0}; double lambda[Sp]={0.0};

	beta[0] = -d - 4*D[0]/(dx*dx);
	lambda_exp[0]= exp( (beta[0])*dt); lambda[0] = 2*(beta[0])/(sigma[0]*sigma[0]*(lambda_exp[0] -1.0));
	for(int s=0; s < Sp; s++)
	{
		diff_coefficient[s] = D[s]/(dx2);
	}
	
	//double rho_rep_avg_var[tot_iter][Sp4_1] ={0.0}; 
	D2Vec_Double rho_rep_avg_var(tot_iter, vector<double> (Sp4_1, 0.0));
	//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
	rho_rep_avg_var[0][0] = 0;
	for(int i =0; i< t_meas.size(); i++)
	{ rho_rep_avg_var[1+i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep

	double alpha_prime = diff_coefficient[0]*4; //Assume each Rho in neighbourhood for alpha estimation averages 1.
	double poiss_ru = lambda[0]*lambda_exp[0]*2.5; //Mean
	double mu_nought = 2.0*alpha_prime/(sigma[0]*sigma[0]) + poiss_ru;  // For Gamma, Beta =1. So Mean = 2*alpha/(sigma^2) + lambda = alpha/beta
	poiss_ru += 5*sqrt(poiss_ru); //Mean of Poisson Sampler + 5 SD.
	mu_nought += 6*sqrt(mu_nought); //Mean of  Mean of Gamma Sampler + 6 SD.

	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m0 << "Parameters:\t a:\t" << a << "\t c:\t " << c << "\t dt,dx:\t " << dt << "," << dx << "\t gmax:\t " << gmax 
	<< "\n Lambda:\t" << lambda[0] << "\t Beta:\t " << beta[0] << "\t dt,dx:\t " << dt << "," << dx << "\t Poisson Cut-Off  For Veg:\t " << poiss_ru 
	<< "\t and  Gamma Cut-Off For Veg:\t " << mu_nought  << "\t and  Alpha  For Veg:\t " << alpha  <<endl; //cout << m1.str();
	cout << m0.str();
	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
	errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	double kappa = c*gmax - d;
  	double p0istar, p0jstar, p0mstar; // Analytic steady state values.
	double astar = rW*d*K[1]/kappa; //Analytic critical point for a.

	if(a < astar)
	{
		p0istar = 0.0; p0jstar = a/rW; p0mstar = a/(alpha*W0);
	}
	else
	{
		p0jstar = d*K[1]/kappa; 
 		p0istar = (c/d)*(a - rW*p0jstar);
  		p0mstar = (a/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
	}

	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
		<< "\t with Total Iterations:\t " << tot_iter  <<endl; cout << m1.str();
		errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close();


		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number



		D2Vec_Double Rho_dt(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Updated through the turn, in turn used to update DRho.
		D2Vec_Double DRho(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.

		D2Vec_Double Rho_M(tot_iter, vector<double> (Sp2, 0.0)); //Stores <rho>x values per time step for a single replicate.
		double perc = 0.015; double c_high[Sp] ={p0istar + dP, p0jstar, p0mstar}; double c_low[Sp] ={p0istar, p0jstar, p0mstar};
		
		init_randconstframe(Rho_dt, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
		
		for(int i=0; i < g*g; i++)
		{
			for(int s=0; s < Sp; s++)
			{
				DRho[s][i] = Rho_dt[s][i]; //Assigning initial conditions
				Rho_tsar[s][i] =Rho_dt[s][i];
			}
		}
		stringstream m1_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1_1 << "Initial Conditions:\t Per:\t" << perc  << "\t C_High[0], C_High[1], C_High[2]:\t " << c_high[0] << "," << c_high[1] << "," << c_high[2]   
		<< "\t C_Lo[0], C_Lo[1], C_Lo[2]:\t " << c_low[0] << "," << c_low[1] << "," << c_low[2] << ",\t R: " << a << endl; //cout << m1.str();
		cout << m1.str();
		//std::ofstream errout; //Save Error Logs
		//std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
		errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close();

		poisson_distribution<int> poisson; gamma_distribution <double> gamma; normal_distribution<double> norm(0.0, 1.0);

		double t=0; int index = 0;  //Initialise t
		int iter_index=0; int po=1; int so =1; int lo=1;

		while( t < t_max + dt)
		{
			// Basic Book-keeping below. Prints frames out to file at given time points. 
			//Also updates number and avg density of active sites at given time points.
			if(t == 0 || t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{
				double rhox_avg, rhox_num; //Stores spatial average and number of occupied sites at given t.
				for(int s=0; s<Sp; s++)
				{
					vector <double> temp= {Rho_dt[s].begin(),Rho_dt[s].end()}; //Rho_dt for species 's'
                    rhox_avg = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
				    rhox_num = occupied_sites_of_vector(temp, g*g); //Finds number of occupied at given t.
				    Rho_M[iter_index][2*s] = rhox_num; Rho_M[iter_index][2*s +1] = rhox_avg; 

					vector<double>().swap(temp); //Flush temp out of memory.0

                    if(rhox_num == 0.0 and s== 0)
					    break; //Finish negotiations as all basal species are dead.
				}
				iter_index+=1;
				if( t >= 2000 && t <= 5000 || t== 0 ||  index >= tot_iter -3 &&  index <= tot_iter-1)
				{
					//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.
					
					stringstream L, tm ,d3, p1, rini, cgm, a1, a2, Dm0, Dm1, alph, w0t, dix, dimitri, sig0, sig1;

  					L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  					rini << j; Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; cgm <<  c*gmax;
					w0t << setprecision(3) << W0; alph <<  alpha; dimitri <<  dP;
					a1 << a_st; a2  << a_end; sig0 << sigma[0]; sig1 << sigma[1];
					// Three replicates are over.

					ofstream frame_dp;
  					// Creating a file instance called output to store output data as CSV
					frame_dp.open("../Data/Rietkerk/Frames/Stochastic/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RAND_Basic_P_c_DP_G_" 
					+ L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_alph_"+ alph.str() + "_cgmax_"+ cgm.str() 
					+ "_a_"+ p1.str() 
					+ "_dx_"+ dix.str() +  "_R_"+ rini.str() + ".csv"); //+ "_DV_"+ Dm0.str() + "_DW_"+ Dm1.str() 

					// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|
					frame_dp << "a_c,  x,  P(x; tmax), W(x; tmax), O(x; tmax) \n";
					for(int i=0; i< g*g; i++)
					{
						frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] <<endl;
				
					}

					frame_dp.close();
				}
				if( t > 0.0)
					index+=1;


			}

			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC. ONLY PLANTS EXPERIENCING DORNIC INTEGRATION.
			for(int s=0; s< 1; s++)
            {
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					double alpha_i = diff_coefficient[s]*(1*(DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
					DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]]));

					if(alpha_i == 0 && DRho[s][i] == 0) 
					{
						continue;  
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					double mu = -1.0 + 2.0*alpha_i/(sigma[s]*sigma[s]);
					double ziggy = lambda[s]*lambda_exp[s]*Rho_dt[s][i]; 
					double gru; // Stores the Gamma random variable.
					if(ziggy == 0.0)
					{
						gru = mu + 1.0;  // The Poisson distribution returns 0 in this case
					}
					else if(ziggy > 100)
					{
						// For large values of lambda, the Poisson distribution is approximated by a Gaussian distribution with mean lambda and variance lambda.
						long gauss = long(norm(rng)*sqrt(ziggy) + ziggy); // mu = 0, sigma = sqrt(lambda) and lambda = ziggy.
						gru = mu + 1.0 + gauss;
					}
					else
					{
						poisson = poisson_distribution<int>(ziggy);
						gru = mu + 1.0 + poisson(rng);
					}
					if(gru > 100)
					{
						// For large shape parameters (alpha = gru), the Gamma distribution is approximated by a Gaussian distribution with mean alpha/beta and variance alpha/(beta^2)
						double gauss = norm(rng)*sqrt(gru) + gru; // mu = 0, sigma = sqrt(gru) and lambda = gru. (as beta = 1)
						Rho_dt[s][i]= gauss/lambda[s];
						/**
						stringstream m6_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6_2 << "LARGE GRU: " << gru << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
						<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and Gamma (Gaussian) Approx Value:\t" << gauss << "\t and VEG RHO*:\t" << Rho_dt[s][i] << endl; 
						cout << m6_2.str(); errout.open(thr, std::ios_base::app); errout << m6_2.str(); errout.close(); */
					}
					else
					{
						gamma = gamma_distribution<double>(gru, 1.0);
						Rho_dt[s][i]= gamma(rng)/lambda[s];
					}
						if(isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
						{
							stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
							m6 << "YOU HIT ROCK BOTTOM WITH: Rho*[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
							<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i << 
							"\t and GAMMA ENTRY:\t" << gru << endl; cout << m6.str();
							errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
						}
						Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
				}
			}

			// Finally RK4 integration of the remaining terms.

			RK4_Integrate_Stochastic(Rho_dt, Rho_tsar, nR2, a, c, gmax, alpha, d, rW, W0, D, K, t, dt, dx, g);

			for(int s=0; s< Sp; s++)
			{
					
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					if( Rho_dt[s][i] < 0 || isinf(Rho_dt[s][i]) == true and so==1)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO FALLS BELOW O:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<<  endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					}

					if( isnan(Rho_dt[s][i] == true))
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO NAAN PARATHA:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<<  endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
						lo == -1;
					}
					DRho[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
					Rho_tsar[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
				}
			}

			t+=dt; //Update timestep

			if(lo == -1)
			{
				exit(3);
			}
		} // End of time loop.

		// Freeing up memory.
		vector<vector <double>>().swap(Rho_dt); vector<vector <double>>().swap(DRho); vector<vector <double>>().swap(Rho_tsar);

		if(int(omp_get_thread_num()) == 5)
		{
			size_t currentSize = getCurrentRSS( ); //Check Memory usage.
			size_t peakSize    = getPeakRSS( );
			stringstream m7;
			m7 << "In replicate " << j << ", current size in MB: " << currentSize/(1024.0*1024.0) 
			<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m7.str();
			errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();
		}

		if(j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "INITIAL UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << tot_iter << endl; cout << m7.str();
			for(int i =0; i< tot_iter; i++)
			{	
				for(int s =0; s< Sp-2; s++)
				{
					rho_rep_avg_var[i][4*s+ 3] = Rho_M[i][2*s + 1]; // Avg density of frame at given time point i.
					rho_rep_avg_var[i][4*s+ 4] = 0.0; // Variance in avg frame density.
					rho_rep_avg_var[i][4*s+ 6] = Rho_M[i][2*s]; // Number of active sites in frame at given time point i.

					if(Rho_M[i][2*s] > 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s+5 ] = 1;  } //One surviving run as of now
					else if(Rho_M[i][2*s] == 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s + 5] = 0;  } //No surviving run as of now
				}
				rho_rep_avg_var[i][1] = Rho_M[i][2*1 + 1]; // Avg density of frame at given time point i.
				rho_rep_avg_var[i][2] = Rho_M[i][2*2 + 1]; // Avg density of frame at given time point i.
				
				

			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{
			// Second or higher replicate, use incremental advances.
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter, j); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}

		// Freeing up memory.
		vector<vector <double>>().swap(Rho_M);

		if((j+1)%3 == 1)
		{	//Start logging at every multiple of 3
			stringstream L, tm ,d3, p1, a1, a2, dimitri, rini, Dm, cgm,  sig0;

  			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dimitri << dP;
  			rini << j; Dm << setprecision(3) << D[0]; cgm << c*gmax; sig0 << sigma[0]; a1 << a_st; a2  << a_end;
			// Three replicates are over.

			ofstream output_dp;
  			// Creating a file instance called output to store output data as CSV.
			output_dp.open("../Data/Rietkerk/Prelims/Stochastic/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/PRELIM_AGGRAND_Basic_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
			"_cgmax_"+ cgm.str() + "_R_"+ rini.str() + ".csv"); //+ "_Sig0_" + sig0.str() //"_D0_"+ Dm.str() + 

			// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
			output_dp << " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs, # Active Sites, \n";
			for(int i=0; i< tot_iter; i++)
			{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |

				output_dp << a << ","<< j << "," << g << "," << rho_rep_avg_var[i][0] << "," << rho_rep_avg_var[i][1] << "," << rho_rep_avg_var[i][2] << ",";
				for(int s=0; s <1; s++)
				{ output_dp << rho_rep_avg_var[i][4*s + 3] << "," <<
				rho_rep_avg_var[i][4*s +4] << "," << rho_rep_avg_var[i][4*s +5] << "," << rho_rep_avg_var[i][4*s +6] << ","; 
				}
				output_dp <<  endl;
				
			}

			output_dp.close();
		}	

	} // End of r loop.


	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m9.str();
	errout.open(thr, std::ios_base::app); errout << m9.str(); errout.close();

	
	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|	r	|    L 		|    t 		|		<<W(x,t)>>x,r	|		<<O(x,t)>>x,r	|
		//     <<Rho(t)>>x,r		|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |

		Rho.push_back({a, static_cast<double>(r), static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], 
					rho_rep_avg_var[i][2]});
		for(int s=0;s < Sp-2; s++)
		{

			
			Rho[i].insert(Rho[i].end(), { rho_rep_avg_var[i][4*s +3], rho_rep_avg_var[i][4*s +4], 
				rho_rep_avg_var[i][4*s +5], rho_rep_avg_var[i][4*s +6]}); 		
			/**if(s == 0)
			{	//Basal species
				Rho.push_back({a, c, static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][4*s + 1], 
					rho_rep_avg_var[i][4*s +2], rho_rep_avg_var[i][4*s + 3], rho_rep_avg_var[i][4*s + 4]});
			}
			else
			{ 	//Extend rows, don't create them.
				Rho[i].insert(Rho[i].end(), { rho_rep_avg_var[i][4*s +1], rho_rep_avg_var[i][4*s +2], 
				rho_rep_avg_var[i][4*s +3], rho_rep_avg_var[i][4*s +4]}); 
			}**/
		}
	}
	

	stringstream m10;
	m10 << "MA PA DHA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m10.str();
	errout.open(thr, std::ios_base::app); errout << m10.str(); errout.close();

}



void first_order_critical_exp_delta_stochastic_OLD(int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double d, double rW, double W0,  double D[], double K[], double sigma[], double dt, double dx, double dP, int r,  int g)
{

	
	//double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2
	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	//double mean[Sp] = {1.0, 0.05}; double sd[Sp] = {0.25, 0.0125};
	//init_randframe(Rho_0, Sp,  g*g); //Returns Rho_0 with a full initial frame filled with 0.2.
	vector<double> a_space = linspace(a_start, a_end, div);
	cout << "NOSTRA" <<endl;
  	// The pspace to iterate over.
    vector <double> t_measure = logarithm10_time_bins(t_max, dt);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

  	cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  	for (int i=0; i< t_measure.size(); i++)
  	{
		cout << t_measure[i] << " ";
  	} 	cout << endl;

  	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
    const int microToSeconds = 1000000;   
    const double delay1 = 5 * microToSeconds;     //5 seconds
    
    cout<<"Delay 1 in progress... ("<<delay1/microToSeconds<<"s)"<<endl;
	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

  	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );

	cout << "On initialisation, current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl;

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 32)
	{
		omp_set_num_threads(32); //Limiting use on Chunk. Don't be greedy.
	}

	//double perc = 0.015; double c_high[Sp] ={dP, p0j, p0m}; double c_low[Sp] ={p0i, p0j, p0m};

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

  	  
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {

		//init_randconstframe(Rh0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho1(t)>>x,r			|    Var[<Rho1(t)>x],r    |
		**/
		rietkerk_Dornic_2D(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K , sigma, a_start, a_end, dt, dx, dP, r, g);
		//RK4_Wrapper_2D(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K , a_start, a_end, dt, dx, dP, r, g);
		//expanded_percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, a_space[i], b, c, D, sigma, dt, dx, r, g);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
		vector<vector <double>>().swap(CExpRho_a);
		// Remove dynamic vectors from memory and free up allocated space.
      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting critical exponent data for each grid size in order.
          stringstream message3;
		  size_t currSize = getCurrentRSS( ); //Check Memory usage.
		  size_t peaksSize    = getPeakRSS( );
			 
          message3 << "Is this happening?\n" << "In thread: " << i << ", current size in MB: " << currSize/(1024.0*1024.0) 
		  << " and Peak Size (in MB): " << peaksSize/(1024.0*1024.0) << endl;
          cout << message3.str();

      }
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, coco, tm ,d3, p1, p2, rini, Dm0, Dm1, cgm, alph, dix, dimitri, Sig0;

	// Creating a string stream instance to store the values of the parameters in the file name.
	// This is done to avoid overwriting of files.

	//Check recursively if a folder exists, if not create one.
	//https://stackoverflow.com/questions/12510874/how-to-check-if-a-directory-exists
	/**
	struct stat info;
	if( stat( "../Data/Rietkerk/Frames", &info ) != 0 )
	{
		cout << "Cannot access ../Data/Rietkerk/Frames. Creating directory." << endl;
		const int dir_err = system("mkdir ../Data/Rietkerk/Frames");
		if (-1 == dir_err)
		{
			printf("Error creating directory!n");
			exit(1);
		}
	}*/

	

  	L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << a_start; p2  << a_end; rini << r; alph << alpha; cgm << c*gmax; coco << setprecision(4) << c;
  	Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; Sig0 << sigma[0]; dix << setprecision(2) << dx;
  	// setprecision() is a stream manipulator that sets the decimal precision of a variable.

	//Input a folder name with user-defined numbers a_st and a_end. Check recursively if this folder exists, if not create one.
	//https://stackoverflow.com/questions/12510874/how-to-check-if-a-directory-exists

	/**
	stringstream foldername;
	foldername << "../Data/Rietkerk/Frames/" << p1.str() << "-" << p2.str() << "/";
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
	*/

	ofstream output_1stdp;
  // Creating a file instance called output to store output data as CSV.
	output_1stdp.open("../Data/Rietkerk/Stochastic/1stOrderCC_Rietkerk_STOC_P_c_Delta_DP_G_" + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + "_cgmax_"+ cgm.str() + "_alpha_"+ alph.str() + "_R_"+ rini.str() + ".csv");
	 //+ "_D0_"+ Dm0.str() + "_D1_"+ Dm1.str() + "_S0_" + Sig0.str() +

	cout << "Save file name: " <<endl;
	cout << "../Data/Rietkerk/Stochastic/1stOrderCC_Rietkerk_STOC_P_c_Delta_DP_G_" + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() + "_D0_"+ Dm0.str() + "_D1_"+ Dm1.str() + "_S0_" + Sig0.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + "_cgmax_"+ cgm.str() + "_alpha_"+ alph.str() + rini.str() + ".csv";

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

	string header = " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs, # Active Sites";
	output_1stdp << header << "\n";
	cout << "The vector elements are: "<< endl;
  	cout << header << "\n";

	for(int i=0; i< vec.size(); i++)
	{
		for(int j=0; j< vec[i].size() -1; j++)
			output_1stdp << setprecision(16) << vec[i][j] << ",";
		output_1stdp << setprecision(16) << vec[i][vec[i].size() -1] <<  endl;
		if( i%(1000) ==1)
    	{
			for(int j=0; j< vec[i].size() -1; j++)
				cout << setprecision(16) << vec[i][j] << ",";
			cout << setprecision(16) << vec[i][vec[i].size() -1] <<  endl;
    	}
	}
	output_1stdp.close();
}



// --------------------- OLD NON-GENERALISED FUNCTIONS FOR RIETKERK 2SP (VEG & GRAZER) MODEL--------------------- //


void rietkerk_Dornic_2D_2Sp_OLD(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double rW, double W0, 
		double D[], double v[], double K[], double sigma[], double a_st, double a_end, double A[Sp][Sp], double H[Sp][Sp], double E[], double M[], double pR[], double dt, double dx, double dP, int r, int g)
{
	//Store max value in pR[] in a variable r_max.
	double r_max_effective = (*max_element(pR, pR + Sp))/dx;
	std::vector<std::pair<int, int>> origin_Neighbourhood = computeNeighboringSitesCentral(r_max_effective);
	auto range_origin = std::ranges::subrange(origin_Neighbourhood.begin(), origin_Neighbourhood.end()); 

	D3Vec_Int nR2 (g*g, D2Vec_Int (2, vector<int> (2, 0)));
    //nR2 is 3D [(L*L)x2x2] vector initialised to 0, which stores neighbours of all sites i (second dim) in a ball of NORM 1 (SQUARE) radius 1 (2ND & 3RD dim).
	// nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 

	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_Sq4(g, nR2); //Assigning neigbours.
	//Append the value 0.0 to the start of t_meas.
	//t_meas.insert(t_meas.begin(), 0.0); //Inserting 0.0 at the start of t_meas.
	int tot_iter = t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events

	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";

	/**
	stringstream m0_1;	//To make cout thread-safe as well as non-garbled due to race conditions.
    m0_1 << "Values in t_measure (which is of total size " << t_meas.size() << ") are :" <<endl;
  	for (int i=0; i< t_meas.size(); i++)
  	{
		m0_1 << t_meas[i] << " ";
  	} 	m0_1 << endl;//cout << m1.str();
	cout << m0_1.str();
	errout.open(thr, std::ios_base::out); errout << m0_1.str(); errout.close();
	**/

	//Defining Dornic variables for integration schema.
	//double bdt = b*dt; double cdt = c*dt; 
	double dx2 = dx*dx; double dx1 = 1/dx;
	double diff_coefficient[Sp]={0.0}; double beta[Sp]={0.0}; double lambda_exp[Sp]={0.0}; double lambda[Sp]={0.0};

	beta[0] = -M[0] - 4*D[0]/(dx2);
	lambda_exp[0]= exp( (beta[0])*dt); lambda[0] = 2*(beta[0])/(sigma[0]*sigma[0]*(lambda_exp[0] -1.0));
	diff_coefficient[0] = D[0]/(dx2);
	for(int s=1; s < Sp-2; s++)
	{
		diff_coefficient[s] = D[s]/(dx2);
		//diff_coefficient[s] = (D[s] -(dx/2.0)*(v[s]*(1 - dt*vdx[s])))/(dx2);
		// Advection leads to excess diffusion, hence the correction term.
		if(diff_coefficient[s] <= 0)
		{
			stringstream m0;	//To make cout thread-safe as well as non-garbled due to race conditions.
    		m0 << "You Fucked Up With Species:\t" << s << "\t with D_eff[s]:\t " << diff_coefficient[s] 
			<< "\t dt,dx:\t " << dt << "," << dx << "\t v:\t " << v[s] <<endl; //cout << m1.str();
			cout << m0.str();
			errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();
		}
	}
	
	//double rho_rep_avg_var[tot_iter][Sp4_1] ={0.0}; 
	D2Vec_Double rho_rep_avg_var(tot_iter, vector<double> (Sp4_1, 0.0));
	//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively. rho_rep_avg_var[0][0] = 0;
	for(int i =0; i< t_meas.size(); i++)
	{ rho_rep_avg_var[i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep

	double alpha_prime = diff_coefficient[0]*4; //Assume each Rho in neighbourhood for alpha estimation averages 1.
	double poiss_ru = lambda[0]*lambda_exp[0]*2.5; //Mean
	double mu_nought = 2.0*alpha_prime/(sigma[0]*sigma[0]) + poiss_ru;  // For Gamma, Beta =1. So Mean = 2*alpha/(sigma^2) + lambda = alpha/beta
	poiss_ru += 5*sqrt(poiss_ru); //Mean of Poisson Sampler + 5 SD.
	mu_nought += 6*sqrt(mu_nought); //Mean of  Mean of Gamma Sampler + 6 SD.

	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m0 << "Parameters:\t a:\t" << a << "\t c:\t " << c << "\t dt,dx:\t " << dt << "," << dx << "\t gmax:\t " << gmax 
	<< "\n Lambda:\t" << lambda[0] << "\t Beta:\t " << beta[0] << "\t dt,dx:\t " << dt << "," << dx << "\t Poisson Cut-Off  For Veg:\t " << poiss_ru 
	<< "\t and  Gamma Cut-Off For Veg:\t " << mu_nought  << "\t and  Alpha  For Veg:\t " << alpha  <<endl; //cout << m1.str();
	cout << m0.str();
	errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	double kappa = c*gmax - M[0];
  	double p0istar, p0jstar, p0mstar, p0nstar; // Analytic steady state values.

	//Represent soil water, surface water, basal vegetation and grazer respectively.
	
	double astar = rW*M[0]*K[1]/kappa; //Analytic critical point for a.
	if(a < 1/24.0)
	{
		p0istar = 0.0; p0mstar = a/rW; p0nstar = a/(alpha*W0); p0jstar =  dP/50000.0; //= 0.4
		//p0istar = a/rW; p0jstar = a/(alpha*W0); p0mstar = 0; p0nstar = dP/50.0;
	}
	else
	{
		//p0jstar = M[0]*K[1]/kappa; p0istar = (c/M[0])*(a - rW*p0jstar); p0mstar = (a/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
		p0mstar = M[0]*K[1]/kappa; //rhoi*= 8000; pojstar = 160;
 		p0istar = (c/M[0])*(a - rW*p0mstar);  p0jstar = p0istar/50.0;
  		p0nstar = (a/alpha)*(p0istar + K[2] )/(p0istar + K[2]*W0);
	}
	
	//const vector<int> initcsv_columns = {2, 3, 4}; //Columns to be read from csv file.
	//D2Vec_Double const_ind_val = {{1, dP/50000.0 }}; 
	// Const value of dp/50000 for Rho indexed 1 (i.e. GRAZER). Other initial values determined by burn-in csv files through init_csvconstframe(...)

	auto start_t = high_resolution_clock::now();
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
		<< "\t with Total Iterations:\t " << tot_iter  <<endl; cout << m1.str();
		errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close();


		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number



		D2Vec_Double Rho_dt(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Updated through the turn, in turn used to update DRho.
		D2Vec_Double DRho(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		//double Rho_M[tot_iter][Sp2] ={0.0}; //Stores <rho>x values per time step for a single replicate.
		D2Vec_Double Rho_M(tot_iter, vector<double> (Sp2, 0.0)); //Stores <rho>x values per time step for a single replicate.
		double perc = 0.015; double c_high[Sp] ={p0istar +dP, p0jstar, p0mstar, p0nstar}; double c_low[Sp] ={p0istar, p0jstar, p0mstar, p0nstar};
		init_randconstframe(Rho_dt, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.


		/**

		stringstream  rain, rep_j, grid; 
		rain  << a; rep_j << int(j%2); grid << g;
    	string initcsv_filename= "Burning_In_Frames/FRAME_RAND_Basic_P_c_DP_G_" + grid.str() +"_T_91201_dt_0.2_alph_0.00833333_cgmax_0.0208333_a_" 
		+ rain.str() +"_dx_0.1_R_"  + rep_j.str() + ".csv";
		const vector<int> initcsv_columns = {2, 3, 4}; //Columns to be read from csv file.
		init_csvconstframe(Rho_dt, const_ind_val, initcsv_filename, initcsv_columns, g*g); //Initialise Rho_dt with burn-in frame.
		
		vector<double> temp_vec= {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species 's'
		double init_veg_num = occupied_sites_of_vector(temp_vec, g*g); //Finds number of occupied at given t.
		vector<double>().swap(temp_vec); //Flush temp_vec out of memory.

		

		if(init_veg_num == 0.0)
		{
			stringstream m1_1;
			m1_1 << "RUN-TIME WARNING: ZERO Active BURN-IN veg sites for Thread Rank:\t " << omp_get_thread_num() 
			<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << " CONSIDER RE-RUNNING OVER A DIFFERENT RANGE." <<endl; 
			cout << m1_1.str(); errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();
		}

		stringstream m1_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1_1 << "Initial Conditions: BURN-IN # Active Veg sites:\t"<< init_veg_num << " for Thread Rank:\t " << omp_get_thread_num() 
		<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << endl; 
		cout << m1_1.str(); errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();

		*/

		
		for(int i=0; i < g*g; i++)
		{
			for(int s=0; s < Sp; s++)
			{
				DRho[s][i] = Rho_dt[s][i]; //Assigning initial conditions
				Rho_tsar[s][i] =Rho_dt[s][i];
			}
		}
		
		stringstream m1_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1_1 << "Initial Conditions:\t Per:\t" << perc  << "\t C_High[0], C_High[1], C_High[2], C_High[3] :\t " << c_high[0] << "," << c_high[1] << "," << c_high[2] << "," << c_high[3]   
		<< "\t C_Lo[0], C_Lo[1], C_Lo[2], C_Lo[3]:\t " << c_low[0] << "," << c_low[1] << "," << c_low[2] << "," << c_low[3] << ",\t R: " << a << endl; //cout << m1.str();
		cout << m1_1.str();
		//std::ofstream errout; //Save Error Logs
		//std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";
		errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();
		

		poisson_distribution<int> poisson; gamma_distribution <double> gamma; 
		uniform_real_distribution<double> unif; normal_distribution<double> norm(0.0, 1.0);
		

		double t=0; int index = 0;  //Initialise t
		int iter_index=0; int po=1; int so =1; int lo=1; int counter =0;

		while( t < t_max + dt)
		{

			// Basic Book-keeping below. Prints frames out to file at given time points. 
			//Also updates number and avg density of active sites at given time points.
			if(t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{

				double rhox_avg, rhox_num; //Stores spatial average and number of occupied sites at given t.
				for(int s=0; s<Sp; s++)
				{
					vector <double> temp= {Rho_dt[s].begin(),Rho_dt[s].end()}; //Rho_dt for species 's'
                    rhox_avg = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
				    rhox_num = occupied_sites_of_vector(temp, g*g); //Finds number of occupied at given t.
				    Rho_M[index][2*s] = rhox_num; Rho_M[index][2*s +1] = rhox_avg; 

					vector<double>().swap(temp); //Flush temp out of memory.0

                    if(rhox_num == 0.0 and s== 0)
					    break; //Finish negotiations as all basal species are dead.
				}

				vector <double> temp_alt= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
				rhox_num = occupied_sites_of_vector(temp_alt, g*g); //Finds number of occupied at given t.
				vector<double>().swap(temp_alt); //Flush temp out of memory.0

				/**
				stringstream m4;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m4 << "STATUS UPDATE AT TIME [t,i]\t" << t_meas[index] << " with # non-zero sites in VEG Rho_dt_0:  " << Rho_M[0][0] 
				<< " AND # Non-Zero in VEG DRho_0:  " << rhox_num <<   "  and # Non-Zero in Grazer Rho_dt_1:  " <<  Rho_M[0][2] << endl; cout << m4.str();
				errout.open(thr, std::ios_base::app); errout << m4.str(); errout.close(); */
				if( t >= 99 && t <= 10000 || t== 0 ||  index >= tot_iter -6 &&  index <= tot_iter-1)
				{
					//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.
					
					stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, alph, w0t, aij, hij, dix, dimitri, sig0, sig1;

  					L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  					rini << j; Dm0 << D[0]*pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; gm << setprecision(3) << gmax;
					w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; dimitri  << dP;
					a1 << a_st; a2  << a_end; sig0 << sigma[0]; sig1 << sigma[1];
					aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1];
					// Three replicates are over.

					ofstream frame_dp;

					if(t < 760)
					{
						//Only save one in three frames here.
						if(index%3 ==2 || t== 0)
						{
							frame_dp.open("../Data/Rietkerk/Frames/Stochastic/2Sp/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RAND_TwoSp_P_c_DP_G_" 
							+ L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str()  
							+ "_D1_"+ Dm1.str() + "_dx_"+ dix.str() + "_R_"+ rini.str() + ".csv");

							// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|   W(x, tmax) 		|    O(x, tmax) 		|
							frame_dp << "a_c,  x,  P(x; t), G(x; t), W(x; t), O(x; t) \n";
							for(int i=0; i< g*g; i++)
							{	frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] << "," << Rho_dt[3][i] << endl;	}

							frame_dp.close();

							stringstream m3;
							m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
							cout << m3.str(); errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
						}
					}
					else if(t >= 760 && t < 6000)
					{
						//Only save one in two frames here.
						if(index%2 ==0)
						{
							frame_dp.open("../Data/Rietkerk/Frames/Stochastic/2Sp/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RAND_TwoSp_P_c_DP_G_" 
							+ L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str()  
							+ "_D1_"+ Dm1.str() + "_dx_"+ dix.str() + "_R_"+ rini.str() + ".csv");

							// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|   W(x, tmax) 		|    O(x, tmax) 		|
							frame_dp << "a_c,  x,  P(x; t), G(x; t), W(x; t), O(x; t) \n";
							for(int i=0; i< g*g; i++)
							{	frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] << "," << Rho_dt[3][i] << endl;	}

							frame_dp.close();
							stringstream m3;
							m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
							cout << m3.str(); errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
						}
					}
					else
					{
						//Save all frames here.
						frame_dp.open("../Data/Rietkerk/Frames/Stochastic/2Sp/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RAND_TwoSp_P_c_DP_G_"
						+ L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str()
						+ "_D1_"+ Dm1.str() + "_dx_"+ dix.str() + "_R_"+ rini.str() + ".csv");

						// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|   W(x, tmax) 		|    O(x, tmax) 		|
						frame_dp << "a_c,  x,  P(x; t), G(x; t), W(x; t), O(x; t) \n";
						for(int i=0; i< g*g; i++)
						{	frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] << "," << Rho_dt[3][i] << endl;	}

						frame_dp.close();
						stringstream m3;
						m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
					}
  					// Creating a file instance called output to store output data as CSV
					/**
					frame_dp.open("../Data/Rietkerk/Frames/Stochastic/2Sp/BURNIN_" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/FRAME_RANDBURNIN_TwoSp_P_c_DP_G_" 
					+ L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str()  
					+ "_D1_"+ Dm1.str() + "_dx_"+ dix.str() + "_R_"+ rini.str() + ".csv");

					// Output =  | 	x		|    Rho0(x, tmax) 		|    Rho1(x, tmax) 		|   W(x, tmax) 		|    O(x, tmax) 		|
					frame_dp << "a_c,  x,  P(x; t), G(x; t), W(x; t), O(x; t) \n";
					for(int i=0; i< g*g; i++)
					{	frame_dp << a << "," << i << ","<< Rho_dt[0][i] << "," << Rho_dt[1][i] << "," << Rho_dt[2][i] << "," << Rho_dt[3][i] << endl;	}

					frame_dp.close();
					


					if(index%2 ==0)
					{
						stringstream m3;
						m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
						cout << m3.str(); errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
					} */
					
				}
				index+=1;
				/**
				stringstream m4_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m4_1 << "NEW INDEX VALUE\t" << index << " , WITH NEW T_MEAS:  " << t_meas[index] << endl; cout << m4_1.str();
				errout.open(thr, std::ios_base::app); errout << m4_1.str(); errout.close(); */
			}
			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC. ONLY LIVING MATTER EXPERIENCES DORNIC INTEGRATION.
			//First up vegetation.

			vector <double> temp_1= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
			double rhox_DR = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.
			vector<double>().swap(temp_1); //Flush temp out of memory.0

			temp_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
			double  rhox_dt = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.q
			vector<double>().swap(temp_1); //Flush temp out of memory.0
			/**
			stringstream m4_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
			m4_1 << "AFTER SAVING FRAME, STATUS UPDATE AT TIME [t,i]\t" << t_meas[index] << " with # non-zero sites in VEG Rho_dt_0:  " << rhox_DR 
			<< " AND # Non-Zero in VEG DRho_0:  " << rhox_dt <<  endl; cout << m4_1.str();
			errout.open(thr, std::ios_base::app); errout << m4_1.str(); errout.close(); **/
			//int counter_nzero_alpha=0;
			for(int s=0; s< 1; s++)
            {
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					double alpha_i = diff_coefficient[s]*((DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
					DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]]));

					if(alpha_i == 0 && DRho[s][i] == 0) 
					{
						/**
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "Much of VEG TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i << " \t with alpha_i, Rho(i), DRho(i):  "
						<< alpha_i << " , " << Rho_dt[s][i] << " , " << DRho[s][i] << "\t and Diff Const (D0/dx2): " << setprecision(16) << diff_coefficient[s] << endl; cout << m6.str();
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close(); */
						continue;  
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					//counter_nzero_alpha +=1;

					double mu = -1.0 + 2.0*alpha_i/(sigma[s]*sigma[s]);
					double ziggy = lambda[s]*lambda_exp[s]*Rho_dt[s][i]; 
					double gru; // Stores the Gamma random variable.
					if(ziggy == 0.0)
					{
						gru = mu + 1.0;  // The Poisson distribution returns 0 in this case
					}
					else if(ziggy > 100)
					{
						// For large values of lambda, the Poisson distribution is approximated by a Gaussian distribution with mean lambda and variance lambda.
						long gauss = long(norm(rng)*sqrt(ziggy) + ziggy); // mu = 0, sigma = sqrt(lambda) and lambda = ziggy.
						gru = mu + 1.0 + gauss;

						/**
						stringstream m6_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6_2 << "LARGE ZIGGY: " << ziggy << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
						<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and Gaussian Approx Value:\t" << gauss << "\t and GRU:\t" << gru << endl; 
						cout << m6_2.str(); errout.open(thr, std::ios_base::app); errout << m6_2.str(); errout.close(); */
					}
					else
					{
						poisson = poisson_distribution<int>(ziggy);
						gru = mu + 1.0 + poisson(rng);
					}
					if(gru > 100)
					{
						// For large shape parameters (alpha = gru), the Gamma distribution is approximated by a Gaussian distribution with mean alpha/beta and variance alpha/(beta^2)
						double gauss = norm(rng)*sqrt(gru) + gru; // mu = 0, sigma = sqrt(gru) and lambda = gru. (as beta = 1)
						Rho_dt[s][i]= gauss/lambda[s];
						/**
						stringstream m6_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6_2 << "LARGE GRU: " << gru << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
						<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and Gamma (Gaussian) Approx Value:\t" << gauss << "\t and VEG RHO*:\t" << Rho_dt[s][i] << endl; 
						cout << m6_2.str(); errout.open(thr, std::ios_base::app); errout << m6_2.str(); errout.close(); */
					}
					else
					{
						gamma = gamma_distribution<double>(gru, 1.0);
						Rho_dt[s][i]= gamma(rng)/lambda[s];
					}
					/**
					stringstream m6_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
					m6_1 << "DORNIC VEG RHO* " << Rho_dt[s][i] << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
					<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i  << "\t and POISSON ENTRY:\t" << ziggy << "\t and GAMMA ENTRY:\t" << gru << endl; 
					cout << m6_1.str(); errout.open(thr, std::ios_base::app); errout << m6_1.str(); errout.close(); **/
					if(isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "VEG HIT ROCK BOTTOM WITH: Rho*"<< s <<"[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i << 
						"\t and GAMMA ENTRY:\t" << gru << endl; cout << m6.str();
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					}
					Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
				}
			}
			// Now for higher order species.
			vector <double> temp_veg= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
			double rhox_num_veg = occupied_sites_of_vector(temp_veg, g*g); //Finds number of occupied sites at given t.
			double rhox_avg_veg = mean_of_vector(temp_veg, g*g); //Finds spatial average of densities at given t.
			vector<double>().swap(temp_veg); //Flush temp out of memory. 
			double nR_fac = 1 - rhox_num_veg/(g*g); //Factor to reduce the number of neighbours for gamma_i estimation
			if (nR_fac < 0.35)
			{	nR_fac = 0.35; }
			/**
			if(counter%50000 == 0)
			{

				temp_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
				rhox_dt = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.q
				vector<double>().swap(temp_1); //Flush temp out of memory.0

				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, thr]\t" << t << " , " << omp_get_thread_num()  << "\t VEG DORNIC DONE with non-zero elements in DRho_0:  " 
				<< rhox_num_veg << "\t and non-zero elements in Rho_dt_0:  " << rhox_dt << "\t and INDEX VAL, T_MEAS:\t" << index <<  " , " << t_meas[index] << endl; cout << m5_1.str();
				errout.open(thr, std::ios_base::app); errout << m5_1.str(); errout.close();
			}
			*/
			for(int s=1; s< Sp-2; s++)
            {
				for(int i=0;i<Rho_dt[0].size();i++)
				{

					//Firstly gamma needs to be estimated for each site.
					double gamma_i = 0.0;
					//To do this, first get neighbouring sites (i.e. perceptual range of species)
					int c_i = int(i/g); int c_j = i%g; //Current x and y coordinates of site i.
					std::vector<int> nR_Perp(origin_Neighbourhood.size(), 0); //Nearest Neighbours at site i
					

					generateNeighboringSitesFromCentral(range_origin, nR_Perp, c_i, c_j, g);
					// If nR_fac < 1, then reduce the number of neighbours to estimate gamma_i.
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i  
						<< "\t NEIGHBOURS FOUND WITH SIZE:" << nR_Perp.size() << endl; cout << m5_2.str();
						errout.open(thr, std::ios_base::app); errout << m5_2.str(); errout.close();
					}
					if(nR_fac < 1)
					{	nR_Perp.resize(int(nR_fac*nR_fac*nR_Perp.size())); } //= getNeighborsInBall(nR_Perp, nR_fac); }

					// Summing over all values in the perceptual range of species s (i.e. nR_Perp)
					
					for(int k=0; k< nR_Perp.size(); k++)
					{
						gamma_i += DRho[s-1][nR_Perp[k]];
					}
					//Finally normalising gamma_i by the number of sites in the perceptual range.
					//vector <double> temp= {Rho_dt[s].begin(),Rho_dt[s].end()}; //Rho_dt for species 's'
                    //double rhox_avg = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
					if(rhox_avg_veg == 0)
						gamma_i = 0;
					else
						gamma_i = gamma_i/(nR_Perp.size()*rhox_avg_veg); //Normalising gamma_i by the number of sites in the perceptual range.

					vector <int>().swap(nR_Perp); // vector <double>().swap(temp);  Remove dynamic vectors from memory and free up allocated space.

					//Finally check if gamma_i > 1, if so set it to 1.
					if(gamma_i > 1)
					{	gamma_i = 1; }
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i  
						<< "\t GAMMA FOUND WITH GAMMA_I:" << gamma_i << endl; cout << m5_2.str();
						errout.open(thr, std::ios_base::app); errout << m5_2.str(); errout.close();
					}

					//Next sampling the advection vector from v[s].
					unif = uniform_real_distribution<double>(0.0, 1.0);
					double ran = unif(rng);
					double vx = v[s]*cos(2*PI*ran); double vy = v[s]*sin(2*PI*ran); //Advection vector.
					double vx_abs = vx*sgn(vx); double vy_abs = vy*sgn(vy); //Absolute value of advection vector.

					// RECALL:	nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 

					double alpha_i = gamma_i*diff_coefficient[s]*(1*(DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
									DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]])) 
									+dx1*(1-gamma_i)*(vx_abs*DRho[s][nR2[i][1][sgn_index(vx)]]+ vy_abs*DRho[s][nR2[i][0][sgn_index(vy)]]);
					// NOTE: dx1 is the inverse of dx. sgn_index(vx) returns 0 if vx is positive, else 1.
					// Recall nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1;
					// ALSO NOTE: diff_coefficient[s] = D[s]/(dx*dx)
					if(alpha_i == 0 && DRho[s][i] == 0) 
					{
						continue;  
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					//Beta depends on the gamma_i value at each lattice site for each species, hence needs to be updated at each iteration.

					beta[s] = -M[s] - 4*gamma_i*diff_coefficient[s] - (1-gamma_i)*(vx_abs + vy_abs)*dx1; // NOTE: diff_coefficient[s] = D[s]/(dx*dx)
					lambda_exp[s]= exp( (beta[s])*dt); lambda[s] = 2*(beta[s])/(sigma[s]*sigma[s]*(lambda_exp[s] -1.0));

					double mu = -1.0 + 2.0*alpha_i/(sigma[s]*sigma[s]);
					double ziggy = lambda[s]*lambda_exp[s]*Rho_dt[s][i]; 
					double gru; // Stores the Gamma random variable.
					if(ziggy == 0.0)
					{
						gru = mu + 1.0;  // The Poisson distribution returns 0 in this case
					}
					else if(ziggy > 100)
					{
						// For large values of lambda, the Poisson distribution is approximated by a Gaussian distribution with mean lambda and variance lambda.
						long gauss = long(norm(rng)*sqrt(ziggy) + ziggy); // mu = 0, sigma = sqrt(lambda) and lambda = ziggy.
						gru = mu + 1.0 + gauss;
					}
					else
					{
						poisson = poisson_distribution<int>(ziggy);
						gru = mu + 1.0 + poisson(rng);
					}
					if(gru > 100)
					{
						// For large shape parameters (alpha = gru), the Gamma distribution is approximated by a Gaussian distribution with mean alpha/beta and variance alpha/(beta^2)
						double gauss = norm(rng)*sqrt(gru) + gru; // mu = 0, sigma = sqrt(gru) and lambda = gru. (as beta = 1)
						Rho_dt[s][i]= gauss/lambda[s];
					}
					else
					{
						gamma = gamma_distribution<double>(gru, 1.0);
						Rho_dt[s][i]= gamma(rng)/lambda[s];
					}
					if(isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "YOU HIT ROCK BOTTOM WITH: Rho*"<< s <<"[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i << 
						"\t and GAMMA ENTRY:\t" << gru << endl; cout << m6.str();
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
						exit(3);
					}
					Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
					/**
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i  
						<< "\t DORNIC VALUE OF GRAZER UPDATED WITH RHO*(I,t): " << Rho_tsar[s][i] << endl; cout << m5_2.str();
						errout.open(thr, std::ios_base::app); errout << m5_2.str(); errout.close();
					}
					*/
				}
			}
			/**
			if(counter%50000 == 0)
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, thr]\t" << t << " , " << omp_get_thread_num()  << "\t GRAZER DORNIC DONE." << endl; cout << m5_1.str();
				errout.open(thr, std::ios_base::app); errout << m5_1.str(); errout.close();
			}
			*/

			// Finally RK4 integration of the remaining terms.

			RK4_Integrate_Stochastic_2Sp(Rho_dt, Rho_tsar, nR2, a, c, gmax, alpha, rW, W0, D, K, A,H,E, t, dt, dx, g);
			/**
			if(counter%50000 == 0)
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, thr]\t" << t << " , " << omp_get_thread_num()  << "\t RK4 INTEGRATION OF REMAINDER TERMS DONE." << endl; cout << m5_1.str();
				errout.open(thr, std::ios_base::app); errout << m5_1.str(); errout.close();
			}
			*/


			for(int s=0; s< Sp; s++)
			{
					
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					if( Rho_dt[s][i] < 0 || isinf(Rho_dt[s][i]) == true and so==1)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO "<<s<<"  FALLS BELOW O:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<<  endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
					}

					if( isnan(Rho_dt[s][i] == true))
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO "<<s<<"  NAAN PARATHA:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<<  endl; cout << m6.str(); //so=0;
						errout.open(thr, std::ios_base::app); errout << m6.str(); errout.close();
						lo == -1;
					}
					DRho[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
					Rho_tsar[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
				}
			}

			t+=dt; //Update timestep
			counter+=1; //Update counter

			if(counter%90000 ==1)
			{
				//Find memory usage and time taken for every 50000 iterations.
				//size_t currentSize = getCurrentRSS( ); //Check Memory usage.
				//size_t peakSize    = getPeakRSS( );
				auto end_t = high_resolution_clock::now();
				auto elapsed_min = duration_cast<minutes>(end_t-start_t);
				// Reset the clock
				start_t = high_resolution_clock::now();
				stringstream m7;
				m7 << "STATUS UPDATE AT [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num() << " , " << j 
				//<< " ,with current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) 
				<< " and time taken for " << counter << " iterations: " << elapsed_min.count() << " min." <<  endl; cout << m7.str();
				errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();

				//Report syst every 50000 iterations.
			

			}

			if(lo == -1)
			{
				exit(3);
			}
		} // End of time loop.
		
		vector<vector <double>>().swap(Rho_dt); vector<vector <double>>().swap(DRho); vector<vector <double>>().swap(Rho_tsar);
		// Freeing up memory.
		if(int(omp_get_thread_num()) == 5)
		{
			size_t currentSize = getCurrentRSS( ); //Check Memory usage.
			size_t peakSize    = getPeakRSS( );
			stringstream m7;
			m7 << "In replicate " << j << ", current size in MB: " << currentSize/(1024.0*1024.0) 
			<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m7.str();
			errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();
		}

		if(j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "INITIAL UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << tot_iter << endl; cout << m7.str();
			for(int i =0; i< tot_iter; i++)
			{	
				for(int s =0; s< Sp-2; s++)
				{
					rho_rep_avg_var[i][4*s+ 3] = Rho_M[i][2*s + 1]; // Avg density of frame at given time point i.
					rho_rep_avg_var[i][4*s+ 4] = 0.0; // Variance in avg frame density.
					rho_rep_avg_var[i][4*s+ 6] = Rho_M[i][2*s]; // Number of active sites in frame at given time point i.

					if(Rho_M[i][2*s] > 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s+5 ] = 1;  } //One surviving run as of now
					else if(Rho_M[i][2*s] == 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s + 5] = 0;  } //No surviving run as of now
				}
				rho_rep_avg_var[i][1] = Rho_M[i][2*(Sp-2) + 1]; // Avg density of frame at given time point i.
				rho_rep_avg_var[i][2] = Rho_M[i][2*(Sp-1) + 1]; // Avg density of frame at given time point i.
				
			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{
			// Second or higher replicate, use incremental advances.
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter, j); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}

		// Freeing up memory.
		vector<vector <double>>().swap(Rho_M); 

		if((j+1)%3 == 1)
		{	//Start logging at every multiple of 3
			stringstream L, tm ,d3, p1, a1, a2, dimitri, rini, Dm, cgm,  sig0;

  			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dimitri << dP;
  			rini << j; Dm << D[0]; cgm << c*gmax; sig0 << sigma[0];
			// Three replicates are over.

			ofstream output_dp;
  			// Creating a file instance called output to store output data as CSV.
			output_dp.open("../Data/Rietkerk/Prelims/Stochastic/2Sp/" + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/PRELIM_AGGRAND_P_c_Delta_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
			"_D0_"+ Dm.str() + "_cgmax_"+ cgm.str() + "_R_"+ rini.str() + ".csv");

			// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
			output_dp << " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs Rho0, # Active Sites Rho0, <<Rho1(x; t)>_x>_r, Var[<Rho1(t)>_x]_r, # Surviving Runs Rho1, # Active Sites Rho1, \n";
			for(int i=0; i< tot_iter; i++)
			{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |

				output_dp << a << ","<< j << "," << g << "," << rho_rep_avg_var[i][0] << "," << rho_rep_avg_var[i][1] << "," << rho_rep_avg_var[i][2] << ",";
				for(int s=0; s <Sp-2; s++)
				{
					output_dp << rho_rep_avg_var[i][4*s + 3] << "," <<
					rho_rep_avg_var[i][4*s +4] << "," << rho_rep_avg_var[i][4*s +5] << "," << rho_rep_avg_var[i][4*s +6] << ","; 
				}
				output_dp <<  endl;
				
			}

			output_dp.close();
		}
		

	} // End of r loop.

	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m9.str();
	errout.open(thr, std::ios_base::app); errout << m9.str(); errout.close();

	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|    r 		|     L 		|    t 		|     <<W(t)>>x,r			|      <<O(t)>>x,r		|
		//..     <<Rho0(t)>>x,r			|    Var0[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |....

		Rho.push_back({a, static_cast<double>(r), static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], 
					rho_rep_avg_var[i][2]});
		for(int s=0;s < Sp-2; s++)
		{
			Rho[i].insert(Rho[i].end(), { rho_rep_avg_var[i][4*s +3], rho_rep_avg_var[i][4*s +4], 
			rho_rep_avg_var[i][4*s +5], rho_rep_avg_var[i][4*s +6]});		
		}
	}
	stringstream m10;
	m10 << "MA PA DHA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m10.str();
	errout.open(thr, std::ios_base::app); errout << m10.str(); errout.close();
}



void first_order_critical_exp_delta_stochastic_2Sp_OLD(int div, double t_max, double a_start, double a_end,  double c, double gmax, double alpha,
	double rW, double W0,  double D[], double v[], double K[], double sigma[], double A[Sp][Sp], double H[Sp][Sp], double E[], double M[], double pR[],
	double dt, double dx, double dP, int r,  int g)
{

	//double p0i = 0.5; double p0j= 2.25; double p0m= 8; // In g/m^2
	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	//double mean[Sp] = {1.0, 0.05}; double sd[Sp] = {0.25, 0.0125};
	//init_randframe(Rho_0, Sp,  g*g); //Returns Rho_0 with a full initial frame filled with 0.2.
	vector<double> a_space = linspace(a_start, a_end, div);
	cout << "NOSTRA" <<endl;
  	// The pspace to iterate over.
    vector <double> t_measure = logarithm10_time_bins(t_max, dt);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

	t_measure.insert(t_measure.begin(), 0.0);

  	cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  	for (int i=0; i< t_measure.size(); i++)
  	{
		cout << t_measure[i] << " ";
  	} 	cout << endl;

  	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
    const int microToSeconds = 1000000;   
    const double delay1 = 5 * microToSeconds;     //5 seconds
    
    cout<<"Delay 1 in progress... ("<<delay1/microToSeconds<<"s)"<<endl;
	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

  	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );

	cout << "On initialisation, current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl;

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 32)
	{
		omp_set_num_threads(32); //Limiting use on Chunk. Don't be greedy.
	}

	//double perc = 0.015; double c_high[Sp] ={dP, p0j, p0m}; double c_low[Sp] ={p0i, p0j, p0m};

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

  	  
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {

		//init_randconstframe(Rh0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho1(t)>>x,r			|    Var[<Rho1(t)>x],r    |
		**/
		rietkerk_Dornic_2D_2Sp(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, rW, W0, D, v, K , sigma, a_start, a_end, A,H,E,M, pR, dt, dx, dP, r, g);
		//RK4_Wrapper_2D(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K , a_start, a_end, dt, dx, dP, r, g);
		//expanded_percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, a_space[i], b, c, D, sigma, dt, dx, r, g);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
		vector<vector <double>>().swap(CExpRho_a);
		// Remove dynamic vectors from memory and free up allocated space.
      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting critical exponent data for each grid size in order.
          stringstream message3;
		  size_t currSize = getCurrentRSS( ); //Check Memory usage.
		  size_t peaksSize    = getPeakRSS( );
			 
          message3 << "Is this happening?\n" << "In thread: " << i << ", current size in MB: " << currSize/(1024.0*1024.0) 
		  << " and Peak Size (in MB): " << peaksSize/(1024.0*1024.0) << endl;
          cout << message3.str();

      }
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, coco, tm ,d3, p1, p2, rini, Dm0, Dm1, aij, hij, gm, alph, dix, dimitri, Sig0;

	

  	L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << a_start; p2  << a_end; rini << r; alph << alpha; gm << gmax; coco << setprecision(4) << c;
  	Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; Sig0 << sigma[0]; dix << setprecision(2) << dx; aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1];

	ofstream output_1stdp;
  // Creating a file instance called output to store output data as CSV.
	output_1stdp.open("../Data/Rietkerk/Stochastic/2Sp/1stOrderCC_Rietkerk_STOC_P_c_Delta_DP_G_" + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() + "_D1_"+ Dm1.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + 
	"_alpha_"+ alph.str() + "_R_"+ rini.str() + ".csv");

	cout << "Save file name: " <<endl;
	cout << "../Data/Rietkerk/Stochastic/2Sp/1stOrderCC_Rietkerk_STOC_P_c_Delta_DP_G_" + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() + "_D1_"+ Dm1.str() + "_S0_" + Sig0.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_aij_"+ aij.str() + "_hij_"+ hij.str() + "_dx_"+ dix.str() + 
	"_gmax_"+ gm.str() + "_alpha_"+ alph.str() + "_R_"+ rini.str() + ".csv" <<endl;

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

	string header = " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<Rho0(x; t)>_x>_r, Var[<Rho0(t)>_x]_r, # Surviving Runs Rho0, # Active Sites Rho0, <<Rho1(x; t)>_x>_r, Var[<Rho1(t)>_x]_r, # Surviving Runs Rho1, # Active Sites Rho1";
	output_1stdp << header << "\n";
	cout << "The vector elements are: "<< endl;
  	cout << header << "\n";

	for(int i=0; i< vec.size(); i++)
	{
		for(int j=0; j< vec[i].size() -1; j++)
			output_1stdp << setprecision(16) << vec[i][j] << ",";
		output_1stdp << setprecision(16) << vec[i][vec[i].size() -1] <<  endl;
		if( i%(1000) ==1)
    	{
			for(int j=0; j< vec[i].size() -1; j++)
				cout << setprecision(16) << vec[i][j] << ",";
			cout << setprecision(16) << vec[i][vec[i].size() -1] <<  endl;
    	}
	}
	output_1stdp.close();
}

// Function to compute the Dornic integration schema for the 2D 2 species model with stochastic terms.
// Superseded by MultiSp version.


void rietkerk_Dornic_2D_2Sp(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double rW, double W0, 
	double (&D)[Sp], double v[], double (&K)[3], double sigma[], double a_st, double a_end, double a_c, double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double pR[], 
	double chigh[], double clow[], double dt, double dx, double dP, int r, int g, double Gstar /* = -1.*/,  double Vstar /* = -1.*/ )
{
	double epsilon = 1.0e-12; //Small number to avoid division by zero.
	double perc =0.015; //Percentage of high density patches.
	//Store max value in pR[] in a variable r_max.
	//NOTE: SpB = Sp - 2; //Number of species excluding water terms.
	double r_max_effective = (*max_element(pR, pR + SpB))/dx;
	// Store pR[]/dx as a fraction of r_max_effective in r_frac[].
	//Initialise vector<std::pair<int, int>> r_frac with initial size SpB and set all values to 0.0.
	vector <std::pair<double, int>> r_frac(SpB, {0.0, 0});
 
	for(int s=0; s< SpB; s++)
	{	r_frac[s] = {pR[s]/(dx*r_max_effective), s};	}

	// Sort the r_frac[] array in descending order and store the indices in r_sort[].
	sort(r_frac.begin(), r_frac.end(), [](std::pair<double, int> &a, std::pair<double, int> &b) 
				{ return a.first > b.first; });
	std::vector<std::pair<int, int>> origin_Neighbourhood = computeNeighboringSitesCentral(int(r_max_effective));

	if(omp_get_thread_num()== 1)
	{
		for(int i=0; i< r_frac.size(); i++)
		{
			stringstream m0;	//To make cout thread-safe as well as non-garbled
			m0 << "r_frac[" << i << "]:\t" << r_frac[i].first << "\t and r_sort[" << i << "]:\t" << r_frac[i].second << endl;
			cout << m0.str();
		}
	}

	//Initialise variables for the RK4 integration.
	D3Vec_Int nR2 (g*g, D2Vec_Int (2, vector<int> (2, 0)));
	//nR2 is 3D [(L*L)x2x2] vector initialised to 0, which stores neighbours of all sites i (second dim) in a ball of NORM 1 (SQUARE) radius 1 (2ND & 3RD dim).
	// nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 

	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_Sq4(g, nR2); //Assigning neigbours.
	//Append the value 0.0 to the start of t_meas.
	//t_meas.insert(t_meas.begin(), 0.0); //Inserting 0.0 at the start of t_meas.
	int tot_iter = t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	t_max = t_meas[t_meas.size() -1]; //Only compute until end of sampling events

	//double rho_rep_avg_var[tot_iter][Sp4_1] ={0.0}; 
	D2Vec_Double rho_rep_avg_var(tot_iter, vector<double> (Sp4_1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively. rho_rep_avg_var[0][0] = 0;
	//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively. rho_rep_avg_var[0][0] = 0;
	for(int i =0; i< t_meas.size(); i++)
	{ rho_rep_avg_var[i][0] = t_meas[i];}  //Setting up time index of rho_t_N_R2_rep


	std::ofstream errout; //Save Error Logs
	std::string thr = "ErrLog_" + std::to_string(omp_get_thread_num()) + ".txt";

	//Defining Dornic variables for integration schema.
	//double bdt = b*dt; double cdt = c*dt; 
	double dx2 = dx*dx; double dx1 = 1/dx; double dx1_2 = 1/(dx2); double dxb2 = dx/2.0; 
	double diff_coefficient[Sp]={0.0}; double beta[Sp]={0.0}; double lambda_exp[Sp]={0.0}; double lambda[Sp]={0.0}; double sigma2_1[Sp]={0.0};
	vector<std::pair<double, double>> diff_eff(SpB, {0.0, 0.0}); //Stores effective diffusion coefficient for each species after correcting for advection.
	//vector <double> gamma_i(SpB, 0.0); //Stores gamma for each species.
	

	beta[0] = -M[0] - 4*D[0]/(dx2);
	lambda_exp[0]= exp( (beta[0])*dt); lambda[0] = 2*(beta[0])/(sigma[0]*sigma[0]*(lambda_exp[0] -1.0));
	diff_coefficient[0] = D[0]/(dx2); sigma2_1[0] = 1/(sigma[0]*sigma[0]);
	for(int s=1; s < SpB; s++)
	{
		diff_coefficient[s] = D[s]/(dx2);
		sigma2_1[s] = 1/(sigma[s]*sigma[s]);
		//diff_coefficient[s] = (D[s] -(dx/2.0)*(v[s]*(1 - dt*vdx[s])))/(dx2);
		// Advection leads to excess diffusion, hence the correction term.
		if(diff_coefficient[s] <= 0)
		{
			stringstream m0;	//To make cout thread-safe as well as non-garbled due to race conditions.
    		m0 << "You Fucked Up With Species:\t" << s << "\t with D_eff[s]:\t " << diff_coefficient[s] 
			<< "\t dt,dx:\t " << dt << "," << dx << "\t v:\t " << v[s] <<endl; //cout << m1.str();
			cout << m0.str();
			//errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();
		}
	}
	for(int s= SpB; s < Sp; s++)
		diff_coefficient[s] = D[s]/(dx2);
	double alpha_prime = diff_coefficient[0]*4; //Assume each Rho in neighbourhood for alpha estimation averages 1.
	double poiss_ru = lambda[0]*lambda_exp[0]*2.5; //Mean
	double mu_nought = 2.0*alpha_prime/(sigma[0]*sigma[0]) + poiss_ru;  // For Gamma, Beta =1. So Mean = 2*alpha/(sigma^2) + lambda = alpha/beta
	poiss_ru += 5*sqrt(poiss_ru); //Mean of Poisson Sampler + 5 SD.
	mu_nought += 6*sqrt(mu_nought); //Mean of  Mean of Gamma Sampler + 6 SD.

	stringstream m0;     //To make cout thread-safe as well as non-garbled due to race conditions.
    m0 << "Parameters:\t a:\t" << a << "\t c:\t " << c << "\t dt,dx:\t " << dt << "," << dx << "\t gmax:\t " << gmax 
	<< "\n Lambda:\t" << lambda[0] << "\t Beta:\t " << beta[0] << "\t dt,dx:\t " << dt << "," << dx << "\t Poisson Cut-Off  For Veg:\t " << poiss_ru 
	<< "\t and  Gamma Cut-Off For Veg:\t " << mu_nought  << "\t and  Alpha  For Veg:\t " << alpha  
	<< "\n with MFT biomass density of Vegetation at Coexistance = " << Vstar << " kg/km^2\n" << endl; //cout << m1.str();
	omp_get_thread_num() == 1 ? cout << m0.str() : cout << "";
	//errout.open(thr, std::ios_base::out); errout << m0.str(); errout.close();

	//const vector<int> initcsv_columns = {2, 3, 4}; //Columns to be read from csv file.
	
	// Const value of dp/50000 for Rho indexed 1 (i.e. GRAZER). Other initial values determined by burn-in csv files through init_csvconstframe(...)

	auto start_t = high_resolution_clock::now(); //Start the clock
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
		<< "\t with Total Iterations:\t " << tot_iter  <<endl; cout << m1.str();
		//errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close();


		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number

		D2Vec_Double Rho_dt(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Updated through the turn, in turn used to update DRho.
		D2Vec_Double DRho(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.
		D2Vec_Double Rho_tsar(Sp, vector<double> (g*g, 0.0)); // 2D vector of dim: Sx(L*l) Dummy variable, Kept constant, only updated at end of turn.

		//	/** MORE 2D VECTORS FOR ERROR -ANALYSIS 
		//D2Vec_Double Rho_dornic(Sp, vector<double> (g*g, 0.0)); //Stores <rho>x values after Dornic update for a single replicate. 2D vector of dim: Sx(L*l)

		//*/

		D2Vec_Double gamma(SpB, vector<double> (g*g, 0.0)); //Stores gamma for each species at each site.

		// 2D vectors of dim: Sx(L*l). Used for RK4 integration (of linear terms) of Rho_dt , after Dornic update.
		D2Vec_Double K1(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K2(Sp, vector<double> (g*g, 0.0));
    	D2Vec_Double K3(Sp, vector<double> (g*g, 0.0)); D2Vec_Double K4(Sp, vector<double> (g*g, 0.0));

		D2Vec_Double Rho_M(tot_iter, vector<double> (Sp2, 0.0)); //Stores <rho>x values per time step for a single replicate.
		//double Rho_M[tot_iter][Sp2] ={0.0}; //Stores <rho>x values per time step for a single replicate.

		double Rhox_avg[Sp] = {0.0}; double Rhox_num[Sp] = {0.0}; //Stores spatial average and number of occupied sites at given t.

		//init_randbistableframe(Rho_dt, g*g, a, a_c,  perc, chigh, clow); // Returns a frame with random speckles of high and low density.

		// MFT PERTURBATION BASED FRAME INITIALISATION
		// NOTE: clow[] is used to store the fractional change from MFT steady state for some species, corresponds to c_spread[] in init_exprtk_randbiMFTframe(...).
		#if defined(INIT) && INIT == 0
			// CORRESPONDS TO HOMOGENEOUS MFT FRAME INITIALISATION
			init_exprtk_homogenousMFTframe(Rho_dt, g*g, a, a_c, clow); // Returns a frame with random speckles of high and low density.
		#elif defined(INIT) && INIT == 2
			// BURN-IN FRAME INITIALISATION (EXPRTK)
		
			stringstream  rain, dPO, Lgrid, t_val; rain  << a; Lgrid << g; dPO << dP; t_val << 91201; 
			map<string, string> local_input_keys = input_keys; // Copy input_keys to local_input_keys.
			local_input_keys["T"] = t_val.str(); local_input_keys["a"] = rain.str(); local_input_keys["dP"] = dPO.str();
			local_input_keys["L"] = Lgrid.str();
			//Update the values of T, a, dP and L in local_input_keys.
			string local_input_frame_subfolder = input_frame_subfolder; //Copy input_frame_subfolder to local_input_frame_subfolder.
			// Formatted string for input_frame_subfolder.
			string local_input_folder = format_str(local_input_frame_subfolder, local_input_keys); //Format input_folder_pattern with local_input_keys.

			string initcsv_parendir= input_frame_parenfolder + local_input_folder; //Parent directory for csv files.

			/**INPUT Filename Patterns are of the form: FRAME_T_{}_a_{}_R_{}.csv (WHEN READING RIETKERK FRAMES)
			string csv_filename_pattern = input_prefix + t_val.str() + "_a_" + rain.str() + "_R_";
			const vector<string> initcsv_columns = {"P(x;t)", "W(x;t)", "O(x;t)"}; //Columns to be read from csv file.
			vector<int> const_species;
			for(int s=1; s< SpB; s++)
				const_species.push_back(s); //Species to be initialised with constant values.
			//*/
			
			///** INPUT WHEN READING ARTIFICIAL HEX FILES (Filename patterns: FRAME_T_0_a_0_R_{}.csv) 
			string csv_filename_pattern = input_prefix + "0_a_0_R_";
			const vector<string> initcsv_columns = {"P(x;t)"};
			vector<int> const_species;
			for(int s=1; s< Sp; s++)
				const_species.push_back(s); //Species to be initialised with constant values (ALL EXCEPT VEGETATION).
			//*/


			//Check first if initcsv_parendir exists, if not throw an error and exit.
			if(!std::filesystem::exists(initcsv_parendir))
			{
				stringstream m1_1;
				m1_1 << "RUN-TIME ERROR: Parent Directory " << initcsv_parendir << 
				"\nfor Initialisation CSV Files does not exist for Thread Rank:\t " << omp_get_thread_num() 
				<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << " EXITING." <<endl; 
				cout << m1_1.str(); cerr << m1_1.str(); errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();
				exit(5);
			}

			
			init_exprtk_readCSVcolumns_frame(Rho_dt, const_species, initcsv_parendir, csv_filename_pattern, initcsv_columns, a, a_c, g, j, clow);

			vector<double> temp_vec= {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species 's'
			double init_veg_num = occupied_sites_of_vector(temp_vec, g*g); //Finds number of occupied at given t.
			vector<double>().swap(temp_vec); //Flush temp_vec out of memory.

			

			if(init_veg_num == 0.0)
			{
				stringstream m1_1;
				m1_1 << "RUN-TIME WARNING: ZERO Active BURN-IN veg sites for Thread Rank:\t " << omp_get_thread_num() 
				<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << " CONSIDER RE-RUNNING OVER A DIFFERENT RANGE." <<endl; 
				cout << m1_1.str(); errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();
			}

			stringstream m1_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
			m1_1 << "Initial Conditions: BURN-IN # Active Veg sites:\t"<< init_veg_num << " for Thread Rank:\t " << omp_get_thread_num() 
			<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << endl; 
			cout << m1_1.str(); errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();
		#else
			// DEFAULT FRAME INITIALISATION (RANDOM MFT BASED SPECKLES)
			init_exprtk_randbiMFTframe(Rho_dt, g*g, a, a_c, dP, perc, clow); // Returns a frame with random speckles of high and low density.
		#endif
		/** // GAUSSIAN FRAME INITIALISATION  
		// Initialise vector <double> amp to elements of clow[].
		vector <double> amp(Sp, 500.0);  //Setting amplitude of gaussian distributions of vegetation to 500.
		for (int s=0; s< Sp; s++)
			amp[s] = clow[s+ Sp];
		
		vector <double> sd{g/8.0, g/16.0}; //Setting standard deviation of gaussian distributions of vegetation to g/8 and g/16.

		init_gaussframe(Rho_dt, g*g, sd, amp); // **/
	
		// The first SpB species stored in Rho_dt are to be initialised on frame as Gaussian distributions.
		// Species 0 is the vegetation, Species 1 is the grazer and Species 2 is the predator.
		// Species 0 should be centered near the top right corner of the grid, Species 1 near the bottom right corner and Species 2 near the bottom left corner. */
		

		/** // BURN-IN FRAME INITIALISATION
		
		stringstream  rain, rep_j, grid; 
		rain  << a; rep_j << int(j%2); grid << g;
    	string initcsv_filename= "Burning_In_Frames/FRAME_RAND_Basic_P_c_DP_G_" + grid.str() +"_T_91201_dt_0.2_alph_0.00833333_cgmax_0.0208333_a_" 
		+ rain.str() +"_dx_0.1_R_"  + rep_j.str() + ".csv";
		const vector<int> initcsv_columns = {2, 3, 4}; //Columns to be read from csv file.
		D2Vec_Double const_ind_val = {{1, dP/1000.0 }, {2, dP/5000.0}}; 
		init_csvconstframe(Rho_dt, const_ind_val, initcsv_filename, initcsv_columns, g*g); //Initialise Rho_dt with burn-in frame.
		
		vector<double> temp_vec= {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species 's'
		double init_veg_num = occupied_sites_of_vector(temp_vec, g*g); //Finds number of occupied at given t.
		vector<double>().swap(temp_vec); //Flush temp_vec out of memory.

		

		if(init_veg_num == 0.0)
		{
			stringstream m1_1;
			m1_1 << "RUN-TIME WARNING: ZERO Active BURN-IN veg sites for Thread Rank:\t " << omp_get_thread_num() 
			<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << " CONSIDER RE-RUNNING OVER A DIFFERENT RANGE." <<endl; 
			cout << m1_1.str(); errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();
		}

		stringstream m1_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1_1 << "Initial Conditions: BURN-IN # Active Veg sites:\t"<< init_veg_num << " for Thread Rank:\t " << omp_get_thread_num() 
		<< "  with a_value:\t" << a << "\t and Replicate:\t" << j << endl; 
		cout << m1_1.str(); errout.open(thr, std::ios_base::app); errout << m1_1.str(); errout.close();
		
		// */

		for(int i=0; i < g*g; i++)
		{
			for(int s=0; s < Sp; s++)
			{
				DRho[s][i] = Rho_dt[s][i]; //Assigning initial conditions
				Rho_tsar[s][i] =Rho_dt[s][i];
			}
		}
		
		stringstream m1_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
    	m1_2 << "Initial Conditions:\t Per:\t" << perc  << "\t C_High[0], C_High[1], C_High[2]:\t " << chigh[0] << "," << chigh[1] << "," << chigh[2]   
		<< "\t C_Lo[0], C_Lo[1], C_Lo[2]:\t " << clow[0] << "," << clow[1] << "," << clow[2] << ",\t R: " << a << "\n"; //cout << m1.str();
		if(omp_get_thread_num()== 1)
			cout << m1_2.str();  
		//errout.open(thr, std::ios_base::app); errout << m1.str(); errout.close(); 
		

		poisson_distribution<int> poisson; gamma_distribution <double> gamma_distr; 
		uniform_real_distribution<double> unif(0.0, 1.0); normal_distribution<double> norm(0.0, 1.0);

		double t=0; int index = 0;  //Initialise t
		int iter_index=0; int po=1; int so =1; int lo=1; int counter =0;

		while( t < t_max + dt)
		{
			// Basic Book-keeping below. Prints frames out to file at given time points. 
			//Also updates number and avg density of active sites at given time points.
			if(t == 0 || t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{
				double rhox_avg, rhox_num; //Stores spatial average and number of occupied sites at given t.
				for(int s=0; s<Sp; s++)
				{
					vector <double> temp= {Rho_dt[s].begin(),Rho_dt[s].end()}; //Rho_dt for species 's'
                    rhox_avg = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
				    rhox_num = occupied_sites_of_vector(temp, g*g); //Finds number of occupied at given t.
				    Rho_M[index][2*s] = rhox_num; Rho_M[index][2*s +1] = rhox_avg; 

					vector<double>().swap(temp); //Flush temp out of memory.0

                    //if(rhox_num == 0.0 and s== 0)
					//    break; //Finish negotiations as all basal species are dead.
				}

				vector <double> temp_alt= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
				rhox_num = occupied_sites_of_vector(temp_alt, g*g); //Finds number of occupied at given t.
				vector<double>().swap(temp_alt); //Flush temp out of memory.0

				// FRAME SAVING
				if(index >= tot_iter -10 &&  index <= tot_iter-1  ||  t >= 60000 && t <= 150000 || t >= 200 && t <= 12000 
					|| t== 0)
				{
					//Saving Rho_dt snapshots to file. This is done at times t= 0, t between 100 and 2500, and at time points near the end of the simulation.
					
					stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, Dm2, alph, w0t, aij, hij, dix, dimitri, sig0, sig1, veq;

  					L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  					rini << j; Dm0 << D[0]*pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2]; 
					a1 << a_st; a2  << a_end; sig0 << sigma[0]; sig1 << sigma[1]; dimitri  << dP; veq << setprecision(5) << Vstar;
					//gm << setprecision(3) << gmax; w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; 
					// Three replicates are over.

					string filenamePattern = frame_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
					+ "_a_" + p1.str()  + "_D1_"+ Dm1.str() + "_D2_"+ Dm2.str() + "_dx_"+ dix.str() + "_R_";
					//Next, designate the naive filename
					string filename = filenamePattern + rini.str() + ".csv";
					string parendir = "";

					// Creating a file instance called output to store output data as CSV.
					if(Vstar != -1)
						parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str();
					else
						parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();
					
					/** // SAVE ALL FRAMES
					save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);

					stringstream m3;
					m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
					cout << m3.str();
					// */
					

					// SAVE SELECTED FRAMES
					/** FRAME SAVING STRATEGY FOR LOW DT (DT ~ 0.01)
					if(t < 250)
					{
						//Only save one in FOUR frames here.
						if(index%4 ==0 || t== 0)
						{
							save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);

							stringstream m3;
							m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
							cout << m3.str(); //errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
						}
					}
					else if(t >= 250 && t < 600)
					{	//Only save one in two frames here.
						if(index%2 ==0)
						{
							save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
							stringstream m3;
							m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << "\n";
							cout << m3.str(); //errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
						}
					}//*/  
					///**  FRAME SAVING STRATEGY (DT ~ 0.1)
					if(t < 760)
					{
						//Only save one in three frames here.
						if(index%3 ==2 || t== 0)
						{
							save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);

							stringstream m3;
							m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << endl;
							cout << m3.str(); //errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
						}
					}
					else if(t >= 760 && t < 6000)
					{	//Only save one in two frames here.
						if(index%2 ==0)
						{
							save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
							stringstream m3;
							m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << "\n";
							cout << m3.str(); //errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
						}
					} //*/
					else
					{
						//Save all frames here.
						save_frame(Rho_dt, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g);
						stringstream m3;
						m3 << "FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << "\n";
						cout << m3.str(); //errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();
					}
					//*/
				}
				//

				// BLOCK FOR CALCULATING AND TEMP PRELIMINARY FRAMES
				if( index == int(tot_iter*0.8) || index == int(tot_iter*0.85) || index == int(tot_iter*0.9) || index == int(tot_iter*0.95) || index == int(tot_iter-1))
				{
					// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
					// and save to file.
					D2Vec_Double rho_rep_avg_var_temp(index+1, vector<double> (Sp4_1, 0.0)); //Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
					//std::copy(rho_rep_avg_var.begin(), rho_rep_avg_var.begin() + index, rho_rep_avg_var_temp.begin());
					var_mean_incremental_surv_runs(rho_rep_avg_var_temp, Rho_M, index+1, 0);

					//Time is stored in the first column.
					for(int i=0; i < index+1; i++)
					   rho_rep_avg_var_temp[i][0] = t_meas[i];

					//Finally save to file.
					stringstream L, tm ,d3, p1, a1, a2, dimitri, rini, Dm, veq, jID; // cgm, sig0;
					a1 << a_st; a2 << a_end;
					L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(5) << a; dimitri << dP;
					rini << j; Dm << setprecision(4) << D[2]; veq << setprecision(5) << Vstar;// cgm << c*gmax; sig0 << sigma[0]; 


					string parendir= "";
					if(Vstar != -1)
						parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str() + "/TimeSeries";
					else
						parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/TimeSeries";

					//double ran_jid = (unif(rng)*1000.0)/1000.0; jID << ran_jid; // Random number between 0 and 1.
					
					string filenamePattern = replicate_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
					"_D2_"+ Dm.str() + "_R_";

					save_prelimframe(rho_rep_avg_var_temp, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);
					
					//stringstream m3;
					//m3 << "TEMP PRELIM FRAME SAVED at time:\t" << t << " for Thread Rank:\t " << omp_get_thread_num() << "  with a_value:\t" << a << " and Replicate:\t" << j << "\n";
					//cout << m3.str(); errout.open(thr, std::ios_base::app); errout << m3.str(); errout.close();

					// Check if no vegetation is left at this current index. If so, break out of the while loop.
					if( Rho_M[index][0] == 0.0)
					{
						stringstream m3_1;
						m3_1 << "RUN-TIME WARNING: ZERO Active VEG sites for TIME [t, a, thr, j]\t" 
						<< t << " , " << a << " , " << omp_get_thread_num()  << " , " << j << " Skipping to next replicate .... \n"; 
						cout << m3_1.str(); cerr << m3_1.str();

												//Next save future preliminary frames with 0.0 values for population densities after this time point.
						int index_prelim_range[4] = {int(tot_iter*0.85), int(tot_iter*0.9), int(tot_iter*0.95), int(tot_iter-1)};
						int iter_index =0;
						for(int iter_index =0; iter_index < 4; iter_index++)
						{	
							int new_index = index_prelim_range[iter_index];
							if( index > new_index)
								continue; //Skip to next iteration if index is greater than the current range.
							// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
							// and save to file.
							D2Vec_Double rho_rep_avg_var_future_temp(new_index+1, vector<double> (Sp4_1, 0.0)); 
							//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
							var_mean_incremental_surv_runs(rho_rep_avg_var_future_temp, Rho_M, new_index+1, 0);

							//Time is stored in the first column.
							for(int i=0; i < new_index+1; i++)
								rho_rep_avg_var_future_temp[i][0] = t_meas[i];

							stringstream tm_future; tm_future << t_meas[new_index]; // Time at which the frame is saved.

							string filenamePattern_future = replicate_prefix + L.str() + "_T_" + tm_future.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
							"_D2_"+ Dm.str() + "_R_";

							//Finally save to file.
							save_prelimframe(rho_rep_avg_var_future_temp, parendir, filenamePattern_future, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);
							vector<vector<double>>().swap(rho_rep_avg_var_future_temp); //Flush temp out of memory.
						}
						break; //Break out of the while (time) loop.
					}
					

					vector<vector<double>>().swap(rho_rep_avg_var_temp); //Flush temp out of memory.
				} 
				// */

				vector <double> temp_1= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
				double rhox_DR = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.
				vector<double>().swap(temp_1); //Flush temp out of memory.0

				temp_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
				double  rhox_dt = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.q
				vector<double>().swap(temp_1); //Flush temp out of memory.0

				index+=1;
				/**
				stringstream m4_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m4_1 << "AFTER SAVING FRAME, STATUS UPDATE AT TIME [t,i]\t" << t_meas[index] << " with # non-zero sites in VEG Rho_dt_0:  " << rhox_DR 
				<< " AND # Non-Zero in VEG DRho_0:  " << rhox_dt <<  "\n"; cout << m4_1.str();
				errout.open(thr, std::ios_base::app); errout << m4_1.str(); errout.close(); **/
				//int counter_nzero_alpha=0;
			}

			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC. ONLY LIVING MATTER EXPERIENCES DORNIC INTEGRATION.
			//First up vegetation.

			for(int s=0; s< 1; s++)
            {
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					double alpha_i = diff_coefficient[s]*((DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
					DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]]));

					if(alpha_i == 0 && DRho[s][i] == 0) 
					{
						continue;  
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					//counter_nzero_alpha +=1;

					double mu = -1.0 + 2.0*alpha_i*sigma2_1[s]; // NOTE: sigma2_1 = 1/(sigma^2)
					double ziggy = lambda[s]*lambda_exp[s]*Rho_dt[s][i]; 
					double gru; // Stores the Gamma random variable.
					if(ziggy == 0.0)
					{
						gru = mu + 1.0;  // The Poisson distribution returns 0 in this case
					}
					else if(ziggy > 100)
					{
						// For large values of lambda, the Poisson distribution is approximated by a Gaussian distribution with mean lambda and variance lambda.
						long gauss = long(norm(rng)*sqrt(ziggy) + ziggy); // mu = 0, sigma = sqrt(lambda) and lambda = ziggy.
						gru = mu + 1.0 + gauss;
					}
					else
					{
						poisson = poisson_distribution<int>(ziggy);
						gru = mu + 1.0 + poisson(rng);
					}
					if(gru > 100)
					{
						// For large shape parameters (alpha = gru), the Gamma distribution is approximated by a Gaussian distribution with mean alpha/beta and variance alpha/(beta^2)
						double gauss = norm(rng)*sqrt(gru) + gru; // mu = 0, sigma = sqrt(gru) and lambda = gru. (as beta = 1)
						Rho_dt[s][i]= gauss/lambda[s];
					}
					else
					{
						gamma_distr = gamma_distribution<double>(gru, 1.0);
						Rho_dt[s][i]= gamma_distr(rng)/lambda[s];
					}
					/**
					stringstream m6_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
					m6_1 << "DORNIC VEG RHO* " << Rho_dt[s][i] << "  AT TIME [t, thr, i]\t" << t << " , " << omp_get_thread_num() << " , " << i 
					<< " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i  << "\t and POISSON ENTRY:\t" << ziggy << "\t and GAMMA ENTRY:\t" << gru << "\n"; 
					cout << m6_1.str(); errout.open(thr, std::ios_base::app); errout << m6_1.str(); errout.close(); **/
					if(isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "VEG HIT ROCK BOTTOM WITH: Rho*"<< s <<"[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and alpha_i:\t" << alpha_i << 
						"\t and GAMMA ENTRY:\t" << gru << "\n"; cout << m6.str(); cerr << m6.str();
					}
					//Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
				}
			} // End of Vegetation Integration
			
			// Book-keeping for determining gamma_i for higher order species. Calculating Rho averages per species at each time step.
			for(int s=0; s < SpB; s++)
			{
				vector <double> temp= {DRho[s].begin(),DRho[s].end()}; //Rho_dt for species '0'
				Rhox_avg[s] = mean_of_vector(temp, g*g); //Finds spatial average of densities at given t.
				vector<double>().swap(temp); //Flush temp out of memory.0
			}
			vector <double> temp_veg= {DRho[0].begin(),DRho[0].end()}; //Rho_dt for species '0'
			double rhox_num_veg = occupied_sites_of_vector(temp_veg, g*g); //Finds number of occupied sites at given t.
			vector<double>().swap(temp_veg); //Flush temp out of memory. 
			double nR_fac = 1 - rhox_num_veg/(g*g); //Factor to reduce the number of neighbours for gamma_i estimation
			if (nR_fac < 0.35)
			{	nR_fac = 0.35; }
			calc_gamma_2Sp_NonRefugia(origin_Neighbourhood, DRho, gamma, Rhox_avg, r_frac, nR_fac, r_max_effective, g); 
			//Calculates gamma for each species at each site.

			/** // BLOCK FOR SAVING GAMMA FRAMES
			if(t == 0 || t >= t_meas[index-1] -dt/2.0 && t < t_meas[index-1] +dt/2.0)// && index >= tot_iter -9 &&  index <= tot_iter)
			{
				// Saving gamma values to file at given time points.
				stringstream L, tm ,d3, p1, rini, a1, a2, Dm1, Dm2, dix, dimitri, veq;

  				L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  				rini << j; Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2];
				a1 << a_st; a2  << a_end; dimitri  << dP; veq << setprecision(5) << Vstar;

				string filenamePattern = gamma_prefix + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_" + p1.str()  + "_D1_"+ Dm1.str() + "_D2_"+ Dm2.str() + "_dx_"+ dix.str() + "_R_";
				
				string parendir = "";
				// Creating a file instance called output to store output data as CSV.
				if(Vstar != -1)
					parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str();
				else
					parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();

				string gammaheader = "a_c,  x,  Gam0(x; t), Gam1(x; t) \n"; //Header for output frame.
				save_frame(gamma, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
				
			}
			// */

			/** STATUS UPDATE
			if(counter%50000 == 0)
			{

				temp_1 = {Rho_dt[0].begin(),Rho_dt[0].end()}; //Rho_dt for species '0'
				rhox_dt = occupied_sites_of_vector(temp_1, g*g); //Finds number of occupied at given t.q
				vector<double>().swap(temp_1); //Flush temp out of memory.0

				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, thr]\t" << t << " , " << omp_get_thread_num()  << "\t VEG DORNIC DONE with non-zero elements in DRho_0:  " 
				<< rhox_num_veg << "\t and non-zero elements in Rho_dt_0:  " << rhox_dt << "\t and INDEX VAL, T_MEAS:\t" << index <<  " , " << t_meas[index] << "\n"; cout << m5_1.str();
				errout.open(thr, std::ios_base::app); errout << m5_1.str(); errout.close();
			}
			*/
			// Now for higher order species.
			for(int s=1; s< SpB; s++)
            {
				if(Rhox_avg[s] == 0.0)
					continue;
				//If the average density of a species is 0, then skip the Dornic integration for that species.
				
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					//Finally check if gamma_i > 1, if so set it to 1.
					if(gamma[s][i] > 1)
					{	gamma[s][i] = 1; }

					/** COUNTER GAMMA MESSAGE
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j 
						<< "\t GAMMA FOUND WITH GAMMA_I:" << gamma[s][i] << "\n"; cout << m5_2.str();
					} */
					
					if( gamma[s][i] < 0.0)
					{
						stringstream m5_3;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_3 << "GAMMA_I FALLS BELOW 0,  AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j  
						<< "\t WITH GAMMA_I:" << gamma[s][i] << " for species: " << s << "\n"; cout << m5_3.str(); cerr << m5_3.str(); 
					}
					/** // Check if gamma_i is NaN or Inf, if so, SAVE GAMMA & exit the program.
					if(isnan(gamma[s][i]) == true || isinf(gamma[s][i]) == true)
					{
						stringstream m5_3;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_3 << "GAMMA_I BLOWS THE LID,  AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j  
						<< "\t WITH GAMMA_I:" << gamma[s][i] << " for species: " << s << " with Rho_avg[0]: " << Rhox_avg[0] << " and Rho_avg[1]: " << Rhox_avg[1] << " and Rho_avg[2]: " << Rhox_avg[2] 
					  	<< " \n and Rho_t[0][i]: " << Rho_dt[0][i] << " and Rho_t[1][i]: " << Rho_dt[1][i] << " and Rho_t[2][i]: " << Rho_dt[2][i]<< "\n"; cout << m5_3.str();
						cerr << m5_3.str();

						stringstream L, tm ,d3, p1, rini, a1, a2, Dm1, Dm2, dix, dimitri, geq;

  						L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
  						rini << j; Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2];
						a1 << a_st; a2  << a_end; dimitri  << dP; geq << setprecision(5) << Gstar;

						string ErrGamfilenamePattern = "/ERROR_GAMMA_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
						+ "_a_" + p1.str()  + "_D1_"+ Dm1.str() + "_D2_"+ Dm2.str() + "_dx_"+ dix.str() + "_R_";
						string parendir = "";

						// Creating a file instance called output to store output data as CSV.
						if(Gstar != -1)
							parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Geq_" + geq.str();
						else
							parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();

						string gammaheader = "a_c,  x,  Gam0(x; t), Gam1(x; t), Gam2(x; t) \n"; //Header for output frame.
						save_frame(Rho_dt, parendir, ErrGamfilenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
						exit(4);
						
					}
					// */

					//Next sampling the advection vector from v[s].
					//unif = uniform_real_distribution<double>(0.0, 1.0);
					double ran = unif(rng);
					double vx = v[s]*cos(2*PI*ran); double vy = v[s]*sin(2*PI*ran); //Advection vector.
					double vx_abs = vx*sgn(vx); double vy_abs = vy*sgn(vy); //Absolute value of advection vector.

					//RECALL: dx1_2 = 1/(dx2); dxb2 = dx/2.0;
					diff_eff[s].first = (D[s] -(dxb2)*(vx_abs*(1 - dt*dx1*vx_abs)))*dx1_2;
					diff_eff[s].second = (D[s] -(dxb2)*(vy_abs*(1 - dt*dx1*vy_abs)))*dx1_2;
					/** CORRECT DsC VERSION 
					diff_eff[s].first = (D[s] -(dxb2)*(gamma[s][i]*vx_abs*(1 - dt*dx1*gamma[s][i]*vx_abs)))*dx1_2;
					diff_eff[s].second = (D[s] -(dxb2)*(gamma[s][i]*vy_abs*(1 - dt*dx1*gamma[s][i]*vy_abs)))*dx1_2;
					*/
					//diff_coefficient[s] = (D[s] -(dx/2.0)*(v[s]*(1 - dt*vdx[s])))/(dx2);
					// Advection leads to excess diffusion, hence the correction term.

					if(diff_eff[s].first < diff_coefficient[0])
						diff_eff[s].first = diff_coefficient[0];
					
					if(diff_eff[s].second < diff_coefficient[0])
						diff_eff[s].second = diff_coefficient[0];
					// If the effective diffusion coefficient is less than D[0]/dx2, set it to D[0]/dx2.
					

					// RECALL:	nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1; 
					// NOTE: dx1 is the inverse of dx. sgn_index(vx) returns 0 if vx is positive, else 1.
					// Recall nR2[i][0][0] = i - g; nR2[i][0][1] = i + g;  nR2[i][1][0] = i - 1; nR2[i][1][1] = i + 1;
					// ALSO NOTE: diff_coefficient[s] = D[s]/(dx*dx); diff_eff[s].first = D_X_eff[s]/(dx*dx); diff_eff[s].second = D_Y_eff[s]/(dx*dx);
					
					// Standard diffusion term.
					//double alpha_i = gamma[s][i]*diff_coefficient[s]*(1*(DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]] +
					//				DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]])) 
					//				+dx1*(1-gamma[s][i])*(vx_abs*DRho[s][nR2[i][1][sgn_index(vx)]]+ vy_abs*DRho[s][nR2[i][0][sgn_index(vy)]]);

					// Advection corrected diffusion term.
					double alpha_i = gamma[s][i]*(diff_eff[s].second*(DRho[s][nR2[i][0][0]] + DRho[s][nR2[i][0][1]]) 
												+ diff_eff[s].first*(DRho[s][nR2[i][1][0]] + DRho[s][nR2[i][1][1]])) 
									+dx1*(1-gamma[s][i])*(vx_abs*DRho[s][nR2[i][1][sgn_index(vx)]]+ vy_abs*DRho[s][nR2[i][0][sgn_index(vy)]]);
					
					if(alpha_i == 0 && DRho[s][i] == 0) 
					{
						continue;  
					} //Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.

					//Beta depends on the gamma_i value at each lattice site for each species, hence needs to be updated at each iteration.
					// Standard diffusion term.
					//beta[s] = -M[s] - 4*gamma[s][i]*diff_coefficient[s] - (1-gamma[s][i])*(vx_abs + vy_abs)*dx1; // NOTE: diff_coefficient[s] = D[s]/(dx*dx)
					// Advection corrected diffusion term.
					beta[s] = -M[s] - 2*gamma[s][i]*(diff_eff[s].first + diff_eff[s].second) - (1-gamma[s][i])*(vx_abs + vy_abs)*dx1;

					lambda_exp[s]= exp( (beta[s])*dt); lambda[s] = 2*(beta[s])*(sigma2_1[s])/(lambda_exp[s] -1.0); // NOTE: sigma2_1 = 1/(sigma^2);

					double mu = -1.0 + 2.0*alpha_i*sigma2_1[s]; // NOTE: sigma2_1 = 1/(sigma^2)
					double ziggy = lambda[s]*lambda_exp[s]*DRho[s][i]; 
					double gru; // Stores the Gamma random variable.
					if(ziggy == 0.0)
					{
						gru = mu + 1.0;  // The Poisson distribution returns 0 in this case
					}
					else if(ziggy > 100)
					{
						// For large values of lambda, the Poisson distribution is approximated by a Gaussian distribution with mean lambda and variance lambda.
						long gauss = long(norm(rng)*sqrt(ziggy) + ziggy); // mu = 0, sigma = sqrt(lambda) and lambda = ziggy.
						gru = mu + 1.0 + gauss;
					}
					else
					{
						poisson = poisson_distribution<int>(ziggy);
						gru = mu + 1.0 + poisson(rng);
					}
					if(gru > 100)
					{
						// For large shape parameters (alpha = gru), the Gamma distribution is approximated by a Gaussian distribution with mean alpha/beta and variance alpha/(beta^2)
						double gauss = norm(rng)*sqrt(gru) + gru; // mu = 0, sigma = sqrt(gru) and lambda = gru. (as beta = 1)
						Rho_dt[s][i]= gauss/lambda[s];
					}
					else
					{
						gamma_distr = gamma_distribution<double>(gru, 1.0);
						Rho_dt[s][i]= gamma_distr(rng)/lambda[s];
					}
					if(isnan(Rho_dt[s][i]) == true || isinf(Rho_dt[s][i]) == true)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "YOU HIT ROCK BOTTOM WITH: Rho*"<< s <<"[t,i]\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho"<< s <<"(t-dt,i):   " << DRho[s][i] << "\t and ( alpha_i, beta_i ):\t" << alpha_i << " , " << beta[s] 
						<< "\t and ( POISSON, GAMMA) ENTRIES:\t" << ziggy << " , " << gru << "\t and (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]" 
						<<  " AND (vx, vy, gamma_i):\t [" << vx << " , " << vy << " , " << gamma[s][i] << " ]" << "\n"; cout << m6.str(); cerr << m6.str();
						cout << m6.str(); cerr << m6.str();

						stringstream m6_1;     //To make cout thread-safe as well as non-garbled
						m6_1 << "Saving error frames and exiting program .... \n"; cout << m6_1.str(); cerr << m6_1.str();

						// SAVING ERROR Rho_dt, DRho, and gamma vectors TO FILE
						stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, Dm2, alph, w0t, aij, hij, dix, dimitri, sig0, sig1, veq;

						L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
						rini << j; Dm0 << D[0]*pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2]; 
						a1 << a_st; a2  << a_end; sig0 << sigma[0]; sig1 << sigma[1]; dimitri  << dP; veq << setprecision(5) << Vstar;
						//gm << setprecision(3) << gmax; w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; 
						// Three replicates are over.

						string frame_parendir = "";

						// Creating a file instance called output to store output data as CSV.
						if(Vstar != -1)
							frame_parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str();
						else
							frame_parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();

						// SAVING RHO_DT FRAME
						string filenamePattern_Rhodt = "/ERROR_RHO_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
						+ "_a_" + p1.str()  + "_dx_"+ dix.str() + "_R_";
						//Next, designate the naive filename
						string filename = filenamePattern_Rhodt + rini.str() + ".csv";
						save_frame(Rho_dt, frame_parendir, filenamePattern_Rhodt, a, a_st, a_end, t, dt, dx, dP, j, g);

						// SAVING DRHO FRAME
						string filenamePattern_DRho = "/ERROR_DRHO_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
						+ "_a_" + p1.str()  + "_dx_"+ dix.str() + "_R_";
						//Next, designate the naive filename
						filename = filenamePattern_DRho + rini.str() + ".csv";
						save_frame(DRho, frame_parendir, filenamePattern_DRho, a, a_st, a_end, t, dt, dx, dP, j, g);

						/** // SAVING GAMMA FRAME
						string filenamePattern_Gamma = "/ERROR_GAMMA_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
						+ "_a_" + p1.str()  + "_dx_"+ dix.str() + "_R_";
						//Next, designate the naive filename
						filename = filenamePattern_Gamma + rini.str() + ".csv";
						 
						string gammaheader = "a_c,  x,  Gam0(x; t), Gam1(x; t) \n"; //Header for output frame.
						save_frame(gamma, frame_parendir, filenamePattern_Gamma, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
						*/
						
						// SAVING PRELIMINARY FRAME
						// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
						// and save to file.
						D2Vec_Double rho_rep_avg_var_temp(index, vector<double> (Sp4_1, 0.0)); 
						
						//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
						//std::copy(rho_rep_avg_var.begin(), rho_rep_avg_var.begin() + index, rho_rep_avg_var_temp.begin());
						var_mean_incremental_surv_runs(rho_rep_avg_var_temp, Rho_M, index, 0);

						string filenamePattern_Prelim = "/ERROR_SERIES_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
						"_D1_"+ Dm1.str() + "_R_";
						//Next, designate the naive filename

						string prelim_parendir= "";
						if(Vstar != -1)
							prelim_parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str() + "/TimeSeries";
						else
							prelim_parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/TimeSeries";


						save_prelimframe(rho_rep_avg_var_temp, prelim_parendir, filenamePattern_Prelim, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);

						exit(3);
					}
					//Rho_tsar[s][i] = Rho_dt[s][i];		// This is the rho* value (refer to Dornic et 2005)
					/** COUNTER STATUS UPDATE
					if(counter%50000 == 0 && counter == i)
					{
						stringstream m5_2;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m5_2 << "STATUS UPDATE AT TIME [t, thr, i, j]\t" << t << " , " << omp_get_thread_num() << " , " << i << " , " << j 
						<< "\t DORNIC VALUE OF SP " << s << " UPDATED WITH RHO*(I,t): " << Rho_tsar[s][i] 
						<< " and (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]" << "\n"; cout << m5_2.str();
						errout.open(thr, std::ios_base::app); errout << m5_2.str(); errout.close();
					}
					*/
				}
			}
			/** // Finally, update Rho_tsar with the new values of Rho_dt for the abiotic species.
			for(int s=SpB; s< Sp; s++)
			{
				for(int i=0;i<Rho_dt[0].size();i++)
					Rho_tsar[s][i] = Rho_dt[s][i];		
					// This is the rho* value (refer to Dornic et 2005)
			} //**/

			
			
			if(counter%50000 == 0)
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num() << " , " << j << "\t GRAZER + PRED DORNIC DONE, FOR " << prefix + "\n"; 
				cout << m5_1.str(); cerr << m5_1.str();
			}
			
			// NOTE: diff_coefficient[s] = D[s]/(dx*dx)
			// Finally RK4 integration of the remaining terms.

			RK4_Integrate_Stochastic_2Sp(Rho_dt, Rho_tsar, K1, K2, K3, K4, nR2, a, c, gmax, alpha, rW, W0, diff_coefficient, K, A,H,E, t, dt, dx, g);
			
			if(counter%50000 == 0)
			{
				stringstream m5_1;     //To make cout thread-safe as well as non-garbled due to race conditions.
				m5_1 << "STATUS UPDATE AT TIME [t, a, thr, j]\t" << t << " , " << a << " , " << omp_get_thread_num()  << " , " << j << "\t RK4 INTEGRATION OF REMAINDER TERMS DONE." 
				<< "\n"; cout << m5_1.str(); cerr << m5_1.str();
			}
			
			for(int s=0; s< Sp; s++)
			{
				so = s;	
				for(int i=0;i<Rho_dt[0].size();i++)
				{
					// CHECK FOR NAN AND INF
					if( Rho_dt[s][i] < 0 || isinf(Rho_dt[s][i]) == true and so==s)
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO "<<s<<"  FALLS BELOW O:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i] << "\t and (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]" <<  "\n"; 
						cout << m6.str(); cerr << m6.str(); so == -1; lo == -1;
					}
					// */ 

					if( isnan(Rho_dt[s][i] == true))
					{
						stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
						m6 << "RHO "<<s<<"  NAAN PARATHA:\t" << Rho_dt[s][i] << "\t at index:  " << i << " and thread:  " << omp_get_thread_num() 
						<< " at time:\t:" << t << " with Rho(t-dt,i):   " << DRho[s][i] << "\t and Analytic Integration Term:\t" 
						<< Rho_tsar[s][i]<< "\t and (DX, DY):  [" << diff_eff[s].first << " , " << diff_eff[s].second <<  " ]" <<  "\n"; 
						cout << m6.str(); cerr << m6.str();
						lo == -1;
					}
					// Set all lattice sites that are less than 10^-12 to 0.0.
					if(Rho_dt[s][i] < epsilon && Rho_dt[s][i] > 0.0)
						Rho_dt[s][i] = 0.0;	
					DRho[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
					//Rho_tsar[s][i] =Rho_dt[s][i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
				}
			}

			t+=dt; //Update timestep
			counter+=1; //Update counter

			if(counter%40000 ==1)
			{
				//Find memory usage and time taken for every 50000 iterations.
				//size_t currentSize = getCurrentRSS( ); //Check Memory usage.
				//size_t peakSize    = getPeakRSS( );
				auto end_t = high_resolution_clock::now();
				auto elapsed_min = duration_cast<minutes>(end_t-start_t);
				// Reset the clock
				start_t = high_resolution_clock::now();
				stringstream m7;
				m7 << "STATUS UPDATE AT [t, a, thr, j]\t " << t << " , "  << a << " , " << omp_get_thread_num() << " , " << j 
				//<< " ,with current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) 
				<< " and time taken for " << counter << " iterations: " << elapsed_min.count() << " min." <<  "\n"; cout << m7.str();
				//errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();

				//Report syst every 50000 iterations.
			

			}

			if(lo == -1)
			{
				stringstream m6_1;     //To make cout thread-safe as well as non-garbled
				m6_1 << "Saving error frames and exiting program .... \n"; cout << m6_1.str(); cerr << m6_1.str();

				// SAVING ERROR Rho_dt, DRho, and gamma vectors TO FILE
				stringstream L, tm ,d3, p1, rini, gm, a1, a2, Dm0, Dm1, Dm2, alph, w0t, aij, hij, dix, dimitri, sig0, sig1, veq;

				L << g; tm << t; d3 << setprecision(3) << dt; p1 << setprecision(4) << a; dix << setprecision(2) << dx;
				rini << j; Dm0 << D[0]*pow(10.0, 7.0); Dm1 << setprecision(3) << D[1]; Dm2 << setprecision(3) << D[2]; 
				a1 << a_st; a2  << a_end; sig0 << sigma[0]; sig1 << sigma[1]; dimitri  << dP; veq << setprecision(5) << Vstar;
				//gm << setprecision(3) << gmax; w0t << setprecision(3) << W0; alph << setprecision(3) << alpha; aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; 
				// Three replicates are over.

				string frame_parendir = "";

				// Creating a file instance called output to store output data as CSV.
				if(Vstar != -1)
					frame_parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str();
				else
					frame_parendir = frame_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();

				// SAVING RHO_DT FRAME
				string filenamePattern_Rhodt = "/ERROR_RHO_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
				+ "_a_" + p1.str()  + "_dx_"+ dix.str() + "_R_";
				//Next, designate the naive filename
				string filename = filenamePattern_Rhodt + rini.str() + ".csv";
				save_frame(Rho_dt, frame_parendir, filenamePattern_Rhodt, a, a_st, a_end, t, dt, dx, dP, j, g);

				// SAVING DRHO FRAME
				string filenamePattern_DRho = "/ERROR_DRHO_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
				+ "_a_" + p1.str()  + "_dx_"+ dix.str() + "_R_";
				//Next, designate the naive filename
				filename = filenamePattern_DRho + rini.str() + ".csv";
				save_frame(DRho, frame_parendir, filenamePattern_DRho, a, a_st, a_end, t, dt, dx, dP, j, g);

				/** // SAVING GAMMA FRAME
				string filenamePattern_Gamma = "/ERROR_GAMMA_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str()
				+ "_a_" + p1.str()  + "_dx_"+ dix.str() + "_R_";
				//Next, designate the naive filename
				filename = filenamePattern_Gamma + rini.str() + ".csv";
					
				string gammaheader = "a_c,  x,  Gam0(x; t), Gam1(x; t) \n"; //Header for output frame.
				save_frame(gamma, frame_parendir, filenamePattern_Gamma, a, a_st, a_end, t, dt, dx, dP, j, g, gammaheader);
				*/
				
				// SAVING PRELIMINARY FRAME
				// In this case, copy the first "index" rows of rho_rep_avg_var to a new 2D vector, update the values using var_mean_incremental_surv_runs()
				// and save to file.
				D2Vec_Double rho_rep_avg_var_temp(index, vector<double> (Sp4_1, 0.0)); 
				
				//Stores time, running avg, var (over replicates) of <rho(t)>x and number of surviving runs (at t) respectively.
				//std::copy(rho_rep_avg_var.begin(), rho_rep_avg_var.begin() + index, rho_rep_avg_var_temp.begin());
				var_mean_incremental_surv_runs(rho_rep_avg_var_temp, Rho_M, index, 0);

				string filenamePattern_Prelim = "/ERROR_SERIES_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
				"_D1_"+ Dm1.str() + "_R_";
				//Next, designate the naive filename

				string prelim_parendir= "";
				if(Vstar != -1)
					prelim_parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str() + "/TimeSeries";
				else
					prelim_parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "/TimeSeries";

				save_prelimframe(rho_rep_avg_var_temp, prelim_parendir, filenamePattern_Prelim, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, false, false);
				exit(3);
			}
		} // End of time loop.

		vector<vector<double>>().swap(gamma); vector<vector<double>>().swap(K1); 
		vector<vector<double>>().swap(K2); vector<vector<double>>().swap(K3); vector<vector<double>>().swap(K4);
		vector<vector <double>>().swap(Rho_dt); vector<vector <double>>().swap(DRho); vector<vector <double>>().swap(Rho_tsar);
		// Freeing up memory.
		if(int(omp_get_thread_num()) == 5)
		{
			size_t currentSize = getCurrentRSS( ); //Check Memory usage.
			size_t peakSize    = getPeakRSS( );
			stringstream m7;
			m7 << "In replicate " << j << ", current size in MB: " << currentSize/(1024.0*1024.0) 
			<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << "\n"; cout << m7.str();
			//errout.open(thr, std::ios_base::app); errout << m7.str(); errout.close();
		}

		if(j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "INITIAL UPDATE FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << tot_iter << "\n"; cout << m7.str();
			for(int i =0; i< tot_iter; i++)
			{	
				for(int s =0; s< SpB; s++)
				{
					rho_rep_avg_var[i][4*s+ 3] = Rho_M[i][2*s + 1]; // Avg density of frame at given time point i.
					rho_rep_avg_var[i][4*s+ 4] = 0.0; // Variance in avg frame density.
					rho_rep_avg_var[i][4*s+ 6] = Rho_M[i][2*s]; // Number of active sites in frame at given time point i.

					if(Rho_M[i][2*s] > 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s+5 ] = 1;  } //One surviving run as of now
					else if(Rho_M[i][2*s] == 0) //Ensuring only surviving runs are considered
					{ rho_rep_avg_var[i][4*s + 5] = 0;  } //No surviving run as of now
				}
				rho_rep_avg_var[i][1] = Rho_M[i][2*(Sp-2) + 1]; // Avg density of frame at given time point i.
				rho_rep_avg_var[i][2] = Rho_M[i][2*(Sp-1) + 1]; // Avg density of frame at given time point i.
				
			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{
			// Second or higher replicate, use incremental advances.
			var_mean_incremental_surv_runs(rho_rep_avg_var, Rho_M, tot_iter, j); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}

		// Freeing up memory.
		vector<vector <double>>().swap(Rho_M); 

		if((j+1)%2 == 1 || (j+1)%2 == 0)
		{	//Start logging at every multiple of 2
			stringstream L, tm ,d3, p1, a1, a2, dimitri, rini, rini_prev, Dm, cgm,  sig0, veq, jID;
			
			a1 << a_st; a2 << a_end;
  			L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(5) << a; dimitri << dP;
  			rini << j; Dm << setprecision(4) << D[2]; cgm << c*gmax; sig0 << sigma[0]; veq << setprecision(5) << Vstar;
			string parendir ="";
			if(Vstar != -1)
				parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str() + "_Veq_" + veq.str();
			else
				parendir = prelim_folder + a1.str() + "-" + a2.str() +  "_dP_" + dimitri.str();

			double ran_jid = (unif(rng)*1000.0)/1000.0; jID << ran_jid; // Random number between 0 and 1.
			
			string filenamePattern = prelim_prefix + jID.str() +"_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a_"+ p1.str() +
			"_D2_"+ Dm.str() + "_R_";

			save_prelimframe(rho_rep_avg_var, parendir, filenamePattern, a, a_st, a_end, t, dt, dx, dP, j, g, prelimheader, true);
		}

	} // End of r loop.

	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );
	stringstream m9;
	m9 << "SA RE GA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m9.str();
	//errout.open(thr, std::ios_base::app); errout << m9.str(); errout.close();

	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|    r 		|     L 		|    t 		|     <<W(t)>>x,r			|      <<O(t)>>x,r		|
		//..     <<Rho0(t)>>x,r			|    Var0[<Rho(t)>x],r    |    #Surviving Runs    |   #Active Sites |....

		Rho.push_back({a, static_cast<double>(r), static_cast<double>(g), rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], 
					rho_rep_avg_var[i][2]});
		for(int s=0;s < Sp-2; s++)
		{
			Rho[i].insert(Rho[i].end(), { rho_rep_avg_var[i][4*s +3], rho_rep_avg_var[i][4*s +4], 
			rho_rep_avg_var[i][4*s +5], rho_rep_avg_var[i][4*s +6]});		
		}
	}
	stringstream m10;
	m10 << "MA PA DHA, for Thread Rank:\t " << omp_get_thread_num() << ", current size in MB: " << currentSize/(1024.0*1024.0) 
	<< " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl; cout << m10.str();
	//errout.open(thr, std::ios_base::app); errout << m10.str(); errout.close();
}

void first_order_critical_exp_delta_stochastic_2Sp(int div, double t_max, double a_start, double a_end, double a_c,  double c, double gmax, double alpha,
	double rW, double W0,  double (&D)[Sp], double v[], double (&K)[3], double sigma[], double (&A)[SpB][SpB], double (&H)[SpB][SpB], double (&E)[SpB], double (&M)[SpB], double pR[], double ch[], double clo[],
	double dt, double dx, double dP, int r,  int g, double Gstar /* = -1.*/,  double Vstar /* = -1.*/ )
{
	//init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	//double mean[Sp] = {1.0, 0.05}; double sd[Sp] = {0.25, 0.0125};
	//init_randframe(Rho_0, Sp,  g*g); //Returns Rho_0 with a full initial frame filled with 0.2.

	//set_density_dependentmortality(M[2]); //Sets the density dependent mortality rate.

	vector<double> a_space = linspace(a_start, a_end, div);
	cout << "NOSTRA" <<endl;
  	// The pspace to iterate over.
    vector <double> t_measure = logarithm10_time_bins(t_max, dt);
	//vector <double> t_measure = linspace(40, 880, 22);
	// Computes and returns ln-distributed points from t= 10^{0} to log10(t_max) (the latter rounded down to 1 decimal place) 
  	// Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

	//vector <double> t_measure = linspace(40.0, 1000.0, 25); // Computes and returns linearly distributed points from t= 20 to 240 (the latter rounded down to 1 decimal place
	t_measure.insert(t_measure.begin(), 0.0);

	//Remove last element of t_measure
	//t_measure.pop_back();
	
	// Remove all elements of t_measure that lie between 750 and 1750.
	//t_measure.erase(std::remove_if(t_measure.begin(), t_measure.end(), [](double x){return (x > 600.0 && x < 1500.0);}), t_measure.end());


	//vector <double> t_linearwindow = linspace(600, 1500, 19); // Computes and returns linearly distributed points from t= 20 to 240 (the latter rounded down to 1 decimal place

	// Insert the linear window to t_measure vector after the element in t_measure that is just less than the first element in t_linearwindow.
	//auto it = std::upper_bound(t_measure.begin(), t_measure.end(), t_linearwindow[0]);
	//t_measure.insert(it, t_linearwindow.begin(), t_linearwindow.end());

	// Sort the t_measure vector in ascending order and remove duplicates (if any).
	sort( t_measure.begin(), t_measure.end() );

	t_measure.erase( unique( t_measure.begin(), t_measure.end() ), t_measure.end() );



  	cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  	for (int i=0; i< t_measure.size(); i++)
  	{
		cout << t_measure[i] << " ";
  	} 	cout << endl;

  	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
    const int microToSeconds = 1000000;   
    const double delay1 = 5 * microToSeconds;     //5 seconds

    
    cout<<"Delay 1 in progress... ("<<delay1/microToSeconds<<"s)"<<endl;
	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

  	size_t currentSize = getCurrentRSS( ); //Check Memory usage.
	size_t peakSize    = getPeakRSS( );

	cout << "On initialisation, current size in MB: " << currentSize/(1024.0*1024.0) << " and Peak Size (in MB): " << peakSize/(1024.0*1024.0) << endl;

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 32)
	{
		//omp_set_num_threads(32); //Limiting use on Chunk. Don't be greedy.
		cout << "WARNING: Number of Processors is: " << nProcessors << ". Don't be greedy." << endl;
	}

	//double perc = 0.015; double c_high[Sp] ={dP, p0j, p0m}; double c_low[Sp] ={p0i, p0j, p0m};

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

  	  
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {

		//init_randconstframe(Rh0, Sp,  g*g, perc, c_high, c_low); // Returns a frame with random speckles of high and low density.
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho1(t)>>x,r			|    Var[<Rho1(t)>x],r    |
		**/
		rietkerk_Dornic_2D_2Sp(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, rW, W0, D, v, K , sigma, a_start, a_end, a_c, A,H,E,M, pR, ch, clo, dt, dx, dP, r, g, Gstar, Vstar);
		//RK4_Wrapper_2D(CExpRho_a, t_measure, t_max, a_space[i], c, gmax, alpha, d, rW, W0, D, K , a_start, a_end, dt, dx, dP, r, g);
		//expanded_percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, a_space[i], b, c, D, sigma, dt, dx, r, g);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);
		//void rietkerk_Dornic_2D_MultiSp(D2Vec_Double &Rho, vector <double> &t_meas, double t_max, double a, double c, double gmax, double alpha, double rW, double W0,  double D[], double v[], double K[], double sigma[], double a_st, double a_end, double a_c, double A[SpB][SpB], double H[SpB][SpB], double E[], double M[], double pR[], chigh[],  clow[],  dt, dx, dP,  r, g)

        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
		vector<vector <double>>().swap(CExpRho_a);
		// Remove dynamic vectors from memory and free up allocated space.
      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting critical exponent data for each grid size in order.
          stringstream message3;
		  size_t currSize = getCurrentRSS( ); //Check Memory usage.
		  size_t peaksSize    = getPeakRSS( );
			 
          message3 << "Is this happening?\n" << "In thread: " << i << ", current size in MB: " << currSize/(1024.0*1024.0) 
		  << " and Peak Size (in MB): " << peaksSize/(1024.0*1024.0) << endl;
          cout << message3.str();

      }
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	int rd = std::random_device{}();
	std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
	rng.seed(rd);
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	//Store a uniformly distributed random number between 0 and 1 (rounded to 3 decimal places).
	double id = round(unif(rng)* 1000.0) / 1000.0; //Random number between 0 and 1.
	//Round to 3 decimal places.
	stringstream L, coco, tm ,d3, p1, p2, rini, Dm0, Dm1, aij, hij, cgm, alph, dix, dimitri, Sig0, veq, ID;

	//double Gstar = M[2]/((E[2] -M[2]*H[1][2])*A[1][2]); // MFT estimate of Grazer density at coexistence.

  	L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << a_start; p2  << a_end; rini << r; 
	alph << alpha; cgm << c*gmax; coco << setprecision(4) << c;
  	Dm0 << setprecision(3) << D[0]; Dm1 << setprecision(3) << D[1]; Sig0 << sigma[0]; dix << setprecision(2) << dx; 
	aij << setprecision(3) << A[0][1]; hij << setprecision(3) << H[0][1]; veq << setprecision(5) << Vstar;
	ID << id;

	ofstream output_1stdp;
  // Creating a file instance called output to store output data as CSV.
	output_1stdp.open(stat_prefix + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() + "_D1_"+ Dm1.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + 
	"_Veq_"+ veq.str() /* + "_ID_" + ID.str() */ + "_R_"+ rini.str() + ".csv");

	cout << "Save file name: " <<endl;
	cout << stat_prefix + L.str() 
	+ "_T_" + tm.str() + "_dt_" + d3.str() + "_D1_"+ Dm1.str() + "_S0_" + Sig0.str() +
	"_a1_"+ p1.str() + "_a2_"+ p2.str() + "_dx_"+ dix.str() + 
	"_cgmax_"+ cgm.str() + "_Veq_"+ veq.str() /* + "_ID_" + ID.str() */ + "_R_"+ rini.str() + ".csv" <<endl;

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

	string header = prelimheader;
	output_1stdp << header;
	cout << "The vector elements are: "<< endl;
  	cout << header << "\n";

	for(int i=0; i< vec.size(); i++)
	{
		for(int j=0; j< vec[i].size() -1; j++)
			output_1stdp << setprecision(16) << vec[i][j] << ",";
		output_1stdp << setprecision(16) << vec[i][vec[i].size() -1] <<  endl;
		if( i%(1000) ==1)
    	{
			for(int j=0; j< vec[i].size() -1; j++)
				cout << setprecision(16) << vec[i][j] << ",";
			cout << setprecision(16) << vec[i][vec[i].size() -1] <<  endl;
    	}
	}
	output_1stdp.close();
}

