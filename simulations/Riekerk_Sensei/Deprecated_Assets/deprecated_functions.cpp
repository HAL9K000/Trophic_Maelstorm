#include "rietkerk_bjork_basic.h"
#include "Debug.h"


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
