#include "SPDP.h"

//--------------------Defining Primitives------------------------------------//

void increase_stack_limit(long stack_size){

	//Source: https://stackoverflow.com/questions/2279052/increase-stack-size-in-linux-with-setrlimit
	//Credit: https://stackoverflow.com/users/253056/paul-r
	//Used with modification

	// Necessary because of complexity of integrations.

	const rlim_t kStackSize = stack_size * 1024L * 1024L;   // min stack size = 128 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    }
}

void init_fullframe(double array[], int size)
{
	//Returns array with all elements initialised to 1
	for(int i=0; i< size; i++)
	{	array[i] =1 ;}
}

void init_solitarytear(double array[], int size)
{
	//Returns array with only central element initalised to 1, all others 0.
	for(int i=0; i< size; i++)
	{	
		array[i] =0;
		if ( i == static_cast<int>(size/2)) 
		{	array[i] =1; }
	}
}

double mean_of_array(double array[],int size){

	double sum = 0.0;

	for (int i =0; i<size; ++i){
		sum += array[i];
	}
	return sum/(double)size;
}

double standard_deviation_of_array(double array[],int size){

	double mean = mean_of_array(array,size);
	double sum = 0.0;

	for (int i = 0; i < size; ++i){
		sum += (array[i]-mean)*(array[i]-mean);
	}

	double variance = sum/(double)size;

	return sqrt(variance);
}

double mean_of_vector(vector<double> array,int size){

	double sum = 0.0;

	for (int i =0; i<size; ++i){
		sum += array[i];
	}
	return sum/(double)size;
}

double standard_deviation_of_vector(vector<double> array,int size){

  double mean = mean_of_vector(array,size);
	float sum = 0.0;

	for (int i =0; i<size; ++i){
		sum += (array[i]-mean)*(array[i]-mean);
	}
  double variance = sum/(double)size;
	return sqrt(variance);
}

void var_mean_incremental(double rep_avg_var[][3], vector<vector <double>> &X_curr, int size, int r)
{
	/** X_curr contains time-series data (indexed as (t, data)) from the current replicate (given by r) with a length "size".
	 *  rep_avg_var stores the running average and variance of X_curr for each time-step from previous replicates ( <= r-1).
	 * In this function, X_curr is used to update the running avg and variance of rep_avg_var to account for current replicate r.
	*/
	double mean_prev[size];
	for(int i=0; i<size; i++)
	{
		rep_avg_var[i][0] = X_curr[i][0]; //Store time points.
		mean_prev[i] = rep_avg_var[i][1]; //Stores X_n-1 mean data as temp variable.
		//Next calculate and store incremental mean (X_n).
		rep_avg_var[i][1] = (X_curr[i][1] + (r-1)*rep_avg_var[i][1])/r; //Note r>= 2.
		//Next calculate and store incremental variance (V^2_n) [For formula derivation, refer to notes.]
		rep_avg_var[i][2] = (rep_avg_var[i][1] - mean_prev[i])*(rep_avg_var[i][1] - mean_prev[i]) +
			(1/(r-1))*( (r-2)*rep_avg_var[i][2] + (X_curr[i][1] - rep_avg_var[i][1])*(X_curr[i][1] - rep_avg_var[i][1]));
	}

}

template<typename T>std::vector<double> linspace(T start_in, T end_in, int num_in)
{
	//Equivalent of numpy's linspace method.

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1)
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end);
  return linspaced;
}

std::vector<double> logarithmic_time_bins(double t_max, double dt)
{	// Computes and returns ln-distributed points from t= e^2 to ln(t_max) (the latter rounded down to 1 decimal place)
	std::vector<double> logspaced;
	if(t_max <= exp(2.0))
	{
		cout << "GOUT!!"<<endl;
		return logspaced; // No points recorded if t_max < e^{2}
	}

	double tend_one = static_cast<double>(floor(log(t_max)* 10.)) / 10.0;
	// Finds the closest power of e (rounded down to one decimal place) from t_max.
	double pow =2.0;
	while(pow <= tend_one)
	{
		double t_raw = exp(pow); // Getting raw time-stamp.
		//Rounding down raw time-stamp to dt.
		double t_fin = static_cast<int>(t_raw/dt)*dt;

		if(t_fin < exp(2))
			t_fin = 7.4; //Ensure t_fin is strictly above e^{2}

		if(logspaced.size() > 0  && t_fin != logspaced[logspaced.size() -1])
		{ 	logspaced.push_back(t_fin);	} // This is to prevent duplicate values from arising in logspaced.
		else if(logspaced.size() == 0)
			logspaced.push_back(t_fin);
		pow += 0.04; //25 time-stamps per decade.
	}
	cout << "TROUT! " << logspaced.size() <<endl;
	return logspaced;


}

//----------------------------- Stochastic SPDP RK4 Integration Machinery --------------------------------------------//

void determine_neighbours_R2(int neighbours_R2[][2][5], int g)
{
	//Stores neighbouring sites (Van-Neumann Radius of 2) of each cell site in a 
	// g^2*2*5 matrix.
	int a=0; int b= 0; //Useful for accounting purposes
	for (int i=0; i<g*g; i++)
	{
		int p = int(i/g); int q = i%g;
		//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.
		for(int j=0; j < 5; j++)
		{
			int k= j-2;
			a= ((p+k)%g + g)%g; // Eq to (p+k)mod g. -2 <= k <= 2
			b= ((q+k)%g + g)%g; // Eq to (q+k)mod g. -2 <= k <= 2

			//Recall, first row stores vertical neighbours, second row stores horizontal neighbours.
			neighbours_R2[i][0][j] = a*g + q; neighbours_R2[i][1][j] = p*g + b;
		}
	}

}

void determine_neighbours_S8(int neighbours_R2[][3][3], int g)
{
	//Stores neighbouring sites in a square (Norm 1) radius of 1
	int a=0; int b= 0; //Useful for accounting purposes
	for (int i=0; i<g*g; i++)
	{
		int p = int(i/g); int q = i%g;
		//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.
		for(int y=0; y < 3; y++)
		{	int ky= y-1; a= ((p+ky)%g + g)%g; // Eq to (p+ky)mod g. -1 <= ky <= 1
			for(int x=0; x < 3; x++)
			{
				int kx = x-1;
				b= ((q+kx)%g + g)%g; // Eq to (q+kx)mod g. -1 <= kx <= 1
				//Recall, first row stores vertical neighbours, second row stores horizontal neighbours.
				neighbours_R2[i][y][x] = a*g + b; 
			}	
		}
	}

}

int twilight_zoneS8(vector<double> &Rho_t, int neighbours_S8[][3][3], int i)
{
	//Checks if Rho_dt[i] is 0 and is surrounded by empty patches (radius of 1). Returns 0 if so, else 1.

	for(int x =1; x<3; x++)
	{
		for(int y=1; y<3; y++) //Recall, S8 resembles Laplacian schematic.
		{		
			if(Rho_t[neighbours_S8[i][x][y]] != 0.0) 
				return 1;
		}
	}
	return 0; //Only occurs all neighbours in a radius of 1 are empty.

}

int twilight_zoneR2(vector<double> &Rho_t, int neighbours_R2[][2][5], int i)
{
	//Checks if Rho_dt[i] is 0 and is surrounded by empty patches (radius of 2). Returns 0 if so, else 1.

	for(int j =0; j<5; j++)
	{
		//Recall, first row stores horizontal neighbours, second row stores vertical neighbours.
		if(Rho_t[neighbours_R2[i][0][j]] != 0)
			return 1;
		else if(Rho_t[neighbours_R2[i][1][j]] != 0)
			return 1;
	}
	return 0; //Only occurs all neighbours in a radius of 2 are empty.

}



void f_2D(double f[], vector<double> &Rho_t, int n[][2][5], double p, double D, double t, double dt, double dx, int g)
{
  //Vector function that evaluates drho/dt = f(p,t) at given p and t.

	for(int i=0; i < g*g; i++)
	{
		f[i] = -Rho_t[i]*Rho_t[i] - D*(1/(12.0*dx*dx))*(1*(Rho_t[n[i][0][0]] + Rho_t[n[i][0][4]] 
		+ Rho_t[n[i][1][0]] + Rho_t[n[i][1][4]]) -16*(Rho_t[n[i][0][1]] + Rho_t[n[i][0][3]] 
		+ Rho_t[n[i][1][1]] + Rho_t[n[i][1][3]]) +60*Rho_t[n[i][0][2]]);

		// RK4 update of Non-linear parts plus diffusion term
		// 4th order Laplacian from 2SHOC.
	}
}

void fDo_2D(double f[], vector<double> &Rho_t, int n[][2][5], double b, double t, double dt, int g)
{
  //RK4 Integration For Old Dornic Methode, where, as you may recall the Diffusion term is treated analytically.

	for(int i=0; i < g*g; i++)
	{
		f[i] = -b*Rho_t[i]*Rho_t[i];
		// RK4 update of Non-linear parts
	}
}

void eulerf_2D(double f[], vector<double> &Rho_t, int n[][2][5], double p, double D, double t, double dt, double dx, int g)
{
  //Performs Euler update of non-linear and diffusion terms and stores it in f
	for(int i=0; i < g*g; i++)
	{
		f[i] = Rho_t[i] + (-Rho_t[i]*Rho_t[i] - D*(1/(12.0*dx*dx))*(1*(Rho_t[n[i][0][0]] + Rho_t[n[i][0][4]] 
		+ Rho_t[n[i][1][0]] + Rho_t[n[i][1][4]]) -16*(Rho_t[n[i][0][1]] + Rho_t[n[i][0][3]] 
		+ Rho_t[n[i][1][1]] + Rho_t[n[i][1][3]]) +60*Rho_t[n[i][0][2]]))*dt;

		// RK4 update of Non-linear parts plus diffusion term
		// 4th order Laplacian from 2SHOC.
	}
}

void fP_2D(double f[], vector<double> &Rho_t, int n[][2][5], double a, double b, double D, double t, double dt, double dx, int g)
{
  //Vector function that evaluates drho/dt = f(p,t) at given p and t.

	for(int i=0; i < g*g; i++)
	{
		f[i] = a*Rho_t[i] -b*Rho_t[i]*Rho_t[i] - D*(1/(12.0*dx*dx))*(1*(Rho_t[n[i][0][0]] + Rho_t[n[i][0][4]] 
		+ Rho_t[n[i][1][0]] + Rho_t[n[i][1][4]]) -16*(Rho_t[n[i][0][1]] + Rho_t[n[i][0][3]] 
		+ Rho_t[n[i][1][1]] + Rho_t[n[i][1][3]]) +60*Rho_t[n[i][0][2]]);

		// RK4 update of Non-linear parts plus diffusion term
		// 4th order Laplacian from 2SHOC.
	}
}


void find_neighbours_R2(double neighbours_R2[2][5], vector<double> &Rho_t, int p, int q, int g)
{
	//Stores neighbour density in a ball of radius 2 around p,q.
	int a=0; int b= 0; //Useful for accounting purposes

	//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.
	for(int j=0; j < 5; j++)
	{
			int k= j-2;
			a= ((p+k)%g + g)%g; // Eq to (p+k)mod g. -2 <= k <= 2
			b= ((q+k)%g + g)%g; // Eq to (q+k)mod g. -2 <= k <= 2

			//Recall, first row stores horizontal neighbours, second row stores vertical neighbours.
			neighbours_R2[0][j] = Rho_t[a*g + q]; neighbours_R2[1][j] = Rho_t[p*g + b];
	}
}

void percNoDiff_2D(vector<vector<double>> &Rho, vector<double> &t_stop, double Rh0[],
	 double t_max, double a, double b, double D, double dt, double dx, int r,  int g)
{
	int n_R2[g*g][2][5] ={0};
	//n_R2 stores neighbours of all sites in a ball of radius 2.
	//Assigning neigbours.
	determine_neighbours_R2(n_R2, g);
	
	for(int j=0; j< 1; j++)
	{
		// Iterating over replicates.
		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		thread_local std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd);

		vector <double> Rho_dt(g*g, 0.0); //Kept constant, only updated at end of turn
		double K1[g*g] ={0.0}; double K2[g*g] ={0.0}; double K3[g*g] ={0.0}; double K4[g*g] ={0.0};	
		//Initialise RK4 terms for future use.
		vector <double> Rho_M(g*g, 0.0); //Stores intermediate values of Rho at dt/2 and dt for RK4 calculations.
		cout<< "Initially: "<<endl;
	  for(int i=0; i< g*g; i++)
	  {
				Rho_dt[i] = Rh0[i]; //Assigning initial conditions
				cout<< Rho_dt[i] << " ";
	  } cout<<endl;

	   	stringstream m2;     //To make cout thread-safe as well as non-garbled due to race conditions.
     	m2 << "POCKET!" <<endl; cout << m2.str(); 

		//int t_max = static_cast<int>(t_stop[t_stop.size() -1]); //Storing last element of the same here.
		double t=0; //Initialise t
		int s=0;
		
		while (t < t_max + dt)
		{
			//End at the end of all time.
			//cout << "SOCKET" << endl;

			stringstream m3;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m3 << "B4RK4" <<endl; cout << m3.str();

			fP_2D(K1, Rho_dt, n_R2, a, b, D, t, dt, dx, g); //K1 updated. 
			cout<< "K1" << endl;
			for(int i = 0; i<g*g; dx, i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K1[i];  
			  //cout << K1[i] << " ";
			} //Recall Runge Kutta-Update Schema for K2 = f(t+dt/2, y+dt/2*K1)
			cout << endl;
			fP_2D(K2, Rho_M, n_R2, a, b, D, t +dt/2.0, dt, dx, g); //K2 updated.
			cout<< "K2" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K2[i]; cout << K2[i] << " ";  } //Recall Runge Kutta-Update Schema for K3
		  	cout << endl;
			fP_2D(K3, Rho_M, n_R2, a, b, D, t +dt/2.0, dt, dx, g); //K3 updated.
			cout<< "K3" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt)*K3[i]; cout << K3[i] << " ";  } //Recall Runge Kutta-Update Schema for K4
			cout << endl;
			fP_2D(K4, Rho_M, n_R2, a, b, D, t + dt, dt, dx, g); //K4 updated.
			cout<< "K4" << endl;
			for(int i = 0; i<g*g; i++)
			{
				cout << K4[i] << " ";
			}
			cout << endl;
			 stringstream m5;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m5 << "A7RK4" <<endl; cout << m5.str();
			//Now that K1, K2, K3, K4 updated for all sites, do the stochastic update.
			for(int i = 0; i<g*g; i++)
			{
			Rho_dt[i]+= (dt/6.0)*( K1[i] + 2*K2[i] + 2*K3[i] + K4[i]);
			}
			//Impproved Dornic  Integration of Linear & Stochastic Term
			//CREDITS: PAULA VILLA MARTIN.
			t+=dt; //Update timestep
			cout << "MOULIN ROUGE Final Rho(t):" <<endl;
			for(int i=0;i<Rho_dt.size();i++)
			{ cout<< Rho_dt[i] << " "; } cout<<endl;

			

			stringstream m4;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m4 << "SOCKET!" <<endl; cout << m4.str();
			

			//Updating Pijm at end of sweep.

			//cout << t << endl;
			if(t > t_stop[s] -dt/2.0 && t <= t_stop[s] + dt/2.0)
			{
				stringstream m4;     //To make cout thread-safe as well as non-garbled due to race conditions.
      			m4 << "LOCKET!:\t" << t << " " << a  <<endl; cout << m4.str();
				for ( int k = 0; k< Rho_dt.size(); k++)
			  {	// Recall Rho = | p (Order Parameter) | Replicate (r) | Time(t)  | a*L + b | Cell Density[a][b] |
					Rho.push_back({a, static_cast<double>(j), t, static_cast<double>(k), Rho_dt[k]});
			  }
				s = (s+1)%t_stop.size(); // To prevent out of bounds exception.
			}

		} // End of t loops



		vector<double>().swap(Rho_M); vector<double>().swap(Rho_dt);
		// Remove dynamic vectors from memory and free up allocated space.

		

	}



}


void percolation_2D(vector<vector<double>> &Rho, vector<double> &t_stop, double Rh0[],
	 double t_max, double p, double D, double sigma, double dt, double dx, int r,  int g)
{
	int n_roll =10000; int n_star = 1000; int n_R2[g*g][2][5] ={0};
	//n_R2 stores neighbours of all sites in a ball of radius 2.
	//Assigning neigbours.
	determine_neighbours_R2(n_R2, g);

	//For debugging.
	int i1 =0 ; int i2 = int(g*g/2); int i3 = g -1;
	vector <vector <double>> KRi1; vector <vector <double>> KRi2;
	vector <vector <double>> KRi3; double diff[2][g*g]; double Rhog[g*g];
	//Stores K1, K2, K3, K4 plus Rho_dt and Diffusion details for each index.

	int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j <<endl; cout << m1.str();
		
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);
		//thread_local std::default_random_engine rng(rd); //engine
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number


		/** Recall in Mean Field Limit (MFT), DP is given by:
		d(rho)/dt = (2p -1)rho - 1*rho^(2) **/
		//Defining the Dornic variables
		double b = 2*p -1.0;
		double lambda_exp = exp( (b)*dt); double lambda = 2*(b)/(sigma*sigma*(lambda_exp -1.0));

		vector <double> Rho_dt(g*g, 0.0); //Kept constant, only updated at end of turn
		double K1[g*g] ={0.0}; double K2[g*g] ={0.0}; double K3[g*g] ={0.0}; double K4[g*g] ={0.0};	
		//Initialise RK4 terms for future use.
		vector <double> Rho_M(g*g, 0.0); //Stores intermediate values of Rho at dt/2 and dt for RK4 calculations.
		cout<< "Initially: "<<endl;
	  for(int i=0; i< g*g; i++)
	  {
				Rho_dt[i] = Rh0[i]; //Assigning initial conditions
				cout<< Rho_dt[i] << " ";
	  } cout<<endl;

	   	stringstream m2;     //To make cout thread-safe as well as non-garbled due to race conditions.
     	m2 << "POCKET!" <<endl; cout << m2.str(); 
		double t=0; //Initialise t
		int s=0;
		
		while (t < t_max + dt)
		{
			//End at the end of all time.
			//cout << "SOCKET" << endl;

			stringstream m3;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m3 << "B4RK4" <<endl; cout << m3.str();

			f_2D(K1, Rho_dt, n_R2, p, D, t, dt, dx, g); //K1 updated. 
			cout<< "K1" << endl;
			for(int i = 0; i<g*g; dx, i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K1[i];  
			  cout << K1[i] << " ";
			} //Recall Runge Kutta-Update Schema for K2 = f(t+dt/2, y+dt/2*K1)
			cout << endl;
			f_2D(K2, Rho_M, n_R2, p, D, t +dt/2.0, dt, dx, g); //K2 updated.
			cout<< "K2" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K2[i]; cout << K2[i] << " ";  } //Recall Runge Kutta-Update Schema for K3
		  	cout << endl;
			f_2D(K3, Rho_M, n_R2, p, D, t +dt/2.0, dt, dx, g); //K3 updated.
			cout<< "K3" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt)*K3[i]; cout << K3[i] << " ";  } //Recall Runge Kutta-Update Schema for K4
			cout << endl;
			f_2D(K4, Rho_M, n_R2, p, D, t + dt, dt, dx, g); //K4 updated.
			cout<< "K4" << endl;
			/**for(int i = 0; i<g*g; i++)
			{
				cout << K4[i] << " ";
			}*/
			cout << endl;

			//Debugging section.
			for(int i = 0; i<g*g; i++)
			{
				diff[0][i] = - D*(1/(12.0*dx*dx))*(1*(Rho_dt[n_R2[i][0][0]] + Rho_dt[n_R2[i][0][4]] 
				+ Rho_dt[n_R2[i][1][0]] + Rho_dt[n_R2[i][1][4]]) -16*(Rho_dt[n_R2[i][0][1]] + Rho_dt[n_R2[i][0][3]] 
				+ Rho_dt[n_R2[i][1][1]] + Rho_dt[n_R2[i][1][3]]) +60*Rho_dt[n_R2[i][0][2]]);
				diff[1][i] = - D*(1/(12.0*dx*dx))*(1*(Rho_M[n_R2[i][0][0]] + Rho_M[n_R2[i][0][4]] 
				+ Rho_M[n_R2[i][1][0]] + Rho_M[n_R2[i][1][4]]) -16*(Rho_M[n_R2[i][0][1]] + Rho_M[n_R2[i][0][3]] 
				+ Rho_M[n_R2[i][1][1]] + Rho_M[n_R2[i][1][3]]) +60*Rho_M[n_R2[i][0][2]]);
			}

			//eulerf_2D(K1,Rho_dt, n_R2, p, D, t, dt, dx, g); //Euler update, stored in K1

			 stringstream m5;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m5 << "A7RK4 Rho at t" << t <<endl; cout << m5.str();
			//Now that K1, K2, K3, K4 updated for all sites, do the stochastic update.

			poisson_distribution<int> poisson; gamma_distribution <double> gamma;

			//Impproved Dornic  Integration of Linear & Stochastic Term
			//CREDITS: PAULA VILLA MARTIN.
			

			for(int i=0;i<Rho_dt.size();i++)
			{
	        Rho_dt[i]+= (dt/6.0)*( K1[i] + 2*K2[i] + 2*K3[i] + K4[i]);
			// Deterministic RK4 complete.
			//Rho_dt[i] = K1[i];
			Rhog[i] = Rho_dt[i]; //Stores Rho_dt RK4 update values (which serves as initial conditions for non-linear terms).
			double mu = lambda*lambda_exp*Rho_dt[i];
			cout << Rho_dt[i] << " ";

	        poisson = poisson_distribution<int>(lambda*lambda_exp*Rho_dt[i]);
			//cout << "Poisson Pill:\t" << poisson(rng) <<endl;
	        gamma = gamma_distribution<double>(poisson(rng), 1.0);
					//Defining shapes of poisson and gamma distributions.
	        Rho_dt[i]= gamma(rng)/lambda;

	    	} cout<<endl;
			cout << "MOULIN ROUGE Final Rho(t) at t:" << t <<endl;
			for(int i=0;i<Rho_dt.size();i++)
			{ cout<< Rho_dt[i] << " "; } cout<<endl;

			

			stringstream m4;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m4 << "SOCKET!" <<endl; cout << m4.str();
			
			t+=dt; //Update timestep

			//Updating Pijm at end of sweep.

			//cout << t << endl;
			if(t > t_stop[s] -dt/2.0 && t <= t_stop[s] + dt/2.0)
			{
				stringstream m4;     //To make cout thread-safe as well as non-garbled due to race conditions.
      			m4 << "LOCKET!:\t" << t << " " << p  <<endl; cout << m4.str();
				for ( int k = 0; k< Rho_dt.size(); k++)
			  {	// Recall Rho = | p (Order Parameter) | Replicate (r) | Time(t)  | a*L + b | Cell Density[a][b] |
					Rho.push_back({p, static_cast<double>(j), t, static_cast<double>(k), Rho_dt[k]});
			  }
				s = (s+1)%t_stop.size(); // To prevent out of bounds exception.
			}

			KRi1.push_back({static_cast<double>(j), t, K1[i1], K2[i1], K3[i1], K4[i1], 
			static_cast<double>(Rho_dt[i1]-Rhog[i1]), Rho_dt[i1], diff[0][i1], diff[1][i1]});

			KRi2.push_back({static_cast<double>(j), t, K1[i2], K2[i2], K3[i2], K4[i2], 
			static_cast<double>(Rho_dt[i2]-Rhog[i2]), Rho_dt[i2], diff[0][i2], diff[1][i2]});
			
			KRi3.push_back({static_cast<double>(j), t, K1[i3], K2[i3], K3[i3], K4[i3], 
			static_cast<double>(Rho_dt[i3]-Rhog[i3]), Rho_dt[i3], diff[0][i3], diff[1][i3]});

		} // End of t loops



		vector<double>().swap(Rho_M); vector<double>().swap(Rho_dt);
		// Remove dynamic vectors from memory and free up allocated space.

		

	}

	stringstream L, dT ,p1, dX, rini, Dm, Sm;

  	L << g; dT << dt; p1 << setprecision(4) << p; dX << setprecision(4) << dx;
  	rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma;
  	// setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_dp1; ofstream output_dp2; ofstream output_dp3;
  	// Creating a file instance called output to store output data as CSV.
	output_dp1.open("New/cDebug_Modulei1_G_" + L.str() + "_dt_" + dT.str() + "_dx_"+ dX.str() +
	"_p_"+ p1.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");
	output_dp2.open("New/cDebug_Modulei2_G_" + L.str() + "_dt_" + dT.str() + "_dx_"+ dX.str() +
	"_p_"+ p1.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");
	output_dp3.open("New/cDebug_Modulei3_G_" + L.str() + "_dt_" + dT.str() + "_dx_"+ dX.str() +
	"_p_"+ p1.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	// Output = | Replicate (r) | Time(t)  | K1 | ... | K4 | Stoc Rho_dt | Stoc - Dem Rho_dt |
	output_dp1 << " # Tr No , t ,  K1 , K2 , K3 , K4 , Stoc - Dem Rho_ij(t), Stoc Rho_ij(t), Diff_K1, Diff_K2, \n";
	output_dp2 << " # Tr No , t ,  K1 , K2 , K3 , K4 , Stoc - Dem Rho_ij(t), Stoc Rho_ij(t), Diff_K1, Diff_K2, \n";
	output_dp3 << " # Tr No , t ,  K1 , K2 , K3 , K4 , Stoc - Dem Rho_ij(t), Stoc Rho_ij(t), Diff_K1, Diff_K2, \n";
	cout << "The vector elements are: "<< endl;
  	cout << " p , # Tr No , t ,  i*L + j ,  Rho_ij(t) \n";

	for(int i=0; i< KRi1.size(); i++)
	{
		output_dp1 << static_cast<int> (KRi1[i][0])  << "," << setprecision(6) << KRi1[i][1] << "," << setprecision(9)
		<< KRi1[i][2] << "," << setprecision(9) << KRi1[i][3] << "," << setprecision(9) << KRi1[i][4] << "," 
		<< setprecision(9) << KRi1[i][5] << ","<< setprecision(16) << KRi1[i][6] << "," << setprecision(16) << KRi1[i][7] 
		<< "," << setprecision(16) << KRi1[i][8] << ","<< setprecision(16) << KRi1[i][9] << endl;

		output_dp2 << static_cast<int> (KRi2[i][0])  << "," << setprecision(6) << KRi2[i][1] << "," << setprecision(9)
		<< KRi2[i][2] << "," << setprecision(9) << KRi2[i][3] << "," << setprecision(9) << KRi2[i][4] << "," 
		<< setprecision(9) << KRi2[i][5] << ","<< setprecision(16) << KRi2[i][6] << ","<< setprecision(16) << KRi2[i][7] 
		<< ","  << setprecision(16) << KRi2[i][8] << ","<< setprecision(16) << KRi2[i][9] << endl;

		output_dp3 << static_cast<int> (KRi3[i][0])  << "," << setprecision(6) << KRi3[i][1] << "," << setprecision(9)
		<< KRi3[i][2] << "," << setprecision(9) << KRi3[i][3] << "," << setprecision(9) << KRi3[i][4] << "," 
		<< setprecision(9) << KRi3[i][5] << ","<< setprecision(16) << KRi3[i][6] << ","<< setprecision(16) << KRi3[i][7] 
		<< "," << setprecision(16) << KRi3[i][8] << "," << setprecision(16) << KRi3[i][9] <<  endl;

	}
	output_dp1.close(); output_dp2.close(); output_dp3.close();



}

void test_gamma_poisson()
{
	int nrolls=10000;  // number of experiments
  	int nstars=600;    // maximum number of stars to distribute
	for(int j=0; j< 4; j++)
	{
		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j <<endl; cout << m1.str();
		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		thread_local std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd);

		poisson_distribution<int> poisson(3.7); gamma_distribution <double> gamma(2.0,2.0);

		int g[24]={};

  		for (int i=0; i<nrolls; ++i) 
		{
    	double number = gamma(rng);
    	if (number<12) ++g[int(2*number)];
  		}

  		std::cout << "gamma_distribution (2.0,2.0):" << std::endl;

  		for (int i=0; i<24; ++i) 
		{
			double vi = i/2.0;
			if(vi < 9.5)
			{
				std::cout << setprecision(3) << vi << "-" << setprecision(3) << (vi+ 0.5) << ": ";
			}
			else if(vi == 9.5)
			{
				std::cout << setprecision(3)  << vi << "-" << setprecision(2) << (vi+ 0.5) << ": ";
			}
			else
			{ std::cout << setprecision(3)  << vi << "-" << setprecision(3) << (vi+ 0.5) << ": "; }
    		
    		std::cout << std::string(g[i]*nstars/nrolls,'*') << std::endl;
  		}
		

		int p[12]={};
		for (int i=0; i<nrolls; ++i) 
		{
    		int number = poisson(rng);
    		if (number<12) ++p[number];
  		}

  		std::cout << "poisson_distribution (mean=3.75):" << std::endl;
  		for (int i=0; i<12; ++i)
    		std::cout << i << ": " << std::string(p[i]*nstars/nrolls,'*') << std::endl;
		
	}
}


//----------------------------- Stochastic SPDP Dornic Integration Machinery --------------------------------------------//



void percolationDornic_2D(vector<vector<double>> &Rho, vector <double> &t_meas, double Rh0[],
	 double t_max, double a, double b, double D, double sigma, double dt, double dx, int r,  int g)
{
	//int nR2[g*g][2][5] ={0}; //n_R2 stores neighbours of all sites in a ball of radius 2.
	int nR2[g*g][3][3] ={0}; //nR2 stores neighbours of all sites in a ball of NORM 1 (SQUARE) radius 1.
	//determine_neighbours_R2(nR2, g); //Assigning neigbours.
	determine_neighbours_S8(nR2, g); //Assigning neigbours.

	//int tot_iter = static_cast<int>(t_max/dt); // Stores total number of iterations (i.e. time steps logged) per replicate.
	int tot_iter = static_cast<int>(exp(2.0)/dt) + 1 + t_meas.size(); // Stores total number of iterations (i.e. time steps logged) per replicate.
	double rho_rep_avg_var[tot_iter][3] ={0.0}; //Stores time, running avg and var (over replicates) of <rho(t)>x respectively.
	double diff_coefficient = D*(1/(6.0*dx*dx)); //Saves us some trouble later.
	double beta = a - (10.0/3.0)*D/(dx*dx); //double beta = a - 4*D/(dx*dx);
	double lambda_exp = exp( (beta)*dt); double lambda = 2*(beta)/(sigma*sigma*(lambda_exp -1.0));

	//stringstream m5;     //To make cout thread-safe as well as non-garbled due to race conditions.
    //m5 << "On Thread:\t" << omp_get_thread_num() << "\t with Total Iterations:\t " << tot_iter <<endl; cout << m5.str();
	
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
		<< "\t with Total Iterations:\t " << tot_iter  <<endl; cout << m1.str();
		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd+j);
		//thread_local std::default_random_engine rng(rd); //engine
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number

		/** Recall in Mean Field Limit (MFT), DP is given by:
		d(rho)/dt = (2p -1)rho - 1*rho^(2) **/
		//Defining the Dornic variables
		
		//double K1[g*g] ={0.0}; double K2[g*g] ={0.0}; double K3[g*g] ={0.0}; double K4[g*g] ={0.0};	
		//Initialise RK4 terms for future use.
		vector <double> Rho_dt(g*g, 0.0); //Updated through the turn, in turn used to update DRho
		vector <double> DRho(g*g, 0.0); //Dummy variable, Kept constant, only updated at end of turn.
		vector <vector <double>> Rho_M; //Stores <rho>x values per time step for a single replicate.

		//cout<< "Initially: "<<endl;
	  	for(int i=0; i< g*g; i++)
	  	{
				Rho_dt[i] = Rh0[i]; //Assigning initial conditions
				DRho[i] = Rh0[i];
	  	} //cout<<endl;

		

		//int t_max = static_cast<int>(t_stop[t_stop.size() -1]); //Storing last element of the same here.
		double t=0; //Initialise t
		int s=1; int po=1; int so =1; int index = 0; 
		int co=1; int eo=1;

		//stringstream m6;     //To make cout thread-safe as well as non-garbled due to race conditions.
        //m6 << "INITIAL CONDITIONS ASSIGNED FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() <<endl; cout << m6.str();
		while (t < t_max + dt)
		{
			double rhox_avg;
			//Start by updating the Rho vector.
			if(t < exp(2.0) || t >= t_meas[index] -dt/2.0 && t < t_meas[index] +dt/2.0)
			{
				rhox_avg = mean_of_vector(Rho_dt, g*g); //Finds spatial average of densities at given t.
				Rho_M.push_back({t, rhox_avg});
				if( t >= exp(2.0))
					index+=1;

			
			}

			/**if(t == 7.4)
			{
				stringstream m10;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m10 << "FML............." << endl; cout << m10.str();
			}
			else if(t >= 7.4 -dt/2.0 && t < 7.4 +dt/2.0)
			{
				stringstream m10;     //To make cout thread-safe as well as non-garbled due to race conditions.
        		m10 << "SUX............." << setprecision(16) << t <<  endl; cout << m10.str();
			}**/

			
			
			// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |

			
			poisson_distribution<int> poisson; gamma_distribution <double> gamma;

			//Basic Dornic  Integration of Linear & Stochastic Term
			//CREDITS: DORNIC.
			for(int i=0;i<Rho_dt.size();i++)
			{
			int check= twilight_zoneS8(DRho, nR2, i);
			//Checks if Rho_dt[i] is 0 and is surrounded by empty patches (Norm 1 radius of 1). Returns 0 if so, else 1.
			if(check == 0) 
			{	continue; } //Rho_dt[i] will remain 0 as neighbours are also at 0.

			/** double alpha_i = D/(dx*dx)*(1*(Rho_dt[nR2[i][0][1]] + Rho_dt[nR2[i][0][3]] +
			Rho_dt[nR2[i][1][1]] + Rho_dt[nR2[i][1][3]]));*/
			/**double alpha_i = - diff_coefficient*(1*(Rho_dt[nR2[i][0][0]] + Rho_dt[nR2[i][0][4]] 
				+ Rho_dt[nR2[i][1][0]] + Rho_dt[nR2[i][1][4]]) -16*(Rho_dt[nR2[i][0][1]] + Rho_dt[nR2[i][0][3]] 
				+ Rho_dt[nR2[i][1][1]] + Rho_dt[nR2[i][1][3]])); **/
			/**double alpha_i = - diff_coefficient*(1*(DRho[nR2[i][0][0]] + DRho[nR2[i][0][4]] 
				+ DRho[nR2[i][1][0]] + DRho[nR2[i][1][4]]) -16*(DRho[nR2[i][0][1]] + DRho[nR2[i][0][3]] 
				+ DRho[nR2[i][1][1]] + DRho[nR2[i][1][3]]));*/
			double alpha_i = diff_coefficient*(1*(DRho[nR2[i][0][0]] + DRho[nR2[i][0][2]] 
				+ DRho[nR2[i][2][0]] + DRho[nR2[i][2][2]]) +4*(DRho[nR2[i][0][1]] + DRho[nR2[i][2][1]] 
				+ DRho[nR2[i][1][0]] + DRho[nR2[i][1][2]]));
			
			double mu = -1 + 2.0*alpha_i/(sigma*sigma);
	        poisson = poisson_distribution<int>(lambda*lambda_exp*Rho_dt[i]);
			
	        gamma = gamma_distribution<double>(mu + 1 + poisson(rng), 1.0);
					
	        Rho_dt[i]= gamma(rng)/lambda;

			}

			//Euler Integration of remaining terms
			for(int i=0;i<Rho_dt.size();i++)
			{
			Rho_dt[i] = (Rho_dt[i])/( 1 + b*Rho_dt[i]*dt); 
			DRho[i] =Rho_dt[i]; //Updating dummy variable (necessary to prevent mixing of time in diffusion)
			}

			t+=dt; //Update timestep

			/**cout << "MOULIN ROUGE Final Rho(t):" <<endl;
			for(int i=0;i<Rho_dt.size();i++)
			{ cout<< Rho_dt[i] << " "; } cout<<endl;*/

		} // End of t loops

		if(j == 0) //Namely, the first replicate, set up initial incremental SD and mean of replicates accordingly.
		{
			stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "INITIAL BLOWHARD FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << Rho_M.size() << endl; cout << m7.str();
			for(int i =0; i< tot_iter; i++)
			{	rho_rep_avg_var[i][0] = Rho_M[i][0]; rho_rep_avg_var[i][1] = Rho_M[i][1]; rho_rep_avg_var[i][2] = 0.0;
				
			} //Namely Var(t,r=1) = 0, Mean_Rho(t, r=1) = Rho_M(t)
		}
		else
		{
			/**stringstream m7;     //To make cout thread-safe as well as non-garbled due to race conditions.
        	m7 << "LATER BLOWHARD FOR Replicate:\t" << j << "\t for Thread Rank:\t " << omp_get_thread_num() 
			<< " THE SIZE OF RHO_M: " << Rho_M.size() << endl; cout << m7.str();**/
			// Second or higher replicate, use incremental advances.
			var_mean_incremental(rho_rep_avg_var, Rho_M, tot_iter, j+1); 
			//Updates SD and Mean values for <Rho>_x across replicates as new replicate data (Rho_M) becomes available.
		}


		vector<double>().swap(Rho_dt); vector<vector <double>>().swap(Rho_M);
		// Remove dynamic vectors from memory and free up allocated space.

		

	} // End of r loop.

	for(int i=0; i< tot_iter; i++)
	{	// Recall Rho is: | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
		Rho.push_back({a, rho_rep_avg_var[i][0], rho_rep_avg_var[i][1], rho_rep_avg_var[i][2]});
	}



}

void critical_exp_delta(double Rho_0[], int div, double t_max, double a_start, double a_end, 
	double b, double D, double sigma, double dt, double dx, int r,  int g)
{

	init_fullframe(Rho_0, g*g); //Returns Rho_0 with a full initial frame filled with ones.
	vector<double> a_space = linspace(a_start, a_end, div);
  // The pspace to iterate over.
    vector <double> t_measure = logarithmic_time_bins(t_max, dt); 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

  cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  for (int i=0; i< t_measure.size(); i++)
  {
	cout << t_measure[i] << " ";
  } cout << endl;

  	//usleep will pause the program in micro-seconds (1000000 micro-seconds is 1 second)
    const int microToSeconds = 1000000;   
    const double delay1 = 5 * microToSeconds;     //5 seconds
    
    cout<<"Delay 1 in progress... ("<<delay1/microToSeconds<<"s)"<<endl;
	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 20)
	{
		omp_set_num_threads(20); //Limiting use on Chunk.
	}

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */
	  
      #pragma omp for nowait schedule(static)
      for (int i=0; i < a_space.size(); i++)
      {
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Value:\t" << a_space[i] << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
		**/
		percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, a_space[i], b, D, sigma, dt, dx, r, g);
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
          message3 << "Is this happening?\n";
          cout << message3.str();

      }
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, tm ,d3, p1, p2, rini, Dm, Sm;

  L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a_start; p2 << setprecision(4) << a_end;
  rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_dp;
  // Creating a file instance called output to store output data as CSV.
	output_dp.open("../Data/SPDP/P_c_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a1_"+ p1.str() +
	"_a2_"+ p2.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
	output_dp << " a , t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>_x]_r \n";
	cout << "The vector elements are: "<< endl;
  	cout << " a , t ,  <<Rho(x; t)>_x>_r, Var[<Rho(t)>x]_r \n";

	for(int i=0; i< vec.size(); i++)
	{
		output_dp << setprecision(8) << vec[i][0] << "," << vec[i][1] << "," << setprecision(16)
		<< vec[i][2] << "," << setprecision(16) << vec[i][3] << endl;
		if( i%(10000) ==1)
    {
			cout << setprecision(5) << vec[i][0] << "," << vec[i][1] << "," << setprecision(9)
			<< vec[i][2] << "," <<  setprecision(9) << vec[i][3] << endl;
    }
	}
	output_dp.close();
}


void critical_exp_theta_delta(double Rho_0[], int div, double t_max, double a_c, double b, 
		double D, double sigma, double dt, double dx, int r,  int g)
{

	init_solitarytear(Rho_0, g*g); //Returns Rho_0 with a single central seed admist a sea of drowned despair.

    vector <double> t_measure = logarithmic_time_bins(t_max, dt); 
  // Returns time-points measured on a natural logarithmic scale from e^{2} to e^ln(t_max) rounded down to one decimal place.

  cout << "Values in t_measure (which is of total size " << t_measure.size() << ") are :" <<endl;
  for (int i=0; i< t_measure.size(); i++)
  {
	cout << t_measure[i] << " ";
  } cout << endl;

	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

	auto start = high_resolution_clock::now();

	int nProcessors=omp_get_max_threads();
	if(nProcessors > 20)
	{
		omp_set_num_threads(20); //Limiting use on Chunk.
	}

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */
	  
      #pragma omp for nowait
      for (int i=0; i < div; i++)
      {
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on a Supra-Replicate:\t" << i*20 << "\t on Thread No:\t" << omp_get_thread_num() <<endl;
        cout << message.str();
        std::vector<vector <double>> CExpRho_a; //Stores relevant details for each time step to compute delta.
		/**
		 * Namely CExpRho_a is structured as:
		 * | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
		**/
		percolationDornic_2D(CExpRho_a, t_measure, Rho_0,  t_max, a_c, b, D, sigma, dt, dx, r, g);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), CExpRho_a.begin(), CExpRho_a.end());
		vector<vector <double>>().swap(CExpRho_a);
		// Remove dynamic vectors from memory and free up allocated space.
      }

      #pragma omp critical
      
      vec.insert(vec.end(), vec_private.begin(), vec_private.end());
      // Inserting critical exponent data for each grid size in order.
      stringstream message3;
      message3 << "Is this happening?\n";
      cout << message3.str();

      
	  vector<vector <double>>().swap(vec_private);
	  // Remove dynamic vectors from memory and free up allocated space.
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Dornic Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, tm ,d3, p1, p2, rini, Dm, Sm;

  L << g; tm << t_max; d3 << setprecision(3) << dt; p1 << setprecision(4) << a_c; //p2 << setprecision(4) << a_end;
  rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_dp;
  // Creating a file instance called output to store output data as CSV.
	output_dp.open("../Data/SPDP/P_c_ThDel_DP_G_" + L.str() + "_T_" + tm.str() + "_dt_" + d3.str() + "_a1_"+ p1.str() +
	"_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	// Output =  | 	a		|    t 		|     <<Rho(t)>>x,r			|    Var[<Rho(t)>x],r    |
	output_dp << " a , t ,  <<Rho(x, t)>_x>_r, Var[<Rho(t)>_x]_r \n";
	cout << "The vector elements are: "<< endl;
  	cout << " a , t ,  <<Rho(x, t)>_x>_r, Var[<Rho(t)>x]_r \n";

	for(int i=0; i< vec.size(); i++)
	{
		output_dp << setprecision(8) << vec[i][0] << "," << vec[i][1] << "," << setprecision(16)
		<< vec[i][2] << "," << setprecision(16) << vec[i][3] << endl;
		if( i%(10000) ==1)
    {
			cout << setprecision(5) << vec[i][0] << "," << vec[i][1] << "," << setprecision(9)
			<< vec[i][2] << "," <<  setprecision(9) << vec[i][3] << endl;
    }
	}
	output_dp.close();
}
