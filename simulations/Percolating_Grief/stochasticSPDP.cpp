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
	return sum/(float)size;
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

//----------------------------- Stochastic SPDP RK4 Integration Machinery --------------------------------------------//

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

void f_2D(double f[], vector<double> &Rho_t, double p, double D, double t, double dt, double dx, int g)
{
  //Vector function that evaluates drho/dt = f(p,t) at given p and t.

	for(int i=0; i < g*g; i++)
	{
		int a = int(i/g); int b = i%g; //Getting 2D coordinates.
		double n_R2[2][5] ={0.0}; //Stores neighbours in a ball of radius 2 around a,b.
		//First row stores horizontal neighbours, second row stores vertical neighbours.

		find_neighbours_R2(n_R2, Rho_t, a, b, g); //Finds these neighbours.

		f[i] = -Rho_t[i]*Rho_t[i] - D*(1/(12.0*dx*dx))*(1*(n_R2[0][0] + n_R2[0][4] + n_R2[1][0] + n_R2[1][4])
		-16*(n_R2[0][1] + n_R2[0][3] + n_R2[1][1] + n_R2[1][3]) +60*n_R2[0][2]);
		// RK4 update of Non-linear parts plus diffusion term
		// 4th order Laplacian from 2SHOC.
	}
}

void percolation_2D(vector<vector<double>> &Rho, vector<double> &t_stop, double Rh0[],
	 double t_max, double p, double D, double sigma, double dt, double dx, int r,  int g)
{
	int n_roll =10000; int n_star = 1000; 
	
	for(int j=0; j< r; j++)
	{
		// Iterating over replicates.

		stringstream m1;     //To make cout thread-safe as well as non-garbled due to race conditions.
        m1 << "Replicate:\t" << j <<endl; cout << m1.str();
		int rd = std::random_device{}(); // random device engine, usually based on /dev/random on UNIX-like systems
		//https://cplusplus.com/forum/beginner/220120/
		thread_local std::mt19937_64 rng; // initialize Mersennes' twister using rd to generate the seed
		rng.seed(rd);
		//thread_local std::default_random_engine rng(rd); //engine
		//Credits: vsoftco, https://stackoverflow.com/questions/29549873/stdmt19937-doesnt-return-random-number


		/** Recall in Mean Field Limit (MFT), DP is given by:
		d(rho)/dt = (2p -1)rho - 1*rho^(2) **/
		//Defining the Dornic variables
		double lambda_exp = exp( (2*p -1.0)*dt); double lambda = 2*(2*p -1)/(sigma*sigma*(lambda_exp -1.0));

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

			f_2D(K1, Rho_dt, p, D, t, dt, dx, g); //K1 updated. 
			cout<< "K1" << endl;
			for(int i = 0; i<g*g; dx, i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K1[i];  
			  cout << K1[i] << " ";
			} //Recall Runge Kutta-Update Schema for K2 = f(t+dt/2, y+dt/2*K1)
			cout << endl;
			f_2D(K2, Rho_M, p, D, t +dt/2.0, dt, dx, g); //K2 updated.
			cout<< "K2" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt/2.0)*K2[i]; cout << K2[i] << " ";  } //Recall Runge Kutta-Update Schema for K3
		  	cout << endl;
			f_2D(K3, Rho_M, p, D, t +dt/2.0, dt, dx, g); //K3 updated.
			cout<< "K3" << endl;
			for(int i = 0; i<g*g; i++)
		  	{ Rho_M[i] = Rho_dt[i] + (dt)*K3[i]; cout << K3[i] << " ";  } //Recall Runge Kutta-Update Schema for K4
			cout << endl;
			f_2D(K4, Rho_M, p, D, t + dt, dt, dx, g); //K4 updated.
			cout<< "K4" << endl;
			for(int i = 0; i<g*g; i++)
			{
				cout << K4[i] << " ";
			}
			cout << endl;
			 stringstream m5;     //To make cout thread-safe as well as non-garbled due to race conditions.
      		m5 << "A7RK4" <<endl; cout << m5.str();
			//Now that K1, K2, K3, K4 updated for all sites, do the stochastic update.

			poisson_distribution<int> poisson; gamma_distribution <double> gamma;

			//Impproved Dornic  Integration of Linear & Stochastic Term
			//CREDITS: PAULA VILLA MARTIN.
			t+=dt; //Update timestep

			for(int i=0;i<Rho_dt.size();i++)
			{
	        Rho_dt[i]+= (dt/6.0)*( K1[i] + 2*K2[i] + 2*K3[i] + K4[i]);
					// Deterministic RK4 complete.
			double mu = lambda*lambda_exp*Rho_dt[i];
			if ( int(t)%50 == 0)
			{ cout << "Debug module: i*g+j: " << i << " Rho_ij: " << Rho_dt[i] << "Mu: "  
			<< mu << endl; }

	        poisson = poisson_distribution<int>(lambda*lambda_exp*Rho_dt[i]);
			//cout << "Poisson Pill:\t" << poisson(rng) <<endl;
	        gamma = gamma_distribution<double>(poisson(rng), 1.0);
					//Defining shapes of poisson and gamma distributions.
	        Rho_dt[i]= gamma(rng)/lambda;
			//cout << "Moulin Rouge:\t" << Rho_dt[i] <<endl;
	    	}
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
      			m4 << "LOCKET!:\t" << t << " " << p  <<endl; cout << m4.str();
				for ( int k = 0; k< Rho_dt.size(); k++)
			  {	// Recall Rho = | p (Order Parameter) | Replicate (r) | Time(t)  | a*L + b | Cell Density[a][b] |
					Rho.push_back({p, static_cast<double>(j), t, static_cast<double>(k), Rho_dt[k]});
			  }
				s = (s+1)%t_stop.size(); // To prevent out of bounds exception.
			}

		} // End of t loops



		vector<double>().swap(Rho_M); vector<double>().swap(Rho_dt);
		// Remove dynamic vectors from memory and free up allocated space.

		

	}



}

void stoc_dem_percolation(double Rho_0[], vector<double> &t_stop, int div, double t_max,
	double p_start, double p_end, double D, double sigma, double dt, double dx, int r,  int g)
{
	vector<double> p_space = linspace(p_start, p_end, div);
  // The pspace to iterate over.
	t_stop = linspace(t_max/3.0, t_max, 2); // 2 Time-stamps at which screenshots will be captured.
	for( int i =0 ; i< t_stop.size(); i++)
		{ cout << t_stop[i] <<endl; }

	std::vector<vector <double>> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.

	auto start = high_resolution_clock::now();

	#pragma omp parallel
  {
      std::vector<vector<double>> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

      #pragma omp for nowait schedule(static)
      for (int i=0; i < p_space.size(); i++)
      {
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on p Value:\t" << p_space[i] <<endl;
        cout << message.str();
        std::vector<vector <double>> Rho_p;
		double Rh0[g*g]={0.0}; //Stores initial density conditions
		for(int j =0; j < g*g; j++)
		{	Rh0[j] = Rho_0[j]; }
		percolation_2D(Rho_p, t_stop, Rh0, t_max, p_space[i], D, sigma, dt, dx, r, g);
        //crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), Rho_p.begin(), Rho_p.end());

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
  } //End of Pragma.

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "RK4 Integration Time: " << duration.count() << " seconds" << endl;

	stringstream L, tm ,p1, p2, rini, Dm, Sm;

  L << g; tm << t_max; p1 << setprecision(4) << p_start; p2 << setprecision(4) << p_end;
  rini << r; Dm << setprecision(3) << D; Sm << setprecision(3) << sigma;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
	ofstream output_dp;
  // Creating a file instance called output to store output data as CSV.
	output_dp.open("../Data/SPDP/P_c_DP_G_" + L.str() + "_T_" + tm.str() + "_p1_"+ p1.str() +
	"_p2_"+ p2.str() + "_Sig_"+ Sm.str() + "_R_"+ rini.str() + ".csv");

	// Output = | p (Order Parameter) | Replicate (r) | Time(t)  | a*L + b | Cell Density[a][b] |
	output_dp << " p , # Tr No , t ,  i*L + j ,  Rho_ij(t) \n";
	cout << "The vector elements are: "<< endl;
  cout << " p , # Tr No , t ,  i*L + j ,  Rho_ij(t) \n";

	for(int i=0; i< vec.size(); i++)
	{
		output_dp << setprecision(8) << vec[i][0] << "," << static_cast<int>(vec[i][1]) << "," << setprecision(7)
		<< vec[i][2] << "," << static_cast<int>(vec[i][3]) << "," << setprecision(16) << vec[i][4] << endl;
		if( i%(5000) ==1)
    {
			cout << setprecision(5) << vec[i][0] << "," << static_cast<int>(vec[i][1]) << "," << setprecision(7)
			<< vec[i][2] << "," << vec[i][3] << "," << setprecision(16) << vec[i][4] << endl;
    }
	}
	output_dp.close();


}
