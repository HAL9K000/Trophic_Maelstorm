#include "trophic_dynamics.h"

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

//----------------------------- Theoretical Constructs --------------------------------------------//

vector<double> rho_eqpk(map<string, double> &c)
{
  //One of two (possible) sets of positive theoretical equilibria denoting coexistance.

  //Recall alpha := mm/em
  //Using https://www.udacity.com/blog/2020/03/c-maps-explained.html
  double rhoj = c["alp"] /(c["ajm"]*(1.0- c["alp"]*c["hjm"]));
  double rhoi= (c["a"]*c["aij"]*c["hij"] - c["b"] + sqrt( pow(c["b"] + c["a"]*c["aij"]*c["hij"], 2.0)
                - 4*pow(c["aij"],2.0)*c["hij"]*c["b"]*rhoj ))/(2*c["b"]*c["aij"]*c["hij"]);

  double rhom= (c["ej"]*(c["a"] - c["b"]*rhoi)*rhoi - c["m"]*rhoj)/c["alp"];

  vector<double> rome={rhoi, rhoj, rhom};

  return rome;
}

vector<double> rhono_eq(map<string, double> &c)
{
  //One of two (possible) sets of positive theoretical equilibria denoting coexistance.

  //Recall alpha := mm/em
  //Using https://www.udacity.com/blog/2020/03/c-maps-explained.html
  double gamma = c["aij"]*(c["ej"] - c["m"]*c["hij"]);
  double rhoi= c["m"]/gamma;

  double rhoj= (c["a"]*gamma - c["b"]*c["m"])*c["ej"]/pow(gamma,2.0);

  vector<double> rome={rhoi, rhoj, 0};

  return rome;
}


int check_stability(map<string, double> &c, string flag)
{
  // Checks stability of equilibria points as per equivalent RH Criterion

  vector <double> th_eq; //Theoretical Equilibria
  if (flag == "+")
  {
    //Checking for positive equilibria
    th_eq = rho_eqpk(c);
    double rhs = (c["a"] - c["b"]*th_eq[0])*th_eq[0]*(c["hij"]*(c["a"] -c["b"]*th_eq[0])
        + c["hjm"]*c["ej"]*c["alp"])/(c["hjm"]*c["m"]*c["alp"] +c["b"]*th_eq[0]);
    if ( th_eq[1] >  rhs)
    {
      cout << "b2, b0 > 0 as per our analysis" << endl;
    }
    else
    {
      cout << "(R*+, G*+, P*+) is UNSTABLE"<< endl;
      return 0;
    }
    double abp = c["a"] -c["b"]*th_eq[0]; //A placeholder to make our life simpler.
    rhs = pow(abp,4)*c["ajm"]*(c["aij"]*c["hij"]*c["hjm"]*c["alp"]*th_eq[0]*( c["ej"]*th_eq[0]
      - c["m"]*th_eq[1]/c["alp"] ) + c["ej"]*th_eq[0] )*( ((c["a"]- c["alp"]*c["hjm"]*c["m"])/abp -1)*th_eq[1])
      + c["alp"]*c["hjm"]*c["ej"]*th_eq[0] -c["hij"]*abp*th_eq[0];
    double lhs = c["aij"]*c["mm"]*c["hjm"]*pow(th_eq[2]*c["alp"], 2.0);
    //Now checking if b2*b1 > b0
    if(lhs < rhs){
      cout << "b2 x b1 > b0 as per our analysis. (R*+, G*+, P*+) is STABLE."  << endl;
      return 1;
    }
    else{
      cout << "(R*+, G*+, P*+) is likely UNSTABLE"  << endl;
      return 0;
    }
    /** f ( (abp**4)*ajm*( aij*hij*hjm*rhoi*alpha*( ej*rhoi - m*rhoj/abp ) +ej*rhoi )*( ((a- alpha*hjm*m)/abp -1)*rhoj + alpha*hjm*ej*rhoi -hij*abp*rhoi  )
      - aij*mm*hjm*(rhom*alpha**2)**2 > 0 ): */
  }  //End of flag ="+"
  else if (flag == "0")
  { //Here rho_m* =0
    th_eq = rhono_eq(c);
    if( (c["em"]*c["ajm"]*th_eq[1]/(1+c["ajm"]*c["hjm"]*th_eq[0])) - c["mm"] > 0){
      cout<< "(R*, G*, 0) is UNSTABLE"  << endl;
      return 0;
    }
    else{
      double ahp2= pow(1 + c["aij"]*c["hij"]*th_eq[0], 2.0); //Defined for our convenience.
      double K2= ( pow(c["b"]*th_eq[0]*ahp2 - c["aij"]*c["aij"]*c["hij"]*th_eq[0]*th_eq[1], 2.0)
      - c["aij"]*c["m"]*th_eq[1]*ahp2 );
      if ( K2 < 0){
        //Imaginary lambdas.
        if ( c["aij"]*c["aij"]*c["hij"]*th_eq[0]*th_eq[1]/ahp2 -c["b"]*th_eq[0] > 0){
          //Real parts positive.
          cout<< "(R*, G*, 0) is UNSTABLE" <<endl;
          return 0;
        }
        else{
          cout << "(R*, G*, 0) is STABLE" <<endl;
          return 1;
        }
      }
      else{
        //Lamdas real.
        double r1 = c["aij"]*c["aij"]*c["hij"]*th_eq[0]*th_eq[1] - c["b"]*th_eq[0]*ahp2 + sqrt(K2);
          //Lambda =r1/(2*abp2)
          if(r1 > 0){
            cout<< "(R*, G*, 0) is UNSTABLE" <<endl;
            return 0;
          }
          else{
            cout<< "(R*, G*, 0) is STABLE" <<endl;
            return 1;
          }
        }

      }
    } //End of flag ="0"
    else if(flag == "00")
    {
      double rhoi= c["a"]/c["b"];
      if( c["a"] >= 0 && c["ej"]*c["aij"]*rhoi/(1+c["aij"]*c["hij"]*rhoi) - c["m"] < 0){
        cout<< "(R*, 0, 0) is STABLE" <<endl;
        return 1;
      }
      else{
        cout<< "(R*, 0, 0) is UNSTABLE" <<endl;
        return 0;
      }

    }  //End of flag ="00"

    return 0;
}

//----------------------------------- Basic Integration Schema ------------------------------------//

void f(double f[], double Pijm[], double t, map<string, double> &c)
{
  //Vector function that updates an array containing (dR/dt, dN1/dt, dN2/dt)

  f[0] = Pijm[0]*(c["a"] - Pijm[0]*c["b"]) - c["aij"]*Pijm[0]*Pijm[1]/(1+c["aij"]*c["hij"]*Pijm[0]);
  f[1] = c["ej"]*c["aij"]*Pijm[0]*Pijm[1]/(1+ c["aij"]*c["hij"]*Pijm[0]) - c["m"]*Pijm[1]
   - c["ajm"]*Pijm[1]*Pijm[2]/(1+c["ajm"]*c["hjm"]*Pijm[1]);
  f[2] = c["em"]*c["ajm"]*Pijm[1]*Pijm[2]/(1+c["ajm"]*c["hjm"]*Pijm[1]) - c["mm"]*Pijm[2];
}


void simulateRK4(vector<vector <double>> &Pijm, vector <double> &t_list, float tmax, double dt, map<string, double> &c)
{
  int n=3; //Number of species
  // Carries out RK4 integration.
  double K1[n] ={0.0}; double K2[n] ={0.0}; double K3[n] ={0.0}; double K4[n] ={0.0};
  //Initialise all frames to 0.
  double Pij_M[n] = {0.0}; double Pij_M2[n] ={0.0}; double Pij_M3[n] ={0.0}; double Pij_New[n] ={0.0};

  double Pijm0[n] ={Pijm[0][0], Pijm[0][1], Pijm[0][2]}; //Turning a new leaf.

  double t=0;
  int count =0;

  while (t < tmax)
  {
    f(K1, Pijm0, t, c); //K1 updated.
    for(int i = 0; i<n; i++)
    { Pij_M[i] = Pijm0[i] + (dt/2.0)*K1[i];  }
    f(K2, Pij_M, t +dt/2.0, c); //K2 updated.
    for(int i = 0; i<n; i++){
      Pij_M2[i] = Pijm0[i] + (dt/2.0)*K2[i]; }
    f(K3, Pij_M2, t +dt/2.0, c); //K3 updated.
    for(int i = 0; i<n; i++){
      Pij_M3[i] = Pijm0[i] + (dt)*K3[i]; }
    f(K4, Pij_M3, t +dt, c); //K4 updated.

    t+=dt; //Updating time step.
    count+=1;
    for(int i = 0; i<n; i++){
      Pij_New[i] = Pijm0[i] + (dt/6.0)*(K1[i]+ 2*K2[i] + 2*K3[i] +K4[i]);
      Pijm0[i] = Pij_New[i]; } //Updating Pijm0

    Pijm.push_back({Pij_New[0], Pij_New[1], Pij_New[2]}); t_list.push_back(t);
    //Updating results.

    //Debugging
    if(count%50000 == 2)
    {
      cout << "K1:\t [ " << setprecision(7) << K1[0] << " , " << setprecision(8) << K1[1] <<  " , " <<
        setprecision(7) << K1[2] << " ]" << endl;
      cout << "K2:\t [ " << setprecision(7) << K2[0] << " , " << setprecision(8) << K2[1] <<  " , " <<
          setprecision(7) << K2[2] << " ]" << endl;
      cout << "K3:\t [ " << setprecision(7) << K3[0] << " , " << setprecision(8) << K3[1] <<  " , " <<
          setprecision(7) << K3[2] << " ]" << endl;
      cout << "K4:\t [ " << setprecision(7) << K4[0] << " , " << setprecision(8) << K4[1] <<  " , " <<
          setprecision(7) << K4[2] << " ]" << endl;
      cout << "Pijm_New:\t [ " << setprecision(7) << Pijm0[0] << " , " << setprecision(8) << Pijm0[1] <<  " , " <<
          setprecision(7) << Pijm0[2] << " ]" << endl;
      cout << "t:\t" << setprecision(7) << t << " hr" << endl;
    }
  }
  cout << "Length of Pijm is:\t" << Pijm.size() << endl;
}

//------------------------- Trophic Chain Simulation (No Diff) --------------------------------------//

void topical_trophics(map<string, double> &c, vector<vector<double>> &Pijm, double dt, float t)
{
  vector <double> t_list; t_list.push_back(0); //t_list stores timestamps.
  ofstream output_troph;
  // Creating a file instance called output_troph to store output data as CSV.

  stringstream aeon, peoff,  iaggo, guano, tallaq;

  aeon << setprecision(6) << c["a"];
  iaggo << setprecision(6) << c["alp"];
  peoff << setprecision(6) << c["b"];
  //p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.

  guano << setprecision(5) << dt; tallaq << t;

  output_troph.open("Data/3Species/Rho_ijm_r" + aeon.str()  + "_q_" + peoff.str() + "_alph_" + iaggo.str()
    + "_dt_" + guano.str() + "_T_" + tallaq.str() + ".csv");

  auto start = high_resolution_clock::now();

  simulateRK4(Pijm, t_list, t, dt, c); //Simulate away to full glory.

  auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

  output_troph << "t, R0(t), R1(t), R2(t) \n";
  cout << "t,\t R0(t),\t R1(t),\t R2(t) \n";

  for ( int i = 0; i< Pijm.size(); i++)
  {
    output_troph << setprecision(10) << t_list[i] << "," << setprecision(16) << Pijm[i][0] << ","
      << setprecision(16) << Pijm[i][1] << "," << setprecision(16) <<  Pijm[i][2] << endl;
    if(i%500 == 0)
    {
      cout << setprecision(7) << t_list[i] << "," << setprecision(7) << Pijm[i][0] << ","
        << setprecision(8) << Pijm[i][1] << "," << setprecision(7) <<  Pijm[i][2] << endl;
    }
  }
  output_troph.close();

  cout << endl << "RK4 Integration Time: " << duration.count() << " seconds" << endl;

}

//------------------------- Trophic Chain Simulation (With Diff) --------------------------------------//

void find_neighbours_R2(double neighbours_R2[][2][5], vector<vector<double>> &Pijm, int n, int p, int q, int g)
{
	//Stores neighbour density in a ball of radius 2 around p,q.
	//int n=g*g;
	int a=0; int b= 0; //Useful for accounting purposes

	//ASSUMING REFLECTIVE BOUNDARY CONDITIONS.
	for(int i=0; i <n; i++)
	{ //Iterating over species.
		for(int j=0; j < 5; j++)
		{
			int k= j-2;
			a= ((p+k)%g + g)%g; // Eq to (p+k)mod g. -2 <= k <= 2
			b= ((q+k)%g + g)%g; // Eq to (q+k)mod g. -2 <= k <= 2

			//Recall, first row stores horizontal neighbours, second row stores vertical neighbours.
			neighbours_R2[i][0][j] = Pijm[a*g + q][i]; neighbours_R2[i][1][j] = Pijm[p*g + b][i];
		}
	}

}


void f_2D(double f[], double Pijm_f[], vector<vector<double>> &Pijm, int n, int p, int q, int g, double t, double dt, map<string, double> &c)
{

	//int n=3; //Number of species
  //Vector function that updates an array containing (dR/dt, dN1/dt, dN2/dt)

	double n_R2[n][2][5] ={0.0}; //Stores neighbours in a ball of radius 2 around p,q.
	//First row stores horizontal neighbours, second row stores vertical neighbours.

	find_neighbours_R2(n_R2, Pijm, n, p, q, g);
	//Finds these neighbours.

  f[0] = Pijm_f[0]*(c["a"] - Pijm_f[0]*c["b"]) - c["aij"]*Pijm_f[0]*Pijm_f[1]/(1+c["aij"]*c["hij"]*Pijm_f[0])
		- c["D0"]*(1/(12.0*dt*dt))*(1*(n_R2[0][0][0] + n_R2[0][0][4] + n_R2[0][1][0] + n_R2[0][1][4])
		-16*(n_R2[0][0][1] + n_R2[0][0][3] + n_R2[0][1][1] + n_R2[0][1][3]) +60*n_R2[0][0][2]);

		/**(1*(N_R[0,-2,0] + N_R[0,2,0] + N_R[0,0,-2] + N_R[0, 0, 2])
        - 16*(N_R[0,-1,0] + N_R[0, 1,0] + N_R[0, 0,-1] + N_R[0, 0, 1]) +60*N_R[0, 0,0]));*/

  f[1] = c["ej"]*c["aij"]*Pijm_f[0]*Pijm_f[1]/(1+ c["aij"]*c["hij"]*Pijm_f[0]) - c["m"]*Pijm_f[1]
   - c["ajm"]*Pijm_f[1]*Pijm_f[2]/(1+c["ajm"]*c["hjm"]*Pijm_f[1]) - c["D1"]*(1/(12.0*dt*dt))*(1*(n_R2[1][0][0] + n_R2[1][0][4]
	 + n_R2[1][1][0] + n_R2[1][1][4]) -16*(n_R2[1][0][1] + n_R2[1][0][3] + n_R2[1][1][1] + n_R2[1][1][3]) +60*n_R2[1][0][2]);


  f[2] = c["em"]*c["ajm"]*Pijm_f[1]*Pijm_f[2]/(1+c["ajm"]*c["hjm"]*Pijm_f[1]) - c["mm"]*Pijm_f[2]
		- c["D2"]*(1/(12.0*dt*dt))*(1*(n_R2[2][0][0] + n_R2[2][0][4] + n_R2[2][1][0] + n_R2[2][1][4])
		-16*(n_R2[2][0][1] + n_R2[2][0][3] + n_R2[2][1][1] + n_R2[2][1][3]) +60*n_R2[2][0][2]);
}


void simulate_RK42D(vector<vector<double>> &Pijm, vector<vector<double>> &Pij, double t, double dt, int p, int q, int g, map<string, double> &c)
{ //Run a single step of RK4 integration.
	int n=3; //Number of species
  // Carries out RK4 integration.
  double K1[n] ={0.0}; double K2[n] ={0.0}; double K3[n] ={0.0}; double K4[n] ={0.0};
  //Initialise all frames to 0.
  double Pij_M[n] = {0.0}; double Pij_M2[n] ={0.0}; double Pij_M3[n] ={0.0}; double Pij_New[n] ={0.0};

  double Pijm0[n] ={Pijm[p*g+q][0], Pijm[p*g+q][1], Pijm[p*g+q][2]}; //Turning a new leaf.


  f_2D(K1, Pijm0, Pijm, n, p, q, g, t, dt, c); //K1 updated.
  for(int i = 0; i<n; i++)
  { Pij_M[i] = Pijm0[i] + (dt/2.0)*K1[i];  }

  f_2D(K2, Pij_M, Pijm, n, p, q, g, t +dt/2.0, dt, c); //K2 updated.

  for(int i = 0; i<n; i++){
  	Pij_M2[i] = Pijm0[i] + (dt/2.0)*K2[i]; }

  f_2D(K3, Pij_M2, Pijm, n, p, q, g, t +dt/2.0, dt, c); //K3 updated.

  for(int i = 0; i<n; i++){
      Pij_M3[i] = Pijm0[i] + (dt)*K3[i]; }

  f_2D(K4, Pij_M3, Pijm, n, p, q, g, t +dt, dt, c); //K4 updated.

  t+=dt; //Updating time step.
  for(int i = 0; i<n; i++){
      Pij_New[i] = Pijm0[i] + (dt/6.0)*(K1[i]+ 2*K2[i] + 2*K3[i] +K4[i]);
      Pijm0[i] = Pij_New[i]; } //Updating Pijm0

  //Pijm.push_back({Pij_New[0], Pij_New[1], Pij_New[2]}); t_list.push_back(t);
	for(int i = 0; i<n; i++){
		Pij[p*g+q][i] = Pij_New[i]; }
    //Updating Pij results.

    //Debugging
    /** if(count%50000 == 2)
    {
      cout << "K1:\t [ " << setprecision(7) << K1[0] << " , " << setprecision(8) << K1[1] <<  " , " <<
        setprecision(7) << K1[2] << " ]" << endl;
      cout << "K2:\t [ " << setprecision(7) << K2[0] << " , " << setprecision(8) << K2[1] <<  " , " <<
          setprecision(7) << K2[2] << " ]" << endl;
      cout << "K3:\t [ " << setprecision(7) << K3[0] << " , " << setprecision(8) << K3[1] <<  " , " <<
          setprecision(7) << K3[2] << " ]" << endl;
      cout << "K4:\t [ " << setprecision(7) << K4[0] << " , " << setprecision(8) << K4[1] <<  " , " <<
          setprecision(7) << K4[2] << " ]" << endl;
      cout << "Pijm_New:\t [ " << setprecision(7) << Pijm0[0] << " , " << setprecision(8) << Pijm0[1] <<  " , " <<
          setprecision(7) << Pijm0[2] << " ]" << endl;
      cout << "t:\t" << setprecision(7) << t << " hr" << endl;
    } */

  //cout << "Length of Pijm is:\t" << Pijm.size() << endl;
}

void topical_trophics_2D(map<string, double> &c, vector<vector<double>> &Pij, vector<double> &t_stop,
	 int div, double t_max, double dt,  int g)
{
	t_stop = linspace(10.0, t_max, div);
	for( int i =0 ; i< t_stop.size(); i++)
		{ cout << t_stop[i] <<endl; }
	//vector <double> t_list; t_list.push_back(0); //t_list stores timestamps.
	auto start = high_resolution_clock::now();

	stringstream aeon, peoff,  iaggo, guano, guinea, tallaq;

  aeon << setprecision(6) << c["a"];
  iaggo << setprecision(6) << c["alp"];
  peoff << setprecision(6) << c["b"];
  //p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
	guinea << g;
  guano << setprecision(5) << dt;

	vector<vector<double>> Pijm; //Kept constant, only updated at end of turn
  for(int i=0; i< g*g; i++)
  {   //
      Pijm.push_back({ Pij[i][0],  Pij[i][1],  Pij[i][2]});
  }

	//int t_max = static_cast<int>(t_stop[t_stop.size() -1]); //Storing last element of the same here.
	double t=0; //Initialise t
	int r=0;
	cout << t_stop[r] << endl;
	while (t < t_max + dt)
	{
		//End at the end of all time.
		//cout << "Cuck" << endl;
		for(int i=0; i < g*g; i++)
		{
			//Going over all the sites.
			int p = int(i/g); int q = i%g; //Getting 2D coordinates.

			simulate_RK42D(Pijm, Pij, t, dt, p, q, g, c);
			//One step of RK4.
			//cout << t << endl;

		}
		t+=dt; //Update timestep

		//Updating Pijm at end of sweep.

		for(int i=0; i< g*g; i++)
	  {   //
	      Pijm[i][0]= Pij[i][0]; Pijm[i][1]= Pij[i][1]; Pijm[i][2]=  Pij[i][2];
	  }

		//cout << t << endl;
		if(t > t_stop[r] -dt/2.0 && t <= t_stop[r] + dt/2.0)
		{
			//Save at this time step:
			cout << "HIT!!!" << t << endl;
			tallaq << t; ofstream output_troph;
			output_troph.open("Data/3Species/Diff/Rho_ijm_L_" + guinea.str() + "_r_" + aeon.str()  + "_q_"
				+ peoff.str() + "_alph_" + iaggo.str() + "_dt_" + guano.str() + "_T_" + tallaq.str() + ".csv");

			// Creating save file.
			for ( int k = 0; k< Pijm.size(); k++)
		  {
		    output_troph << setprecision(8) << t_stop[r] << "," << k << "," << setprecision(16) << Pijm[k][0] << ","
		      << setprecision(16) << Pijm[k][1] << "," << setprecision(16) <<  Pijm[k][2] << endl;
		    if(k%(g) == 0)
		    {		//Debugging
		      cout << setprecision(7) << t_stop[r] << "," << k << "," << setprecision(7) << Pijm[k][0] << ","
		        << setprecision(8) << Pijm[k][1] << "," << setprecision(7) <<  Pijm[k][2] << endl;
		    }
		  }
		  output_troph.close(); tallaq.str(""); //Flush stringstream variable content and output file from memory.
			r = (r+1)%t_stop.size(); // To prevent out of bounds exception.
		}

	} // End of t loops

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "RK4 Integration Time: " << duration.count() << " seconds" << endl;

}
