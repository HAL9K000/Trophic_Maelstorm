#include "trophic_dynamics.h"

int main()
{
  increase_stack_limit(128L); //Increase stack limit to 128 MB.

  double p0i, p0j, p0m;
  double dt, t; int g, div;
  double mR = 1; double mC = 10; double mP=100; //Mass of creatures.
  double a = (1)*1.71*pow(10.0,-6.0)*pow(mR,-0.25)*3600;
  double b = a*1*pow(mR, 0.25); double aij =pow(10,-3.08)*pow(mC,-0.37)*3600;
  double hij = 1*pow(mC,-0.02)/3600.0; double m= 4.15*pow(10,-8.0)*pow(mC,-0.25)*3600;
  double ajm = pow(10,-3.08)*pow(mP,-0.37)*3600; double hjm = 1*pow(mP,-0.02)/3600.0;
  double mm= 4.15*pow(10,-8.0)*pow(mP,-0.25)*3600; double ej =0.45; double em = 0.85;
  double d0 = pow(10,-5.0); double d1 = pow(10,-3.0); double d2 = pow(10,-3.0);

  cout << "Enter desired time-step: ";
  cin >> dt;

  cout << "Enter maximum duration of simulation (in hr): ";
  cin >> t;

  cout << "Enter number of frames to be captured (set <=10 for efficiency): ";
  cin >> div;

  cout << "Enter the desired grid-size (L): ";
  cin >> g;

  //cout << "Enter 1 to read from csv or 0 to input standardised initial condions"

  /* cout << "Enter initial density of Primary Producer (kg/m^2): ";
  cin >> p0i;

  cout << "Enter initial density of Primary Consumer (kg/m^2): ";
  cin >> p0j;

	cout << "Enter initial density of Secondary Producer (kg/m^2): ";
  cin >> p0m; */

  std::map<string, double> c = {{"a", a}, {"b", b}, {"aij", aij}, {"hij", hij},
          {"ajm", ajm}, {"hjm", hjm}, {"ej", ej}, {"em", em},
          {"m", m}, {"mm", mm}, {"alp", mm/em}, {"D0", d0}, {"D1", d1},{"D2", d2} };

 float data[g*g][4];
 std::ifstream file("Data/3Species/Diff/Init_Con/L8_Tier1_Random.csv");
 // Read initial densities from a csv file and save the same to data.

 file.ignore(140, '\n'); //ignore the first 140 characters, or until first \n, whichever is met first

 for(int row = 0; row < g*g; ++row)
 {
       std::string line;
       std::getline(file, line);
       if ( !file.good() )
           break;

       std::stringstream iss(line);

       for (int col = 0; col < 4; ++col)
       {
           std::string val;
           std::getline(iss, val, ',');
           if ( !iss.good() )
               break;

           std::stringstream convertor(val);
           convertor >> data[row][col];
       }
 }

  stringstream aeon, peoff,  iaggo, guano, guinea;

  aeon << setprecision(6) << c["a"];
  iaggo << setprecision(6) << c["alp"];
  peoff << setprecision(6) << c["b"];
  guinea << g;
  guano << setprecision(5) << dt;
  ofstream parm_troph;
  parm_troph.open("Data/3Species/Diff/ParameterLog_L_" + guinea.str() + "_r_" + aeon.str()  + "_q_"
    + peoff.str() + "_alph_" + iaggo.str() + "_dt_" + guano.str() +  ".txt");

  parm_troph << "L:\t" << g << endl;
  parm_troph << "dt:\t" << dt << endl;
  parm_troph << "mR:\t" << mR << endl;
  parm_troph << "mC Escher:\t" << mC << endl;
  parm_troph << "mP:\t" << mP << "\n" << endl;

  for (auto const& [key, val] : c)
  {
    parm_troph << key        // string (key)
               << ":\t" << setprecision(14) << val << endl;       // string's value
  }

  parm_troph << "Initial Densities (kg/m^2):\t [" << p0i << " , " << p0j
             << " , " << p0m << " ]" <<endl;

  parm_troph << "\n" << endl; //Blank space.


  cout << "Positive Coexistance Equilibria:\n";
  vector <double> pos_eq = rho_eqpk(c);
  cout << "[ " << pos_eq[0] << " , " << pos_eq[1] << " , " << pos_eq[2] << " ]" <<endl;

  parm_troph << "Positive Coexistance Equilibria:\t [" << pos_eq[0] << " , " << pos_eq[1]
             << " , " << pos_eq[2] << " ]" <<endl;

  bool flag = check_stability(c, "+");
  if( flag == true)
  { cout << "Above equilibria is STABLE" << '\n';}
  else { cout << "Above equilibria is UNSTABLE" << '\n';}

  cout << "Extinct Predator Equilibria:\n";
  vector <double> nopm_eq = rhono_eq(c);
  cout << "[ " << nopm_eq[0] << " , " << nopm_eq[1] << " , " << nopm_eq[2] << " ]" <<endl;

  parm_troph << "Extinct Predator Equilibria:\t [" << nopm_eq[0] << " , " << nopm_eq[1]
             << " , " << nopm_eq[2] << " ]" <<endl;
  parm_troph.close();

  flag = check_stability(c, "0");
  if( flag == true)
  { cout << "Above equilibria is STABLE" << '\n';}
  else { cout << "Above equilibria is UNSTABLE" << '\n';}

  flag = check_stability(c, "00");

  vector<vector<double>> Pijm;
  for(int i=0; i< g*g; i++)
  {   //
      p0i = data[i][1]; p0j = data[i][2]; p0m = data[i][3];
      Pijm.push_back({p0i, p0j, p0m});
      cout << i << "\t[ " << p0i << " , " << p0j << " , " << p0m << " ]" <<endl;
  }

  t = static_cast<double>(t); div = static_cast<int>(div);
  std::vector<double> t_stop;// = linspace(0.0, t, div);

  topical_trophics_2D(c, Pijm, t_stop, div, t, dt, g);

  return 0;
}
