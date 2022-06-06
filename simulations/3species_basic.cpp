#include "trophic_dynamics.h"

int main()
{
  increase_stack_limit(256);

  double p0i, p0j, p0m;
  double dt, t;
  double mR = 1; double mC = 10; double mP=100; //Mass of creatures.
  double a = (1)*1.71*pow(10.0,-6.0)*pow(mR,-0.25)*3600;
  double b = a*1*pow(mR, 0.25); double aij =pow(10,-3.08)*pow(mC,-0.37)*3600;
  double hij = 1*pow(mC,-0.02)/3600.0; double m= 4.15*pow(10,-8.0)*pow(mC,-0.25)*3600;
  double ajm = pow(10,-3.08)*pow(mP,-0.37)*3600; double hjm = 1*pow(mP,-0.02)/3600.0;
  double mm= 4.15*pow(10,-8.0)*pow(mP,-0.25)*3600; double ej =0.45; double em = 0.85;



  cout << "Enter desired time-step: ";
  cin >> dt;

  cout << "Enter total duration of simulation (in hr): ";
  cin >> t;

  cout << "Enter initial density of Primary Producer (kg/m^2): ";
  cin >> p0i;

  cout << "Enter initial density of Primary Consumer (kg/m^2): ";
  cin >> p0j;

	cout << "Enter initial density of Secondary Producer (kg/m^2): ";
  cin >> p0m;

  std::map<string, double> c = {{"a", a}, {"b", b}, {"aij", aij}, {"hij", hij},
          {"ajm", ajm}, {"hjm", hjm}, {"ej", ej}, {"em", em},
          {"m", m}, {"mm", mm}, {"alp", mm/em}};

  cout << "Positive Coexistance Equilibria:\n";
  vector <double> pos_eq = rho_eqpk(c);
  cout << "[ " << pos_eq[0] << " , " << pos_eq[1] << " , " << pos_eq[2] << " ]" <<endl;

  bool flag = check_stability(c, "+");
  if( flag == true)
  { cout << "Above equilibria is STABLE" << '\n';}
  else { cout << "Above equilibria is UNSTABLE" << '\n';}

  cout << "Extinct Predator Equilibria:\n";
  vector <double> nopm_eq = rhono_eq(c);
  cout << "[ " << nopm_eq[0] << " , " << nopm_eq[1] << " , " << nopm_eq[2] << " ]" <<endl;

  flag = check_stability(c, "0");
  if( flag == true)
  { cout << "Above equilibria is STABLE" << '\n';}
  else { cout << "Above equilibria is UNSTABLE" << '\n';}

  flag = check_stability(c, "00");

  vector<vector<double>> Pijm; Pijm.push_back({p0i, p0j, p0m});

  topical_trophics(c, Pijm,  dt,  t);

  return 0;
}
