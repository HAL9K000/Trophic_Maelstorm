#include "2species_stochastic.h"

int main()
{
    increase_stack_limit(512LL);
    cout << maxis(3.14, 1.78) <<endl;

    int c =8;
    add_three(1,2,c); int g =3; 

    D2Vec_Double Rho_0(Sp, vector<double> (g*g, 1.0));
    double p0i = 1.0; double p0j= 0.05;
    double mean[Sp] = {p0i, p0j}; double sd[Sp] = {p0i/4.0, p0j/4.0};
	init_randframe(Rho_0, Sp,  g*g, mean, sd); //Returns Rho_0 with a full initial frame filled with 0.2.

    for(int i =0; i <Sp; i++)
    {
        for(int j=0; j<g*g; j++)
            cout << Rho_0[i][j] << " ";
        cout << endl;
    }

}

