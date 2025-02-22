#include <string>
#include <stdlib.h>
#include <stdio.h>

inline const int Sp = 5; //Total number of species in the system.
inline const std::string frame_header = "a_c,  x,  P(x; t), G(x; t), Pr(x; t), W(x; t), O(x; t) \n"; 
//Header for frame files.
inline std::string prelimheader = " a , r, L, t , <<W(x; t)>_x>_r, <<O(x; t)>_x>_r,  <<P(x; t)>_x>_r, Var[<P(x; t)>_x]_r, # Surviving Runs P(x; t),"
					" # Active Sites P(x; t), <<G(x; t)>_x>_r, Var[<G(x; t)>_x]_r, # Surviving Runs G(x; t), # Active Sites G(x; t)," 
					"<<Pr(x; t)>_x>_r, Var[<Pr(x; t)>_x]_r, # Surviving Runs Pr(x; t), # Active Sites Pr(x; t), \n";
//Header for preliminary files.