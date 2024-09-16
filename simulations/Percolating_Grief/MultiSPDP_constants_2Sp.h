#include <string>
#include <stdlib.h>
#include <stdio.h>

inline const int Sp = 2; //Total number of species in the system.
inline const std::string frame_header = "a_c,  x,  P(x; t), G(x; t) \n"; 
//Header for frame files.
inline std::string prelimheader = " a , r, L, t , <<P(x; t)>_x>_r, Var[<P(x; t)>_x]_r, # Surviving Runs P(x; t),"
					" # Active Sites P(x; t), <<G(x; t)>_x>_r, Var[<G(x; t)>_x]_r, # Surviving Runs G(x; t), # Active Sites G(x; t), \n";
//Header for preliminary files.