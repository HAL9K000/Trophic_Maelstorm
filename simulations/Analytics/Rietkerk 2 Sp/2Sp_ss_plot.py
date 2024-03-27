import numpy as np
import matplotlib.pyplot as plt
import pandas as pan
import seaborn as sea
import cmath
import math
import numpy.linalg as la



# Define global variables.

c = 10000; gmax = 0.05*pow(10, -3.0)/24.0; d = 0.25/24.0; alpha =0.2/24.0; W0 = 0.2; rW = 0.2/24.0; # From Bonachela et al 2015
dx= 0.1 ; dt = 0.1; #From Bonachela et al 2015 (in km)
d0 = 0.00025/24.0; d1=0.0298; d2 = 0.00025/24.0; d3= 0.025/24.0; #From Bonachela et al 2015 (in km^2/hr)
k0= 0; k1 = 5; k2 =5000;

# Parmeters for grazer (Pawar & Co)
aij = 3.6*pow(10.0, -6.08)*pow(20.0, -0.37); # in km^2/(hr kg)
hij = 1; #Handling time in hrs
ej =0.45; mj = 0.061609*pow(20.0, -0.25)/8760.0; # Mortality rate in hr^{-1}

# Parameters for grazer (Kefi and Brose)
#Assuming mass of producer (mi) = 1 kg, mass of grazer (mj) = 10 kg
#ej = 0.45; mj = 0.314*pow(10, -0.25); y =8; B0 =0.5
#aij = (mj*y)/B0; hij = 1/(mj*y);

print("Parameters:\t")
print("c = %g, gmax = %g, d = %g, alpha = %g, W0 = %g, rW = %g" %(c, gmax, d, alpha, W0, rW))
print("cgmax = %g, d = %g" %(c*gmax, d))
print("K1 = %g, K2 = %g" %(k1, k2))
print("aij = %g, hij = %g, ej = %g, mj = %g" %(aij, hij, ej, mj))



# Define derived parameters/terms.

kappa = c*gmax - d;
Vstar_coex = mj/((ej- mj*hij)*aij);
Wstar_veg = d*k1/kappa;

print("Derived parameters:\t")
print("Kappa = %g, Vstar_coex = %g, Wstar_veg = %g" %(kappa, Vstar_coex, Wstar_veg))


R_trans = d*rW*k1/kappa;
omega= rW*k1 + gmax*Vstar_coex;
R_plus_min = 2*rW*k1 + omega + 2*rW*k1*math.sqrt(1 +omega/(rW*k1));
R_minus_min = 2*rW*k1 + omega - 2*rW*k1*math.sqrt(1 +omega/(rW*k1));


''' 
This script will plot certain functions of a variable "R" with defined functional form, across a range of values of R (specified by user).
Each of these functions will be plotted on the same graph, with the x-axis being the range of R values, 
and the y-axis being the value of the function.
Pretty colour palette will be used to distinguish between the functions (using Seaborn).
Some of the functions can return imaginary values for certain values of R, in which case they will not be plotted for those values of R.
'''

# Define the functions to be plotted.

# First functions for equilbirum values with no vegetation and grazer.

def Vstar_noveg(R):
    return np.zeros_like(R)

def Wstar_noveg(R):
    W = R/rW
    return W

def Ostar_noveg(R):
    O = R/(alpha*W0)
    return O

def Gstar_noveg(R): 
    return np.zeros_like(R)

# Next functions for equilibria with vegetation but no grazer.

def Vstar_veg(R):
    V = R - rW*Wstar_veg
    V *= (c/d)
    return V

def Wastar_veg(R):
    return Wstar_veg*np.ones_like(R)

def Ostar_veg(R):
    V = Vstar_veg(R)
    #num = V + k2; denom = V  + k2*W0;
    O = R/(alpha)*(V + k2)/(V+ k2*W0);# O *= R/(alpha)
    return O

def Gstar_veg(R): 
    return np.zeros_like(R)

# Next functions for equilibria with grazer and vegetation (coexistence).

def Vegstar_coex(R):
    return Vstar_coex*np.ones_like(R)

def Ostar_coex(R):
    V = Vegstar_coex(R)
    #num = V + k2; denom = V  + k2*W0;
    #O = num/denom; O *= R/(alpha)
    O = R/(alpha)*(V + k2)/(V+ k2*W0); #O *= R/(alpha)
    return O

# Wstar_coex and Gstar_coex can be calculated from the other equilibria.
# They can also return imaginary values for certain R, so we will not plot them for those values.

def Wstar_coex(R):
    V = Vegstar_coex(R)
    b = R  - gmax*V -rW*k1; #b.astype(complex);
    riyal = b/(2*rW); #riyal.astype(complex);
    Wstar_plus = riyal + np.emath.sqrt(b**2 - 4*rW*k1*R)/(2*rW)
    Wstar_minus = riyal - np.emath.sqrt(b**2 - 4*rW*k1*R)/(2*rW)
    return Wstar_plus, Wstar_minus

def Gstar_coex(R, Wstar_plus, Wstar_minus):
    V = Vegstar_coex(R)
    common= (ej*c/mj)*R - (ej*d/mj)*V; common.astype(complex);
    Gstar_plus = common- (ej*rW*c/mj)*Wstar_plus
    Gstar_minus = common - (ej*rW*c/mj)*Wstar_minus
    return Gstar_plus, Gstar_minus

def eigenval_check(R, Vstar, Wstar, Ostar, Gstar):
    # Check the eigenvalues of the Jacobian matrix at the equilibria.
    # If the real part of the eigenvalues are negative, the equilibria are stable.
    #Input: R, Vstar, Wstar, Ostar, Gstar are arrays of equilibrium values of the variables.
    #Output: Returns a boolean array of the same length as R, where True indicates stability and False indicates instability.
    #Output: Also returns the eigenvalues of the Jacobian matrix at the equilibria as a 2D array of size (len(R), 4).

    # Define the Jacobian matrix.
    J = np.zeros((len(R), 4, 4), dtype=complex)
    J[:, 0, 0] = -d + c*gmax*Wstar/(Wstar + k1) -aij*Gstar/np.emath.power(1 + aij*hij*Vstar, 2)
    J[:, 0, 1] = c*gmax*k1*Vstar/np.emath.power(Wstar + k1, 2)
    J[:, 0, 2] = 0
    J[:, 0, 3] = -aij*Vstar/(1 + aij*hij*Vstar)

    J[:, 1, 0] = alpha*k2*(1-W0)*Ostar/np.emath.power(Vstar + k2, 2) -gmax*Wstar/(Wstar + k1)
    J[:, 1, 1] = - gmax*k1*Vstar/np.emath.power(Wstar + k1, 2) - rW
    J[:, 1, 2] = alpha*(Vstar + k2*W0)/(Vstar + k2)
    J[:, 1, 3] = 0

    J[:, 2, 0] = -alpha*k2*(1-W0)*Ostar/np.emath.power(Vstar + k2, 2)
    J[:, 2, 1] = 0
    J[:, 2, 2] = -alpha*(Vstar + k2*W0)/(Vstar + k2)
    J[:, 2, 3] = 0

    J[:, 3, 0] = ej*aij*Gstar/np.emath.power(1 + aij*hij*Vstar, 2)
    J[:, 3, 1] = 0
    J[:, 3, 2] = 0
    J[:, 3, 3] = ej*aij*Vstar/(1 + aij*hij*Vstar) - mj

    print("At R = ", R[0], "Vstar, Wstar, Ostar, Gstar: ", 
          Vstar[0], "  ", Wstar[0], "  ", Ostar[0], "  ", Gstar[0])
    print("Jacobian matrix at R =  " +str(R[0]))
    print(J[0, 0, :]); print(J[0, 1, :]); print(J[0, 2, :]); print(J[0, 3, :]); print("\n")
    print("At R = ", R[len(R)//2], "Vstar, Wstar, Ostar, Gstar: ",
            Vstar[len(R)//2], "  ", Wstar[len(R)//2], "  ", Ostar[len(R)//2], "  ", Gstar[len(R)//2])
    print("Jacobian matrix at R = " +str(R[len(R)//2]))
    print(J[len(R)//2, 0, :]); print(J[len(R)//2, 1, :]); print(J[len(R)//2, 2, :]); 
    print(J[len(R)//2, 3, :]); print("\n")

    # Calculate the eigenvalues of the Jacobian matrix.
    eigvals = np.zeros((len(R), 4), dtype=complex)
    for i in range(len(R)):
        try:
            eigvals[i] = la.eigvals(J[i])
        except la.LinAlgError:
            print("LinAlgError at R = ", R[i])
            print("Jacobian matrix: ", J[i, :, :])
            print("Vstar, Wstar, Ostar, Gstar: ", Vstar[i], Wstar[i], Ostar[i], Gstar[i])
            eigvals[i] = np.array([np.nan, np.nan, np.nan, np.nan], dtype=complex)

        eigvals[i] = la.eigvals(J[i])
    
    # Check the stability of the equilibria.
    # If the real part of the eigenvalues are negative, the equilibria are stable.
    # Return a boolean array of the same length as R, where True indicates stability and False indicates instability.
    return np.all(eigvals.real < 0, axis=1), eigvals




v =1;
# Define the range of R values to be plotted.
R = np.linspace(0, 0.3, 1000)
#R = np.array([0.1, 0.15, 0.2, 0.25])



def ss_plot():
    # Plot the various functions defined above over the range of R values.
    # Use Seaborn to make the plot look pretty.
    # Use a different colour for each function.
    # Use a legend and axis labels to make the plot easy to understand.
    # Use ax.text to label the functions on the plot.
    # Wstar_coex and Gstar_coex can return imaginary values for certain R, so we will not plot them for those values.

    print("R_min_+, R_min_-, R_c = (%g, %g, %g)" %(R_plus_min, R_minus_min, R_trans))

    print("Testing the functions...")
    R_test =np.array([0.04, 0.06])
    print("At R = ", R_test)
    print("Vstar_noveg: ", Vstar_noveg(R_test))
    print("Wstar_noveg: ", Wstar_noveg(R_test))
    print("Ostar_noveg: ", Ostar_noveg(R_test))
    print("Gstar_noveg: ", Gstar_noveg(R_test))
    print("\n=====================================\n")
    print("Vstar_veg: ", Vstar_veg(R_test))
    print("Wstar_veg: ", Wastar_veg(R_test))
    print("Ostar_veg: ", Ostar_veg(R_test))
    print("Gstar_veg: ", Gstar_veg(R_test))
    print("\n=====================================\n")
    print("Vegstar_coex: ", Vegstar_coex(R_test))
    print("Ostar_coex: ", Ostar_coex(R_test))
    Wstar_plus, Wstar_minus = Wstar_coex(R_test)
    print("Wstar_plus: ", Wstar_plus)
    print("Wstar_minus: ", Wstar_minus)
    Gstar_plus, Gstar_minus = Gstar_coex(R_test, Wstar_plus, Wstar_minus)
    print("Gstar_plus: ", Gstar_plus)
    print("Gstar_minus: ", Gstar_minus)
    print("\n=====================================\n")
    


    fig, ax = plt.subplots(2,2, sharex=True, figsize=(15,10))
    sea.set(style='whitegrid')
    sea.set_palette('husl')
    ax[0,0].plot(R, Vstar_noveg(R), label='$V^*_\mathrm{Ext}$')
    ax[0,1].plot(R, Wstar_noveg(R), label='$W^*_\mathrm{Ext}$')
    ax[1,0].plot(R, Ostar_noveg(R), label='$O^*_\mathrm{Ext}$')
    ax[1,1].plot(R, Gstar_noveg(R), label='$G^*_\mathrm{Ext}$')

    ax[0,0].plot(R, Vstar_veg(R), label='$V^*_\mathrm{Veg}$')
    ax[0,1].plot(R, Wastar_veg(R), label='$W^*_\mathrm{Veg}$')
    ax[1,0].plot(R, Ostar_veg(R), label='$O^*_\mathrm{Veg}$')
    ax[1,1].plot(R, Gstar_veg(R), label='$G^*_\mathrm{Veg}$')

    ax[0,0].plot(R, Vegstar_coex(R), label='$V^*_\mathrm{Coexist}$')
    ax[1,0].plot(R, Ostar_coex(R), label='$O^*_\mathrm{Coexist}$')

    Wstar_plus, Wstar_minus = Wstar_coex(R)
    Gstar_plus, Gstar_minus = Gstar_coex(R, Wstar_plus, Wstar_minus)
    
    #Plot these complex arrays only for real element values.
    # Use numpy splicing.
    mask1 = np.isreal(Wstar_plus)

    ax[0,1].plot(R[mask1], Wstar_plus[mask1], label='$W^*_\mathrm{Coexist+}$')
    
    mask2 = np.isreal(Wstar_minus)
    ax[0,1].plot(R[mask2], Wstar_minus[mask2], label='$W^*_\mathrm{Coexist-}$')

    ax[1,1].plot(R[mask1], Gstar_plus[mask1], label='$G^*_\mathrm{Coexist+}$')
    ax[1,1].plot(R[mask2], Gstar_minus[mask2], label='$G^*_\mathrm{Coexist-}$')

    # Finally add three dashed horizontal lines to show values of R_plus_min, R_minus_min and R_trans ON ALL SUBPLOTS with appropriate axes.text labels (rotated 90)
    #axi = [ax[0], ax[1], ax[2], ax[3]]

    ax[0,0].set_ylim(-200, 20000); ax[0,0].set_title('$V^*$ vs R'); ax[0,1].set_title('$W^*$ vs R')
    ax[1,0].set_ylim(-20, 200); ax[1,0].set_title('$O^*$ vs R'); ax[1,1].set_title('$G^*$ vs R')

    for i in range(2):
      for j in range(2):
        axes = ax[i,j]
        axes.axvline(x=R_plus_min, color='tomato', linestyle='-.' )#, label='$R^+_{min}$')
        axes.axvline(x=R_minus_min, color='tomato', linestyle='-.')# , label='$R^-_{min}$')
        axes.axvline(x = R_trans, color='gray', linestyle=':', label= '$R_c$')
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axes.set_xlabel('R $(mm/hr)$')
        axes.set_ylabel('Equilibrium values $kg m^{-2}$')
        axes.tick_params(axis='both')

        # Find max and min values of y for each subplot.
        ymax = max(axes.get_ylim()); xmax = max(axes.get_xlim())
        axes.text(R_plus_min - xmax/50.0, ymax/5.0, '$R^{+}_{min}$', rotation=90, color='tomato')
        axes.text( R_minus_min - xmax/50.0, ymax/5.0, '$R^{-}_{min}$', rotation=90, color='tomato')
        axes.text(R_trans - xmax/50.0, ymax/10.0, '$R_{c}$', rotation=90, color='grey')

    plt.tight_layout()
    figure = plt.gcf() # get current figure
    figure.set_size_inches(15, 10)
    plt.savefig('2Sp_ss_plot Rst --- %g Rend --- %g.png' %(R[0], R[-1]), dpi=1000)
    plt.show()
    plt.close()


def stable_ss_plot():
    # Plot the various functions defined above over the range of R values.
    # Use Seaborn to make the plot look pretty.
    # Use a different colour for each function.
    # Use a legend and axis labels to make the plot easy to understand.
    # Use ax.text to label the functions on the plot.
    # Wstar_coex and Gstar_coex can return imaginary values for certain R, so we will not plot them for those values.

    print("R_min_+, R_min_-, R_c = (%g, %g, %g)" %(R_plus_min, R_minus_min, R_trans))

    #R= np.array([0.04, 0.25])
    # Extinction equilibria
    Veq_noveg = Vstar_noveg(R); Weq_noveg = Wstar_noveg(R); Oeq_noveg = Ostar_noveg(R); Geq_noveg = Gstar_noveg(R)
    # Vegetation only equilibria
    Veq_veg = Vstar_veg(R); Weq_veg = Wastar_veg(R); Oeq_veg = Ostar_veg(R); Geq_veg = Gstar_veg(R)
    # Coexistence equilibria
    Veq_coex = Vegstar_coex(R); Oeq_coex = Ostar_coex(R)
    Wstar_plus, Wstar_minus = Wstar_coex(R)
    Gstar_plus, Gstar_minus = Gstar_coex(R, Wstar_plus, Wstar_minus)

    # Check the stability of the equilibria using the eigenvalue check function.
    # If the real part of the eigenvalues are negative, the equilibria are stable.
    # We will plot the stable equilibria only.
    checkstable_noveg, eigenval_noveg = eigenval_check(R, Veq_noveg, Weq_noveg, Oeq_noveg, Geq_noveg)
    print("Stability profile of extinction equilibria: ")
    for i in range(0, len(R), len(R)//10):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_noveg[i], "with stability: ", checkstable_noveg[i], "\n\n")
    
    checkstable_veg, eigenval_veg = eigenval_check(R, Veq_veg, Weq_veg, Oeq_veg, Geq_veg)
    print("Stability profile of vegetation equilibria: ")
    for i in range(0, len(R), len(R)//10):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_veg[i], "with stability: ", checkstable_veg[i], "\n\n")

    checkstable_coex_plus, eigenval_coex_plus = eigenval_check(R, Veq_coex, Wstar_plus, Oeq_coex, Gstar_plus)
    print("Stability profile of coexistence equilibria (Wstar_plus): ")
    for i in range(0, len(R), len(R)//10):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_coex_plus[i], "with stability: ", checkstable_coex_plus[i], "\n\n")
    
    checkstable_coex_minus, eigenval_coex_minus = eigenval_check(R, Veq_coex, Wstar_minus, Oeq_coex, Gstar_minus)
    print("Stability profile of coexistence equilibria (Wstar_minus): ")
    for i in range(0, len(R), len(R)//10):
        print("Eigenvalues at R = %g: " %(R[i]), eigenval_coex_minus[i], "with stability: ", checkstable_coex_minus[i], "\n\n")
    
    checkstable_coex = np.logical_or(checkstable_coex_plus, checkstable_coex_minus) 
    # Coexistence equilibria for V and O are stable if either Wstar_plus or Wstar_minus are stable.
    print("Number of true values in checkstable_coex_plus: ", np.sum(checkstable_coex_plus))
    print("Number of true values in checkstable_coex_minus: ", np.sum(checkstable_coex_minus))
    print("Number of true values in checkstable_coex: ", np.sum(checkstable_coex))

    
    fig, ax = plt.subplots(2,2, sharex=True, figsize=(15,10))
    sea.set(style='whitegrid')
    sea.set_palette('husl')
    ax[0,0].plot(R[checkstable_noveg], Veq_noveg[checkstable_noveg], label='$V^*_\mathrm{Ext}$', alpha=0.75)
    ax[0,1].plot(R[checkstable_noveg], Weq_noveg[checkstable_noveg], label='$W^*_\mathrm{Ext}$', alpha=0.75)
    ax[1,0].plot(R[checkstable_noveg], Oeq_noveg[checkstable_noveg], label='$O^*_\mathrm{Ext}$', alpha=0.75)
    ax[1,1].plot(R[checkstable_noveg], Geq_noveg[checkstable_noveg], label='$G^*_\mathrm{Ext}$', alpha=0.75)

    ax[0,0].plot(R[checkstable_veg], Veq_veg[checkstable_veg], label='$V^*_\mathrm{Veg}$', alpha=0.75)
    ax[0,1].plot(R[checkstable_veg], Weq_veg[checkstable_veg], label='$W^*_\mathrm{Veg}$', alpha=0.75)
    ax[1,0].plot(R[checkstable_veg], Oeq_veg[checkstable_veg], label='$O^*_\mathrm{Veg}$', alpha=0.75)
    ax[1,1].plot(R[checkstable_veg], Geq_veg[checkstable_veg], label='$G^*_\mathrm{Veg}$', alpha=0.75)

    ax[0,0].plot(R[checkstable_coex], Veq_coex[checkstable_coex], label='$V^*_\mathrm{Coexist}$', alpha=0.75)
    ax[1,0].plot(R[checkstable_coex], Oeq_coex[checkstable_coex], label='$O^*_\mathrm{Coexist}$', alpha=0.75)

    ax[0,1].plot(R[checkstable_coex_plus], Wstar_plus[checkstable_coex_plus], label='$W^*_\mathrm{Coexist+}$', alpha=0.75)
    ax[0,1].plot(R[checkstable_coex_minus], Wstar_minus[checkstable_coex_minus], label='$W^*_\mathrm{Coexist-}$', alpha=0.75)

    ax[1,1].plot(R[checkstable_coex_plus], Gstar_plus[checkstable_coex_plus], label='$G^*_\mathrm{Coexist+}$', alpha=0.75)
    ax[1,1].plot(R[checkstable_coex_minus], Gstar_minus[checkstable_coex_minus], label='$G^*_\mathrm{Coexist-}$', alpha=0.75)

    ax[0,0].set_ylim(-100, 100); ax[0,0].set_title('Locally Stable $V^*$ vs R'); ax[0,1].set_title('Locally Stable $W^*$ vs R')
    ax[1,0].set_ylim(-20, 200); ax[1,0].set_title('Locally Stable $O^*$ vs R'); ax[1,1].set_title('Locally Stable $G^*$ vs R')

    # Set title for overall plot
    fig.suptitle('Stable Equilibria for Rietkerk + Grazer Model, Allometric Parameters', fontsize=16)
    
    # Finally add three dashed horizontal lines to show values of R_plus_min, R_minus_min and R_c ON ALL SUBPLOTS with appropriate axes.text labels (rotated 90)

    for i in range(2):
      for j in range(2):
        axes = ax[i,j]
        axes.axvline(x=R_plus_min, color='tomato', linestyle='-.' )#, label='$R^+_{min}$')
        axes.axvline(x=R_minus_min, color='tomato', linestyle='-.')# , label='$R^-_{min}$')
        axes.axvline(x = R_trans, color='gray', linestyle=':', label= '$R_c$')
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axes.set_xlabel('R $(mm/hr)$')
        axes.set_ylabel('Stable Equilibrium values $kg m^{-2}$')
        axes.tick_params(axis='both')

        # Find max and min values of y for each subplot.
        ymax = max(axes.get_ylim()); xmax = max(axes.get_xlim())
        axes.text(R_plus_min - xmax/50.0, ymax/5.0, '$R^{+}_{min}$', rotation=90, color='tomato')
        axes.text( R_minus_min - xmax/50.0, ymax/5.0, '$R^{-}_{min}$', rotation=90, color='tomato')
        axes.text(R_trans - xmax/50.0, ymax/10.0, '$R_{c}$', rotation=90, color='grey')

    plt.tight_layout()
    figure = plt.gcf() # get current figure
    figure.set_size_inches(15, 10)
    plt.savefig('STANDARD Stable_ALLO_2Sp_ss_plot Rst --- %g Rend --- %g.png' %(R[0], R[-1]), dpi=1000)
    plt.show()
    plt.close()

    





stable_ss_plot()
#ss_plot()

'''
print("Testing the functions...")
print("At R = ", R)
print("Vstar_noveg: ", Veq_noveg)
print("Wstar_noveg: ", Weq_noveg)
print("Ostar_noveg: ", Oeq_noveg)
print("Gstar_noveg: ", Geq_noveg)
print("\n=====================================\n")
print("Vstar_veg: ", Veq_veg)
print("Wstar_veg: ", Weq_veg)
print("Ostar_veg: ", Oeq_veg)
print("Gstar_veg: ", Geq_veg)
print("\n=====================================\n")
print("Vegstar_coex: ", Veq_coex)
print("Ostar_coex: ", Oeq_coex)
print("Wstar_plus: ", Wstar_plus)
print("Gstar_plus: ", Gstar_plus)
print("\n=====================================\n")
print("Vegstar_coex: ", Veq_coex)
print("Ostar_coex: ", Oeq_coex)
print("Wstar_minus: ", Wstar_minus)
print("Gstar_minus: ", Gstar_minus)
print("\n=====================================\n")

print(c*gmax*Weq_noveg[0]/(Weq_noveg[0] + k1)-d)
print(alpha*k2*(1-W0)*Oeq_noveg[0]/((Veq_noveg[0] + k2)**2) - gmax*Weq_noveg[0]/(Weq_noveg[0] + k1))
print(-gmax*k1*Veq_noveg[0]/(Weq_noveg[0] + k1)**2 - rW)
print(alpha*(Veq_noveg[0] + k2*W0)/(Veq_noveg[0] + k2))
print(ej*aij*Geq_noveg[0]/(1 + aij*hij*Veq_noveg[0])**2)
print(ej*aij*Veq_noveg[0]/(1 + aij*hij*Veq_noveg[0]) - mj)

print(c*gmax*Weq_veg[0]/(Weq_veg[0] + k1)-d)
print(alpha*k2*(1-W0)*Oeq_veg[0]/((Veq_veg[0] + k2)**2) - gmax*Weq_veg[0]/(Weq_veg[0] + k1))
print(-gmax*k1*Veq_veg[0]/(Weq_veg[0] + k1)**2 - rW)
print(alpha*(Veq_veg[0] + k2*W0)/(Veq_veg[0] + k2))
print(ej*aij*Geq_veg[0]/(1 + aij*hij*Veq_veg[0])**2)
print(ej*aij*Veq_veg[0]/(1 + aij*hij*Veq_veg[0]) - mj)

print(c*gmax*Wstar_plus[1]/(Wstar_plus[1] + k1)-d -aij*Gstar_plus[1]/(1 + aij*hij*Veq_coex[1])**2)
print(alpha*k2*(1-W0)*Oeq_coex[1]/((Veq_coex[1] + k2)**2) - gmax*Wstar_plus[1]/(Wstar_plus[1] + k1))
print(-gmax*k1*Veq_coex[1]/(Wstar_plus[1] + k1)**2 - rW)
print(alpha*(Veq_coex[1] + k2*W0)/(Veq_coex[1] + k2))
print(ej*aij*Gstar_plus[1]/(1 + aij*hij*Veq_coex[1])**2)
print(ej*aij*Veq_coex[1]/(1 + aij*hij*Veq_coex[1]) - mj)
'''