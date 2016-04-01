from __future__ import division
from math import pi
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.optimize import minimize, curve_fit
import sys

# Python script to fit (nuclear) curves to data. Requires Scipy.
# Made by Joergen E. Midtboe, University of Oslo
# This version 9/11/2015

## Parameters for the nucleus
A = 240          # (residual nucleus)



# == Import data ==
# Import data as matrices using built-in functionality.
data_ocl = np.loadtxt('AA_gSF_240Pu_strength.dat', skiprows=1)
data_experimental1 = np.loadtxt('AA_gSF_239Pu_gurevich_1976_g_abs.txt', skiprows=1)
data_experimental2 = np.loadtxt('AA_gSF_239Pu_moraes_1993_g_abs.txt', skiprows=1)

#### preparation
# raw_input returns the empty string for "enter"
yes = set(['yes','y', 'ye', ''])
no = set(['no','n'])
###  /end preparation 

# Potential HACK:Skip part in the middle
print "Shall the (middel/otherwise defined) part of the OCL_data be skipped? (y/n)"
choice_skip = raw_input().lower()
if choice_skip in yes:
   # HACK: Skip part in the middle
    n_points_skip      = 2   # number of points that should be skipped
    n_point_start_skip = 2   # first point to be skipped, start counting from 1
    n_point_continue   = n_point_start_skip + n_points_skip + 1
    data_ocl_sliced = np.vstack((data_ocl[0:n_point_start_skip,:], data_ocl[n_point_continue:,:]))
    print "choice: skip some data"
    #print data_ocl[0:n_point_start_skip,:]
    #print data_ocl_sliced
elif choice_skip in no:
    print "choice: take all data"
    data_ocl_sliced = data_ocl
else:
   sys.stdout.write("Please respond with 'yes' or 'no'")

# Combine OCL and experimental1 data for this analysis (skip last two points from OCL - outliers)
data_exp_ocl = np.vstack((data_ocl_sliced, data_experimental1, data_experimental2))

# Define some parameters
factor = 8.6737E-08	# commonly used const. factor in mb^(-1) MeV^(-2)

# == Define fitting functions ==
# Define the different types of functions
# that we want to fit the curve to. We typically
# sum together several instances of each.
#
# Note the possibility to add "hacks" to force the function
# away from unphysical parameter values, such as negative 
# values, by adding a huge contribution (in a continuous way)
# when it tries to go in the unphysical directions.
# However, if this is done, we must take care not to accept the
# best-fit value if it lands on the edge of the accepted region,
# since this is no real minimum.

Tmin = 0.0001  # minimum accepted temperature (see "hack")
Tmax = 0.6   # maximum accepted temperature (see "hack")

def SLO(E, E0, Gamma0, sigma0):
	# Special Lorentzian,
	# adapted from Kopecky & Uhl (1989) eq. (2.1)
	f = factor * sigma0 * E * Gamma0**2 / ( (E**2 - E0**2)**2 + E**2 * Gamma0**2 )
	
	# "the hack"
	if sigma0 < 0:
		return f + sigma0**2*1e10
	elif E0 < 0:
		return f + E0**2*1e10
	elif Gamma0 < 0:
		return f + Gamma0**2*1e10

	else:
		return f

def GLO(E, E0, Gamma0, sigma0, T):
	# Generalized Lorentzian,
	# adapted from Kopecky & Uhl (1989) eq. (2.3-2.4)
	Gamma = Gamma0 * (E**2 + 4* pi**2 * T**2) / E0**2
	f1    = (E * Gamma)/( (E**2 - E0**2)**2 + E**2 * Gamma**2 )
	f2    = 0.7 * Gamma0 *  4* pi**2 * T**2 / E0**5

	f = factor * sigma0 * Gamma0 * ( f1 + f2 )
	
    # "the hack"
	if sigma0 < 0:
		return f + sigma0**2*1e10
	elif E0 < 0:
		return f + E0**2*1e10
	elif Gamma0 < 0:
		return f + Gamma0**2*1e10
	elif T < 0:
		return f + T**2*1e10
	elif T > Tmax:
		return f + (T-Tmax)**2*1e10
	elif T < Tmin:
		return f + (T-Tmin)**2*1e10

	else:
		return f

def EGLO(E, E0, Gamma0, sigma0, T):
    # Enhanced Generalized Lorentzian,
    # adapted from
    # [92] S.G. Kadmenskii, V.P. Markushev, and V.I. Furman. Radiative width of neutron reson-
    # ances. Giant dipole resonances. Sov. J. Nucl. Phys., 37:165, 1983.
    # [91] J. Kopecky and R.E. Chrien. Observation of the M 1 giant resonance by resonance
    # averaging in 106 P d. Nuclear Physics A, 468(2):285-300, 1987. ISSN 0375-9474. doi:
    # 10.1016/0375-9474(87)90518-5.
    # and RIPL3 documentation
    # -- but constant temperature dependece!
    epsilon_0 = 4.5    # (MeV); adopted from RIPL: However, potentially depends on model for state density
    
    # also k is adopted from RIPL: However,"depends on model for state density! (they assume Fermi gas!)
    if A<148:
       k = 1.0
    if(A>=148): 
       k = 1. + 0.09*(A-148)*(A-148)*math.exp(-0.18*(A-148))

    Kappa =    k + (1.0-k)   * (E-epsilon_0)/(E0-epsilon_0)
    Kappa_0 =  k + (k-1.)    * (epsilon_0)  /(E0-epsilon_0)
    Gamma_k =  Kappa *   Gamma0 * (E**2 + (2.0 * pi * T)**2) / E0**2
    Gamma_k0 = Kappa_0 * Gamma0 * (2. * pi * T)**2 / E0**2
    denominator = ( E**2 - E0**2 )**2  + E**2 * E0**2
    
    f = factor * sigma0 * Gamma0 * ( (E*Gamma_k) / denominator + 0.7*Gamma_k0 / E0**3 )

    # "the hack"
    if sigma0 < 0:
        return f + sigma0**2*1e10
    elif E0 < 0:
        return f + E0**2*1e10
    elif Gamma0 < 0:
        return f + Gamma0**2*1e10
    elif T < 0:
        return f + T**2*1e10
    elif T > Tmax:
        return f + (T-Tmax)**2*1e10
    elif T < Tmin:
        return f + (T-Tmin)**2*1e10

    else:
        return f

# Define the Function wot all the required parameters

# The E1 strength - in our case the sum of the GDRs
def f_E1(E, E01, Gamma01, sigma01, T, E02, Gamma02, sigma02):
    f_E1 = GLO(E, E01, Gamma01, sigma01, T) + GLO(E, E02, Gamma02, sigma02, T)
    return f_E1

# Define potentially two scissors
def f_scissors(E, E03, Gamma03, sigma03, E04, Gamma04, sigma04):
    return SLO(E, E03, Gamma03, sigma03) + SLO(E, E04, Gamma04, sigma04)

# The M1 strength - Scissors + 1 extra
def f_M1(E, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05):
    f_M1 = f_scissors(E, E03, Gamma03, sigma03, E04, Gamma04, sigma04) + SLO(E, E05, Gamma05, sigma05)
    return f_M1

# Define the combined function that we want to fit to data, with all required parameters.
def f(E, E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05):
	return f_E1(E, E01, Gamma01, sigma01, T, E02, Gamma02, sigma02) + f_M1(E, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05)

# == Do the curve fit ==
# Starting parameters (by-eye fit) 
p0=[
 # omega, Gamma, sigma
 # MeV,   MeV,   mb
	11.3, 3.2, 290, 	# (E)GLO number 1
	14.15, 3, 300,  	# (E)GLO number 2
	0.34, 				# Common (E)GLO temperature in MeV
	1.9, 0.5, 0.68,	    # SLO number 1 (scissors 1)
	2.7, 0.5, 0.35, 	# SLO number 2 (scissors 2)
	7.5, 5.45, 20		# SLO number 3
	]

parameter_names = [
"(E)GLO1_E", "(E)GLO1_gamma", "(E)GLO1_sigma",
"(E)GLO2_E", "(E)GLO2_gamma", "(E)GLO2_sigma",
"(E)GLO_T",
"SLO1_E", "SLO1_gamma", "SLO1_sigma",
"SLO2_E", "SLO2_gamma", "SLO2_sigma",
"SLO3_E", "SLO3_gamma", "SLO3_sigma"
]

  
# Known value and definition of the function that should run though it
known_value_E1 = [6.534, 1.99e-7, 0.65e-7] # E value, y axis value, uncertainty y value
# for 240Pu we know from Chrien1985: F_E1 = (1.99+-0.65) * 1e-7 at Bn
# the error on the y-axsis is chosen deliberatly small such as to give a big weight to the point!
weight_known_value_E1 = 1/known_value_E1[2]**2

#known_value_M1 = [5.671, 0.31e-7, 0.65e-8] # E value, y axis value, uncertainty y value
# for 240Pu we know from Chrien1985: F_E1 = (1.99+-0.65) * 1e-7 at Bn
# the error on the y-axsis is chosen deliberatly small such as to give a big weight to the point!
#weight_known_value_E1 = 1/known_value_E1[2]**2

def error(p):
	E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05 = p

	weight = 1/np.power(data_exp_ocl[:,2],2)

	sum = np.sum( np.power( f(data_exp_ocl[:,0], E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05) - data_exp_ocl[:,1], 2) * weight )
	# here the contraint is added as a function 
	sum += np.power( (f_E1(known_value_E1[0], E01, Gamma01, sigma01, T, E02, Gamma02, sigma02)) - known_value_E1[1], 2) * weight_known_value_E1
	return sum

print error(p0)

result = minimize(error, p0, method="BFGS")
print result


# Run curve fitting algorithm! (Taking uncertainties into account through the argument "sigma")
#popt, pcov = curve_fit(f, data_exp_ocl[:,0], data_exp_ocl[:,1], p0=p0,
#						sigma=data_exp_ocl[:,2], 
#						maxfev=100000)
#print popt 	# Print resulting best-fit parameters.
		   	# Note that we can also investigate the 
		   	# covariance matrix of the fit if desired.

popt = result.x
print popt 	# Print resulting best-fit parameters.
		   	# Note that we can also investigate the 
		   	# covariance matrix of the fit if desired.
pcov = result.hess_inv


print       "Fit results: \t\t Starting value \t Optimal value \t Standard deviation"
for i in range(len(p0)):
    print   (  "Parameter %i %13s: \t\t %.2f \t\t %.2f \t\t %.2f" %(i, parameter_names[i], p0[i], 
popt[i], np.sqrt(pcov[i,i]))   ).expandtabs(6)

n_temp=6
print "\n (E)GLO temp: \t Start \t Tmin \t Tmax \t Opt \t Std"
print (  "\t\t %.2f \t %.2f \t %.2f \t %.2f +- %.2f" %(p0[n_temp], 
Tmin,  Tmax, popt[n_temp], np.sqrt(pcov[n_temp,n_temp]) )   )


# Calculate the resulting B(M1) values

# some constant
#pi     = 3.14159265358979323846;
alpha            = 7.2973525664E-3 # approx. 1/137
value_protonmass = 938.272046 # in MeV/c^2
value_hbarc      = 197.32697 # in MeV fm
factor_BM1       = 9./(32. * pi**2 ) * 1/alpha * (2.* value_protonmass / value_hbarc)**2. *1/10

# (sigma1, Gamma1, omega2, sigma2, Gamma2, omega2) = popt[7:12]
(omega1, Gamma1, sigma1, omega2, Gamma2, sigma2) = popt[7:13]

B1 = factor_BM1 *(sigma1 * Gamma1 / omega1)
B2 = factor_BM1 *(sigma2 * Gamma2 / omega2)
B = B1 + B2
Omega = (omega1 * B1 + omega2 * B2) / (B1+B2)

#print popt[8]
print "Estiamtion of the B(M1) value resulting from the scissors mode"
print "\n Summed B(M1) strength \t B(M1) \t Centroid"
print " \t\t\t mu_n**2 \t MeV"
print (  "\t\t\t %.2f \t %.2f" %(B , Omega) )

# == Plotting ==

#creating subplots of different size
gs = gridspec.GridSpec(2, 2,
                       width_ratios=[4,1.5]
                       #,height_ratios=[1,1]
                       )

# Make x-axis array to plot from
Earray = np.linspace(0,20,800)

# Initialize figure
#plt.figure()
ax = plt.subplot(gs[0:2,0],)

# Plot data points with error bars
# experimental1
plt.semilogy(data_experimental1[:,0], data_experimental1[:,1], 'o', color="blue")
plt.ylim([1e-10, 1e-5]) # Set y-axis limits
ax.errorbar(data_experimental1[:,0], data_experimental1[:,1], yerr=data_experimental1[:,2], fmt='o', color="blue")
# experimental2
plt.semilogy(data_experimental2[:,0], data_experimental2[:,1], 'o', color="green")
plt.ylim([1e-10, 1e-5]) # Set y-axis limits
ax.errorbar(data_experimental2[:,0], data_experimental2[:,1], yerr=data_experimental2[:,2], fmt='o', color="green")

#Plot constraint (known value), for 240Pu this was the f_E1 value by Chrien
plt.semilogy(known_value_E1[0], known_value_E1[1], 'o', color="black")
plt.ylim([1e-10, 1e-5]) # Set y-axis limits
ax.errorbar(known_value_E1[0], known_value_E1[1], yerr=known_value_E1[2], fmt='o', color="black")


# OCL
# two different version for the plotting; 1) if data is skipped, 2) if all data is taken
if choice_skip in no:
   ax.errorbar(data_ocl[:,0], data_ocl[:,1], yerr=data_ocl[:,2], fmt='o', color='red')
elif choice_skip in yes:
   ax.errorbar(data_ocl_sliced[:,0], data_ocl_sliced[:,1], yerr=data_ocl_sliced[:,2], fmt='o', color='red')
   ax.errorbar(data_ocl[n_point_start_skip+1:n_point_continue,0], data_ocl[n_point_start_skip+1:n_point_continue,1], yerr=data_ocl[n_point_start_skip+1:n_point_continue,2], fmt='o', color='#FFC0CB')

# Make the plot stick around to wait for more
plt.hold('on')

# Plot fit curve(s):
# Initial by-eye fit
E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05 = p0
farray_0 = f(Earray, E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05)
line1, = plt.plot(Earray, farray_0, color="magenta", label="line 1")
# Actual best-fit curve
E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05 = popt
farray_opt = f(Earray, E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05)
line2, = plt.plot(Earray, farray_opt, color="cyan", label="line 2")

# Plot fitted partial spectra
line3, = plt.plot(Earray, GLO(Earray, E01, Gamma01, sigma01, T), '--', color="grey")
plt.plot(Earray, GLO(Earray, E02, Gamma02, sigma02, T), '--', color="grey")
plt.plot(Earray, SLO(Earray, E03, Gamma03, sigma03), '--', color="grey")
plt.plot(Earray, SLO(Earray, E04, Gamma04, sigma04), '--', color="grey")
plt.plot(Earray, SLO(Earray, E05, Gamma05, sigma05), '--', color="grey")


# ==LEGEND==
# Create a legend for the first line.
plt.legend([line1,line2,line3], ["input estimate","optimized","optimized SLO/(E)GLOs"], loc=4)

# Show the whole thing
#plt.show()

#Plotting the covariance matrixces
plt.subplot(gs[0,1])

plt.imshow(pcov,interpolation="nearest")
plt.colorbar()
#plt.show()


### Scissors####
# Calculate the "Scissors more plot" ("background substracted")
def f_bg(E, E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05):
    f_bg = f(Energy, E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05) \
           - f_scissors(Energy, E03, Gamma03, sigma03, E04, Gamma04, sigma04)
    return f_bg


for i in range(len(data_ocl[:,0])):
    E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05 = popt
    #substract the "bg" by f_M1 and f_E1 apart from scissors
    Energy = data_ocl[i,0]
    #print (  "Test before %.2f \t %.2f" %(data_ocl[i,0] , data_ocl[i,1]*10**7) )
    data_ocl[i,1] = data_ocl[i,1] \
                    - f_bg(Energy, E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05)
    #calculate the y-axis staring value

for i in range(len(data_ocl_sliced[:,0])):
    E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05 = popt
    #substract the "bg" by f_M1 and f_E1 apart from scissors
    Energy = data_ocl_sliced[i,0]
    data_ocl_sliced[i,1] = data_ocl_sliced[i,1] \
                           - f_bg(Energy, E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05)

    #print (  "Test after %.2f \t %.2f" %(data_ocl[i,0] , data_ocl[i,1]*10**7) )

#### Plot #######

ax = plt.subplot(gs[1,1])

# OCL
# two different version for the plotting; 1) if data is skipped, 2) if all data is taken
#plt.semilogy(data_ocl[:,0], data_ocl[:,1], 'o', color="grey")
ax.set_yscale("log", nonposy='clip')

plt.ylim([1e-9, 1e-6]) # Set y-axis limits
if choice_skip in no:
   ax.errorbar(data_ocl[:,0], data_ocl[:,1], yerr=data_ocl[:,2], fmt='o', color='red')
elif choice_skip in yes:
   ax.errorbar(data_ocl_sliced[:,0], data_ocl_sliced[:,1], yerr=data_ocl_sliced[:,2], fmt='o', color='red')
   ax.errorbar(data_ocl[n_point_start_skip+1:n_point_continue,0], data_ocl[n_point_start_skip+1:n_point_continue,1], yerr=data_ocl[n_point_start_skip+1:n_point_continue,2], fmt='o', color='#FFC0CB')

# Make x-axis array to plot from
Earray = np.linspace(0,5,800)

# Plot fit curve(s):
# Actual best-fit curve
E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05 = popt
f_scissors_opt = f_scissors(Earray, E03, Gamma03, sigma03, E04, Gamma04, sigma04)
line1, = plt.plot(Earray, f_scissors_opt, color="cyan", label="line 1")

# Plot fitted partial spectra
# line3, = plt.plot(Earray, GLO(Earray, E01, Gamma01, sigma01, T), '--', color="grey")
# plt.plot(Earray, GLO(Earray, E02, Gamma02, sigma02, T), '--', color="grey")
plt.plot(Earray, SLO(Earray, E03, Gamma03, sigma03), '--', color="grey")
plt.plot(Earray, SLO(Earray, E04, Gamma04, sigma04), '--', color="grey")
# plt.plot(Earray, SLO(Earray, E05, Gamma05, sigma05), '--', color="grey")

plt.show()