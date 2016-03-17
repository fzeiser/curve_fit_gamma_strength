from __future__ import division
from math import pi
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import sys

# Python script to fit (nuclear) curves to data. Requires Scipy.
# Made by Joergen E. Midtboe, University of Oslo
# This version 9/11/2015


# == Import data ==
# Import data as matrices using built-in functionality.
data_ocl = np.loadtxt('Output_data_OCL_strength_U234.txt', skiprows=1)
data_berman = np.loadtxt('rsf_gx_233U_Berman_1986_unw.txt', skiprows=1)

# HACK: Skip part in the middle
data_ocl_sliced = np.vstack((data_ocl[0:17,:], data_ocl[24:-3,:]))
print data_ocl[0:17,:]
print data_ocl_sliced

# Combine OCL and Berman data for this analysis (skip last two points from OCL - outliers)
data_berm_ocl = np.vstack((data_ocl_sliced, data_berman))


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
def SLO(E, E0, Gamma0, sigma0):
	# Special Lorentzian,
	# adapted from Kopecky & Uhl (1989) eq. (2.1)
	f = 8.68e-8 * sigma0 * E * Gamma0**2 / ( (E**2 - E0**2)**2 + E**2 * Gamma0**2 )
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
	# adapted from Kopecky & Uhl (1989) eq. (2.2-2.3)
	Gamma = Gamma0 * (E**2 + 4* pi**2 * T**2) / E0**2
	f = 8.68e-8 * sigma0 * E * Gamma0 * Gamma / ( (E**2 - E0**2)**2 + E**2 * Gamma**2 )
	if sigma0 < 0:
		return f + sigma0**2*1e10
	elif E0 < 0:
		return f + E0**2*1e10
	elif Gamma0 < 0:
		return f + Gamma0**2*1e10
	elif T < 0:
		return f + T**2*1e10
	elif T > 0.6:
		return f + (T-0.6)**2*1e10
	elif T < 0.01:
		return f + (T-0.01)**2*1e10
	else:
		return f

# Define the combined function that we want to fit to data, with all required parameters.
def f(E, E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05):
	return GLO(E, E01, Gamma01, sigma01, T) + GLO(E, E02, Gamma02, sigma02, T) + SLO(E, E03, Gamma03, sigma03) + SLO(E, E04, Gamma04, sigma04) + SLO(E, E05, Gamma05, sigma05)


# == Do the curve fit ==
# Starting parameters (by-eye fit) 
p0=[
	11.5, 3.5, 371, 	# GLO number 1
	14, 4.46, 300,  	# GLO number 2
	0.2, 				# Common GLO temperature
	2.55, 1.2, 1.55,	# SLO number 1 (scissors?)
	5.0, 2.5, 7, 		# SLO number 2
	7.5, 2, 10			# SLO number 3
	]
# Run curve fitting algorithm! (Taking uncertainties into account through the argument "sigma")
popt, pcov = curve_fit(f, data_berm_ocl[:,0], data_berm_ocl[:,1], p0=p0,
						sigma=data_berm_ocl[:,2], 
						maxfev=100000)
print popt 	# Print resulting best-fit parameters.
		   	# Note that we can also investigate the 
		   	# covariance matrix of the fit if desired.


# == Plotting ==
# Make x-axis array to plot from
Earray = np.linspace(0,20,800)

# Initialize figure
plt.figure()
ax = plt.subplot(111)

# Plot data points with error bars
# Berman
plt.semilogy(data_berman[:,0], data_berman[:,1], 'o', color="blue")
plt.ylim([1e-10, 1e-5]) # Set y-axis limits
ax.errorbar(data_berman[:,0], data_berman[:,1], yerr=data_berman[:,2], fmt='o', color="blue")
# OCL
ax.errorbar(data_ocl[:,0], data_ocl[:,1], yerr=data_ocl[:,2], fmt='o', color='red')

# Make the plot stick around to wait for more
plt.hold('on')

# Plot fit curve(s):
# Initial by-eye fit
E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05 = p0
farray_0 = f(Earray, E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05)
plt.plot(Earray, farray_0, color="magenta")
# Actual best-fit curve
E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05 = popt
farray_opt = f(Earray, E01, Gamma01, sigma01, E02, Gamma02, sigma02, T, E03, Gamma03, sigma03, E04, Gamma04, sigma04, E05, Gamma05, sigma05)
plt.plot(Earray, farray_opt, color="cyan")

# Plot fitted partial spectra
plt.plot(Earray, GLO(Earray, E01, Gamma01, sigma01, T), '--', color="grey")
plt.plot(Earray, GLO(Earray, E02, Gamma02, sigma02, T), '--', color="grey")
plt.plot(Earray, SLO(Earray, E03, Gamma03, sigma03), '--', color="grey")
plt.plot(Earray, SLO(Earray, E04, Gamma04, sigma04), '--', color="grey")
plt.plot(Earray, SLO(Earray, E05, Gamma05, sigma05), '--', color="grey")

# Show the whole thing
plt.show()