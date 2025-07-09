"""
===============================================================================
Purpose of this code:

Curve fitting CME data provided by `ht4.txt` from Dr. Zhang, which contains 
data for all three phases of the CME:
    - Slow rise phase
    - Exponential rise phase
    - Propagation phase

References:
- https://www.youtube.com/watch?v=peBOquJ3fDo&t=609s
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
- ChatGPT

Names: Omar, Sam
===============================================================================
"""

# =============================================================================
# GLOBAL VARIABLES & MODULE IMPORTS
# -----------------------------------------------------------------------------

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import timedelta
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from scipy.optimize import least_squares

# -----------------------------------------------------------------------------
# END GLOBAL VARIABLES & MODULE IMPORTS
# =============================================================================

# =============================================================================
# Position vs Time Functions
# =============================================================================

#phase 1 and 2 
def linear_func(x, a, h0):
    return a * x + h0

def exp_func(x, a, b):
    return a * np.exp(b * x)
def exp_linear_func(x, a, b, d, h0):
    return exp_func(x, a, b) + linear_func(x, d, h0)

#phase 2 and 3 

def logstic_func(x, a1, a2, a3):
    return a1 / (1 + np.exp(-a2 * (x - a3)))

def quadratic_func(x, b1, b2, b3):
    return b1 * (x)**2 + b2 * (x) + b3

def logstic_quadratic_func(x, a1, a2, a3, b1, b2, b3):
    return logstic_func(x, a1, a2, a3) + quadratic_func(x, b1, b2, b3)


# ----Function for interaction point between two functions----

#phase 1 and 2 
def exp_m_linear(x, a, b, d, h0):
    return linear_func(x, d, h0) - exp_func(x, a, b)

#phase 2 and 3
def logstic_m_quadratic(x, a, b, x0, b1, b2, b3):
    return logstic_func(x, a, b, x0) - quadratic_func(x, b1, b2, b3)

# =============================================================================
# Velocity vs Time Functions
# =============================================================================

#phase 1 and 2 
def exp_linear_func_dx(x, a, b, d, h0):
    return a * b * np.exp(b * x) + d

def linear_func_dx(x, a, h0):
    return a * ((1 + x) - x)

def exp_func_dx(x, a, b):
    return a * b * np.exp(b * x)

#pashe 2 and 3
def quadratic_func_dx(x, a, b, h0):
    return 2 * a * x + b

def logstic_func_dx(x,a,b,x0):
    return a* ((b * np.exp(-b * (x - x0))) / (1 + np.exp(-b * (x - x0))))

def softplus_quadratic_func(x, p0, p1, p2, p3, p4, p5):
    return p0 * np.log(1 + np.exp(p1 * (x - p2))) + p3 * x**2 + p4 * x + p5

def softplus_quadratic_residuals(p, x, y, yerr):
    model = p[0] * np.log(1 + np.exp(p[1] * (x - p[2]))) + p[3] * x**2 + p[4] * x + p[5]
    residuals = (y - model) / yerr
    return residuals

def logstic_quadratic_func_dx(x, a, b, x0, c, d, h0):
    return a* ((b * np.exp(-b * (x - x0))) / (1 + np.exp(-b * (x - x0)))) -  2 * c * x + d



# --Function for interaction point between two functions--

#pashe 1 and 2 
def exp_m_linear_dx(x, a, b, d, h0):
    return linear_func_dx(x, d, h0) - exp_func_dx(x, a, b)

#pashe 2 and 3 
def logstic_m_quadratic_dx(x, a, b, x0, c, d, h0):
    return logstic_func_dx(x, a, b, x0) - quadratic_func_dx(x, c, d, h0 )

def softplus_quadratic_func_dx(x, p0, p1, p2, p3, p4, p5):
    softplus_derivative = (p0 * p1 * np.exp(p1 * (x - p2))) / (1 + np.exp(p1 * (x - p2)))
    quadratic_derivative = 2 * p3 * x + p4
    return softplus_derivative + quadratic_derivative
def softplus_func_dx(x, p0, p1, p2):
    z = p1 * (x - p2)
    return p0 * p1 * np.exp(z) / (1 + np.exp(z))




# -----------------------------------------------------------------------------
# END FUNCTIONS
# =============================================================================

# =============================================================================
# MAIN SCRIPT
# -----------------------------------------------------------------------------

#------------------------------------------------------------------
# Importing the Data
#-----------------------------------------------------------------
# Define file path
dname = '../samcraig'  # One directory up
fname = 'ht4.txt'  # File name with extension
dfname = f'{dname}/{fname}'  # Full file path

# Load data, assuming the first two columns are numbers and the third is a date-time string
data = np.genfromtxt(dfname, skip_header=1, dtype=None, encoding=None)

# Extract height, height error, and date-time
heights = data['f0']           # First column (heights in km)
height_errors = data['f1']     # Second column (height errors in km)
date_time_strings = data['f2']  # Third column (date-time strings)

nstart = 0
nend = len(heights)

heights = heights[nstart:nend]
height_errors = height_errors[nstart:nend]
date_time_strings = date_time_strings[nstart:nend]

# Convert date-time strings to pandas datetime objects in UTC
date_times = pd.to_datetime(date_time_strings)

# Convert UTC to TAI by adding 37 seconds
tai_times = date_times + timedelta(seconds=37)

# Convert TAI times to seconds since the epoch
tai_seconds = tai_times.astype('int64') // 10**9

# Shift to start at zero
xtime = tai_seconds - tai_seconds[0]
heights = heights - heights[0]

# Organize variables for function inputs
x = xtime/60  #converting seconds into minutes 
y = heights/659700 +1.0  #convert km to solar raduis and plus one to shift to surtafec of the sun
yerr = height_errors/659700 #same ocnevartin for the err as well 

# =============================================================================
# Cutting Up the Data
# =============================================================================

# phase one and pashe 2  region
x12 = x[0:35]    
y12 = y[0:35]
yerr12 = yerr[0:35]

# phase 2 and phase 3  region
x23 = x[1:90]    
y23 = y[1:90]
yerr23 = yerr[1:90] 




# =============================================================================
# Finding Velocity of the Data
# =============================================================================

# Compute velocity using the gradient
vel = np.gradient(y, x)

# =============================================================================
# Curve Fitting
# =============================================================================

# Linear + Exponential Fit (Phase 1-2)
#popt12, pcov12 = curve_fit(exp_linear_func, x12, y12, 
    #                        p0=[1.0, 0.001, 1.0, 1.0], 
 #                           sigma=yerr12, maxfev=5000)

# Logstic + Quadratic Fit (Phase 2-3)
#note here i do the WHOLE pahses togather but will only use points starting in phase 2 
popt23, pcov23 = curve_fit(logstic_quadratic_func, x, y, 
                            p0=[100.0, 10.0, 15.0, 0.001, 0.001, 1.0], 
                            sigma=yerr, maxfev=5000)

# Example loose but realistic parameter bounds (adjust as needed):
lower_bounds = [0.0, 0.0001, 0.0, -1.0, -1.0, -10.0]  # Avoid negative softplus scale or slope
upper_bounds = [100.0, 0.1, 300.0, 0.01, 1.0, 10.0]    # Prevent runaway large values

p0 = [10.0, 0.005, 100.0, 0.0001, 0.01, 1.0]

result = least_squares(softplus_quadratic_residuals, p0, args=(x23, y23, yerr23),
                       bounds=(lower_bounds, upper_bounds), max_nfev=50000)


# =============================================================================
# Finding Heights of the Functions
# =============================================================================
# Generate fine grid for smoother function representation
xf = np.linspace(min(x23), max(x23), 1000)

 

# Phase 1-2
#fit12_e = exp_func(xf, popt12[0], popt12[1])
#yfit12_l = linear_func(xf, popt12[2], popt12[3])

#yfit12 = yfit12_e + yfit12_l

# Phase 2-3
 

def softplus_part(x, p0, p1, p2):
    return p0 * np.log(1 + np.exp(p1 * (x - p2)))

def quadratic_part(x, p3, p4, p5):
    return p3 * x**2 + p4 * x + p5

popt23_softplus = result.x

xf = np.linspace(np.min(x), np.max(x), 1000)

yfit23_softplus = softplus_quadratic_func(xf, *popt23_softplus)
softplus_component = softplus_part(xf, *popt23_softplus[:3])
quadratic_component = quadratic_part(xf, *popt23_softplus[3:])

yfit23_velocity = softplus_quadratic_func_dx(xf, *popt23_softplus)

# =============================================================================
# Finding Velocity of the Functions
# =============================================================================

# Phase 1-2
#vfit12_e = exp_func_dx(xf, popt12[0], popt12[1])
#vfit12_l = linear_func_dx(xf, popt12[2], popt12[3])
#vfit12 = vfit12_e + vfit12_l

# # Phase 2-3
vfit23_lo = logstic_func_dx(xf, popt23[0], popt23[1], popt23[2])
vfit23_q = quadratic_func_dx(xf, popt23[3], popt23[4], popt23[5]) 
softplus_velocity = softplus_func_dx(xf, *popt23_softplus[:3])


vfit23 = vfit23_lo + vfit23_q
    
# =============================================================================
# Finding Transition Times
# =============================================================================

# Phase 1-2 Transition
#x_inter12 = fsolve(lambda x: exp_m_linear_dx(x, *popt12), [12])
#y_inter12 = exp_linear_func(x_inter12, *popt12)
#v_inter12 = exp_linear_func_dx(x_inter12, *popt12)

# Phase 2-3 Transition (need fixing)
x_inter23 = fsolve(lambda x: logstic_m_quadratic_dx(x, *popt23), [1.0])
y_inter23 = logstic_quadratic_func(x_inter23, *popt23)
v_inter23 = logstic_quadratic_func_dx(x_inter23, *popt23)



# =============================================================================
# PLOTTING RESULTS
# =============================================================================

# =============================================================================
# PHASE 1-2: HEIGHT VS. TIME
# =============================================================================

#plt.figure(1)

#plt.title("Phase 1-2 Height", fontsize=14, fontweight="bold")

# True data
#plt.plot(x12, y12, '.', label='True Data')

# Exponential fit
#plt.plot(xf, yfit12_e+1, label='Exponential Fit')

# Linear fit
#plt.plot(xf, yfit12_l, label='Linear Fit')

# Combined fit (Exponential + Linear)
#plt.plot(xf, yfit12, label='Exponential + Linear Fit')

# Transition point
#plt.axvline(x=x_inter12, color='r', linestyle='--', label='Transition Point (1-2)')
#plt.scatter(x_inter12, y_inter12, color='red', zorder=3)

# Customization
#plt.xlim(0, 40)
#plt.ylim(1, 1.3)
#plt.grid()
#plt.legend(loc='upper right')
#plt.xlabel("Time (min)")
#plt.ylabel("Height (R☉)")
#plt.savefig("p12_height_fitting.jpg")

#plt.show()


# =============================================================================
# PHASE 1-2: VELOCITY VS. TIME
# =============================================================================

#plt.figure(2)
#plt.title("Phase 1-2 Velocity", fontsize=14, fontweight="bold")

# True velocity data
#plt.plot(x, vel, '.', label='True Data')

# Exponential fit
#plt.plot(xf, vfit12_e, label='Exponential Fit')

# Linear fit
#plt.plot(xf, vfit12_l, label='Linear Fit')

# Combined fit (Exponential + Linear)
#plt.plot(xf, vfit12, label='Exponential + Linear Fit')

# Transition point
#plt.axvline(x=x_inter12, color='r', linestyle='--', label='Transition Point (1-2)')
#plt.scatter(x_inter12, v_inter12, color='red', zorder=3)

# Customization
#plt.xlim(0, 25)
#plt.ylim(0, .1)
#plt.grid()
#plt.legend(loc='upper right')
#plt.xlabel("Time (min)")
#plt.ylabel("Velocity from Sun Surface (R☉/min)")
# plt.savefig("p12_velocity_fitting.jpg")

#plt.show()


# =============================================================================
# PHASE 2-3: HEIGHT VS. TIME
# =============================================================================

plt.figure(3)
plt.title("Phase 2-3 Height", fontsize=14, fontweight="bold")

# True data (x23 vs y23)
plt.plot(x23, y23, '.', label='True Data')

# Total Softplus + Quadratic Fit (smooth)
plt.plot(xf, yfit23_softplus, label='Softplus + Quadratic Fit​', color='red')

# Softplus component
plt.plot(xf, softplus_component, label='Softplus Component', color='green', linestyle='--')

# Quadratic component
plt.plot(xf, quadratic_component, label='Quadratic Component', color='orange', linestyle='--')

plt.tight_layout()
plt.grid()
plt.xlabel("Time (min)")
plt.ylabel("Height from Sun Surface (R☉)")
plt.legend(loc='upper right')
plt.show()

# Combined fit (Exponential + Quadratic)
#plt.plot(xf, yfit23, label='logstic + Quadratic Fit')

# logistic fit
#plt.plot(xf, yfit23_logs,label='Logtsic Fit')

# Quadratic fit
#plt.plot(xf, yfit23_q, label='Quadratic Fit')




# Transition point (did not get yet need to make sure fitting is right)
#plt.axvline(x=x_inter23, color='r', linestyle='--', label='Transition Point (2-3)')
#plt.scatter(x_inter23, y_inter23, color='red', zorder=3)

# Customization
# plt.xlim(22, 300)
# plt.ylim(0.7, 8)
#plt.legend(loc='upper right')


# plt.savefig("p23_height_fitting.jpg")




# =============================================================================
# PHASE 2-3: VELOCITY VS. TIME
# =============================================================================

plt.figure(4)
plt.title("Phase 2-3 Velocity", fontsize=14, fontweight="bold")

# # True velocity data
plt.plot(x, vel, '.', label='True Data')

plt.plot(xf, yfit23_velocity, color='purple', linestyle='-', label='Velocity (Softplus + Quadratic)')

plt.plot(xf, softplus_velocity, color='green', linestyle='-', label='Velocity (Softplus Only)')

# Quadratic fit
#plt.plot(xf, vfit23_q, label='Quadratic Fit')

# Logstic fit
#plt.plot(xf, vfit23_lo, label='logstic Fit')

# Combined fit (logstic + Quadratic)
#plt.plot(xf, vfit23, label='logstic + Quadratic Fit')

# # Transition point
#plt.axvline(x=x_inter23, color='r', linestyle='--', label='Transition Point (2-3)')
#plt.scatter(x_inter23, v_inter23, color='red', zorder=3)

# # Customization
plt.xlim(0, 100)
plt.ylim(-0.1, .2)
plt.grid()
#plt.legend(loc='upper right')
plt.xlabel("Time (min)")
plt.ylabel("Velocity (R☉/min)")

plt.show()


# =============================================================================
# END OF SCRIPT
# =============================================================================
