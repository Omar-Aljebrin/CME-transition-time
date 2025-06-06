"""
===============================================================================
-what you are doing:

cure fitting CME data provied by ht.txt from Dr.Zhang which has data of the
CME in the slow rising fasts to quick phase

the gaol is to show an undertanding of this code and that requries
-converting the provide IDL code to python
-undertanding the code 80/20 up to 80% at least (alsos fro it to work)
-able to purpose funcattion to fit that are better at fitting that is 
mimzizing the kai squred value

anchor:
    find diffrent filtering methods to see what are the best
        but all use gradent to find velcoty justs filtering chnages 
  
    

-Refrances
https://www.youtube.com/watch?v=peBOquJ3fDo&t=609s
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
ChatGpt

-------------------------------------------------------------------------------
-Date 
10/27/2024
-Name 
Omar Aljebrin
===============================================================================
"""
# =============================================================================
# GLOBAL VAR - GLOBAL MODS
# -----------------------------------------------------------------------------

from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import timedelta
from scipy.optimize import curve_fit
from scipy.stats import chi2
import scipy.ndimage as ndimage
# -----------------------------------------------------------------------------
# END GLOBAL VAR - GLOBAL MODS
# =============================================================================

# =============================================================================
# FUNCTIONS
# -----------------------------------------------------------------------------
# Define model functions (the inies provided by the orginal code)


def exp_func(x, a, b, h0):
    return a * np.exp(b * x) + h0


def exp_linear_func(x, a, b, d, h0):
    return a * np.exp(b * (x)) + d * x + h0


def power_func(x, a, b, h0):
    return a * (x ** b) + h0


def power_linear_func(x, a, b, c, d, h0):
    return a * ((x - b) ** c) + d * x + h0
# -----------------------------------------------------------------------------
# END FUNCTIONS
# =============================================================================

# =============================================================================
# MAIN SCRIPT
# -----------------------------------------------------------------------------

# importing the data and normailzing it
# 3nneed i do this
# if not keyword_set(file) then file='ht.sav'
# restore,file,/ver
# nstart=0
# nend=n_elements(height)-5              ;;adjustable data range, here 5 is an arbitrary value

# xtime=utc2tai(str2utc(anytim(time[nstart:nend],/vms)))
# xtime=xtime-xtime[0]
# time=time[nstart:nend]
# height=height[nstart:nend]
# height_err=height_err[nstart:nend]
# xtime_err=fltarr(n_elements(height)


# Define file path
dname = '../data'   # one up
fname = 'ht.txt'    # file name with extension
dfname = '{0:s}/{1:s}'.format(dname, fname)  # full file path


# Load data, assuming the first two columns are numbers and the third is a date-time string
data = np.genfromtxt(dfname, skip_header=1, dtype=None, encoding=None)

# Extract height, height error, and date-time
heights = data['f0']           # First column (heights in Km)
height_errors = data['f1']     # Second column (height errors in Km)
date_time_strings = data['f2']  # Third column (date-time strings)


# narrow less it by 5 since that is done in IDL code (but why?)
nstart = 0
nend = len(heights)

heights = heights[nstart:nend]
# Second column (height errors in Km)
height_errors = height_errors[nstart:nend]
date_time_strings = date_time_strings[nstart:nend]


# Convert date-time strings to pandas datetime objects in UTC
date_times = pd.to_datetime(date_time_strings)  # correct date time

# Convert UTC to TAI by adding 37 seconds
tai_times = date_times + timedelta(seconds=37)

# Convert TAI times to seconds since the epoch
# Convert nanoseconds to seconds
tai_seconds = tai_times.astype('int64') // 10**9


# shifting to the starting point

xtime = tai_seconds-tai_seconds[0]
heights = heights-heights[0]


# complie all into one and irnized to input intoi functai easier

x = xtime
y = heights
yerr = height_errors

# smoothing data
y_smooth2 = ndimage.median_filter(y, size=(nend,))
y_smooth3 = ndimage.gaussian_filter(y, sigma=2400)  # the worst
ns = 3
nw = (np.ones(ns))/ns
y_smooth4 = np.convolve(y, nw, mode='vaild')
nw = len(y_smooth4)
xw = x[nstart:nw]

# finding velcity
vel = np.gradient(y, x)
# velg = np.gradient(y_smooth2, x)
# velu = np.gradient(y_smooth3, x)
# velw = np.gradient(y_smooth4, xw)
# ===========================================================================----
# curve fitting the data using the funcations and extarctiinf it's infocmation
# --------------------------------------------------------------------
# popt is the optimal aparemtars that is a b c.... that can bets fit to the data
# while pcov is the error but the thing we want is the diginal elenmyts those
# are error on those vlaues like how confdent taht those values are really the true vlaues
# tested by if vary it how prone is it ti changing the fitting if not much then
# we arent as condfcdent its taht values like form  2-5 i and remain the same idk if 2 or 5 or annythign in between


# 1------------------------------------- fitting exp funcation


# popt1,pcov1=curve_fit(exp_func,
#                     x, #x vaivle indpent
#                     y, #y varibvnel depdnt on x
#                     p0=[1.0,1.0,1.0], #intail gusses this case none so place 1
#                     sigma=yerr , #the erroe for height is sigma
#                     maxfev=5000)  #maxuim intration ot must converge to

# yfit1 = exp_func(x,popt1[0],popt1[1],popt1[2])

# finding kai squred vlaue for the expotianl funcation


# 2 ------------------curve fitting  expontail linear functaioin


popt2, pcov2 = curve_fit(exp_linear_func,  # popt is the best values for the paremtsr a and b and so on
                         #
                         x,  # x vaivle indpent
                         y,  # y varibvnel depdnt on x
                         # intail gusses this case none so place 1
                         p0=[1.0, 0.001, 1.0, 100.0],
                         sigma=yerr,  # the erroe for height is sigma
                         maxfev=5000)  # maxuim intration ot must converge to

yfit2 = exp_linear_func(x, popt2[0], popt2[1], popt2[2], popt2[3])


# finding kai squred vlaue for the expintaila linear funcation


# kai=bestnirm/dof in idl


# 3------------------------------ curve fittung power fubcation


popt3, pcov3 = curve_fit(power_func,  # popt is the best values for the paremtsr a and b and so on
                         #
                         x,  # x vaivle indpent
                         y,  # y varibvnel depdnt on x
                         p0=[1, 1, 1],  # intail gusses this case none so place 1
                         sigma=yerr,  # the erroe for height is sigma
                         maxfev=5000)  # maxuim intration ot must converge to

yfit3 = power_func(x, popt3[0], popt3[1], popt3[2])


# finding kai squred vlaue for the power funcation


# 4    -----------------------  curve fitting power linear fucnation


# popt4,pcov4=curve_fit(power_linear_func,   #popt is the best values for the paremtsr a and b and so on
#                                     #
#                     x, #x vaivle indpent
#                     y, #y varibvnel depdnt on x
#                     p0=[1,1,1,1,1], #intail gusses this case none so place 1
#                     sigma=yerr , #the erroe for height is sigma
#                     maxfev=5000)  #maxuim intration ot must converge to

# yfit4 = power_linear_func(x,popt4[0],popt4[1],popt4[2])


# finding kai squred vlaue for thepower linear funcation

# perror=np.sqrt(np.diag(pcov))
# chi_squared = np.sum(((heights - yfit2) / height_errors) ** 2)
# chi_squared = np.sum(((heights - yfit3) / height_errors) ** 2)

# i thnk thos i wring


# ==============================================================================
# FINDING THE VLEOCITY
# ==============================================================================================


# ============================================================================
# GRAPHING
# ==============================================================================


# plotting velcoty

plt.figure(1)
plt.plot(x, vel, '.')
plt.plot(x, velg, '.')
# plt.plot(x,velu)
plt.plot(xw, velw, '.')
plt.show()


# 0------------------------------plotting ture data


plt.plot(x, y, '.')
# plt.errorbar(x, y, yerr=yerr, fmt='.', capsize=3, capthick=1, label="Data with error bars")
plt.plot(x, y_smooth2, '.')  # meh
# plt.plot(x,y_smooth3, '.')   (THE WORST)
plt.plot(xw, y_smooth4, '.')  # (the best)

# 1------------------------------------- plotting exp funcation


# plt.plot(x,yfit1)


# 2-------------------------------------  expontail linear functaioin


plt.plot(x, yfit2)


# 3------------------------------------- plotting power fubcation


plt.plot(x, yfit3)


# 4------------------------------------- plotting pwoer linear funcation


# plt.plot(x,yfit4)


# 00-------------------------------plot cossutmzation


# 01--------------------------------- outputting teh plot

plt.show()


# -----------------------------------------------------------------------------
# END MAIN SCRIPT
# =============================================================================
