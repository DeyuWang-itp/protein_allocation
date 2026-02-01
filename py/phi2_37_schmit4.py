# Using in Fig 5(a).
# Fitting to other temperature growth rate data
# ingraham1958growth.csv mohr1980temperature.csv
# generate file
# fit2_tem_gr_ingraham1958growth_fixed2.csv fit2_tem_gr_mohr1980temperature_fixed2.csv
# When fitting, some parameters, e2, n, will be fixed to the values obtained from table 1 in the main text.
# Due to the varying experimental cultivation conditions, it is necessary to fit the parameters k10, the value of k1 at 37 degree, related to nutritional conditions.
# Due to the need to compare fitting results, this file is primarily used to output the calibretion of our model for other growth rate datasets when e2, n is fixed. For a comparison of the fitting results of the model on different datasets under various constraints, please refer to the 'fit.py' file
import matplotlib.pyplot as plt
import pickle
def readfile(filename):
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    return data

def savefile(filename, data):
    with open(filename, 'wb') as f:
        pickle.dump(data, f)

import numpy as np
from scipy.optimize import curve_fit
def savecsv(datatuple, dataname, filename):
    inter = np.column_stack(datatuple)
    np.savetxt(filename, inter, delimiter = ",", header=dataname, comments="")

import pandas as pd

kr0 = 6
phi8 = 0.55
phi0 = 0.066
phi88 = phi8 - phi0
phim = 0.0061
km0 = phim / phi88

mu = lambda k1,kr,k:k1*kr*phi88/(k1+kr) * 1/(((k)**0.5 + (1 + k) ** 0.5) ** 2)
phir = lambda k1,kr,k: k1*phi88/(k1+kr) * (1 + k)**0.5/((k) ** 0.5 + (1 + k) ** 0.5) + phi0
phi1 = lambda k1,kr,k: kr*phi88/(k1+kr) * (1 + k)**0.5/((k) ** 0.5 + (1 + k) ** 0.5)
phi2 = lambda k1,kr,k: phi88* (k) ** 0.5/((k) ** 0.5 + (1 + k) ** 0.5)

dh = lambda n:4.0 * n + 143
ds = lambda n:0.01327 * n + 0.448
dcp = lambda n: 0.048 * n + 0.85
Gascon1 = 8.314 * 1e-3
Gascon = 8.314
th = 373.5
to = 310.15
dg = lambda t, n, e2: - (e2 + dh(n) - dcp(n) * th)/Gascon1 *(1/t - 1/to) + dcp(n)/Gascon1 * np.log(t / to)

k = lambda t, n, e2, km: np.exp(dg(t,n,e2)) * km
k1 = lambda t, er, k10: k10 * np.exp(- er/Gascon * (1/ t- 1/310.15))
kr = lambda t, er: kr0* np.exp( - er/Gascon * (1/(t) - 1/310.15))

e20 = -367.07147305150454
n0 = 558.3790783842443
muparam = lambda t, er, k10: mu(k1(t, er, k10), kr(t, er), k(t, n0, e20, km0))

mutem = pd.read_csv('ingraham1958growth.csv')
tem = np.array(list(mutem['tem']))
grt = list(mutem['grt'])
popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [60000, 20], maxfev = 100000)
er0, k10= popt1

temt = np.linspace(10, 50, 100) + 273.15
grtt = [muparam(i, er0, k10) for i in temt]

savecsv(([i - 273.15 for i in temt], grtt), 'tem,gr', 'fit2_tem_gr_ingraham1958growth_fixed2.csv')


mutem = pd.read_csv('mohr1980temperature.csv')
tem = np.array(list(mutem['tem']))
grt = list(mutem['grt'])
popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [60000, 20], maxfev = 100000)
er0, k10= popt1

temt = np.linspace(10, 50, 100) + 273.15
grtt = [muparam(i, er0, k10) for i in temt]

savecsv(([i - 273.15 for i in temt], grtt), 'tem,gr', 'fit2_tem_gr_mohr1980temperature_fixed2.csv')




