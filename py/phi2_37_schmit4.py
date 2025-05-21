# Using in S1 Figure.
# Fitting to other temperature growth rate data
# ingraham1958growth.csv mohr1980temperature.csv
# When fitting, some parameters will be fixed to the values obtained from fitting in the main text.
# km1 is the value of km at 37 celsius degree.
# Due to the varying experimental cultivation conditions, it is necessary to fit the parameters k10, the value of k1 at 37 celsius degree, related to nutritional conditions.
# Parameters to be determined in fitting process are k10, er, e2, n, km1
# If saved file end with '_fix', km1 will be fixed in fitting process
# If saved file end with '_fix2', km1, n0, e20 will be fixed in fitting process
# If not any of above situations, no parameter will be fixed.
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

# phi20 = popt[1]
# phi8 = 0.55 + phi20
# phi0 = 0.066
#
# phi88 = phi8 - phi0 - phi20
# Gascon = 8.314
kr0 = 6
# k20 = 1/popt[0]/phi88
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
dg = lambda t, n, e2: - (e2 + dh(n) - dcp(n) * th)/Gascon1 *(1/t - 1/to) + dcp(n)/Gascon1 * np.log(t / to) # - dg(n,T)/RT + dg(n,To)/RTo

k = lambda t, n, e2, km: np.exp(dg(t,n,e2)) * km
# k2 = lambda t, a, e2: k20 * np.exp(- e2/Gascon*(1/(t) - 1/310.15) + a*(t - 310.15))
k1 = lambda t, er, k10: k10 * np.exp(- er/Gascon * (1/ t- 1/310.15))
kr = lambda t, er: kr0* np.exp( - er/Gascon * (1/(t) - 1/310.15))
muparam = lambda t, er, e2, n, k10, km: mu(k1(t, er, k10), kr(t, er), k(t, n, e2, km))

mutem = pd.read_csv('ingraham1958growth.csv')
tem = np.array(list(mutem['tem']))
grt = list(mutem['grt'])
popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [60000, 0.0, 325, 20, km0], maxfev = 100000)
er0, e20, n0, k10, km1= popt1

temt = np.linspace(10, 50, 100) + 273.15
grtt = [muparam(i, er0, e20, n0, k10, km1) for i in temt]

savecsv(([i - 273.15 for i in temt], grtt), 'tem,gr', 'fit2_tem_gr_ingraham1958growth.csv')

mutem = pd.read_csv('mohr1980temperature.csv')
tem = np.array(list(mutem['tem']))
grt = list(mutem['grt'])
popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [60000, 0.0, 325, 20, km0], maxfev = 100000)
er0, e20, n0, k10, km1= popt1

temt = np.linspace(10, 50, 100) + 273.15
grtt = [muparam(i, er0, e20, n0, k10, km1) for i in temt]

savecsv(([i - 273.15 for i in temt], grtt), 'tem,gr', 'fit2_tem_gr_mohr1980temperature.csv')

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


muparam = lambda t, er, k10, n, e2: mu(k1(t, er, k10), kr(t, er), k(t, n, e2, km0))

mutem = pd.read_csv('ingraham1958growth.csv')
tem = np.array(list(mutem['tem']))
grt = list(mutem['grt'])
popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [60000, 20, 325, 0], maxfev = 100000)
er0, k10, n0, e20= popt1

temt = np.linspace(10, 50, 100) + 273.15
grtt = [muparam(i, er0, k10, n0, e20) for i in temt]

savecsv(([i - 273.15 for i in temt], grtt), 'tem,gr', 'fit2_tem_gr_ingraham1958growth_fixed.csv')

mutem = pd.read_csv('mohr1980temperature.csv')
tem = np.array(list(mutem['tem']))
grt = list(mutem['grt'])
popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [60000, 20, 325, 0], maxfev = 100000)
er0, k10, n0, e20= popt1

temt = np.linspace(10, 50, 100) + 273.15
grtt = [muparam(i, er0, k10, n0, e20) for i in temt]

savecsv(([i - 273.15 for i in temt], grtt), 'tem,gr', 'fit2_tem_gr_mohr1980temperature_fixed.csv')



