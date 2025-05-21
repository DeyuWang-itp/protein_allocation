# Used in Figure S2.
# used file
# mu_tem.csv
# generated file
# fit4_tem_gr.csv fit4_tem_phir.csv fit4_tem_phi1.csv fit4_tem_phi2.csv fit4_tem_phir_phi0_17.csv
import matplotlib.pyplot as plt

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
km = phim / phi88



dh = lambda n:4.0 * n + 143
ds = lambda n:0.01327 * n + 0.448
dcp = lambda n: 0.048 * n + 0.85
Gascon1 = 8.314 * 1e-3
Gascon = 8.314
th = 373.5
to = 310.15
dg = lambda t, n, e2: - (e2 + dh(n) - dcp(n) * th)/Gascon1 *(1/t - 1/to) + dcp(n)/Gascon1 * np.log(t / to) # - dg(n,T)/RT + dg(n,To)/RTo

k = lambda t, n, e2: np.exp(dg(t,n,e2)) * km
k1 = lambda t, er, k10: k10 * np.exp(- er/Gascon * (1/ t- 1/310.15))
kr = lambda t, er: kr0* np.exp( - er/Gascon * (1/(t) - 1/310.15))


# limit for nr is (21,37)
# such that phi_0 increase from 21 degree, and reach phi0 * 2 as t reach 15 degree
def phi0t(t):
    cond = (t - 273.15 <= 21) * (21 + 273.15 - t) / 6 * phi0 * 1 + phi0
    return cond

mu = lambda k1,kr,k,phi0tt:k1*kr*(phi8 - phi0tt)/(k1+kr) * 1/(((k)**0.5 + (1 + k) ** 0.5) ** 2)
phir = lambda k1,kr,k,phi0tt: k1*(phi8 - phi0tt)/(k1+kr) * (1 + k)**0.5/((k) ** 0.5 + (1 + k) ** 0.5) + phi0tt
phi1 = lambda k1,kr,k,phi0tt: kr*(phi8 - phi0tt)/(k1+kr) * (1 + k)**0.5/((k) ** 0.5 + (1 + k) ** 0.5)
phi2 = lambda k1,kr,k,phi0tt: (phi8 - phi0tt)* (k) ** 0.5/((k) ** 0.5 + (1 + k) ** 0.5)
muparam = lambda t, er, e2, n: mu(k1(t, er, 20), kr(t, er), k(t, n, e2), phi0t(t))

mutem = pd.read_csv('mu_tem.csv')

tem = np.array(list(mutem['tem']))
grt = list(mutem['growth rate'])
popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [60000, 0, 325], maxfev = 100000)
er0, e20, n0 = popt1


rho = 0.76
lacz = lambda t, c: c * phi1(k1(t, er0, 2), kr(t, er0), k(t, n0, e20))/(rho + phir(k1(t, er0, 2), kr(t, er0), k(t, n0, e20)))
temt = np.linspace(10, 50, 100) + 273.15
grtt = [muparam(i, er0, e20, n0) for i in temt]


phirt = [phir(k1(i, er0, 20), kr(i, er0), k(i, n0, e20), phi0t(i)) for i in temt]
phi1t = [phi1(k1(i, er0, 20), kr(i, er0), k(i, n0, e20), phi0t(i)) for i in temt]
phi2t = [phi2(k1(i, er0, 20), kr(i, er0), k(i, n0, e20), phi0t(i)) for i in temt]

savecsv(([i - 273.15 for i in temt],grtt), ('tem,gr'), 'fit4_tem_gr.csv')
savecsv(([i - 273.15 for i in temt],phirt), ('tem,phirt'),'fit4_tem_phir.csv')
savecsv(([i - 273.15 for i in temt],phi1t), ('tem,phi1t'),'fit4_tem_phi1.csv')
savecsv(([i - 273.15 for i in temt],phi2t), ('tem,phi2t'),'fit4_tem_phi2.csv')
savecsv(([i - 273.15 for i in temt],[i + 0.17 - phi0 for i in phirt]), 'tem,phirt','fit4_tem_phir_phi0_17.csv')
