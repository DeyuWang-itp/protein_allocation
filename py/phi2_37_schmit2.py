# Used in figure S4
# used file
# schmit2016_37_mu_phir_phi2.pkl (gr phir phi2 data from schmidt2016quantitative with no stress) mu_tem.csv
# generated file
# schmit2016_37_gr_phi2.csv predict_gr_phi2_37.csv fit3_tem_gr.csv fit3_tem_phi2.csv fit3_tem_phir.csv fit3_tem_phir_phi0_17.csv
import matplotlib.pyplot as plt
import pickle
def readfile(filename):
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    return data

def savecsv(datatuple, dataname, filename):
    inter = np.column_stack(datatuple)
    np.savetxt(filename, inter, delimiter = ",", header=dataname, comments="")

import numpy as np
from scipy.optimize import curve_fit
data = readfile('schmit2016_mu_phir_phi2.pkl')
gr = data['mu']
phi2 = data['phi2']
phi2 = [phi2[i] for i in gr if i != '42°C glucose']
gr = [gr[i] for i in gr if i != '42°C glucose']
a,b,c = np.polyfit(phi2, gr, 2)

def model_func(x,a,b):
    return a* (x - b) ** 2

popt, pcov = curve_fit(model_func, phi2, gr, p0 = [a, -b /(2*a)])
data = np.column_stack((np.array(gr), np.array(phi2)))
np.savetxt('schmit2016_37_gr_phi2.csv', data, delimiter=',',header='gr,phi2', comments="")
grn = np.linspace(0, 2, 100)
phi2n = [(i / popt[0]) ** (1/2) + popt[1] for i in grn]

gr_phi2_37 = np.column_stack((grn, phi2n))
np.savetxt("predict_gr_phi2_37.csv", gr_phi2_37, delimiter=",", header="gr,phi2", comments="")

plt.ion()
lin = plt.plot(grn, phi2n, color = 'red', linestyle=':', linewidth=3, label = "$phi_2$ predicted")
poi = plt.scatter(gr, phi2, facecolors='none', edgecolors='black', label = "$phi_2$ from [24]")
plt.xlabel(r'$\mu / h^{-1}$')
plt.ylabel(r'$\phi_2$')
plt.title('37$^\circ C$, nutrient changed')
plt.legend()
plt.xlim([0,2])
plt.ylim([0,0.1])
plt.savefig("phi_2_nutrient_changed_37_degree.png")
plt.savefig("phi_2_nutrient_changed_37_degree.eps",format='eps')

import pandas as pd

phi20 = popt[1]
phi8 = 0.55 + phi20
phi0 = 0.066

phi88 = phi8 - phi0 - phi20
Gascon = 8.314
kr0 = 6
k20 = 1/popt[0]/phi88

mu=lambda k1, k0, k: phi88 * k0 / ((1 + k0 / k1) ** 0.5 + (k) ** 0.5) ** 2
phir = lambda k1, k0, k: phi88 / ((1 + k0/k1) ** 0.5 + (k) ** 0.5) / ((1 + k0 / k1) ** 0.5) + phi0
phi2 = lambda k1, k0, k: phi88 / (1 + ((1 + k0/k1)/(k)) ** 0.5) + phi20
phi1 = lambda k1, k0, k: k0 / k1 * phi88 / ((1 + k0/k1) ** 0.5 + (k) ** 0.5) / ((1 + k0 / k1) ** 0.5)

dh = lambda n:4.0 * n + 143
ds = lambda n:0.01327 * n + 0.448
dcp = lambda n: 0.048 * n + 0.85
Gascon1 = 8.314 * 1e-3
th = 373.5
to = 310.15
dg = lambda t, n, e2: - (e2 + dh(n) - dcp(n) * th)/Gascon1 *(1/t - 1/to) + dcp(n)/Gascon1 * np.log(t / to) # - dg(n,T)/RT + dg(n,To)/RTo

k2 = lambda t, n, e2: np.exp(dg(t,n,e2)) * k20
k1 = lambda t, er, k10: k10 * np.exp(- er/Gascon * (1/ t- 1/310.15))
kr = lambda t, er: kr0* np.exp( - er/Gascon * (1/(t) - 1/310.15))
muparam = lambda t, er, n, k10, e2: mu(k1(t, er, k10), kr(t, er), k2(t, n, e2))

mutem = pd.read_csv('mu_tem.csv')

tem = np.array(list(mutem['tem']))
grt = list(mutem['growth rate'])
popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [17 * 4.132 * 1000, 325, 11.23, 0], maxfev = 100000)
# popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [17 * 4.132 * 1000, 0.0, 11.23, 0], maxfev = 100000)
er0, a0, k100, b0 = popt1

rho = 0.76
lacz = lambda t, c: c * phi1(k1(t, er0, 2), kr(t, er0), k2(t, a0, b0))/(rho + phir(k1(t, er0, 2), kr(t, er0), k2(t, a0, b0)))
temt = np.linspace(10, 50, 100) + 273.15
grtt = [muparam(i, er0, a0, k100, b0) for i in temt]


phirt = [phir(k1(i, er0, k100), kr(i, er0), k2(i, a0, b0)) for i in temt]
phi1t = [phi1(k1(i, er0, k100), kr(i, er0), k2(i, a0, b0)) for i in temt]
phi2t = [phi2(k1(i, er0, k100), kr(i, er0), k2(i, a0, b0)) for i in temt]

savecsv(([i - 273.15 for i in temt],grtt), ('tem,gr'), 'fit3_tem_gr.csv')
savecsv(([i - 273.15 for i in temt],phirt), ('tem,phirt'),'fit3_tem_phir.csv')
# savecsv(([i - 273.15 for i in temt],phi1t), ('tem,phi1t'),'fit3_tem_phi1.csv')
savecsv(([i - 273.15 for i in temt],phi2t), ('tem,phi2t'),'fit3_tem_phi2.csv')
savecsv(([i - 273.15 for i in temt],[i + 0.17 - phi0 for i in phirt]), 'tem,phirt','fit3_tem_phir_phi0_17.csv')

plt.clf()
index = 15
lin1 = plt.plot(temt-273.15, grtt, color='red', linestyle=':', linewidth=3, label = "$\mu$ predicted")
poi1 = plt.scatter(tem[index:], grt[index:], facecolors='none', edgecolors='green', label="$\mu$ from [33]")
poi2 = plt.scatter(tem[:index], grt[:index], facecolors='none', marker = '^', edgecolors='blue',label="$\mu$ from [34]")
plt.xlabel(r'$T / ^\circ C$')
plt.ylabel(r'$\mu / h^{-1}$')
plt.title('$mu-T$')
plt.legend()
plt.xlim([9,50])
plt.ylim([0, 2])
plt.savefig("tem-mu-neidhardt-rich.png")
plt.savefig("tem-mu-neidhardt-rich.eps",format='eps')


