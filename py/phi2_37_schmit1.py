# Used in Figure 4,5(b,c), 6
# Utilize growth rate data to calibrate the model, and simultaneously use it to predict protein allocation, as well as the relative volume and relative constitutive succ expression level of bacteria
# used file
# mu_tem.csv marr1963_succ.csv
# mu_tem.csv is summary of herendeen1979_mu_tem.csv and farewell1998_mu_tem.csv
# generated file
# fit2_tem_k.csv fit2_tem_gr.csv fit2_tem_phir.csv fit2_tem_phi1.csv fit2_tem_phi2.csv fit2_tem_phir_phi0_17.csv fit2_tem_volume.csv fit2_tem_apf.csv
# fit2_tem_gr_nok.csv marr1963_succ_rela.csv fit2_tem_lacz_rela.csv
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

k = lambda t, n, e2: np.exp(dg(t,n,e2)) * km
k1 = lambda t, er, k10: k10 * np.exp(- er/Gascon * (1/ t- 1/310.15))
kr = lambda t, er: kr0* np.exp( - er/Gascon * (1/(t) - 1/310.15))
# Assuming  piecewise linear fitting about elongate rate, with activation energy, er1, er2
# der1 = np.exp(-er1/Gascon * (1 / (294.15) - 1 / 310.15))
# kr = lambda t: kr0 * (np.exp( -er1 / Gascon * (1/t - 1/310.15)) * (t > (273.15 + 21))
#                       + der1 * np.exp(-er2/Gascon * (1/(t) - 1/(294.15))) * (t <= (273.15 + 21)))
# muparam = lambda t, er, e2, n: mu(k1(t, er, 20), kr(t), k(t, n, e2))
muparam = lambda t, er, e2, n: mu(k1(t, er, 20), kr(t, er), k(t, n, e2))

mutem = pd.read_csv('mu_tem.csv')

tem = np.array(list(mutem['tem']))
grt = list(mutem['growth rate'])
popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [60000, 0.0, 325], maxfev = 100000)
er0, e20, n0 = popt1


rho = 0.76
k11 = 1.404 # calcu k1 for minimal + succ from k1=2h^-1 for minimal + glc
lacz = lambda t, c: c * phi1(k1(t, er0, k11), kr(t, er0), k(t, n0, e20))/(rho + phir(k1(t, er0, k11), kr(t, er0), k(t, n0, e20)))
temt = np.linspace(10, 50, 100) + 273.15
grtt = [muparam(i, er0, e20, n0) for i in temt]
grnok = [mu(k1(i, er0, 20), kr(i, er0), km) for i in temt]
ktt = [k(i, n0, e20) for i in temt]
phirt = [phir(k1(i, er0, 20), kr(i, er0), k(i, n0, e20)) for i in temt]
phi1t = [phi1(k1(i, er0, 20), kr(i, er0), k(i, n0, e20)) for i in temt]
phi2t = [phi2(k1(i, er0, 20), kr(i, er0), k(i, n0, e20)) for i in temt]

savecsv(([i - 273.15 for i in temt],grnok), ('tem,gr'), 'fit2_tem_gr_nok.csv')
savecsv(([i - 273.15 for i in temt],ktt), ('tem,k'), 'fit2_tem_k.csv')
savecsv(([i - 273.15 for i in temt],grtt), ('tem,gr'), 'fit2_tem_gr.csv')
savecsv(([i - 273.15 for i in temt],phirt), ('tem,phirt'),'fit2_tem_phir.csv')
savecsv(([i - 273.15 for i in temt],phi1t), ('tem,phi1t'),'fit2_tem_phi1.csv')
savecsv(([i - 273.15 for i in temt],phi2t), ('tem,phi2t'),'fit2_tem_phi2.csv')
savecsv(([i - 273.15 for i in temt],[i + 0.17 - phi0 for i in phirt]), 'tem,phirt','fit2_tem_phir_phi0_17.csv')


laztem = pd.read_csv('marr1963_succ.csv')
ltem = np.array(laztem['tem'])
llaz = laztem['lacz']
popt2, pcov2 = curve_fit(lacz, ltem + 273.15, llaz, p0 = [2000])
c = popt2[0]
lazt = [lacz(i, c) for i in temt]
rlazt = [lacz(i, 1)/lacz(30 + 273.15, 1) for i in temt]
rlaze = [i/llaz[3] for i in llaz]
savecsv((ltem, rlaze), ('tem,lacz'), 'marr1963_succ_rela.csv')
savecsv(([i - 273.15 for i in temt],rlazt), ('tem,lacz'), 'fit2_tem_lacz_rela.csv')



temp = np.linspace(14,46,100)
cphirt = np.array([phir(k1(i, er0, 20), kr(i, er0), k(i, n0, e20)) for i in (temp + 273.15)])
cphi1t = np.array([phi1(k1(i, er0, 20), kr(i, er0), k(i, n0, e20)) for i in (temp + 273.15)])
t35 = 308.15
cphir35 = phir(k1(t35, er0, 20), kr(t35, er0), k(t35, n0, e20))
cphi135 = phi1(k1(t35, er0, 20), kr(t35, er0), k(t35, n0, e20))
v35 = 1 / cphi135
v =  1 / cphi1t

savecsv(([i - 273.15 for i in temt],v/v35), ('tem,volume'), 'fit2_tem_volume.csv')


apf = phi88 * np.array(ktt) / (np.array(phi2t) + np.array(ktt) * phi88)
savecsv(([i - 273.15 for i in temt], apf ), ('temt,apf'), 'fit2_tem_apf.csv')
