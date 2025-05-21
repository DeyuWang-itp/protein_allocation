# Used in Figure 3,4 and Fig s1 (a,b)
# used file
# mu_tem.csv herendeen1979_tem_phi2_phir.csv benjamin2024_lb_tem_phi2_phir.csv lacz_tem.csv shehata1975_lb_tem_size.csv
# mu_tem.csv is summary of herendeen1979_mu_tem.csv and farewell1998_mu_tem.csv
# generated file
# fit2_tem_k.csv fit2_tem_gr.csv fit2_tem_phir.csv fit2_tem_phi1.csv fit2_tem_phi2.csv fit2_tem_phir_phi0_17.csv fit2_tem_lacz.csv fit2_tem_volume.csv fit2_tem_apf.csv
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
muparam = lambda t, er, e2, n: mu(k1(t, er, 20), kr(t, er), k(t, n, e2))

mutem = pd.read_csv('mu_tem.csv')

tem = np.array(list(mutem['tem']))
grt = list(mutem['growth rate'])
popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [60000, 0.0, 325], maxfev = 100000)
er0, e20, n0 = popt1


rho = 0.76
lacz = lambda t, c: c * phi1(k1(t, er0, 2), kr(t, er0), k(t, n0, e20))/(rho + phir(k1(t, er0, 2), kr(t, er0), k(t, n0, e20)))
temt = np.linspace(10, 50, 100) + 273.15
grtt = [muparam(i, er0, e20, n0) for i in temt]

ktt = [k(i, n0, e20) for i in temt]
phirt = [phir(k1(i, er0, 20), kr(i, er0), k(i, n0, e20)) for i in temt]
phi1t = [phi1(k1(i, er0, 20), kr(i, er0), k(i, n0, e20)) for i in temt]
phi2t = [phi2(k1(i, er0, 20), kr(i, er0), k(i, n0, e20)) for i in temt]

savecsv(([i - 273.15 for i in temt],ktt), ('tem,k'), 'fit2_tem_k.csv')
savecsv(([i - 273.15 for i in temt],grtt), ('tem,gr'), 'fit2_tem_gr.csv')
savecsv(([i - 273.15 for i in temt],phirt), ('tem,phirt'),'fit2_tem_phir.csv')
savecsv(([i - 273.15 for i in temt],phi1t), ('tem,phi1t'),'fit2_tem_phi1.csv')
savecsv(([i - 273.15 for i in temt],phi2t), ('tem,phi2t'),'fit2_tem_phi2.csv')
savecsv(([i - 273.15 for i in temt],[i + 0.17 - phi0 for i in phirt]), 'tem,phirt','fit2_tem_phir_phi0_17.csv')

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
plt.savefig("tem-mu-neidhardt-rich_km.png")
plt.savefig("tem-mu-neidhardt-rich_km.eps",format='eps')

sherrie1979 = pd.read_csv('herendeen1979_tem_phi2_phir.csv')
benjamin2024lb = pd.read_csv('benjamin2024_lb_tem_phi2_phir.csv')

stem = list(sherrie1979['tem'])
sphir = list(sherrie1979['phir'])
sphi2 = list(sherrie1979['phi2'])

btem = list(benjamin2024lb['tem'])
bphir1 = list(benjamin2024lb['phir_1'])
bphir2 = list(benjamin2024lb['phir_2'])
bphi21 = list(benjamin2024lb['phi2_1'])
bphi22 = list(benjamin2024lb['phi2_2'])

plt.clf()
line1 = plt.plot(temt-273.15, phirt, color = 'red', linestyle=':', linewidth=3, label = "$\phi_r$ predicted")
line2 = plt.plot(temt-273.15, phi1t, color = 'blue', linestyle='-.', linewidth=3, label = "$\phi_1$ predicted")
line3 = plt.plot(temt-273.15, np.array(phirt) + 0.17 - phi0, color = 'green', linestyle=':', linewidth=3, label = r"$\phi_r$ predicted with $phi_0$ $\approx$ 0.17")
poin1 = plt.scatter(stem, sphir, facecolors='none', edgecolors = 'red', label="$\phi_r$ scaled from [34]")
poin2 = plt.scatter(btem, bphir1, facecolors = 'none', edgecolors = 'green', label="$\phi_r$ LB 1 in [35]")
poin3 = plt.scatter(btem, bphir2, facecolors = 'none', edgecolors = 'green', marker = '^', label = "$\phi_r$ LB 2 in [35]")
plt.xlabel(r'$T / ^\circ C$')
plt.ylabel(r'$\phi_1,\phi_r$')
plt.legend()
plt.xlim([10,50])
plt.ylim([0, 0.55])
plt.legend(loc='upper left', bbox_to_anchor=(1,1))
plt.savefig("phir_phi1_predict_data.png",bbox_inches='tight')
plt.savefig("phir_phi1_predict_data.eps",format='eps',bbox_inches='tight')


plt.clf()
lind1 = plt.plot(temt - 273.15, phi2t, color='blue', linestyle=':', linewidth = 3, label = r"$\phi_2$ predicted")
poiq1 = plt.scatter(stem, sphi2, facecolors = 'none', edgecolors = 'black', label="$\phi_2$ scaled from [34]")
poiq2 = plt.scatter(btem, bphi21, facecolors = 'none', edgecolors = 'green', label = "$\phi_2$ LB 1 in [35]")
poiq3 = plt.scatter(btem, bphi22, facecolors = 'none', edgecolors = 'green', marker = 's', label = "$\phi_r$ LB 2 in [35]")
plt.xlabel(r'$T / ^\circ C$')
plt.ylabel(r'$\phi_2$')
plt.legend()
plt.xlim([10,50])
plt.ylim([0, 0.35])
plt.legend(loc='upper left', bbox_to_anchor=(1,1))
plt.savefig("phi2_predict_data.png",bbox_inches='tight')
plt.savefig("phi2_predict_data.eps",format='eps')


laztem = pd.read_csv('lacz_tem.csv')
ltem = np.array(laztem['tem'])
llaz = laztem[' lacz']
popt2, pcov2 = curve_fit(lacz, ltem + 273.15, llaz, p0 = [2000])
c = popt2[0]
lazt = [lacz(i, c) for i in temt]

savecsv(([i - 273.15 for i in temt],lazt), ('tem,lacz'), 'fit2_tem_lacz.csv')

plt.clf()
linf1 = plt.plot(temt-273.15, lazt, color='blue', linestyle=':', linewidth = 3, label=r"predicted with $k_1=2$")
poim2 = plt.scatter(ltem, llaz, facecolors='none', edgecolors='black', marker='s', label="data from [38]")
plt.xlabel(r'$T/^\circ C$')
plt.ylabel(r'$\beta$-galactosidase activity Z unit/($OD \times$ ml)')
plt.legend(loc='upper right')
plt.xlim([10, 50])
plt.ylim([0, 3100])
plt.savefig("laz_activate_data.png")
plt.savefig("laz_activate_data.eps",format="eps")

cellsize = pd.read_csv('shehata1975_lb_tem_size.csv')
ctem = cellsize['tem']
csize = cellsize[' mean']
# temp = np.linspace(14,40,100)
temp = np.linspace(14,46,100)
cphirt = np.array([phir(k1(i, er0, 20), kr(i, er0), k(i, n0, e20)) for i in (temp + 273.15)])
cphi1t = np.array([phi1(k1(i, er0, 20), kr(i, er0), k(i, n0, e20)) for i in (temp + 273.15)])
t35 = 308.15
cphir35 = phir(k1(t35, er0, 20), kr(t35, er0), k(t35, n0, e20))
cphi135 = phi1(k1(t35, er0, 20), kr(t35, er0), k(t35, n0, e20))
v35 = 1 / cphi135
v =  1 / cphi1t

savecsv(([i - 273.15 for i in temt],v/v35), ('tem,volume'), 'fit2_tem_volume.csv')
plt.clf()
ling1 = plt.plot(temp, v/v35, color="blue", linestyle=":", linewidth = 3, label=r"predicted scaled volume")
poil2 = plt.scatter(ctem, csize/csize[1], facecolors='none', edgecolors='black', marker='s', label='scaled volume from [39]')
plt.xlabel(r'$T/^\circ C$')
plt.ylabel(r'Volume of cell/V(35 $^\circ C$)')
plt.legend(loc='upper right')
plt.xlim([13, 42])
plt.ylim([0.9,1.5])
plt.savefig("scaled_volume_data.png")
plt.savefig("scaled_volume_data.eps",format="eps")


apf = phi88 * np.array(ktt) / (np.array(phi2t) + np.array(ktt) * phi88)
savecsv(([i - 273.15 for i in temt], apf ), ('temt','apf'), 'fit2_tem_apf.csv')
