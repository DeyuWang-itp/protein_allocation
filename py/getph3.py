# This file is used to generate data in fig s6
# generated file about fraction of P3-sector in different proteomic data
# mori2021_phi3.csv
# schmit2016_phi3.csv
# ph3_knapp_glc.csv
# knapp2024_phi3.csv
# schmit_phi23.csv
# ph3_zhu2025.csv
# prediction about P2-sector and P3-sector when P3 is valid
# model_prefit_phi23.csv
# model_prefit_phi23_offset.csv
# fitresult_with_ph3.csv
# glc_fitresult_with_ph3.csv
# used file in diffProteom.py

from diffProteom import *
prosec = {'phi3':pep}

ph3_mori = mori2021(prosec)
ph3_chenhao = chenhao2023(prosec)
ph3_zhu = zhu2025(prosec)
ph3_knapp = knapp2024(prosec)
ph3_schmit = schmit2016(prosec)

prosec['phi2'] = allchape
ph23_schmit = schmit2016(prosec, mode=2)

import matplotlib.pyplot as plt

schmit = {42:[], 37:[]}
schmit[42].append(ph3_schmit['42°C glucose']['phi3'])
schmit[37].append(ph3_schmit['Glucose']['phi3'])

knapp = {16:[], 25:[], 30:[], 37:[], 43:[]}
knapp[16].append(ph3_knapp['LB16_1_norm']['phi3'])
knapp[16].append(ph3_knapp['LB16_2_norm']['phi3'])
knapp[25].append(ph3_knapp['LB25_1_norm']['phi3'])
knapp[25].append(ph3_knapp['LB25_2_norm']['phi3'])
knapp[30].append(ph3_knapp['LB30_1_norm']['phi3'])
knapp[30].append(ph3_knapp['LB30_2_norm']['phi3'])
knapp[37].append(ph3_knapp['LB37_1_norm']['phi3'])
knapp[37].append(ph3_knapp['LB37_2_norm']['phi3'])
knapp[43].append(ph3_knapp['LB43_1_norm']['phi3'])
knapp[43].append(ph3_knapp['LB43_2_norm']['phi3'])

knapp_glu = {25:[], 30:[], 37:[]}
knapp_glu[25].append(ph3_knapp['Glucose25_1_norm']['phi3'])
knapp_glu[25].append(ph3_knapp['Glucose25_2_norm']['phi3'])
knapp_glu[30].append(ph3_knapp['Glucose30_1_norm']['phi3'])
knapp_glu[30].append(ph3_knapp['Glucose30_2_norm']['phi3'])
knapp_glu[37].append(ph3_knapp['Glucose37_1_norm']['phi3'])
knapp_glu[37].append(ph3_knapp['Glucose37_2_norm']['phi3'])

mori = {42:[], 37:[]}
mori[42].append(ph3_mori['42°C']['Lib-01']['phi3'])
mori[37].append(ph3_mori['R-lim']['A2']['phi3'])
mori[37].append(ph3_mori['R-lim']['H1']['phi3'])
mori[37].append(ph3_mori['R-lim']['H5']['phi3'])


def savecsv(phn, csvname):
    result = {'temp':[], 'phi3':[]}
    for i0, i1 in phn.items():
        for i11 in i1:
            result['temp'].append(i0)
            result['phi3'].append(i11)
    df = pd.DataFrame(result)
    df.to_csv(csvname, index=False, encoding='utf-8')

savecsv(mori, 'mori2021_phi3.csv')
savecsv(schmit, 'schmit2016_phi3.csv')
savecsv(knapp_glu, 'ph3_knapp_glc.csv')
savecsv(knapp, 'knapp2024_phi3.csv')

def mycine(k0, k1, k3, phi8, phim):
    q = (1 + k0/k3)/(1 + k0/k1)
    qr1 = (1 + k0/k1)
    qr3 = (1 + k0/k3)
    q = qr3 / qr1
    pm = phim/phi8
    sqpm = (q * pm) ** 0.5
    sqpm1 = (1 + q * pm) ** 0.5
    p2 = phi8 * (sqpm1 - sqpm) * sqpm
    p1 = k0 / (k1 + k0) * phi8 * (sqpm1 - sqpm) ** 2
    p3 = k0 / (k3 + k0) * phi8 * (sqpm1 - sqpm) * sqpm
    pr = phi8 * 1 / qr1 * (sqpm1 - sqpm) * (sqpm1 - sqpm + sqpm / q)
    return (pr, p1, p2, p3, k1 * p1)

kr0 = 6
phi8 = 0.55
phi0 = 0.066
phi88 = phi8 - phi0
phim = 0.0061
# km = phim / phi88


def k3_plot(k3, phim, offset_2 = 0, offset_3 = 0):
    plt.clf()
    k1 = np.linspace(0.1, 50, 1000)
    res = [mycine(kr0, i, k3, phi88, phim) for i in k1]
    plt.scatter([i['gr'] for i in ph23_schmit.values()],[i['phi3'] for i in ph23_schmit.values()])
    plt.scatter([i['gr'] for i in ph23_schmit.values()],[i['phi2'] for i in ph23_schmit.values()])
    plt.plot([i[-1] for i in res], np.array([i[3] for i in res]) + offset_2)
    plt.plot([i[-1] for i in res], np.array([i[2] for i in res]) + offset_3)

inter = {'gr':[], 'phi2':[], 'phi3':[]}
for i in ph23_schmit.values():
    inter['gr'].append(i['gr'])
    inter['phi2'].append(i['phi2'])
    inter['phi3'].append(i['phi3'])

df = pd.DataFrame(inter)
df.to_csv('schmit_phi23.csv', index=False, encoding='utf-8')

k1 = np.linspace(0.1, 50, 100)
res = [mycine(kr0, i, 12.5, phi88, phim) for i in k1]
inter = {'gr':[], 'phi2':[], 'phi3':[]}
for i in res:
    inter['gr'].append(i[-1])
    inter['phi2'].append(i[2])
    inter['phi3'].append(i[3])

df = pd.DataFrame(inter)
df.to_csv('model_prefit_phi23.csv', index=False, encoding='utf-8')

inter = {'gr':[], 'phi2':[], 'phi3':[]}
for i in res:
    inter['gr'].append(i[-1])
    inter['phi2'].append(i[2] + 0.022)
    inter['phi3'].append(i[3] + 0.008)

df = pd.DataFrame(inter)
df.to_csv('model_prefit_phi23_offset.csv', index=False, encoding='utf-8')

dh = lambda n:4.0 * n + 143
ds = lambda n:0.01327 * n + 0.448
dcp = lambda n: 0.048 * n + 0.85
Gascon1 = 8.314 * 1e-3
Gascon = 8.314
th = 373.5
to = 310.15
dg = lambda t, n, e2: - (e2 + dh(n) - dcp(n) * th)/Gascon1 *(1/t - 1/to) + dcp(n)/Gascon1 * np.log(t / to)

km = phim/phi88
k30 = 12.5
k = lambda t, n, e2, km: np.exp(dg(t,n,e2)) * km
k1 = lambda t, er, k10: k10 * np.exp(- er/Gascon * (1/ t- 1/310.15))
kr = lambda t, er: kr0* np.exp( - er/Gascon * (1/(t) - 1/310.15))
k3 = lambda t, er, k30: k30* np.exp( - er/Gascon * (1/(t) - 1/310.15))
mu = lambda k, kr, k1, k3: mycine(kr, k1, k3, phi88, k * phi88)[-1]
muparam = lambda t, e1, er, e2, e3, n, k30, km, k10: mu(k(t, n, e2, km), k1(t, e1, k10), kr(t, er), k3(t, e3, k30))
# muparam = lambda t, er, e2, n, k30, km: mu(k(t, n, e2, km), k1(t, er, 20), kr(t, er), k3(t, er, k30))
mutem = pd.read_csv('mu_tem.csv')

tem = np.array(list(mutem['tem']))
grt = list(mutem['growth rate'])

from scipy.optimize import curve_fit

def fit_with_r2(func, x_data, y_data, p0 = None, bounds=None, maxfev=10000):
    popt, pcov = curve_fit(func, x_data, y_data, p0, bounds = bounds, maxfev=maxfev)
    y_fit = [func(i, *popt) for i in x_data]
    y_exp = np.array(y_data)
    y_fit = np.array(y_fit)
    ym = np.mean(y_exp)
    s_tot = sum((y_exp - ym) ** 2)
    s_res = sum((y_fit - y_exp) ** 2)
    return popt, pcov, 1 - s_res/s_tot

popt1, pcov1, r2 = fit_with_r2(lambda t, er, e2, n: muparam(t, er, er, e2, er, n, k30, km, 20), tem + 273.15, grt, p0 = [6e4, 0.0, 325], maxfev = 100000, bounds =([0,-np.inf, 0],[+np.inf, + np.inf, + np.inf]) )

er0, e20, n0 = popt1
glc_res = [mycine(kr(i, er0), k1(i, er0, 2), k3(i, er0, k30), phi88, k(i, n0, e20, km) * phi88) for i in np.linspace(10,50,100) + 273.15]
ric_res = [mycine(kr(i, er0), k1(i, er0, 20), k3(i, er0, k30), phi88, k(i, n0, e20, km) * phi88) for i in np.linspace(10,50,100) + 273.15]
# er0 77618 +- 3492 j/mol
# e20 -33.117 +- 229.38558 kJ/mol
# n0 326.8695 +- 193.98 kJ/mol

inter = {'gr':[], 'phi2':[], 'phi3':[], 'phi1':[], 'phir':[], 'tem':[]}
temp = np.linspace(10,50,100)
for i in range(100):
    inter['gr'].append(ric_res[i][-1])
    inter['phir'].append(ric_res[i][0] + 0.066)
    inter['phi1'].append(ric_res[i][1])
    inter['phi2'].append(ric_res[i][2] + 0.022)
    inter['phi3'].append(ric_res[i][3] + 0.008)
    inter['tem'].append(temp[i])

df = pd.DataFrame(inter)
df.to_csv('fitresult_with_ph3.csv', index=False, encoding='utf-8')

inter = {'gr':[], 'phi2':[], 'phi3':[], 'phi1':[], 'phir':[], 'tem':[]}
temp = np.linspace(10,50,100)
for i in range(100):
    inter['gr'].append(glc_res[i][-1])
    inter['phir'].append(glc_res[i][0] + 0.066)
    inter['phi1'].append(glc_res[i][1])
    inter['phi2'].append(glc_res[i][2] + 0.022)
    inter['phi3'].append(glc_res[i][3] + 0.008)
    inter['tem'].append(temp[i])

df = pd.DataFrame(inter)
df.to_csv('glc_fitresult_with_ph3.csv', index=False, encoding='utf-8')

inter = {'tem':[], 'phi3':[]}
inter['tem'].append(30)
inter['phi3'].append(ph3_zhu['K1']['phi3'])
inter['tem'].append(37)
inter['phi3'].append(ph3_zhu['F12']['phi3'])
inter['tem'].append(37)
inter['phi3'].append(ph3_zhu['K6']['phi3'])
df = pd.DataFrame(inter)
df.to_csv('ph3_zhu2025.csv', index=False, encoding='utf-8')

