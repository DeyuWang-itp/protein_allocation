# When there is a positive correlation between the folding pressure and the flux rate of the translation flow, calibration of our model and the prediction results for protein allocation
# generated data used in fig s2 and fig3(a)
# generated file
# fit_mu_phi2_3_dataset_bw25113.csv used in fig3(a)
# fit_mu_phi2_3_dataset_ncm3722.csv used in fig3(a)
# fit_press_k_b_pr_with_temp.csv used in figs2

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
phi88 = 0.484
phi2 = lambda alpha, k5, k31: phi88 * np.sqrt(k5 * alpha + k31) / (np.sqrt(k31 + 1) + np.sqrt(k31 + k5 * alpha))
phir = lambda kr, alpha, k5, k31: phi88 /(kr/alpha) * np.sqrt(k31 + 1) / (np.sqrt(k31 + 1) + np.sqrt(k31 + k5 * alpha))
mu = lambda alpha, k5, k31: alpha * phi88 /(np.sqrt(k31 + 1) + np.sqrt(k31 + k5 * alpha))**2

kr0 = 6
k1 = np.linspace(0,20,100)
alpha = k1 * kr0/(kr0 + k1)

data = pd.read_csv('schmit2016_37_gr_phi2.csv')
phi2e = np.array(data['phi2'])
mue = np.array(data['gr'])
# phi2e mue
def loss(params):
    mu = lambda alpha, k5, k31: alpha * phi88 /(np.sqrt(k31 + 1) + np.sqrt(k31 + k5 * alpha))**2
    phi2 = lambda alpha, k5, k31: phi88 * np.sqrt(k5 * alpha + k31) / (np.sqrt(k31 + 1) + np.sqrt(k31 + k5 * alpha))
    k310 = params[0]
    k50  = params[1]
    alpha_pred = params[2:]
    phi2t = phi2(alpha_pred, k50, k310)
    mut = mu(alpha_pred, k50, k310)
    return sum((phi2t - phi2e) ** 2 + (mut - mue) ** 2)

phi88 = 0.484
phim = 0.0061
alphai = mue/phi88 * (np.sqrt(phim) + np.sqrt(phim + 1)) ** 2
alllen = len(alphai) + 2
init_params = np.zeros(alllen)
init_params[2:] = alphai
init_params[0] = phim/phi88
bounds = [(0,phim/phi88),(0,phim/phi88)]
for i in range(alllen - 2):
    bounds.append((0, 5))

result = minimize(loss, init_params, method='L-BFGS-B',bounds=bounds)

k31 = result.x[0]
k5 = result.x[1]

def emf(f1, f2, y1, y2, xlim, init_para, M = 1000, iter_num = 100):
    x_grid = np.linspace(xlim[0], xlim[1], M)
    t_now = init_para
    sigma1_sq = np.var(y1)/2
    sigma2_sq = np.var(y2)/2
    max_iter = iter_num
    tolerance = 1e-5
    leny = len(y1)
    for numi in range(max_iter):
        weights = np.zeros((leny, M))
        for i in range(leny):
            unnorm_weights = np.exp(
                -0.5 * (y1[i] - f1(x_grid, t_now)) **2 /sigma1_sq
                -0.5 * (y2[i] - f2(x_grid, t_now)) **2 /sigma2_sq
                )
            weights[i] = unnorm_weights/np.sum(unnorm_weights)
        def loss_function(t):
            loss = 0.0
            for i in range(leny):
                loss += np.sum(weights[i] * (
                    (y1[i] - f1(x_grid, t)) ** 2 / sigma1_sq
                    + (y2[i] - f2(x_grid,t)) ** 2/ sigma2_sq
                    ))
            return loss
        constraints = (
            {'type':'ineq', 'fun': lambda t:t[0]},
            {'type':'ineq', 'fun': lambda t:t[1]},
            )
        result = minimize(loss_function, t_now, constraints = constraints, method='SLSQP')
        t_new = result.x
        sigma1_sq_new = 0
        sigma2_sq_new = 0
        for i in range(leny):
            weighti = weights[i]
            sigma1_sq_new += np.sum(weighti * (y1[i] - f1(x_grid, t_new)) ** 2)
            sigma2_sq_new += np.sum(weighti * (y2[i] - f2(x_grid, t_new)) ** 2)
        sigma1_sq = sigma1_sq_new / leny
        sigma2_sq = sigma2_sq_new / leny
        if np.linalg.norm(t_new - t_now) < tolerance:
            break
        t_now = t_new
    return t_now, sigma1_sq, sigma2_sq

phi2f = lambda alpha, params: phi2(alpha, params[0], params[1])
muf = lambda alpha, params: mu(alpha, params[0], params[1])

data2 = pd.read_csv('zhu2025distantly_ecoli_37.csv')
mu2 = np.array(data2['grt'])
phi22 = np.array(data2['phi2t'])

data3 = pd.read_csv('zhu2025distantly_ecoli_30.csv')
mu3 = np.array(data3['grt'])
phi23 = np.array(data3['phi2t'])

data4 = pd.read_csv('chenhao2023enzyme.csv')
mu4 = np.array(data4['gr'])
phi24 = np.array(data4['phi2t'])

alphan = np.linspace(0,6,100)

k50 = 0.006523299013662104
k310 = 0.005717987459553773
(k5,k31),_,_ = emf(muf, phi2f, mue, phi2e, (0, 6), (k50,k310), iter_num=1000)

fit_res = {'gr':mu(alphan, k5,k31),'phi2':phi2(alphan, k5, k31)}
fit_df = pd.DataFrame(fit_res)
fit_df.to_csv('fit_mu_phi2_3_dataset_bw25113.csv', index=False, encoding='utf-8')

k50 = k5 *0.36
k310 = k31
(k5,k31),_,_ = emf(muf, phi2f, np.hstack((mu2,mu3,mu4)), np.hstack((phi22,phi23,phi24)), (0, 6), (k50,k310), iter_num=1000)

fit_res = {'gr':mu(alphan, k5,k31),'phi2':phi2(alphan, k5, k31)}
fit_df = pd.DataFrame(fit_res)
fit_df.to_csv('fit_mu_phi2_3_dataset_ncm3722.csv', index=False, encoding='utf-8')

k50 = 0.006523299013662104
k310 = 0.005717987459553773
(k5,k31),_,_ = emf(muf, phi2f, mue, phi2e, (0, 6), (k50,k310), iter_num=1000)


kr0 = 6
phi8 = 0.55
phi0 = 0.066
phi88 = phi8 - phi0
k50,k310 = k5,k31
# new mu phir phi1 phi2
mu = lambda k1,kr,k5,k31: phi88 * k1 * kr/(k1 + kr) / ((k31 + 1) ** 0.5 + (k31 + k5 * (k1 * kr)/(k1 + kr)) ** 0.5) **2
phi2 = lambda k1,kr,k5,k31: phi88 * (k31 + k5 * (k1 * kr)/(k1 + kr)) ** 0.5 / ((k31 + 1) ** 0.5 + (k31 + k5 * (k1 * kr)/(k1 + kr)) ** 0.5)
phir = lambda k1,kr,k5,k31: k1/(kr+k1)* phi88 * (k31 + 1) ** 0.5 / ((k31 + 1) ** 0.5 + (k31 + k5 * (k1 * kr)/(k1 + kr)) ** 0.5) + phi0
phi1 = lambda k1,kr,k5,k31: kr/(kr+k1)* phi88 * (k31 + 1) ** 0.5 / ((k31 + 1) ** 0.5 + (k31 + k5 * (k1 * kr)/(k1 + kr)) ** 0.5)

dh = lambda n:4.0 * n + 143
ds = lambda n:0.01327 * n + 0.448
dcp = lambda n: 0.048 * n + 0.85
Gascon1 = 8.314 * 1e-3
Gascon = 8.314
th = 373.5
to = 310.15
dg = lambda t, n, e2: - (e2 + dh(n) - dcp(n) * th)/Gascon1 *(1/t - 1/to) + dcp(n)/Gascon1 * np.log(t / to)

k5 = lambda t, n, e2: np.exp(dg(t,n,e2)) * k50
k31 = lambda t, n, e2: np.exp(dg(t,n,e2)) * k310
k1 = lambda t, er, k10: k10 * np.exp(- er/Gascon * (1/ t- 1/310.15))
kr = lambda t, er: kr0* np.exp( - er/Gascon * (1/(t) - 1/310.15))
muparam = lambda t,er, n, e2: mu(k1(t,er,20),kr(t,er),k5(t,n,e2),k31(t,n,e2))

mutem = pd.read_csv('mu_tem.csv')
tem = np.array(list(mutem['tem']))
grt = list(mutem['growth rate'])
from scipy.optimize import curve_fit
popt1, pcov1 = curve_fit(muparam, tem + 273.15, grt, p0 = [60000, 325, 0], maxfev = 100000)
er0, n0, e20 = popt1
temt = np.linspace(10,50,100)+273
phirt = [phir(k1(i, er0, 20), kr(i, er0), k5(i, n0, e20),k31(i, n0, e20)) for i in temt]
phi1t = [phi1(k1(i, er0, 20), kr(i, er0), k5(i, n0, e20),k31(i, n0, e20)) for i in temt]
phi2t = [phi2(k1(i, er0, 20), kr(i, er0), k5(i, n0, e20),k31(i, n0, e20)) for i in temt]
mut = [mu(k1(i, er0, 20), kr(i, er0), k5(i, n0, e20),k31(i, n0, e20)) for i in temt]
phirt2 = [phir(k1(i, er0, 20), kr(i, er0), k5(i, n0, e20),k31(i, n0, e20)) + 0.17 - phi0 for i in temt]

res = {
    'tem':temt - 273.15,
    'gr':mut,
    'phirt':phirt,
    'phirt2':phirt2,
    'phi2t':phi2t,
    'phi1t':phi1t,
    }

df = pd.DataFrame(res)
df.to_csv('fit_press_k_b_pr_with_temp.csv', index=False, encoding='utf-8')
