# used to analyse residual and calculate the r2 (coefficient of determination) of the fitting result
# These data is used in fig s5 (a,b) and supplementary table 1
# generate file
# residual_mu_qq.csv
import matplotlib.pyplot as plt

import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
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


def fit_with_r2(func, x_data, y_data, p0 = None, bounds=(-np.inf, +np.inf), maxfev=10000):
    popt, pcov = curve_fit(func, x_data, y_data, p0, bounds = bounds, maxfev=maxfev)
    y_fit = [func(i, *popt) for i in x_data]
    y_exp = np.array(y_data)
    y_fit = np.array(y_fit)
    ym = np.mean(y_exp)
    s_tot = sum((y_exp - ym) ** 2)
    s_res = sum((y_fit - y_exp) ** 2)
    n = len(x_data)
    return popt, pcov, 1 - s_res/s_tot, n* np.log(s_res/n) + 2 * len(popt)

def corr(x_data, y_data):
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    x_mean = np.mean(x_data)
    y_mean = np.mean(y_data)
    return np.mean((x_data - x_mean) * (y_data - y_mean))/np.std(x_data)/np.std(y_data)

def initial_value_feasibility(func, x_data, y_data, p_init, p_bound, n_initializations = 50, rr = 10):
    # used rr-fold of init value as new init value, and limited by p_bound,
    # test the fitting result dependence about init value
    # rr = 10
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    results = []
    p = len(p_init)
    for i in range(n_initializations):
        p_init0 = []
        for pj in range(p):
            pjk = np.random.uniform(low = - rr * p_init[pj], high= rr * p_init[pj], size = 1)[0]
            while pjk < p_bound[0][pj] or pjk > p_bound[1][pj]:
                pjk = np.random.uniform(low = -rr * p_init[pj], high= rr * p_init[pj], size = 1)[0]
            p_init0.append(pjk)
        # breakpoint()
        try:
            popt, pcov = curve_fit(func, x_data, y_data, p0 = p_init0, bounds = p_bound, maxfev=100000)
            results.append({'init':p_init0.copy(), 'final': popt.copy(), 'state': 'success', 'sse': np.sum((y_data - func(x_data, *popt))**2)})
        except Exception as e:
            results.append({'init':p_init0.copy(), 'final': None, 'state': 'failed', 'error': str(e)})
    res_succ = [res for res in results if res['state'] == 'success']
    num = len(res_succ)
    success_rate = num/ n_initializations
    if num > 0:
        final_params = np.array([r['final'] for r in res_succ])
        sse_list = np.array([r['sse'] for r in res_succ])
        param_std = np.std(final_params, axis=0)
        param_cv = param_std / np.abs(np.mean(final_params, axis=0))
        best_idx = np.argmin(sse_list)
        best_params = final_params[best_idx]
        best_sse = sse_list[best_idx]
        sse_diff = np.max(sse_list) - np.min(sse_list)
        relative_se_range = sse_diff / best_sse
        report = {
            'success_rate': success_rate,
            'n_different_solutions': len(np.unique(np.round(final_params, 4), axis=0)),
            'relative_se_range': relative_se_range,
            'best_sse': best_sse,
            'parameter_variation': param_cv,
            'best_param': best_params,
            'all_solutions': final_params
            }
    else:
        report = {'success_rate': 0, 'error': 'no convergence initimal'}
    return report

def normality_tests(residuals, alpha=0.05):
    results = {}
    sw_stat, sw_p = stats.shapiro(residuals)
    results['Shapiro-Wilk'] = {
        'statistic': sw_stat,
        'p_value': sw_p,
        'normality': sw_p > alpha
        }
    ks_stat, ks_p = stats.kstest((residuals - np.mean(residuals))/np.std(residuals),'norm')
    results['Kolmogorov-Smirnov'] = {
        'statistic': ks_stat,
        'p_value': ks_p,
        'normality': ks_p > alpha
        }
    ad_result = stats.anderson(residuals, dist='norm')
    results['Anderson-Darling'] = {
        'statistic': ad_result.statistic,
        'critical_values': ad_result.critical_values,
        'significance_levels': ad_result.significance_level,
        'normality': ad_result.statistic < ad_result.critical_values[2]  # 5%水平
    }
    skewness = stats.skew(residuals)
    kurtosis = stats.kurtosis(residuals)
    skew_stat, skew_p = stats.skewtest(residuals)
    kurt_stat, kurt_p = stats.kurtosistest(residuals)
    results['Moments'] = {
        'skewness': skewness,
        'skewness_p': skew_p,
        'kurtosis': kurtosis,
        'kurtosis_p': kurt_p,
        'is_symmetric': skew_p > alpha,
        'has_normal_tails': kurt_p > alpha
    }
    return results


popt1, pcov1, r21, aic1 = fit_with_r2(muparam, tem + 273.15, grt, p0 = [60000, 0.0, 325], maxfev = 100000)
er1, e21, n1 = popt1
residual1 = np.array(grt) - muparam(np.array(tem) + 273.15, *popt1)
# r21 0.9851499703346154
# aic1 -162.04087805011488
# corr(tem, residual1) = -0.17725004853106385
# initial_value_feasibility(muparam, tem + 273.15, grt, [60000, 0.0, 325], ([0, -np.inf, 0], [+np.inf, + np.inf, + np.inf]))
# normality_tests(residual1)
# ks_test passed, sw-test passed, ad-test not.
# skewness test passed, kurtosis test not.

muparam1 = lambda t, er, e1, e2, n: mu(k1(t, e1, 20), kr(t, er), k(t, n, e2))
popt2, pcov2, r22, aic2 = fit_with_r2(muparam1, tem + 273.15, grt, p0 = [60000, 60000, 0.0, 325], maxfev = 100000)
er2, e12, e22, n2 = popt2
residual2 = np.array(grt) - muparam1(np.array(tem) + 273.15, *popt2)
# r22 0.9851499701178575
# aic2 -160.04087758302822
# corr(tem, residual2) = -0.17730928849455763
# normality_tests(residual2)
# same as above

mutem = pd.read_csv('ingraham1958growth.csv')
tem1 = np.array(list(mutem['tem']))
grt1 = np.array(mutem['grt'])

mutem = pd.read_csv('mohr1980temperature.csv')
tem2 = np.array(list(mutem['tem']))
grt2 = np.array(list(mutem['grt']))
muparam2 = lambda t, er, e1, k10, e2, n: mu(k1(t, e1, k10), kr(t, er), k(t, n, e2))
muparam3 = lambda t, er, k10, e2, n: mu(k1(t, er, k10), kr(t, er), k(t, n, e2))
muparam4 = lambda t, er, k10: mu(k1(t, er, k10), kr(t, er), k(t, n1, e21))
popt31, pcov31, r231, aic31 = fit_with_r2(muparam2, tem1 + 273.15, grt1, p0 = [60000, 60000, 20, 0.0, 325], maxfev = 100000)
popt32, pcov32, r232, aic32 = fit_with_r2(muparam3, tem1 + 273.15, grt1, p0 = [60000, 20, 0.0, 325], maxfev = 100000)
popt33, pcov33, r233, aic33 = fit_with_r2(muparam4, tem1 + 273.15, grt1, p0 = [60000, 20], maxfev=100000)
residual31 = np.array(grt1) - muparam2(np.array(tem1) + 273.15, *popt31)
residual32 = np.array(grt1) - muparam3(np.array(tem1) + 273.15, *popt32)
residual33 = np.array(grt1) - muparam4(np.array(tem1) + 273.15, *popt33)

popt41, pcov41, r241, aic41 = fit_with_r2(muparam2, tem2 + 273.15, grt2, p0 = [60000, 60000, 20, 0.0, 325], maxfev = 100000)
popt42, pcov42, r242, aic42 = fit_with_r2(muparam3, tem2 + 273.15, grt2, p0 = [60000, 20, 0.0, 325], maxfev = 100000)
popt43, pcov43, r243, aic43 = fit_with_r2(muparam4, tem2 + 273.15, grt2, p0 = [60000, 20], maxfev=100000)
residual41 = np.array(grt2) - muparam2(np.array(tem2) + 273.15, *popt41)
residual42 = np.array(grt2) - muparam3(np.array(tem2) + 273.15, *popt42)

popt5, pcov5, r25, aic5 = fit_with_r2(muparam3, tem + 273.15, grt, p0 = [60000, 20, 0, 325], maxfev=100000)
residual5 = np.array(grt) - muparam3(np.array(tem) + 273.15, *popt5)
# variance of fitting results, used in fig s5 (a)
popt1[0], pcov1[0,0] ** 0.5
popt2[0], pcov2[0,0] ** 0.5
popt2[1], pcov2[1,1] ** 0.5

popt31[0], pcov31[0,0] ** 0.5
popt31[1], pcov31[1,1] ** 0.5
popt32[0], pcov32[0,0] ** 0.5

popt41[0], pcov41[0,0] ** 0.5
popt41[1], pcov41[1,1] ** 0.5
popt42[0], pcov42[0,0] ** 0.5


def generate_qq_data(data, alpha=0.05):
    sample_quantiles = np.sort(data)
    n = len(sample_quantiles)
    positions = (np.arange(1, n+1) - 0.5) / n
    mu, sigma = np.mean(data), np.std(data)
    theoretical_quantiles = stats.norm.ppf(positions, loc=mu, scale=sigma)
    se = sigma * stats.norm.ppf((1 - alpha/2)) * np.sqrt(positions * (1 - positions) / n)
    upper = theoretical_quantiles + se
    lower = theoretical_quantiles - se
    return {'theoretical':theoretical_quantiles,
            'sample': sample_quantiles,
            'upper':upper,
            'lower':lower,
        }

# used data in fig s5(b)
daaa = generate_qq_data(residual1)
qq_df = pd.DataFrame(daaa)
qq_df.to_csv('residual_mu_qq.csv', index=False)

import numpy as np
from scipy.optimize import curve_fit

kr0 = 6
phi8 = 0.55
phi0 = 0.066
phi88 = phi8 - phi0

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
dg = lambda t, e2, n: - (e2 + dh(n) - dcp(n) * th)/Gascon1 *(1/t - 1/to) + dcp(n)/Gascon1 * np.log(t / to)

k31_0 = 0.0054225643379801815
k5_0 = 0.006166521010584875
k5 = lambda t, e2, n: np.exp(dg(t,n,e2)) * k5_0
k31 = lambda t, e2, n: np.exp(dg(t,n,e2)) * k31_0
k1 = lambda t, er, k10: k10 * np.exp(- er/Gascon * (1/ t- 1/310.15))
kr = lambda t, er: kr0* np.exp( - er/Gascon * (1/(t) - 1/310.15))
muparam7 = lambda t,er, e2, n: mu(k1(t,er,20),kr(t,er),k5(t,n,e2),k31(t,n,e2))


popt6, pcov6, r26, aic6 = fit_with_r2(muparam7, tem + 273.15, grt, p0 = [60000, 0.0, 325], maxfev = 100000)

mutem = pd.read_csv('ingraham1958growth.csv')
tem1 = np.array(list(mutem['tem']))
grt1 = np.array(mutem['grt'])

mutem = pd.read_csv('mohr1980temperature.csv')
tem2 = np.array(list(mutem['tem']))
grt2 = np.array(list(mutem['grt']))

popt7, pcov7, r27, aic7 = fit_with_r2(muparam7, tem1 + 273.15, grt1, p0 = [60000, 0.0, 325], maxfev = 100000)

popt8, pcov8, r28, aic8 = fit_with_r2(muparam7, tem2 + 273.15, grt2, p0 = [60000, 0.0, 325], maxfev = 100000)
