# Fit the data in  farewell1998_elongate_rate.csv
# and generate data file farewell1998_elongate_rate_re.csv and fit_el_r.csv
import pandas as pd
import numpy as np
def savecsv(datatuple, dataname, filename):
    inter = np.column_stack(datatuple)
    np.savetxt(filename, inter, delimiter = ",", header=dataname, comments="")

elrr = pd.read_csv('farewell1998_elongate_rate.csv')

Gascon = 8.314
def elrrf(ti, e1, e2):
    cond = (1/ti - 273.15 <= 21) * (-1 * e2/Gascon *(ti - 1/(273.15 + 21)) + -1 * e1/Gascon * (1/(273.15 + 21) - 1/(273.15 + 37))) + (1/ti - 273.15 > 21) * (-1 * e1/Gascon * (ti - 1/ (273.15 + 37)))
    return cond

ti = list(elrr['ti'])
elr = np.array(list(elrr['elr']))

i37 = -3
logelr = np.log(elr)
logelr -= logelr[i37]


from scipy.optimize import curve_fit
popt, pcov = curve_fit(elrrf, ti, logelr, p0= [60000, 60000], maxfev=100000)

ti21 = 1/(273.15 + 21)
ti50 = 1/(273.15 + 50)
ti0 = 1/(273.15 + 0)
elfit = [(elrrf(ti0, popt[0],popt[1])), (elrrf(ti21,popt[0],popt[1])), (elrrf(ti50, popt[0],popt[1]))]
tifit = [ti0, ti21, ti50]

savecsv((tifit, elfit), 'ti,elfit', 'fit_el_r.csv')
savecsv((ti, np.log(elr/elr[i37])), 'ti,elrr', 'farewell1998_elongate_rate_re.csv')
