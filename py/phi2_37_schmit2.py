# Used in figure 3(a) and figure s1(c)
# used file
# schmit2016_37_mu_phir_phi2.pkl (gr phir phi2 data from schmidt2016quantitative with no stress)
# generated file
# schmit2016_37_gr_phi2.csv

# Due to changes in the revision, the primary function of this file has been replaced by the 'fold_pressure_with_translation_flux.py' file. Now, this file primarily outputs growth rate and phi2 data for different media at 37 degree
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

