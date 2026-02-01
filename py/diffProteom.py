# Used to calculated fraction of P2-sector and P3-sector
# used file (source)
# mori2021absolutecoli.xlsx (Mori, Matteo, et al. "From coarse to fine: the absolute Escherichia coli proteome under diverse growth conditions." Molecular systems biology 17.5 (2021): 1-23.)
# chenhao2023enzyme.xlsx (Wu, Chenhao, et al. "Enzyme expression kinetics by Escherichia coli during transition from rich to minimal media depends on proteome reserves." Nature Microbiology 8.2 (2023): 347-359.)
# pnas.2427091122.sd01.xlsx (Zhu, Manlu, et al. "Distantly related bacteria share a rigid proteome allocation strategy with flexible enzyme kinetics." Proceedings of the National Academy of Sciences 122.18 (2025): e2427091122.)
# tem_proteomic.csv (Knapp, Benjamin D., et al. "Metabolomic rearrangement controls the intrinsic microbial response to temperature changes." bioRxiv (2023).)
# proteomic_diff_medium.csv (Schmidt, Alexander, et al. "The quantitative and condition-dependent Escherichia coli proteome." Nature biotechnology 34.1 (2016): 104-110.)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

allchape = 'nfua  hdeb  tig  skp  htpg  cspf  cspb  cspa  cspe  cspg  cspd  cspc  csph  cspi  clpa  clpb  clpx  hslu  dnak  hsca  hscc  yegd  grol  dnaj  cbpa  djla  hscb  hslo  ibpa  ibpb  grpe  hslr  gros  degp  ftsh  fimc  secb  ppia  ppib  ppid  fkpa  fklb  fkpb  slyd  ppic  sura  trxa  trxc  ybbn  grxa  grxb  grxc  grxd  dsba  dsbb  dsbc  dsbd  ccmg  dsbg'.split('  ')

acsp = 'cspa cspb cspg csda rbfa nusa pnp tig dead reca infb gyra hns'.split(' ')


pep = 'lon  clpp  pppa  gspo  ompt  hyad  hybd  hyci  guaa  puud  spr  nlpc  purf  glms  asnb  yhbo  yajl  pepn  prlc  dcp  ddpx  ptra  pqql  pepb  pept  pepd  allc  gcp  map  pepq  pepp  ypdf  iap  iada  htpx  yggg  ycal  rsep  mepa  ydgd  degq  degs  ptrb  daca  dacd  dacc  pbpg  dacb  ycbz  lexa  umud  prc  sppa  sohb  pepe  ldca  hslv  iaaa  ggt  ydcp  yegq  yhbu  pmba  tldd'.split('  ')

def isinlist(geneL, genes):
    inter = []
    for i in range(len(genes)):
        if type(genes[i]) != float:
            inter.append(genes[i].lower() in geneL)
        else:
            inter.append(False)
    return inter

def mori2021(prosec):
    xlsxfilename = 'mori2021absolutecoli.xlsx'
    ecoli0 = (
        'EV8-AbsoluteMassFractions-1',{
            'ethanol': [('Lib-18', 1/(59/60) * np.log(2))],
            'kanamycin': [('Lib-07', 1/(56/60) * np.log(2))],
            '42°C': [('Lib-01', 1/(47/60) * np.log(2))],
            'Carbon':[
                ('Lib-06', 1/(83/60) * np.log(2)),
                ('Lib-24', 1/(46/60) * np.log(2)),
                ('Lib-25', 1/(61/60) * np.log(2)),
                ('Lib-26', 1/(37/60) * np.log(2)),
                ('Lib-27', 1/(85/60) * np.log(2)),
                ('Lib-28', 1/(50/60) * np.log(2)),
                ('Lib-29', 1/(54/60) * np.log(2)),
                ('Lib-30', 1/(49/60) * np.log(2)),
                ],
            }
        )
    ecoli1 = (
        'EV9-AbsoluteMassFractions-2',{
            'C-lim':[
                ],
            'A-lim':[
                ],
            'R-lim':[
                ]
            })
    sheetname1 = 'EV3-Samples-2'
    datacond = pd.read_excel(xlsxfilename, sheet_name=sheetname1)
    limn = None
    for i in range(len(datacond)):
        datacgi = datacond['Group'][i]
        if type(datacgi) is str and 'lim' in datacgi:
            limn = datacgi[:5]
        if limn is not None:
            ecoli1[1][limn].append((datacond['Sample ID'][i], datacond['Growth rate (1/h)'][i]))
    def readdata(ecolii):
        sheeti, readn = ecolii
        df = pd.read_excel(xlsxfilename, sheet_name=sheeti)
        genei = df['Gene name']
        inter = {}
        for condj in readn:
            subj = readn[condj]
            inter2 = {}
            for tagk, grk in subj:
                inter3 = {}
                for l in prosec:
                    setl = prosec[l]
                    inter3[l] = df[tagk][isinlist(setl, genes=genei)].sum()
                inter3['gr'] = grk
                inter2[tagk] = inter3
            inter[condj] = inter2
        return inter
    data0 = readdata(ecoli0)
    data1 = readdata(ecoli1)
    data = {}
    for i0, i1 in data0.items():
        data[i0] = i1
    for i0, i1 in data1.items():
        data[i0] = i1
    return data

_ncm3722 = ['M1', 'N5', 'N6', 'N7', 'N8', 'P1', 'A26', 'A27', 'A28', 'A29', 'S1', 'A30', 'A31', 'O8', 'A1', 'A2', 'A3_4', 'A3_5', 'A3_6', 'A3_7', 'A3_8',]
def chenhao2023(prosec, strain = 'ncm3722'):
    xlsxfilename = 'chenhao2023enzyme.xlsx'
    expID = ['A32','A33', 'A34', 'A35', 'S2', 'S3', 'S4', 'S5', 'A36', 'A37', 'A38', 'A39',  # 这些是shift
           'O9', 'O10', 'P3', 'A7', 'A8', 'A9']
    if strain == 'ncm3722':
        idgr = {
            'M1':0.96,
            'M2':0.75,
            'M3':0.54,
            'M4':0.33,
            'N5':1.77,
            'N6':1.42,
            'N7':1.58,
            'N8':1.32,
            'P1':0.95,
            'P8':0.73,
            'P9':0.53,
            'P10':0.35,
            'A26':0.67,
            'A27':0.67,
            'A28':0.67,
            'A29':1.45,
            'S1' :1.45,
            'A30':1.17,
            'A31':0.98,
            'O8' :2.19,
            'A1' :0.96,
            'A2' :0.95,
            'A3_4':0.96,
            'A3_5':0.96,
            'A3_6':0.96,
            'A3_7':0.96,
            'A3_8':0.96,}
    elif strain == 'mg1655':
        idgr = {
            'O9':1.85,
            'O10': 0.66,
            }
    masssheet = ['Group #1 protein mass fractions',
             'Group #2 protein mass fractions',
             'Group #3 protein mass fractions']
    rrsheet = ['Group #1 relative error',
               'Group #2 relative error',
               'Group #3 relative error',
        ]
    dataset = [pd.read_excel(xlsxfilename, sheet_name=i) for i in masssheet]
    rrsheet = [pd.read_excel(xlsxfilename, sheet_name=i) for i in rrsheet]
    # assert idgr.keys() & set(expID) == set()
    def readdata(datas, rdata):
        df = datas
        df1 = rdata
        genei = df['Gene name']
        # assert genei == df1['Gene name']
        inter = {}
        remKey = df.keys() & idgr.keys()
        for keyi in remKey:
            gri = idgr[keyi]
            dfkeyi = df[keyi]
            df1keyi = df1[keyi] * dfkeyi
            inter2 = {}
            # inter3 = {}
            for l in prosec:
                setl = prosec[l]
                inter2[l] = dfkeyi[isinlist(setl, genes=genei)].sum()
                err = (df1keyi[isinlist(setl, genes=genei)] ** 2).sum() ** 0.5
                inter2[l] = (dfkeyi[isinlist(setl, genes=genei)].sum(), err)
            inter2['gr'] = gri
            inter[keyi] = inter2
        return inter
    return {j0:j1 for i in range(3) for j0,j1 in readdata(dataset[i], rrsheet[i]).items()}

def zhu2025(prosec):
    xlsxname = 'pnas.2427091122.sd01.xlsx'
    ecoli = {
        'Ecoli batch F':{
            'F12'  :(37,0.92,'glc'),
            'F13'  :(37,1,'glc+glu'),
            'F13-1':(37,1.9,'LB'),
            'F13-2':(37,1.2,'glc+cAA'),
            'F13-3':(37,0.67,'gly'),
            'F13-4':(37,0.42,'man'),
            },
        'Ecoli batch K':{
            'K1':(30,1.205,'LB'),
            'K2':(30,0.915,'glc+cAA'),
            'K3':(30,0.644,'glc'),
            'K4':(30,0.418,'gly'),
            'K5':(30,0.395,'man'),
            'K6':(37,0.949,'glc'),
            }
        }
    dataset = {}
    for i in ecoli:
        dataset[i] = pd.read_excel(xlsxname, sheet_name=i)
    def readdata(datas, ecol):
        df = datas
        genei = df['gene name']
        inter = {}
        for keyi, (temi, gri, medi) in ecol.items():
            inter2 = {}
            for l in prosec:
                setl = prosec[l]
                inter2[l] = df[keyi][isinlist(setl, genes=genei)].sum()
            inter2['gr'] = gri
            inter2['tem'] = temi
            inter2['med'] = medi
            inter[keyi] = inter2
        return inter
    return {j0:j1 for i0,i1 in dataset.items() for j0,j1 in readdata(i1, ecoli[i0]).items()}


def knapp2024(prosec):
    csvname = 'tem_proteomic.csv'
    data = pd.read_csv(csvname)
    allkey = ['LB16_1_norm', 'LB16_2_norm', 'LB25_1_norm','LB25_2_norm','LB30_1_norm','LB30_2_norm','LB37_1_norm','LB37_2_norm','LB43_1_norm','LB43_2_norm',
              'Glucose25_1_norm', 'Glucose25_2_norm', 'Glucose30_1_norm', 'Glucose30_2_norm', 'Glucose37_1_norm', 'Glucose37_2_norm',
              'Glycerol25_1_norm', 'Glycerol25_2_norm', 'Glycerol30_1_norm', 'Glycerol30_2_norm', 'Glycerol37_1_norm', 'Glycerol37_2_norm']
    def readdata():
        df = data
        genei = df['genename']
        inter = {}
        for keyi in allkey:
            inter2 = {}
            for l in prosec:
                setl = prosec[l]
                inter2[l] = df[keyi][isinlist(setl, genes=genei)].sum() / 1e7
            inter[keyi] = inter2
        return inter
    return readdata()

from copy import deepcopy
def schmit2016(prosec, mode=0, errn=False):
    # mode 0 return results for all condition
    # mode 1 only return results for conditions which mu is accessible. mu dict
    # mode 2 only return results for conditions which only carbon source alter, i.e., conditions in keyL
    csvname = 'proteomic_diff_medium.csv'
    cvname = 'schmit2016_cv.xlsx'
    data = pd.read_csv('proteomic_diff_medium.csv',header=2)
    rr = pd.read_excel(cvname, header=2)
    keyL = ['Glucose', 'LB', 'Glycerol + AA','Acetate', 'Fumarate', 'Glucosamine', 'Glycerol', 'Pyruvate','Chemostat µ=0.5', 'Chemostat µ=0.35', 'Chemostat µ=0.20','Chemostat µ=0.12', 'Xylose','Mannose', 'Galactose ', 'Succinate', 'Fructose']
    stress = ['Stationary phase 1 day', 'Stationary phase 3 days','Osmotic-stress glucose', '42°C glucose', 'pH6 glucose']
    mu = {'LB':1.9, 'Glycerol + AA':1.27, 'Acetate': 0.3, 'Fumarate': 0.42, 'Galactose ':0.26, 'Glucose':0.58, 'Glucosamine':0.46, 'Glycerol':0.47, 'Pyruvate':0.4,'Succinate':0.44,'Fructose':0.65, 'Mannose':0.47, 'Xylose':0.55, '42°C glucose':0.66, 'Chemostat µ=0.5':0.5, 'Chemostat µ=0.35':0.35,'Chemostat µ=0.20':0.2, 'Chemostat µ=0.12':0.12}
    mustress = {'pH6 glucose': 0.63, 'Osmotic-stress glucose':0.55, 'Stationary phase 1 day':0.0, 'Stationary phase 3 days':0}
    allkey = keyL + stress
    def readdata():
        if mode == 0:
            almu = deepcopy(mu)
            for medi, gri in mu.items():
                almu[medi] = gri
        elif mode == 1:
            almu = mu
        elif mode == 2:
            almu = deepcopy(mu)
            del almu['42°C glucose']
        else:
            print("error not consider")
        df = data
        dr = rr
        genei = df['Gene']
        genei1 = df['Gene']
        inter = {}
        for keyi, gri in almu.items():
            dfkeyi = df[keyi]
            drkeyi = dr[keyi]
            dfrkeyi = df[keyi] * dr[keyi] * 0.01
            inter2 = {}
            inter3 = {}
            for l in prosec:
                setl = prosec[l]
                setin = isinlist(setl, genes=genei)
                if errn:
                    err = (dfrkeyi[setin] ** 2).sum() ** 0.5 / df[keyi].sum()
                    inter2[l] = (dfkeyi[setin].sum()/df[keyi].sum(), err)
                else:
                    inter2[l] = dfkeyi[setin].sum()/df[keyi].sum()
            inter2['gr'] = gri
            inter[keyi] = inter2
        return inter
    return readdata()






