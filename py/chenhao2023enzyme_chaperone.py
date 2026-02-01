# used data from  supplementary file of article "Enzyme expression kinetics by Escherichia coli during transition from rich to minimal media depends"
# show relationship between phi2 and mut
# used mainly in fig 3(a)
# generated file
# chenhao2023enzyme.csv

import pandas as pd
xlsxfilename = 'chenhao2023enzyme.xlsx'

expID = ['A32','A33', 'A34', 'A35', 'S2', 'S3', 'S4', 'S5', 'A36', 'A37', 'A38', 'A39',  # shift
           'O9', 'O10', 'P3', 'A7', 'A8', 'A9'] # mg1655
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
    'A3_8':0.96,
    }
masssheet = ['Group #1 protein mass fractions',
             'Group #2 protein mass fractions',
             'Group #3 protein mass fractions']

allchape = 'nfua  hdeb  tig  skp  htpg  cspf  cspb  cspa  cspe  cspg  cspd  cspc  csph  cspi  clpa  clpb  clpx  hslu  dnak  hsca  hscc  yegd  grol  dnaj  cbpa  djla  hscb  hslo  ibpa  ibpb  grpe  hslr  gros  degp  ftsh  fimc  secb  ppia  ppib  ppid  fkpa  fklb  fkpb  slyd  ppic  sura  trxa  trxc  ybbn  grxa  grxb  grxc  grxd  dsba  dsbb  dsbc  dsbd  ccmg  dsbg'.split('  ')

acsp = 'cspa cspb cspg csda rbfa nusa pnp tig dead reca infb gyra hns'.split(' ')

def returnTotFra(genefraL, allgene):
    tot = 0
    for geni, frai in genefraL:
        if geni.lower() in allgene:
            tot += frai
    return tot

def isinlist(geneL, genes):
    inter = []
    for i in range(len(genes)):
        if type(genes[i]) != float:
            inter.append(genes[i].lower() in geneL)
        else:
            inter.append(False)
    return inter

dataset = [pd.read_excel(xlsxfilename, sheet_name=i) for i in masssheet]

assert idgr.keys() & set(expID) == set()

prosec = {'phi2t':allchape,'phict':acsp}
def readdata(datas):
    df = datas
    genei = df['Gene name']
    inter = {}
    remKey = df.keys() & idgr.keys()
    for keyi in remKey:
        gri = idgr[keyi]
        inter2 = {}
        for l in prosec:
            setl = prosec[l]
            inter2[l] = df[keyi][isinlist(setl, genes=genei)].sum()
        inter2['gr'] = gri
        inter[keyi] = inter2
    return inter

resd = [readdata(i) for i in dataset]
resd = {j:i[j] for i in resd for j in i}


import matplotlib.pyplot as plt

res = {'gr':[],'phi2t':[],'phict':[]}
for i in resd:
    resdi = resd[i]
    for j in resdi:
        res[j].append(resdi[j])

df0 = pd.DataFrame(res)
df0.to_csv('chenhao2023enzyme.csv', index=False, encoding='utf-8')
