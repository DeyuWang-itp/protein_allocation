# data from pnas.2427091122.sd01.xlsx
# get the relation between phi2 and mu under 30 and 37 degree
# generated file
# zhu2025distantly_ecoli_30.csv
# zhu2025distantly_ecoli_37.csv
# used in fig s1, fig 3


import numpy as np
import pandas as pd
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

data = {}
for i in ecoli:
    data[i] = pd.read_excel(xlsxname, sheet_name=i)

allchape = 'nfua  hdeb  tig  skp  htpg  cspf  cspb  cspa  cspe  cspg  cspd  cspc  csph  cspi  clpa  clpb  clpx  hslu  dnak  hsca  hscc  yegd  grol  dnaj  cbpa  djla  hscb  hslo  ibpa  ibpb  grpe  hslr  gros  degp  ftsh  fimc  secb  ppia  ppib  ppid  fkpa  fklb  fkpb  slyd  ppic  sura  trxa  trxc  ybbn  grxa  grxb  grxc  grxd  dsba  dsbb  dsbc  dsbd  ccmg  dsbg'.split('  ')


def isinlist(geneL, genes):
    inter = []
    for i in range(len(genes)):
        if type(genes[i]) != float:
            inter.append(genes[i].lower() in geneL)
        else:
            inter.append(False)
    return inter


def returnTotFra(genefraL, allgene):
    tot = 0
    for geni, frai in genefraL:
        if geni.lower() in allgene:
            tot += frai
    return tot

datae = {}
prosec = {'phi2t':allchape}
for i in data:
    datai = data[i]
    ecolii = ecoli[i]
    genei = datai['gene name']
    inter = {}
    for j in ecolii:
        temj, grj, medj = ecolii[j]
        inter2 = {}
        for k in prosec:
            setk = prosec[k]
            inter2[k] = datai[j][isinlist(setk, genes = genei)].sum()
        inter2['medium'] = medj
        if temj in inter:
            inter[temj][grj] = inter2
        else:
            inter[temj] = {grj:inter2}
    datae[i] = inter

datag = {30:{},37:{}}
for i in datae:
    dataei = datae[i]
    for temj in dataei:
        for key,val in dataei[temj].items():
            datag[temj][key] = val


# two file with different temperature, 30 and 37 degree
dag30 = datag[30]
dag37 = datag[37]

grt0 = 0.92
dag37c = {i:[] for i in dag37[grt0]}
dag37c['grt'] = []
dag30c = {i:[] for i in dag37[grt0]}
dag30c['grt'] = []
for grti in dag37:
    dag37c['grt'].append(grti)
    dag37i = dag37[grti]
    for key, val in dag37i.items():
        dag37c[key].append(val)

for grti in dag30:
    dag30c['grt'].append(grti)
    dag30i = dag30[grti]
    for key, val in dag30i.items():
        dag30c[key].append(val)

df30 = pd.DataFrame(dag30c)
df37 = pd.DataFrame(dag37c)
df30.to_csv('zhu2025distantly_ecoli_30.csv', index=False, encoding='utf-8')
df37.to_csv('zhu2025distantly_ecoli_37.csv', index=False, encoding='utf-8')



