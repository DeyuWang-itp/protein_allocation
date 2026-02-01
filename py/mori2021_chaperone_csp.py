# Using proteomic data from mori 2021, we aim to investigate the changes in the mass fraction of chaperone and cold shock protein in the proteome under different conditions at room temperature
# mori 2021, From coarse to fine: the absolute Escherichia coli proteome under diverse growth conditions
# data is mainly used in fig s1
xlsxfilename = 'mori2021absolutecoli.xlsx'
# generated file
# mori2021absolutecoli_calim_phi2_gr_fit.csv
# mori2021absolutecoli_37_42.csv
# mori2021absolutecoli_kanamycin.csv
# mori2021absolutecoli_ethanol.csv
# mori2021absolutecoli_A-lim.csv
# mori2021absolutecoli_R-lim.csv
# mori2021absolutecoli_C-lim.csv

import numpy as np
import pandas as pd
ecoli0 = (
    'EV8-AbsoluteMassFractions-1',{
        'ethanol': [('Lib-18', 1/(59/60) * np.log(2))],
        'kanamycin': [('Lib-07', 1/(56/60) * np.log(2))],
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
        }
    )

sheetname1 = 'EV3-Samples-2'
datacond = pd.read_excel(xlsxfilename, sheet_name=sheetname1)
limn = None
for i in range(len(datacond)):
    datacgi = datacond['Group'][i]
    if type(datacgi) is str and 'lim' in datacgi:
        limn = datacgi[:5]
    if limn is not None:
        ecoli1[1][limn].append((datacond['Sample ID'][i], datacond['Growth rate (1/h)'][i]))

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

prosec = {'phi2t':allchape,'phict':acsp}
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

dataprefix = 'mori2021absolutecoli_'
def savedata(datai,result=False):
    if result is True:
        ress = {}
    for condj in datai:
        dataij = datai[condj]
        res = {
            'tag':[],
            'phi2':[],
            'phic':[],
            'gr':[]
               }
        for tagk in dataij:
            res['tag'].append(tagk)
            dataijk = dataij[tagk]
            res['phi2'].append(dataijk['phi2t'])
            res['phic'].append(dataijk['phict'])
            res['gr'].append(dataijk['gr'])
        df = pd.DataFrame(res)
        df.to_csv(dataprefix + condj + '.csv', index=False, encoding='utf-8')
        if result:
            ress[condj] = res
    if result:
        return ress

resu0 = savedata(data0)
resu1 = savedata(data1,result=True)
k, b=np.polyfit(resu1['A-lim']['gr'] + resu1['C-lim']['gr'], resu1['A-lim']['phi2'] + resu1['C-lim']['phi2'], 1)
df = pd.DataFrame({'gr':[0,1],'phi2':[b,k+b]})
df.to_csv('mori2021absolutecoli_calim_phi2_gr_fit.csv',index=False, encoding='utf-8')

tempglc = [('EV8-AbsoluteMassFractions-1', {42:[('Lib-01',1/(47/60) * np.log(2))]}), ('EV9-AbsoluteMassFractions-2', {37:[('A2',0.98)]})]

data2 = [readdata(i) for i in tempglc]
tempres = {'temp':[],
           'phi2t':[],
           'phict':[],
           'grt':[]}

for tempi in data2:
    for temp in tempi:
        tempdata = tempi[temp]
        for idx in tempdata:
            tempres['temp'].append(temp)
            tempres['phi2t'].append(tempdata[idx]['phi2t'])
            tempres['phict'].append(tempdata[idx]['phict'])
            tempres['grt'].append(tempdata[idx]['gr'])

df = pd.DataFrame(tempres)
df.to_csv('mori2021absolutecoli_37_42.csv',index=False, encoding='utf-8')

