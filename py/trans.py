# Used to process various data to obtain the proportion of protein mass for different sector
# used file
# proteomic_diff_medium.csv slice from schimit2016
# herendeen1979.csv  from herendeen1979 table
# tem_proteomic.csv slice from knapp2024metabolite
# mainly used in fig 5(b,c), fig s4
# generate file
# schmit2016groePhi2.csv. The proportion of total mass of groL+groS and the proportion of mass of phi2 in multiple environments
# schmit2016usedRPPhir.csv. Fraction of R-class protein measured  in herendeen1979 and total fraction of all R-class protein in multiple environments
# usedRPPhirLine.csv. Linearship fitted from above data file.
# knapp2024_phi2_groe.csv. Like schmit2016groePhi2.csv but for proteomic data from knapp2024metabolite
# phi2_groe_line.csv. Linearship fitted from datafiles, schmit2016groePhi2.csv and knapp2024_phi2_groe.csv
# knapp2024_usedRP.csv. Like schmit2016usedRPPhir.csv but for proteomic data from knapp2024metabolite
# knapp2024_usedRP_Line.csv. Linearship fitted from knapp2024_usedRP.csv
# herendeen1979_tem_phi2_phir.csv. Transformed data from herendeen1979 where phi2 used linearship in phi2_groe_line.csv，and phir used linearship in usedRPPhirLine.csv and usedRP15PhirLine, which not save
# additional file
# knapp2024_37_gr_phi2.csv, mu and phi2 at 37 degree from knapp2024metabolite
# knapp2024_30_gr_phi2.csv
# knapp2024_25_gr_phi2.csv

# used proteomap functional clusters from https://www.proteomaps.net/download.html
# allchape, trna, tran, ribo, rnap, pep all from these file
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
data = pd.read_csv('proteomic_diff_medium.csv',header=2)

gene = data['Gene']
groelIdx = data.index[gene == 'groL'].tolist()[0]
groesIdx = data.index[gene == 'groS'].tolist()[0]
keyL = ['Glucose', 'LB', 'Glycerol + AA','Acetate', 'Fumarate', 'Glucosamine', 'Glycerol', 'Pyruvate','Chemostat µ=0.5', 'Chemostat µ=0.35', 'Chemostat µ=0.20','Chemostat µ=0.12', 'Xylose','Mannose', 'Galactose ', 'Succinate', 'Fructose']
stress = ['Stationary phase 1 day', 'Stationary phase 3 days','Osmotic-stress glucose', '42°C glucose', 'pH6 glucose']

totalkey = keyL + stress
totalmass = {i: data[i].sum() for i in totalkey}
groefrac = {i:(data[i][groelIdx] + data[i][groesIdx])/totalmass[i] for i in totalkey}

# mu from si table s23 of schmit2016
mu = {'LB':1.9, 'Glycerol + AA':1.27, 'Acetate': 0.3, 'Fumarate': 0.42, 'Galactose ':0.26, 'Glucose':0.58, 'Glucosamine':0.46, 'Glycerol':0.47, 'Pyruvate':0.4,'Succinate':0.44,'Fructose':0.65, 'Mannose':0.47, 'Xylose':0.55, '42°C glucose':0.66, 'Chemostat µ=0.5':0.5, 'Chemostat µ=0.35':0.35,'Chemostat µ=0.20':0.2, 'Chemostat µ=0.12':0.12}

allchape = 'nfua  hdeb  tig  skp  htpg  cspf  cspb  cspa  cspe  cspg  cspd  cspc  csph  cspi  clpa  clpb  clpx  hslu  dnak  hsca  hscc  yegd  grol  dnaj  cbpa  djla  hscb  hslo  ibpa  ibpb  grpe  hslr  gros  degp  ftsh  fimc  secb  ppia  ppib  ppid  fkpa  fklb  fkpb  slyd  ppic  sura  trxa  trxc  ybbn  grxa  grxb  grxc  grxd  dsba  dsbb  dsbc  dsbd  ccmg  dsbg'.split('  ')
def isinlist(geneL, genes = gene):
    inter = []
    for i in range(len(genes)):
        if type(genes[i]) != float:
            inter.append(genes[i].lower() in geneL)
        else:
            inter.append(False)
    return inter

ischape = isinlist(allchape)

chape = {i: data[i][ischape].sum() for i in totalkey}
chapefrac = {i: chape[i]/totalmass[i] for i in totalkey}
coeff = np.polyfit([groefrac[i] for i in totalkey], [chapefrac[i] for i in totalkey], 1)
# coeff = 2.58274216, 0.01984251
groePhi2 = np.column_stack(([groefrac[i] for i in totalkey],[chapefrac[i] for i in totalkey]))
np.savetxt("schmit2016groePhi2.csv", groePhi2, delimiter=",", header="groe,phi2", comments="")


trna = 'trps  epma  args  lyss  lysu  vals  gluq  hiss  tyrs  thrs  phet  pros  leus  glns  cyss  iles  alas  sers  asps  gltx  phes  metg  alax  alat  alav  alau  alaw  argw  argv  argu  argq  argy  argx  argz  asnu  asnv  asnt  asnw  aspv  aspt  aspu  cyst  glnw  glnv  glnu  glnx  gltu  gltw  gltt  gltv  glyt  glyu  glyw  glyy  glyv  glyx  hisr  ilev  ilet  ileu  leux  leut  leuw  leuz  leup  leuq  leuv  leuu  lysz  lyst  lysq  lysv  lysw  lysy  mett  mety  metv  ilex  metz  metw  iley  metu  pheu  phev  prom  prok  prol  sert  serv  serx  serw  seru  thrw  thrt  thru  thrv  trpt  valx  valz  valv  valy  valw  valu  valt  selc  asns  glyq  glys  sela  fmt'.split('  ')
tran = 'tufb  tufa  rmf  ycih  infb  infc  tsf  fusa  prfa  frr  infa  efp  yeip  prfb  prfc  prfh'.split('  ')
ribo = 'rpld  rpsu  rpma  rplx  rplw  rpse  rpll  rpli  rpsk  rpsj  rplu  rpls  rplt  rplb  rplc  rpsl  rpmb  rpsp  rpsg  rpsf  rpsi  rpln  rplo  rplp  rplq  rpmd  rpso  rpsn  rpmg  rpmf  rpmi  rpmh  rpsr  rpsq  rpla  rplk  rplj  rpsa  rpsb  rpsc  rpsd  rpsh  rpsm  rpss  rpst  sra  rple  rplf  rplm  rplr  rplv  rply  rpmc  ykgm  rpme  rpmj  ykgo  rrsd  rrsa  rrsg  rrsb  rrse  rrsc  rrsh  rrlg  rrlc  rrlb  rrla  rrle  rrld  rrlh'.split('  ')
rnap = 'rpoa  rpob  rpoc  rpod  rpos  flia  rpoe  feci  rpoh  rpon  rpoz'.split('  ')
allr = trna + tran + ribo + rnap
isribo = isinlist(ribo)
isallr = isinlist(allr)
allrfrac = {i:data[i][isallr].sum()/totalmass[i] for i in totalkey}
ribofrac = {i:data[i][isribo].sum()/totalmass[i] for i in totalkey}

usedP = 'rpll  rpoa  groe  rpsa  tsf  lyss  fusa  phet  leus  rpob  tufa  args  glys  vals  asps  iles  glns  gltx'.split('  ')
from copy import deepcopy
usedRP = deepcopy(usedP)
usedRP.remove('groe')
for i in usedRP:
    assert i in allr


usedRP15 = deepcopy(usedRP)
usedRP15.remove('rpll')
isusedRP = isinlist(usedRP)
isusedRP15 = isinlist(usedRP15)
usedRPfrac = {i:data[i][isusedRP].sum()/totalmass[i] for i in totalkey}
usedRP15frac = {i:data[i][isusedRP15].sum()/totalmass[i] for i in totalkey}
# usedRPfrac and usedRP15frac both can show significant linear relationship to allrfrac
coeffrp = np.polyfit([usedRPfrac[i] for i in totalkey],[allrfrac[i] for i in totalkey],1)
coeffrp15 = np.polyfit([usedRP15frac[i] for i in totalkey],[allrfrac[i] for i in totalkey],1)
usedRPPhir = np.column_stack(([usedRPfrac[i] for i in totalkey],[allrfrac[i] for i in totalkey]))
np.savetxt("schmit2016usedRPPhir.csv", usedRPPhir, delimiter=",", header="usedRP,phir",comments="")

usedRPPhirLine = np.column_stack(([0,1],[coeffrp[1], coeffrp[0] + coeffrp[1]]))
usedRP15PhirLine = np.column_stack(([0,1],[coeffrp15[1],coeffrp15[0] + coeffrp15[1]]))
np.savetxt("usedRPPhirLine.csv", usedRPPhirLine, delimiter=",", header="usedRP,phir",comments="")

herendeen1979 = pd.read_csv('herendeen1979.csv')
gene1 = herendeen1979['gene']
groeIdx = herendeen1979.index[gene1 == 'groE'].tolist()[0]
htempkey = ['13.5','15','23','30','42','46']
for keyi in htempkey:
    herendeen1979[keyi] = herendeen1979['37'] * herendeen1979[keyi]

htempkey.append('37')
groTemp = {i:herendeen1979[i][groeIdx]/1000 for i in htempkey}

htempkey.remove('15')
isusedRPH = isinlist(usedRP, gene1)
isusedRP15H = isinlist(usedRP15, gene1)
usedRPtemp = {i:herendeen1979[i][isusedRPH].sum()/1000 for i in htempkey}
phirRPtemp = {i:usedRPtemp[i]*coeffrp[0] + coeffrp[1] for i in htempkey}
usedRPtemp['15'] = herendeen1979['15'][isusedRP15H].sum()/1000
phirRPtemp['15'] = usedRPtemp['15'] * coeffrp15[0] + coeffrp15[1]
htempkey.append('15')
# or
htempkey = ['13.5','15','23','30','37','42','46']

temlabel = ['Glucose', '42°C glucose']
np.savetxt('schmit2016_tem.csv', np.column_stack(([37,42],[mu[i] for i in temlabel],[allrfrac[i] for i in temlabel],[chapefrac[i] for i in temlabel])), delimiter=",", header = "tem,gr,phir,phi2",comments="")

tempP = pd.read_csv('tem_proteomic.csv')
gene2 = tempP['genename']
istrna3 = isinlist(trna, gene2)
istran3 = isinlist(tran, gene2)
isribo3 = isinlist(ribo, gene2)
isrnap3 = isinlist(rnap, gene2)
isgroe3 = isinlist(['gros', 'grol'], gene2)
ischape3 = isinlist(allchape, gene2)
usedRP3 = isinlist(usedRP, gene2)
usedRP153 = isinlist(usedRP15, gene2)
isallr3 = isinlist(allr, gene2)
lbkey = ['LB16_1_norm', 'LB16_2_norm', 'LB25_1_norm','LB25_2_norm','LB30_1_norm','LB30_2_norm','LB37_1_norm','LB37_2_norm','LB43_1_norm','LB43_2_norm']
trnafrac = {i:tempP[i][istrna3].sum()/10000000 for i in lbkey}
tranfrac = {i:tempP[i][istran3].sum()/10000000 for i in lbkey}
ribofrac = {i:tempP[i][isribo3].sum()/10000000 for i in lbkey}
rnapfrac = {i:tempP[i][isrnap3].sum()/10000000 for i in lbkey}
phir3frac = {i:(trnafrac[i] + tranfrac[i] + ribofrac[i] + rnapfrac[i]) for i in lbkey}
lbtem = [16,25,30,37,43]
groe3frac = {i: sum(tempP[f'LB{i}_{j}_norm'][isgroe3].sum() for j in range(1,3))/2/10000000 for i in lbtem}
chape3frac = {i: sum(tempP[f'LB{i}_{j}_norm'][ischape3].sum() for j in range(1,3))/2/10000000 for i in lbtem}
glctem = [25,30,37]
alltem = ['LB16','LB25','LB30','LB37','LB43','Glucose25','Glucose30', 'Glucose37', 'Glycerol25', 'Glycerol30', 'Glycerol37']
groe3frac = {i: sum(tempP[i + f'_{j}_norm'][isgroe3].sum() for j in range(1,3))/2/1e7 for i in alltem}
chaper3frac ={i: sum(tempP[i + f'_{j}_norm'][ischape3].sum() for j in range(1,3))/2/1e7 for i in alltem}
usedRP3frac = {i: sum(tempP[i + f'_{j}_norm'][usedRP3].sum() for j in range(1,3))/2/1e7 for i in alltem}
usedRP153frac = {i: sum(tempP[i + f'_{j}_norm'][usedRP153].sum() for j in range(1,3))/2/1e7 for i in alltem}
allr3frac = {i: sum(tempP[i + f'_{j}_norm'][isallr3].sum() for j in range(1,3))/2/1e7 for i in alltem}

knapp2024_phi2_groe = np.column_stack(([groe3frac[i] for i in alltem],[chaper3frac[i] for i in alltem]))
coeffgroe = np.polyfit([groe3frac[i] for i in alltem] + [groefrac[i] for i in totalkey], [chaper3frac[i] for i in alltem] + [chapefrac[i] for i in totalkey], 1)
phi2groeLine = np.column_stack(([0,1],[coeffgroe[1], coeffgroe[1] + coeffgroe[0]]))
np.savetxt("knapp2024_phi2_groe.csv",knapp2024_phi2_groe, delimiter=",", header="groe,phi2", comments="")
np.savetxt("phi2_groe_line.csv",phi2groeLine, delimiter=",", header="groe,phi2", comments="")

knapp2024_usedRP = np.column_stack(([usedRP3frac[i] for i in alltem],[allr3frac[i] for i in alltem]))
np.savetxt("knapp2024_usedRP.csv", knapp2024_usedRP, delimiter=",", header="usedRP,phir", comments="")
coeff3usedRP = np.polyfit([usedRP3frac[i] for i in alltem],[allr3frac[i] for i in alltem],1)
knapp2024_usedRP_Line = np.column_stack(([0,1],[coeff3usedRP[1], coeff3usedRP[1] + coeff3usedRP[0]]))
np.savetxt("knapp2024_usedRP_Line.csv", knapp2024_usedRP_Line, delimiter=",", header="usedRP,phir", comments="")

phi2Temp = {i:groTemp[i]*coeffgroe[0] + coeffgroe[1] for i in htempkey}
herendeen1979_phi2_phir = np.column_stack(([float(i) for i in htempkey],[phi2Temp[i] for i in htempkey],[phirRPtemp[i] for i in htempkey]))
np.savetxt("herendeen1979_tem_phi2_phir.csv", herendeen1979_phi2_phir, delimiter=",", header="tem,phi2,phir")



data0 = {'mu':mu, 'phir':{i: allrfrac[i] for i in mu},'phi2':{i:chapefrac[i] for i in mu}}
import pickle
with open('schmit2016_mu_phir_phi2.pkl','wb') as f:
    pickle.dump(data0, f)

lbphi2_1 = [tempP[f'LB{i}_1_norm'][ischape3].sum()/1e7 for i in [16,25,30,37,43]]
lbphi2_2 = [tempP[f'LB{i}_2_norm'][ischape3].sum()/1e7 for i in [16,25,30,37,43]]
lbphir_1 = [tempP[f'LB{i}_1_norm'][isallr3].sum()/1e7 for i in [16,25,30,37,43]]
lbphir_2 = [tempP[f'LB{i}_2_norm'][isallr3].sum()/1e7 for i in [16,25,30,37,43]]
# glcphi2_1 = [tempP[f'Glucose{i}_1_norm'][ischape3].sum()/1e7 for i in [25,30,37]]
# glcphi2_2 = [tempP[f'Glucose{i}_2_norm'][ischape3].sum()/1e7 for i in [25,30,37]]
# glcphir_1 = [tempP[f'Glucose{i}_1_norm'][isallr3].sum()/1e7 for i in [25,30,37]]
# glcphir_2 = [tempP[f'Glucose{i}_2_norm'][isallr3].sum()/1e7 for i in [25,30,37]]
# glyphi2_1 = [tempP[f'Glycerol{i}_1_norm'][ischape3].sum()/1e7 for i in [25,30,37]]
# glyphi2_2 = [tempP[f'Glycerol{i}_2_norm'][ischape3].sum()/1e7 for i in [25,30,37]]
# glyphir_1 = [tempP[f'Glycerol{i}_1_norm'][isallr3].sum()/1e7 for i in [25,30,37]]
# glyphir_2 = [tempP[f'Glycerol{i}_2_norm'][isallr3].sum()/1e7 for i in [25,30,37]]
knapp2024_lb = np.column_stack(([16,25,30,37,43], lbphir_1, lbphir_2, lbphi2_1, lbphi2_2))
np.savetxt("knapp2024_lb_tem_phi2_phir.csv", knapp2024_lb, delimiter=",",header="tem,phir_1,phir_2,phi2_1,phi2_2",comments="")
# knapp2024_glc = np.column_stack(([25,30,37], glcphir_1, glcphir_2, glcphi2_1, glcphi2_2))
# np.savetxt("knapp2024_glc_tem_phi2_phir.csv", knapp2024_glc, delimiter=",",header="tem,phir_1,phir_2,phi2_1,phi2_2",comments="")
# knapp2024_gly = np.column_stack(([25,30,37], glyphir_1, glyphir_2, glyphi2_1, glyphi2_2))
# np.savetxt("knapp2024_gly_tem_phi2_phir.csv", knapp2024_gly, delimiter=",",header="tem,phir_1,phir_2,phi2_1,phi2_2",comments="")

# Based on the above content, data on phi2 at different temperatures can be provided. Meanwhile, the growth rate can be obtained from fig 2b or fig 2c in the original text (knapp et.al., 10.1038/s41564-024-01841-4). From this, we can provide
# knapp2024_37_gr_phi2.csv
# knapp2024_30_gr_phi2.csv
# knapp2024_25_gr_phi2.csv
