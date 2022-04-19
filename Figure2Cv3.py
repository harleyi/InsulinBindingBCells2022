#FIgure2Cv3.py

#Analysis for Figure2C and supporting figures
#"B cell receptor repertoire and affinities of high-affinity insulin-binding B cells in transgenic mouse models of autoimmune diabetes"
#Maureen Banach1, Isaac T. W. Harley1,2,3, Andrew Getahun1, John C. Cambier1*

'''
22-01-24 clusters gene UP-DOWN.xlsx
6 tabs as follows: 5 columns X 51 rows – up or down in clusters 1-3

FeatureID	FeatureName	cluster 1 Average	cluster 1 Log2 Fold Change	cluster 1 P-Value
ENSMUSG00000029603	Dtx1	1.237531684	2.4952666	4.35E-26
ENSMUSG00000026616	Cr2	2.473895885	1.614150194	3.65E-12
ENSMUSG00000029322	Plac8	3.732734128	1.556736983	1.80E-11
ENSMUSG00000045092	S1pr1	1.506928321	1.357118117	1.02E-08
'''

'''
GSE109125_Normalized_Gene_count_table.csv
206 columns X 55012 rows (including header) – Normalized gene expression from
GSE109125 – Immgen RNA-Seq

gene_symbol	B.Fem.Sp#1	B.Fo.Sp#1	B.Fo.Sp#2	B.Fo.Sp#3	B.Fo.Sp#4	B.FrE.BM#1	B.FrE.BM#2	B.GC.CB.Sp#1	B.GC.CB.Sp#2	B.GC.CB.Sp#3
0610005C13Rik	1.551979483	2.130245484	7.371546308	7.766246314	1	5.024681072	2.323716869	8.381132403	1	1
0610006L08Rik	1	1	1	1	1	1	1	1	1	1
0610009B22Rik	50.12617397	38.29810096	35.58839424	58.02979036	34.49665979	74.01921373	75.7900031	70.38264459	85.05688861	83.02371675
...
NKT.19-8-TCRb+CD1daGalCerTet+.Th#2	NKT.19-8-TCRb+CD1daGalCerTet+.Th#3	NKT.19-8-TCRb+CD1daGalCerTet+.Lv#1	NKT.19-8-TCRb+CD1daGalCerTet+.Lv#2	NKT.19-8-TCRb+CD1daGalCerTet+.Lv#3	T.4.19-8-TCRb+CD4+.Sp#1	T.4.19-8-TCRb+CD4+.Sp#2
2.046885851	2.997083777	4.083510379	4.179456157	4.315315856	2.802346958	5.879549841
2.046885851	1	1	1	1	7.308214353	1
33.45346139	28.95917288	19.50106227	21.34851941	21.55495831	28.93637785	22.4700193

'''

'''
Plan for main figure: heatmap
left 3 columns:
– gene rank in clusters
middle 3 columns:
– average expression (tpm?)
Right 3 columns:
– normalized average (median) expression across each of 3 top splenic B cell sample types identified as most likely based on W-plot
'''


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


#read in files
path = "~/Desktop/MB InsulinBindingBcells/2022_04_13/"

Immgen = pd.read_csv(path+"GSE109125_Normalized_Gene_count_table.csv")
ImmgenColumns = Immgen.columns.to_numpy()
'''
array(['gene_symbol', 'B.Fem.Sp#1', 'B.Fo.Sp#1', 'B.Fo.Sp#2', 'B.Fo.Sp#3',
       'B.Fo.Sp#4', 'B.FrE.BM#1', 'B.FrE.BM#2', 'B.GC.CB.Sp#1',
       'B.GC.CB.Sp#2', 'B.GC.CB.Sp#3', 'B.GC.CC.Sp#1', 'B.GC.CC.Sp#2',
       'B.GC.CC.Sp#3', 'B.mem.Sp#1', 'B.MZ.Sp#1', 'B.MZ.Sp#2',
       'B.PB.Sp#1', 'B.PB.Sp#2', 'B.PC.BM#1', 'B.PC.Sp#2', 'B.Sp#4',
       'B.T1.Sp#2', 'B.T2.Sp#1', 'B.T2.Sp#2', 'B.T3.Sp#2', 'B1b.PC#1',
       'B1b.PC#2', 'BEC.SLN#1', 'BEC.SLN#2', 'BEC.SLN#3', 'DC.4+.Sp#1',
       'DC.4+.Sp#2', 'DC.4+.Sp#3', 'DC.8+.Sp#1', 'DC.8+.Sp#2',
       'DC.8+.Sp#3', 'DC.pDC.Sp#1', 'DC.pDC.Sp#2', 'Ep.MECHi.Th#1',
       'Ep.MECHi.Th#2', 'FRC.CD140a+.Madcam-.CD35-.SLN#1',
       'FRC.CD140a+.Madcam-.CD35-.SLN#2',
       'FRC.CD140a+.Madcam-.CD35-.SLN#3', 'GN.BM#1', 'GN.BM#2', 'GN.Sp#3',
       'GN.Sp#4', 'GN.Thio.PC#1', 'GN.Thio.PC#2', 'IAP.SLN#1',
       'IAP.SLN#2', 'ILC2.SI#1', 'ILC2.SI#2', 'ILC2.ST2-.SI#2',
       'ILC2.ST2-.SI#1', 'ILC3.CCR6+.SI#1', 'ILC3.CCR6+.SI#2',
       'ILC3.NKp46+.SI#1', 'ILC3.NKp46+.SI#2', 'ILC3.NKp46-CCR6-.SI#2',
       'LEC.SLN#2', 'LEC.SLN#3', 'LTHSC.34-.BM#1', 'LTHSC.34-.BM#2',
       'LTHSC.34+.BM#1', 'LTHSC.34+.BM#2', 'MC.heparinase.PC#3',
       'MF.Alv.Lu#1', 'MF.Alv.Lu#2', 'MF.Fem.PC#1', 'MF.Fem.PC#2',
       'MF.PC#3', 'MF.PC#4', 'MF.pIC.Alv.Lu#2', 'MMP2.150+48+.BM#1',
       'MMP3.48+.BM#1', 'MMP3.48+.BM#2', 'MMP4.135+.BM#1',
       'MMP4.135+.BM#2', 'NK.27+11b-.BM#1', 'NK.27+11b-.BM#2',
       'NK.27+11b-.Sp#1', 'NK.27+11b-.Sp#2', 'NK.27+11b+.BM#1',
       'NK.27+11b+.BM#2', 'NK.27+11b+.Sp#1', 'NK.27+11b+.Sp#2',
       'NK.27-11b+.BM#1', 'NK.27-11b+.BM#2', 'NK.27-11b+.Sp#2',
       'NKT.Sp#3', 'NKT.Sp.LPS.18hr#1', 'NKT.Sp.LPS.18hr#2',
       'NKT.Sp.LPS.3d#2', 'NKT.Sp.LPS.3hr#1', 'NKT.Sp.LPS.3hr#2',
       'preT.DN1.Th#1', 'preT.DN1.Th#2', 'preT.DN2a.Th#2',
       'preT.DN2b.Th#1', 'preT.DN2b.Th#2', 'preT.DN3.Th#1',
       'preT.DN3.Th#2', 'proB.CLP.BM#1', 'proB.CLP.BM#2', 'proB.FrA.BM#1',
       'proB.FrA.BM#2', 'proB.FrBC.BM#1', 'proB.FrBC.BM#2',
       'STHSC.150-.BM#1', 'STHSC.150-.BM#2', 'T.4.Nve.Fem.Sp#1',
       'T.4.Nve.Fem.Sp#2', 'T.4.Nve.Sp#1', 'T.4.Nve.Sp#2',
       'T.4.Sp.aCD3+CD40.18hr#1', 'T.4.Sp.aCD3+CD40.18hr#2', 'T.4.Th#1',
       'T.4.Th#2', 'T.8.Nve.Sp#1', 'T.8.Nve.Sp#2', 'T.8.Th#1', 'T.8.Th#2',
       'T.DN4.Th#1', 'T.DN4.Th#2', 'T.DP.Th#1', 'T.DP.Th#2', 'T.ISP.Th#1',
       'T.ISP.Th#2', 'T8.IEL.LCMV.d7.Gut#1', 'T8.IEL.LCMV.d7.Gut#2',
       'T8.MP.LCMV.d7.Sp#1', 'T8.MP.LCMV.d7.Sp#2',
       'T8.Tcm.LCMV.d180.Sp#1', 'T8.Tcm.LCMV.d180.Sp#2',
       'T8.TE.LCMV.d7.Sp#1', 'T8.TE.LCMV.d7.Sp#2',
       'T8.Tem.LCMV.d180.Sp#2', 'T8.TN.P14.Sp#1', 'T8.TN.P14.Sp#2',
       'Tgd.g1.1+d1.24a+.Th#1', 'Tgd.g1.1+d1.24a+.Th#2',
       'Tgd.g1.1+d1.LN#1', 'Tgd.g1.1+d1.LN#2', 'Tgd.g2+d1.24a+.Th#1',
       'Tgd.g2+d1.24a+.Th#2', 'Tgd.g2+d1.LN#1', 'Tgd.g2+d1.LN#2',
       'Tgd.g2+d17.24a+.Th#2', 'Tgd.g2+d17.LN#1', 'Tgd.g2+d17.LN#2',
       'Tgd.Sp#3', 'Tgd.Sp#4', 'Treg.4.25hi.Sp#1', 'Treg.4.25hi.Sp#2',
       'Treg.4.FP3+.Nrplo.Co#1', 'Treg.4.FP3+.Nrplo.Co#2',
       'MF.102+480+.PC#1', 'MF.102+480+.PC#2', 'MF.microglia.CNS#1',
       'MF.microglia.CNS#2', 'MF.RP.Sp#1', 'MF.RP.Sp#2', 'MF.AT#1',
       'MF.AT#2', 'MF.226+II+480lo.PC#1', 'MF.226+II+480lo.PC#2',
       'Mo.6C-II-.Bl#1', 'Mo.6C-II-.Bl#2', 'Mo.6C+II-.Bl#1',
       'Mo.6C+II-.Bl#2', 'Eo.PC#1', 'Eo.PC#2', 'Eo.Sp#2', 'Eo.Sp#3',
       'Ba.Sp#1', 'Ba.Sp#2', 'Ba.Sp#3', 'MC.Ht#1', 'MC.Ht#2',
       'MC.SCA1hi.CX3CR1hi.PC#1', 'MC.SCA1hi.CX3CR1hi.PC#2',
       'MC.SCA1hi.CX3CR1hi.PC#3', 'MC.SCA1hi.CX3CR1lo.PC#1',
       'MC.SCA1hi.CX3CR1lo.PC#2', 'MC.SCA1hi.CX3CR1lo.PC#3',
       'MC.SCA1lo.CX3CR1hi.PC#1', 'MC.SCA1lo.CX3CR1hi.PC#2',
       'MC.SCA1lo.CX3CR1hi.PC#3', 'MC.SCA1lo.CX3CR1lo.PC#1',
       'MC.SCA1lo.CX3CR1lo.PC#2', 'MC.SCA1lo.CX3CR1lo.PC#3',
       'NKT.19-8-TCRb+CD1daGalCerTet+.Lu#1',
       'NKT.19-8-TCRb+CD1daGalCerTet+.Lu#2',
       'NKT.19-8-TCRb+CD1daGalCerTet+.Sp#1',
       'NKT.19-8-TCRb+CD1daGalCerTet+.Sp#2',
       'NKT.19-8-TCRb+CD1daGalCerTet+.Sp#3',
       'NKT.19-8-TCRb+CD1daGalCerTet+.Th#1',
       'NKT.19-8-TCRb+CD1daGalCerTet+.Th#2',
       'NKT.19-8-TCRb+CD1daGalCerTet+.Th#3',
       'NKT.19-8-TCRb+CD1daGalCerTet+.Lv#1',
       'NKT.19-8-TCRb+CD1daGalCerTet+.Lv#2',
       'NKT.19-8-TCRb+CD1daGalCerTet+.Lv#3', 'T.4.19-8-TCRb+CD4+.Sp#1',
       'T.4.19-8-TCRb+CD4+.Sp#2'], dtype=object)
'''
ImmgenBCells = Immgen[['gene_symbol', 'B.Fem.Sp#1', 'B.Fo.Sp#1', 'B.Fo.Sp#2', 'B.Fo.Sp#3',
       'B.Fo.Sp#4', 'B.FrE.BM#1', 'B.FrE.BM#2', 'B.GC.CB.Sp#1',
       'B.GC.CB.Sp#2', 'B.GC.CB.Sp#3', 'B.GC.CC.Sp#1', 'B.GC.CC.Sp#2',
       'B.GC.CC.Sp#3', 'B.mem.Sp#1', 'B.MZ.Sp#1', 'B.MZ.Sp#2',
       'B.PB.Sp#1', 'B.PB.Sp#2', 'B.PC.BM#1', 'B.PC.Sp#2', 'B.Sp#4',
       'B.T1.Sp#2', 'B.T2.Sp#1', 'B.T2.Sp#2', 'B.T3.Sp#2', 'B1b.PC#1',
       'B1b.PC#2', 'proB.CLP.BM#1', 'proB.CLP.BM#2', 'proB.FrA.BM#1',
       'proB.FrA.BM#2', 'proB.FrBC.BM#1', 'proB.FrBC.BM#2']]

#This is the one that we went with – 
CombinedClusters=pd.read_csv(path+"CombinedClusters123.csv")

SigGenesMask = ((CombinedClusters['cluster 1 P-Value'] <= 0.05) | (CombinedClusters['cluster 2 P-Value'] <= 0.05) | (CombinedClusters['cluster 3 P-Value'] <= 0.05))
SignificantGenes = CombinedClusters[SigGenesMask]
'''>>> SignificantGenes.columns
Index(['FeatureID', 'FeatureName', 'cluster 1 Average',
       'cluster 1 Log2 Fold Change', 'cluster 1 P-Value', 'cluster 2 Average',
       'cluster 2 Log2 Fold Change', 'cluster 2 P-Value', 'cluster 3 Average',
       'cluster 3 Log2 Fold Change', 'cluster 3 P-Value'],
      dtype='object')
'''

SignificantAverages = SignificantGenes[['FeatureID', 'FeatureName', 'cluster 1 Average','cluster 2 Average','cluster 3 Average']].copy(deep=True)
SignificantAverages.fillna(0, inplace=True)

Genes = SignificantAverages.FeatureName
BCellGeneMask = ImmgenBCells.gene_symbol.isin(Genes)
#>>>BCellGeneMask.sum()
#50
#len(Genes)
#>>>52
# 2 gene names not found...
BcellGenes = ImmgenBCells[BCellGeneMask]
BCG = BcellGenes.gene_symbol.tolist()
#>>>set(BCG)^set(Genes)
#{'Hmha1', 'AY036118'}

'''
Now take median across the B cell subsets
'''

#BcellGenes.assign(B_FO_Sp=lambda x: np.median(x[['B.Fo.Sp#1', 'B.Fo.Sp#2', 'B.Fo.Sp#3', 'B.Fo.Sp#4']]))
#median over the whole slice ...
'''
only assign
ImmgenBCells = Immgen[['gene_symbol',
B_Fo_Sp     'B.Fo.Sp#1', 'B.Fo.Sp#2', 'B.Fo.Sp#3', 'B.Fo.Sp#4',
B_mem_Sp    'B.mem.Sp#1',
B_MZ_Sp     'B.MZ.Sp#1', 'B.MZ.Sp#2',
'''
BcellGenes = BcellGenes.assign(B_Fo_Sp=lambda x: np.median(x[['B.Fo.Sp#1', 'B.Fo.Sp#2', 'B.Fo.Sp#3', 'B.Fo.Sp#4']], axis=1))
BcellGenes = BcellGenes.assign(B_mem_Sp=lambda x: np.median(x[['B.mem.Sp#1']], axis=1))
BcellGenes = BcellGenes.assign(B_MZ_Sp=lambda x: np.median(x[['B.MZ.Sp#1', 'B.MZ.Sp#2']], axis=1))

BcellGenes = BcellGenes.assign(B_FrE_BM=lambda x: np.median(x[['B.FrE.BM#1', 'B.FrE.BM#2']], axis=1))
BcellGenes = BcellGenes.assign(B_GC_CB=lambda x: np.median(x[['B.GC.CB.Sp#1','B.GC.CB.Sp#2', 'B.GC.CB.Sp#3']], axis=1))
BcellGenes = BcellGenes.assign(B_GC_CC=lambda x: np.median(x[['B.GC.CC.Sp#1', 'B.GC.CC.Sp#2', 'B.GC.CC.Sp#3']], axis=1))
BcellGenes = BcellGenes.assign(B_PB_Sp=lambda x: np.median(x[['B.PB.Sp#1', 'B.PB.Sp#2']], axis=1))
BcellGenes = BcellGenes.assign(B_PC_Sp=lambda x: np.median(x[['B.PC.Sp#2']], axis=1))
BcellGenes = BcellGenes.assign(B_PC_BM=lambda x: np.median(x[['B.PC.BM#1']], axis=1))
BcellGenes = BcellGenes.assign(B_T1_Sp=lambda x: np.median(x[['B.T1.Sp#2']], axis=1))
BcellGenes = BcellGenes.assign(B_T2_Sp=lambda x: np.median(x[['B.T2.Sp#1', 'B.T2.Sp#2']], axis=1))
BcellGenes = BcellGenes.assign(B_T3_Sp=lambda x: np.median(x[['B.T3.Sp#2']], axis=1))
BcellGenes = BcellGenes.assign(B_B1b_PC=lambda x: np.median(x[['B1b.PC#1','B1b.PC#2']], axis=1))
BcellGenes = BcellGenes.assign(proB_CLP_BM=lambda x: np.median(x[['proB.CLP.BM#1', 'proB.CLP.BM#2']], axis=1))
BcellGenes = BcellGenes.assign(proB_FrA_BM=lambda x: np.median(x[['proB.FrA.BM#1','proB.FrA.BM#2']], axis=1))
BcellGenes = BcellGenes.assign(proB_FrBC_BM=lambda x: np.median(x[['proB.FrBC.BM#1', 'proB.FrBC.BM#2']], axis=1))


#Subset only columns with measures of central tendency.
#Order according to Gene skyline/developmental progression
#http://rstats.immgen.org/Skyline/skyline.html
#BCG_Medians = BcellGenes[['gene_symbol', 'B_Fo_Sp', 'B_MZ_Sp', 'B_mem_Sp']]
BCG_Medians = BcellGenes[['gene_symbol', 'proB_CLP_BM', 'proB_FrA_BM', 'proB_FrBC_BM', 'B_FrE_BM', 'B_T1_Sp', 'B_T2_Sp', 'B_T3_Sp', 'B_Fo_Sp', 'B_MZ_Sp', 'B_mem_Sp', 'B_GC_CB', 'B_GC_CC', 'B_PB_Sp', 'B_PC_Sp', 'B_PC_BM', 'B_B1b_PC']]

#Log normalize Maureen Data
L10NSA = np.log10(SignificantAverages[['cluster 1 Average','cluster 2 Average','cluster 3 Average']]+1)

#Rank the averages across cluster
Ranks = SignificantAverages[['cluster 1 Average','cluster 2 Average','cluster 3 Average']].T.rank(method="first").T
Ranks.columns = ['cluster 1 Rank', 'cluster 2 Rank', 'cluster 3 Rank']

GeneIDs = SignificantAverages[['FeatureID', 'FeatureName']]
RankSorted =pd.DataFrame([])

Up1 = GeneIDs.join(Ranks[Ranks['cluster 1 Rank']==3].join(L10NSA), how='inner').sort_values(by=['cluster 1 Average'])
Up2 = GeneIDs.join(Ranks[Ranks['cluster 2 Rank']==3].join(L10NSA), how='inner').sort_values(by=['cluster 2 Average'])
Up3 = GeneIDs.join(Ranks[Ranks['cluster 3 Rank']==3].join(L10NSA), how='inner').sort_values(by=['cluster 3 Average'])
RankSorted = pd.concat([Up1,Up2,Up3], axis=0)
#RankSorted.concat(GeneIDs.join(Ranks[Ranks['cluster 1 Rank']==3]].join(L10NSA.sort_values(by=['cluster 1 Average']))) ,axis=0)
RankSorted.columns
'''
>>> RankSorted.columns
Index(['FeatureID', 'FeatureName', 'cluster 1 Rank', 'cluster 2 Rank',
       'cluster 3 Rank', 'cluster 1 Average', 'cluster 2 Average',
       'cluster 3 Average'],
      dtype='object')
'''

'''
Now merge the BCG_Medians from before with the RankSorted dataframe on gene_symbol & FeatureName
Index(['FeatureID', 'FeatureName', 'cluster 1 Rank', 'cluster 2 Rank',
       'cluster 3 Rank', 'cluster 1 Average', 'cluster 2 Average',
       'cluster 3 Average', 'gene_symbol', 'B_Fo_Sp',
       'B_MZ_Sp', 'B_mem_Sp'],
      dtype='object')
'''

MaureenImmgen = pd.merge(RankSorted, BCG_Medians, how='right', left_on='FeatureName', right_on='gene_symbol')

MaureenImmgen1 = MaureenImmgen[MaureenImmgen['cluster 1 Rank'] ==3].sort_values(by=['cluster 1 Average'])
MaureenImmgen2 = MaureenImmgen[MaureenImmgen['cluster 2 Rank'] ==3].sort_values(by=['cluster 2 Average'])
MaureenImmgen3 = MaureenImmgen[MaureenImmgen['cluster 3 Rank'] ==3].sort_values(by=['cluster 3 Average'])

MaureenImmgen = pd.concat([MaureenImmgen3,MaureenImmgen2,MaureenImmgen1], axis=0)

SpleenPoplist = ['B_T1_Sp', 'B_T2_Sp', 'B_T3_Sp', 'B_Fo_Sp','B_MZ_Sp', 'B_mem_Sp', 'B_GC_CB', 'B_GC_CC', 'B_PB_Sp', 'B_PC_Sp']
BcellSlice = MaureenImmgen[SpleenPoplist]
rownorm = BcellSlice.apply(lambda row: (row - np.average(row)) / np.std(row), axis=1)

SelectPoplist = ['B_MZ_Sp', 'B_Fo_Sp','B_mem_Sp']
NormBcellMaureenImmgen = MaureenImmgen[['FeatureName', 'cluster 1 Rank', 'cluster 2 Rank','cluster 3 Rank','cluster 1 Average', 'cluster 2 Average','cluster 3 Average']].join(rownorm[SelectPoplist])

'''
SelectPoplist = ['B_Fo_Sp','B_MZ_Sp', 'B_mem_Sp']
BcellSlice = MaureenImmgen[SelectPoplist]
rownorm = BcellSlice.apply(lambda row: (row - np.average(row)) / np.std(row), axis=1)
'''


plt.clf()

fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(25, 15))

d = ax0.matshow(NormBcellMaureenImmgen[['cluster 1 Rank', 'cluster 2 Rank','cluster 3 Rank']], cmap='binary')
e = ax1.matshow(NormBcellMaureenImmgen[['cluster 1 Average', 'cluster 2 Average','cluster 3 Average']], cmap='YlOrRd', alpha=0.75)
#BcellPoplist = ['B_T1_Sp', 'B_T2_Sp', 'B_T3_Sp', 'B_Fo_Sp','B_MZ_Sp', 'B_mem_Sp', 'B_GC_CB', 'B_GC_CC', 'B_PB_Sp', 'B_PC_Sp']
f =  ax2.matshow(NormBcellMaureenImmgen[SelectPoplist], cmap='YlOrRd', alpha=0.75)

cbar = fig.colorbar(f, ax=fig.gca())
cbar.set_ticks([])

y_axis = np.arange(len(MaureenImmgen.FeatureName))
ax0.set_yticks(y_axis)
ax0.tick_params(axis='y', labelsize='large')
ax0.set_yticklabels(MaureenImmgen.FeatureName, fontweight='bold', rotation=45, ha='right')

clusters = ["1","2","3"]
ax0.set_xticks([0,1,2])
ax0.tick_params(axis='x', labelsize='large', rotation=90)
ax0.set_xticklabels(clusters, fontweight='bold')

ax1.set_xticks([0,1,2])
ax1.tick_params(axis='x', labelsize='large', rotation=90)
ax1.set_xticklabels(clusters, fontweight='bold')
ax1.axes.yaxis.set_visible(False)

SelectPopNames = ['Marginal Zone B', 'Follicular B', 'Memory B']
ax2_xaxis = np.arange(len(SelectPopNames))
ax2.set_xticks(ax2_xaxis)
ax2.tick_params(axis='x', labelsize='large')
ax2.set_xticklabels(SelectPopNames, fontweight='bold', rotation=90, ha='left')
ax2.axes.yaxis.set_visible(False)

plt.tight_layout()
plt.savefig('MaureenClustersHeatmap3subsetv5.png')
plt.savefig('MaureenClustersHeatmap3subsetv5.svg', dpi=600)
