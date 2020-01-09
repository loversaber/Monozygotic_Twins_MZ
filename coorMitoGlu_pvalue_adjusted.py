from itertools import combinations
import pandas as pd
import numpy as np
from sys import argv
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import mannwhitneyu
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests

ids=pd.read_table("ids.txt",sep=" ",names=["sample1","sample2"])
d_ids=dict(zip(ids['sample1'],ids['sample2']))

info0=pd.read_table("phy1.txt", )#InfoMaf.xls
lchange=(info0['glu2h']-info0['glu0h'])/info0['glu0h']
info0.insert(loc=6,column='gluChange',value=lchange)
info=info0.sort_index()

#mafmatrix=pd.read_csv("Smarker_Matrix_MZ.xls",sep="\t",header=0,index_col=0)#samrker sep should be careful
mafmatrix=pd.read_csv("maf0.csv",sep=",",header=0,index_col=0)
mafmatrix.columns=mafmatrix.columns.astype(int)######very inportant

subrefandalt=pd.read_csv("sub_altandfeq_MZ.xls",sep=",",header=0,index_col=0)
subrefandalt
subrefandalt=subrefandalt.loc[:,((subrefandalt!='0').sum()>=5) & ((subrefandalt!='0').sum()<=45)]
subrefandalt.columns=subrefandalt.columns.astype(int)

dfsubinfo=pd.merge(subrefandalt,info,left_index=True,right_index=True)

listglu=['glu0h', 'glu2h', 'g2-g1','gluChange','diff1','nheter','avgmaf', 'summaf']
listglu_combination=list(combinations(listglu,2))
listglu2=['glu0h', 'glu2h', 'g2-g1','gluChange','diff1','avgmaf', 'summaf']
listglu_combination2=list(combinations(listglu2,2))


def corrSubGlu(glu,filname):
    d={}
    for col in dfsubinfo.columns:
        if type(col)==int:
            refandsub=dfsubinfo[col].str.split("/",expend=True)

            ref=list(set(refandsub[0]))[0]
            #dalt=dict(refandsub[1].value_counts())
            alt=dict(refandsub[1].value_counts()[[0]]).keys()
            print(col,alt)
            #alt=list(set(refandsub[1]))
            refglu=dfsubinfo[glu][refandsub[1].isna()]
            altglu=dfsubinfo[glu][refandsub[1]==alt]

            if refglu.count()>=5 and altglu.count()>=5:
                d[col] = {}
                u, prob = mannwhitneyu(refglu,altglu)
                statistic,pvalue=stats.ttest_ind(refglu,altglu)
                d[col]['t_statistic']=statistic
                d[col]['t_pvalue']=pvalue
                d[col]['u_test']=u
                d[col]['u_pvalue']=prob
    DF = pd.DataFrame(d)
    DF = DF.T
    DF = DF.sort_values(by='t_pvalue', ascending=True)  # ascending=False means down
    DF = DF.sort_index(axis=1)
    DF.to_csv(glu + "_" + filname + ".xls", sep="\t", index=True)

def mainsub():
    for glu in listglu2:
        corrSubGlu(glu,"MZ50_sub_0427")
        #corrMitoGlu(dfall50, glu, 5, "MZ50_sub")
        ##smarker subtract all zero
mainsub()

def coorInfo50():#info file correlation each other
    for i in listglu_combination:
        rho, pvalue = spearmanr(info[i[0]],info[i[1]])
        #rho, pvalue = pearsonr(info[i[0]], info[i[1]])
        print("%s\t%s\trho:%.4f\tpvalue:%.4f" % (i[0], i[1], rho, pvalue))
#coorInfo50()

def corrMitoGlu(df,glu,nfam,filname):
    d={}
    for col in df.columns:
        #df.loc[df[col]<0.005,col]=0
        nfamily=df[col].astype(bool).sum(axis=0)
        if type(col)==int and nfamily >= nfam:
            #rho,pvalue=pearsonr(df[col],df[glu])
            rho,pvalue=spearmanr(df[col],df[glu])
            if rho != np.nan:
                #print(type(rho), pvalue)
                d[col]={}
                d[col]['rho']=rho
                d[col]['abs_rho']=abs(rho)
                d[col]['pvalue']=pvalue
                d[col]['n_family']=nfamily

    spearmanDF=pd.DataFrame(d)
    spearmanDF=spearmanDF.T

    spearmanDF['BH']=multipletests(spearmanDF['pvalue'],method='fdr_bh')[1]
    spearmanDF['BY']=multipletests(spearmanDF['pvalue'],method='fdr_by')[1]
    spearmanDF['bonferroni']=multipletests(spearmanDF['pvalue'],method='bonferroni')[1]

    spearmanDF=spearmanDF.sort_values(by='pvalue',ascending=True)#ascending=False means down
    spearmanDF=spearmanDF.sort_index(axis=1)
    spearmanDF.to_csv(glu+"_"+filname+".xls",sep="\t",index=True)

def MatchGroup25Delete50(matrixdf,matrixType):
    for k in d_ids.keys():
        group_name=k+"_"+d_ids[k]
        matrixdf.loc[group_name]=matrixdf.loc[k]-matrixdf.loc[d_ids[k]]
        matrixdf.drop(labels=[k,d_ids[k]],axis=0,inplace=True)
    if matrixType=="mito":
        matrixdf=matrixdf.sort_index()
        matrixdf.columns=matrixdf.columns.astype(int)
    else:
        pass
    return matrixdf
dfall50=pd.merge(mafmatrix,info,left_index=True,right_index=True)
mafmatrix_new=MatchGroup25Delete50(mafmatrix,"mito")

#mafmatrix_new.to_csv("smarker25.xls",sep="\t",index=True)
#mafmatrix_new.to_csv("hmarker25.xls",sep="\t",index=True)

info_new=MatchGroup25Delete50(info,"info")

dfall25=pd.merge(mafmatrix_new,info_new,left_index=True,right_index=True)
info_new.to_csv("info25.xls",sep="\t",index=True)
def main():
    for glu in listglu2:
        corrMitoGlu(dfall50, glu, 5, "MZ50_maf_0426")
        #corrMitoGlu(dfall50, glu, 5, "MZ50_sub")
        corrMitoGlu(dfall25,glu,3,"GP25_maf_0426")#smarker subtract all zero
#main()

def InfoPdCorr(dfall):
    def corrInfo(info1,info2):
        rho,pvalue=spearmanr(dfall[info1],dfall[info2])
        #rho, pvalue = pearsonr(dfall[info1], dfall[info2])
        #if pvalue<=0.1:
        print("%s\t%s\trho:%.4f\tpvalue:%.4f"%(info1,info2,rho,pvalue))
    for two_info in listglu_combination2:
        info1,info2=two_info[0],two_info[1]
        corrInfo(info1,info2)
#InfoPdCorr(dfall25)#not well result
def corr3Pos(dfall,glu):#three pos in family and glu
    for pos in [248,8279,14872]:
        #rho, pvalue = spearmanr(dfall[pos], dfall[glu])
        rho, pvalue = pearsonr(dfall[pos], dfall[glu])
        nfamily=dfall[pos].astype(bool).sum(axis=0)
        print("%s\t%s\t%d\trho:%.4f\tpvalue:%.4f" % (pos, glu,nfamily, rho, pvalue))

def coor3():
    for glu in listglu2:
        corr3Pos(dfall25,glu)
#coor3()


