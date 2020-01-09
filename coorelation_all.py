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

info=pd.read_table("all0503.xls",sep="\t",header=0,index_col=0 )#InfoMaf.xls
mafmatrix=pd.read_csv("maf0.csv",sep=",",header=0,index_col=0)
mafmatrix.columns=mafmatrix.columns.astype(int)######very inportant
print(len(mafmatrix.columns))
subrefandalt=pd.read_csv("sub_altandfeq_MZ.xls",sep="\t",header=0,index_col=0)
print(len(subrefandalt.columns))
subrefandalt=subrefandalt.loc[:,((subrefandalt!='0').sum()>=5) & ((subrefandalt!='0').sum()<=45)]
subrefandalt.columns=subrefandalt.columns.astype(int)
print(len(subrefandalt.columns))

dfsubinfo=pd.merge(subrefandalt,info,left_index=True,right_index=True)
#print(dfsubinfo)
listglu=['glu0h', 'glu2h', 'g2-g1','gluChange','diff1','cpy','nheter','avgmaf', 'summaf']
listglu_combination=list(combinations(listglu,2))
listglu2=['glu0h', 'glu2h', 'g2-g1','gluChange','diff1','avgmaf', 'summaf']
listglu_combination2=list(combinations(listglu2,2))
#listglu3=['glu0h']

def adjustTestPvalue(d):
    DF = pd.DataFrame(d)
    DF = DF.T
    # print(DF)
    DF['t_BH'] = multipletests(DF['t_pvalue'], method='fdr_bh')[1]
    DF['t_BY'] = multipletests(DF['t_pvalue'], method='fdr_by')[1]
    DF['t_bonferroni'] = multipletests(DF['t_pvalue'], method='bonferroni')[1]

    DF['u_BH'] = multipletests(DF['u_pvalue'], method='fdr_bh')[1]
    DF['u_BY'] = multipletests(DF['u_pvalue'], method='fdr_by')[1]
    DF['u_bonferroni'] = multipletests(DF['u_pvalue'], method='bonferroni')[1]

    DF = DF.sort_values(by='u_BH', ascending=True)  # ascending=False means down
    DF = DF.sort_index(axis=1)
    return(DF)

def adjustCorrPvalue(d):
    DF = pd.DataFrame(d)
    DF = DF.T
    DF['BH'] = multipletests(DF['pvalue'], method='fdr_bh')[1]
    DF['BY'] = multipletests(DF['pvalue'], method='fdr_by')[1]
    DF['Bonferroni'] = multipletests(DF['pvalue'], method='bonferroni')[1]
    # rho, pvalue = pearsonr(info[i[0]], info[i[1]])

    DF = DF.sort_values(by='BH', ascending=True)  # ascending=False means down
    DF = DF.sort_index(axis=1)
    return(DF)

def UTTest(df,glu,filname,symbol):
    d={}
    for col in df.columns:
        if type(col)==int:
            if symbol=="smarker":
                refandsub=df[col].str.split("/",expand=True)
                #print(col,type(col))
            #ref=list(set(refandsub[0]))[0]
            #dalt=dict(refandsub[1].value_counts())
                alt = [s for s in dict(refandsub[1].value_counts()[[0]]).keys()][0]
            #print(col,alt)
            #alt=list(set(refandsub[1]))
                unmut=df[glu][refandsub[1].isna()]
                mut=df[glu][refandsub[1]==alt]
            #print(refglu.count(),altglu.count())
            elif symbol=="hmarker":
                unmut = df[glu][df[col] == 0]
                mut = df[glu][df[col] != 0]

            if unmut.count()>=5 and mut.count()>=5:
                #print(col,unmut.count(),mut.count())
                d[col] = {}
                u, prob = mannwhitneyu(unmut,mut)
                statistic,pvalue=stats.ttest_ind(unmut,mut)
                d[col]['t_statistic']=statistic
                d[col]['t_pvalue']=pvalue
                d[col]['u_test']=u
                d[col]['u_pvalue']=prob
    dfend=adjustTestPvalue(d)
    dfend.to_csv(glu + "_" + filname + ".xls", sep="\t", index=True)

dfall50=pd.merge(mafmatrix,info,left_index=True,right_index=True)
mafmatrix_5=mafmatrix.loc[:,((mafmatrix!=0).sum()>=5)]#subrefandalt is !='0' but there is 0 not a str
print(len(mafmatrix_5.columns))
mafmatrix_5.to_csv("maf50_big5.xls",sep="\t",index=True)


def mainTest():
    for glu in listglu:
        UTTest(dfsubinfo,glu,"MZ50_sub_U_T_test_0606","smarker")
        UTTest(dfall50,glu,"MZ50_maf_U_T_test_0606","hmarker")
#mainTest()

def coorInfo(info,otfilname):#info file correlation each other
    d={}
    for i in listglu_combination:
        name=i[0]+"_"+i[1]
        rho, pvalue = spearmanr(info[i[0]],info[i[1]])
        print("%s\t%s\trho:%.4f\tpvalue:%.4f" % (i[0], i[1], rho, pvalue))
        d[name]={}
        d[name]['rho']=rho
        d[name]['pvalue']=pvalue
    dfend=adjustCorrPvalue(d)
    dfend.to_csv(otfilname+"_corr.xls", sep="\t", index=True)
#coorInfo(info,"info50")

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
    spearmanDF=adjustCorrPvalue(d)
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

mafmatrix_new=MatchGroup25Delete50(mafmatrix,"mito")
mafmatrix_new_3=mafmatrix_new.loc[:,((mafmatrix_new!=0).sum()>=3)]
print(len(mafmatrix_new_3.columns))
mafmatrix_new_3.to_csv("maf25_big3.xls",sep="\t",index=True)

#mafmatrix_new.to_csv("smarker25.xls",sep="\t",index=True)
#mafmatrix_new.to_csv("hmarker25.xls",sep="\t",index=True)
info_new=MatchGroup25Delete50(info,"info")
#coorInfo(info_new,"info25")
dfall25=pd.merge(mafmatrix_new,info_new,left_index=True,right_index=True)
#info_new.to_csv("info25_0507.xls",sep="\t",index=True)
#dfall25.to_csv("maf25_info.xls",sep="\t",index=True)#not run to get file yet
def mainCorrMaf():
    for glu in listglu2:
        corrMitoGlu(dfall50, glu, 5, "MZ50_maf_0606")
        #corrMitoGlu(dfall50, glu, 5, "MZ50_sub")
        corrMitoGlu(dfall25,glu,3,"GP25_maf_0606")#smarker subtract all zero
#mainCorrMaf()

def corr3Pos(dfall,otfilname):#three pos in family and glu
    dpvalue={}
    for pos in [248,8279,14872]:
        for glu in listglu:
        #rho, pvalue = spearmanr(dfall[pos], dfall[glu])
            rho, pvalue = pearsonr(dfall[pos], dfall[glu])
            nfamily=dfall[pos].astype(bool).sum(axis=0)
            print("%s\t%s\t%d\trho:%.4f\tpvalue:%.4f" % (pos, glu,nfamily, rho, pvalue))
            name=str(pos)+"_"+glu
            dpvalue[name] = {}
            #dpvalue[name]['Phenotype'] = glu
            dpvalue[name]['nfamily']=nfamily
            dpvalue[name]['rho'] = rho
            dpvalue[name]['pvalue'] = pvalue
    df_pvalue = adjustCorrPvalue(dpvalue)
    df_pvalue.to_csv(otfilname+'_3pos_pvalue_adjusted.xls', sep="\t", index=True, float_format='%g')
def coor3():
    corr3Pos(dfall25,"df25")
    corr3Pos(dfall50,"df50")
#coor3()
