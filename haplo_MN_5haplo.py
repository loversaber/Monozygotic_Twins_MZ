from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon#paired samples number of data should be same
from scipy.stats import ranksums
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from statsmodels.sandbox.stats.multicomp import multipletests

haplo=pd.read_table("HaploGroup_2.txt",sep="\t",header=0,index_col=0)
#haplo=haplo.sort_index()
info0=pd.read_table("all0503.xls",sep="\t",header=0,index_col=0)#information include copy number
#lchange=(info0['glu2h']-info0['glu0h'])/info0['glu0h']
#info0.insert(loc=6,column='gluChange',value=lchange)
info=info0.sort_index()
#print(haplo.index)
dfall=pd.merge(haplo,info,left_index=True,right_index=True)
#print(dfall)
listglu=['glu0h', 'glu2h', 'g2-g1','gluChange','diff1','cpy','nheter','avgmaf', 'summaf']


def cprUtest(d,group1, group2):
    u, prob = mannwhitneyu(d[group1], d[group2])
    # T,pvalue=wilcoxon(MgroupGlu,NgroupGlu)
    z_statistic, pvalue = ranksums(d[group1], d[group2])

    #print("%s\t%s\tmanUtest:\t%s\t%.4f\t%.4f" % (group1, group2, glu, u, prob))  # one-sided pvalue
    # print("wilcoxn:\t",T,pvalue)#two-sided pvalue
    #print("%s\t%s\tranksums:\t%s\t%.4f\t%.4f" % (group1, group2, glu, z_statistic, pvalue))
    return (u, prob, z_statistic, pvalue)

def HaploGlu(glu):
    d={}
    M7groupGlu=dfall[dfall.haplogroup=='M7'][glu]
    DgroupGlu=dfall[dfall.haplogroup=='D'][glu]
    CZgroupGlu = dfall[dfall.haplogroup == 'CZ'][glu]
    FgroupGlu = dfall[dfall.haplogroup == 'F'][glu]
    BgroupGlu = dfall[dfall.haplogroup == 'B'][glu]

    d['M7']=list(M7groupGlu)
    d['D']=list(DgroupGlu)
    d['CZ'] = list(CZgroupGlu)
    d['F'] = list(FgroupGlu)
    d['B'] = list(BgroupGlu)
    return d

def HaploGP_MN(glu):
    d = {}
    M = dfall[dfall.haploGP == 'M'][glu]
    N = dfall[dfall.haploGP == 'N'][glu]

    d['M'] = list(M)
    d['N'] = list(N)

    return d

def func_MN_5haplo(func,filname):
    dpvalue = {}
    fig = plt.figure()
    for i, glu in enumerate(listglu):
        # draw picture
        d = func(glu)
        plt.style.use("ggplot")
        ax = fig.add_subplot(3, 3, i + 1)
        plt.boxplot(x=d.values(), labels=d.keys(), showmeans=True)
        ax.set_title(glu)
        # pvalue adjusted
        lgroup = list(combinations(d.keys(), 2))

        for groups in lgroup:
            g1, g2 = groups[0], groups[1]
            U_statistic, pu, z_statistic, pr = cprUtest(d, g1, g2)
            name = g1 + "_" + g2 + "_" + glu
            dpvalue[name] = {}
            # dpvalue[name]['Phenotype'] = glu
            dpvalue[name]['U_statistic'] = U_statistic
            dpvalue[name]['U_pvalue'] = pu
            dpvalue[name]['Z_statistic'] = z_statistic
            dpvalue[name]['Z_pvalue'] = pr

    plt.show()

    df_pvalue = pd.DataFrame(dpvalue)
    df_pvalue = df_pvalue.T
    df_pvalue['U_BH'] = multipletests(df_pvalue['U_pvalue'], method='fdr_bh')[1]
    df_pvalue['U_BY'] = multipletests(df_pvalue['U_pvalue'], method='fdr_by')[1]
    df_pvalue['U_bonferroni'] = multipletests(df_pvalue['U_pvalue'], method='bonferroni')[1]

    df_pvalue['Z_BH'] = multipletests(df_pvalue['Z_pvalue'], method='fdr_bh')[1]
    df_pvalue['Z_BY'] = multipletests(df_pvalue['Z_pvalue'], method='fdr_by')[1]
    df_pvalue['Z_bonferroni'] = multipletests(df_pvalue['Z_pvalue'], method='bonferroni')[1]

    df_pvalue = df_pvalue.sort_values(by='U_BH', ascending=True)
    df_pvalue = df_pvalue.sort_index(axis=1)
    df_pvalue.to_csv(filname+'_haplogroup_0606_adjusted.xls', sep="\t", index=True, float_format='%g')

func_MN_5haplo(HaploGlu,"5haplo")
func_MN_5haplo(HaploGP_MN,"MN")
