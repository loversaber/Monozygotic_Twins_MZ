from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon#paired samples number of data should be same
from scipy.stats import ranksums
import pandas as pd
import matplotlib.pyplot as plt

dfall=pd.read_table("sub16362.txt",sep="\t",header=0,index_col=0)

def cprHaploGlu(glu):
    d={}
    group1=dfall[dfall['16362']!='0'][glu]
    group2=dfall[dfall['16362']=='0'][glu]
    d['C']=list(group1)
    d['T']=list(group2)
    def cprUtest(group1,group2):
        u,prob=mannwhitneyu(d[group1],d[group2])
        #T,pvalue=wilcoxon(MgroupGlu,NgroupGlu)
        z_statistic,pvalue=ranksums(d[group1],d[group2] )
        print("%s\t%s\tmanUtest:\t%s\t%.4f\t%.4f"%(group1,group2,glu,u,prob))#one-sided pvalue
        #print("wilcoxn:\t",T,pvalue)#two-sided pvalue
        print("%s\t%s\tranksums:\t%s\t%.4f\t%.4f"%(group1,group2,glu, z_statistic, pvalue))

    cprUtest('C','T')
    return d

fig=plt.figure()
d = cprHaploGlu('glu0h')
plt.style.use("ggplot")
plt.boxplot(x=d.values(),labels=d.keys(),showmeans=True)
#fig.set_title('glu0h')
plt.ylabel('glu0h',fontsize=25)
plt.title('loc_16362_glu0h',fontsize=25)
plt.tick_params(axis="both",labelsize=25)
plt.show()
