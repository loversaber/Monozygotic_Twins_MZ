#usage:python3 get_samrker_for_cadd.py ids_file1_file2
from sys import argv
from itertools import islice
import os

l1=list(range(302,317))
l2=list(range(513,527))
l3=list(range(566,574))
l4=list(range(16181,16195))
l5=l1+l2+l3+l4
l5.extend([3106,3107,16519])
l6=[str(x) for x in l5]

#sfile=open(argv[2],"w")
#hfile=open(argv[3],"w")

def op(f):
 name=f.split(".")[0]
 smarker,heter=name+".s",name+".h"
 sfile=open(smarker,"w")
 hfile=open(heter,"w")
 with open(f,"r") as f:
  for i in islice(f,1,None):
   i=i.strip().split("\t")
   if i[1] not in l6 and i[26]=='1':
    sfile.write("%s\t%s\t%s\t%s\t%s\n"%("MT",i[1],".",i[0],i[11]))#MT 1747 . G A
  
   if len(i)>27 and i[1] not in l6:
    if i[0]!=i[20]:
     hfile.write("%s\t%s\t%s\t%s\t%s\n"%("MT",i[1],".",i[0],i[20]))
    elif i[0]==i[20] and i[0]!=i[11]:
     hfile.write("%s\t%s\t%s\t%s\t%s\n"%("MT",i[1],".",i[0],i[11]))
    elif i[0]==i[20] and i[0]==i[11]:
     hfile.write("%s\t%s\t%s\t%s\t%s\n"%("MT",i[1],".",i[0],i[16]))
    
 sfile.close()
 hfile.close()
 return(smarker,heter)

def cadd_pro(f):#use cadd to get the tsv file
 cmd1="gzip %s"%(f)
 os.system(cmd1)
 cmd2="/home/bioinfo/liuqi/tool/CADD_v1.3/bin/score_anno.sh %s.gz %s.tsv.gz"%(f,f)
 os.system(cmd2)
 cmd3="gunzip %s.tsv.gz"%(f)
 os.system(cmd3)
 tsvfile=f+".tsv"
 return(tsvfile) 

def op_tsv(f):#change tsv file into dict
 d={}
 with open(f,"r")as fil:
  for line in islice(fil,2,None):
   i=line.strip().split("\t")
   d[i[1]]=list(i[9],i[10],i[95],i[114],i[115])
 return(d)  

def cpr_d(d1,d2):
 d1u,d2u,dsame={},{},{}
 for k in d1:
  if k not in d2:
   d1u[k]=d1[k]
  if k in d2:
   #d1[k].extend(d2[k])
   dsame[k]=list(d1[k],d2[k])
 for j in d2:
  if j not in d1:
   d2u[j]=d2[j]
 return(d1u,d2u,dsame) 

su1file=open("S#uniq_1.txt","w")
su2file=open("S#uniq_2.txt","w")
ssamef= open("s#same.txt","w")

hu1file=open("H#uniq_1.txt","w")
hu2file=open("H#uniq_2.txt","w")
hsamef= open("H#same.txt","w")

def wrt_u(f,d,name):
 for k in d:
  f.write("%s\t%s\t"%(name,k))
  f.write("\t".join(d[k])+"\n")

def wrt_same(f,d,name):
 for k in d:
  f.write("%s\t%s\t"%(name,k))
  f.write("\t".join(d[k][0])+"\t")
  f.write("\t".join(d[k][1])+"\n")
  
with open(argv[1],"r") as ids:
 for line in ids.readlines():
  l0=line.strip().split()
  f1,f2=l0[0],l0[1]
  u1name=f1+"_"+f2
  u2name=f2+"_"+f1
  f1s,f1h=op(f1)
  f2s,f2h=op(f2)
  tsvf1s,tsvf1h,tsvf2s,tsvf2h=cadd_pro(f1s),cadd_pro(f1h),cadd_pro(f2s),cadd_pro(f2h)
 
  t1s,t2s=op_tsv(tsvf1s),op_tsv(tsvf2s)
  d1,d2,d3=cpr_d(t1s,t2s)
  wrt_u(su1file,d1,u1name)
  wrt_u(su2file,d2,u2name)
  wrt_same(ssamef,d3,u1name)
 
  t1h,t2h=op_tsv(tsvf1h),op_tsv(tsvf2h)
  d4,d5,d6=cpr_d(t1h,t2h)
  wrt_u(hu1file,d4,u1name)
  wrt_u(hu2file,d5,u2name)
  wrt_same(hsamef,d6,u1name)

su1file.close()
su2file.close()
ssamef.close()

hu1file.close()
hu2file.close()
hsamef.close()
