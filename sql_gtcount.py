## Get mutation count from trios and genotype files

import sys
import pandas as pd
import numpy as np
import sqlite3

triofile = sys.argv[1]
db1 = sys.argv[2]
db2 = sys.argv[3]
ibdfile = sys.argv[4]
mapfile = sys.argv[5]
idfile = sys.argv[6]
use.phase = sys.argv[7]

connex1 = sqlite3.connect(db1)
connex2 = sqlite3.connect(db2)

trios = pd.read_csv(triofile, sep=" ", header=None)
trioid = len(trios)//4 
maps = pd.read_csv(mapfile, sep="\t", header=None)
mapmin = maps.iloc[:,3].min()
mapmax = maps.iloc[:,3].max()

ibd = pd.read_csv(ibdfile, sep="\t", header=None)
ibd.columns = ['id1','hap1','id2','hap2','chr','start','end','cM']
ibd = ibd[(ibd['start']>=mapmin) & (ibd['end']<=mapmax)]

mapid = dict()
idindex = 1
for line in open(idfile):
  bits = line.split()
  mapid[bits[0]] = idindex
  idindex = idindex+1

def reid(id):
  newid = mapid[id]
  return (newid)

reid = np.vectorize(reid)

def otherallele(id):
  if id%2==0:
    newid = id-1
  else:
    newid = id+1
  return (newid) 

otherallele = np.vectorize(otherallele)

def reorder(m):
  dip = np.append(m*2-1,m*2)
  return (dip)

ibd.iloc[:,0] = reid(ibd.iloc[:,0])
ibd.iloc[:,2] = reid(ibd.iloc[:,2])
newibd = ibd.iloc[:,0:8]
newibd.columns = ['hap1','hap2','chr','start','end','cM','starttrim','endtrim']
newibd.iloc[:,0] = ibd.iloc[:,0]*2+ibd.iloc[:,1]-2
newibd.iloc[:,1] = ibd.iloc[:,2]*2+ibd.iloc[:,3]-2
newibd.iloc[:,2:6] = ibd.iloc[:,4:8]
## trim ibd 0.5cM
starttrim = np.interp(ibd.iloc[:,5],maps.iloc[:,3],maps.iloc[:,2])+0.5
newibd.iloc[:,6] = np.interp(starttrim,maps.iloc[:,2],maps.iloc[:,3]).astype(int)
endtrim = np.interp(ibd.iloc[:,6],maps.iloc[:,3],maps.iloc[:,2])-0.5
newibd.iloc[:,7] = np.interp(endtrim,maps.iloc[:,2],maps.iloc[:,3]).astype(int)
ibd = newibd


for i in range(1,trioid+1):
  trio = trios.iloc[(i*4-4):(i*4-1),:]
  ids = trio.iloc[:,0].append(trio.iloc[:,1]).sort_values(ascending=True).unique()
  otherids = otherallele(ids)
  pos = trios.iloc[i*4-1,0:5]
  start = pos.iloc[1]
  starttrimcm = np.interp(start,maps.iloc[:,3],maps.iloc[:,2])+0.5
  start = np.round(np.interp(starttrimcm,maps.iloc[:,2],maps.iloc[:,3]))
  end = pos.iloc[2]
  endtrimcm = np.interp(end,maps.iloc[:,3],maps.iloc[:,2])-0.5
  end = np.round(np.interp(endtrimcm,maps.iloc[:,2],maps.iloc[:,3]))
  newover = end-start
  newover = newover.astype(int)
  sepcount = np.array([[0,0], [0,0], [0,0]])
  ks = [0,1,2]
  for k in ks:
    patt = [3,2,1][k]
    pair = [1,2,3]
    del pair[patt-1]
    hapindex = np.append(pair,patt)
    hapids = ids[hapindex-1]
    dipindex = np.concatenate([reorder(hapindex[0]),reorder(hapindex[1]),reorder(hapindex[2])])
    dipids = np.append(ids,otherids) 
    dipids.sort()
    dipids = dipids[dipindex-1]
    vcfid = -(-hapids//1000)
    vcfid = sorted(set(vcfid))
    sqlhap = 'SELECT HAP'+str(hapids[0])+', HAP'+str(hapids[1])+', HAP'+str(hapids[2]) 
    sqldip = 'SELECT HAP'+str(dipids[0])+', HAP'+str(dipids[1])+', HAP'+str(dipids[2])+', HAP'+str(dipids[3])+', HAP'+str(dipids[4])+', HAP'+str(dipids[5])
    if len(vcfid)==1:
      sqlpart = ', POS FROM vcf'+str(vcfid[0])+', vcfinfo WHERE vcf'+str(vcfid[0])+'."index" = vcfinfo."index" AND vcfinfo.POS BETWEEN '+str(start)+' AND '+str(end)
    if len(vcfid)==2:
      sqlpart = ', POS FROM vcf'+str(vcfid[0])+', vcf'+str(vcfid[1])+', vcfinfo WHERE vcf'+str(vcfid[0])+'."index" = vcfinfo."index" AND vcf'+str(vcfid[1])+'."index" = vcfinfo."index" AND vcfinfo.POS BETWEEN '+str(start)+' AND '+str(end)
    if len(vcfid)==3:
      sqlpart = ', POS FROM vcf'+str(vcfid[0])+', vcf'+str(vcfid[1])+', vcf'+str(vcfid[2])+', vcfinfo WHERE vcf'+str(vcfid[0])+'."index" = vcfinfo."index" AND vcf'+str(vcfid[1])+'."index" = vcfinfo."index" AND vcf'+str(vcfid[2])+'."index" = vcfinfo."index" AND vcfinfo.POS BETWEEN '+str(start)+' AND '+str(end)
      
    sql2 = sqlhap + sqlpart
    sql1 = sqldip + sqlpart
    ## common variants
    vcfi = pd.read_sql_query(sql2, connex2)
    commontf = ((vcfi.iloc[:,0]==1)&(vcfi.iloc[:,1]==1)&(vcfi.iloc[:,2]==0))
    commoncount = commontf.sum() 
    if commoncount>0:
      index = np.where(commontf==True)[0]
      for m in range(0,commoncount):
        countpos = vcfi.loc[index[m],'POS']
        maf = pd.read_sql_query('SELECT MAF FROM maf WHERE POS = ' + str(countpos), connex2)
        maf = maf.loc[0, 'MAF']
        if maf <= 0.1:
          sepcount[k][0] = sepcount[k][0]+1
        elif maf > 0.1 and maf < 0.3:
          sepcount[k][1] = sepcount[k][1]+1
        else:
          sepcount[k][1] = sepcount[k][1]+0
    ## rare variants
    if use.phase.lower() == 'true':
      vcfi = pd.read_sql_query(sql2, connex1)
      raretf = ((vcfi.iloc[:,0]==1)&(vcfi.iloc[:,1]==1)&(vcfi.iloc[:,2]==0))
      sepcount[k][0] = sepcount[k][0] + raretf.sum() 
    else:
      vcfi = pd.read_sql_query(sql1, connex1)
      vcfgt = vcfi.iloc[:,0].map(str)+vcfi.iloc[:,1].map(str)+vcfi.iloc[:,2].map(str)+vcfi.iloc[:,3].map(str)+vcfi.iloc[:,4].map(str)+vcfi.iloc[:,5].map(str)
      c34 = (vcfgt=='111100') | (vcfgt=='111000') | (vcfgt=='110100') | (vcfgt=='011100')
      c34 = sum(c34)
      c20tf = (vcfgt=='101000') | (vcfgt=='100100') | (vcfgt=='011000') | (vcfgt=='010100')
      c20 = sum(c20tf)
      sepcount[k][0]=sepcount[k][0]+c34+c20
      if c20>0:
        index = np.where(c20tf==True)[0]
        id1 = ids[hapindex[0]-1]
        id2 = ids[hapindex[1]-1]
        def carrier(ibdid):
          otherid = otherallele(ibdid)
          ex = np.zeros(c20).astype(int)
          for m in range(c20):
            countpos = vcfi.loc[index[m],'POS']
            posibd = ibd[(ibd['hap1']==otherid) | (ibd['hap2']==otherid)]
            posibd = posibd[(posibd['starttrim']<=countpos) & (posibd['endtrim']>=countpos)]
            if len(posibd)>0:
              posid = posibd.iloc[:,0].append(posibd.iloc[:,1]).sort_values().unique()
              posid = posid[posid!=otherid]
              ## check it and it's other haplotype are not ibd with IBD haplotype
              checkibd = ibd[(ibd['hap1']==ibdid) | (ibd['hap2']==ibdid)] 
              checkibd = checkibd[(checkibd['starttrim']<=countpos) & (checkibd['endtrim']>=countpos)]
              checkid = checkibd.iloc[:,0].append(checkibd.iloc[:,1]).sort_values().unique() 
              otherposid = otherallele(posid)
              posid = posid[(np.isin(posid,checkid)+np.isin(otherposid,checkid))==0]
              if len(posid)>0:
                otherposid = otherallele(posid)  
                for j in range(len(posid)):
                  vcfidj = -(-posid[j]//1000)
                  sqlj = 'SELECT HAP'+str(posid[j])+', HAP'+str(otherposid[j])+', POS FROM vcf'+str(vcfidj)+', vcfinfo WHERE vcf'+str(vcfidj)+'."index" = vcfinfo."index" AND vcfinfo.POS = '+str(countpos)
                  vcfj = pd.read_sql_query(sqlj, connex1)
                  ex[m] = ex[m]+(vcfj.iloc[:,0]==1)|(vcfj.iloc[:,1]==1)
                  # exclude if either is 1
          return(ex)
        exclude = (carrier(id1)+carrier(id2))>0
        sepcount[k][0]=sepcount[k][0]-sum(exclude)


  mugc = pd.DataFrame(data = {'type':["12-3","13-2","23-1"],'rare':sepcount[:,0].tolist(),'common':sepcount[:,1].tolist(),'len':[newover,newover,newover]})

  print(trio.to_csv(index=False, header=False,sep=" ",float_format='%.4f'),end='')
  print(pd.DataFrame(pos).T.to_csv(index=False, header=False,sep=" ",float_format='%.0f'),end='')
  print(mugc.to_csv(index=False, header=False,sep=" "),end='')

connex1.commit()
connex1.close()
connex2.commit()
connex2.close()

