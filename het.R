s<-0
maf<-c(0.1,0.249999)
thetas<-matrix(NA,nrow=length(maf),ncol=22)
for (i in 1:22){
  plink<-read.table(paste0("chr",i,"/chr",i,".plink.hwe"),header=T)
  colid<-rep('NULL',15)
  colid[2]<-'numeric'
  colid[7]<-'character'
  colid[14]<-'numeric'
  gtstats<-read.table(paste0("chr",i,"/chr",i,".gtstats"),colClasses=colid)
  plink<-cbind(plink,gtstats) 
  dup.pos<-names(which(table(plink$V2)>1))
  dup.index<-which(plink$V2%in%dup.pos) 
  if(length(dup.index)>0){
     plink<-plink[-dup.index,]
  }
  plink<-plink[plink$V7=='PASS',]
  plink<-plink[plink$A1=='A'|plink$A1=='T'|plink$A1=='C'|plink$A1=='G',] 
  plink<-plink[plink$A2=='A'|plink$A2=='T'|plink$A2=='C'|plink$A2=='G',]
  
  s<-s+plink[nrow(plink),'V2']-plink[1,'V2']
  
  het<-numeric(length(maf))
  for(j in 1:length(maf)){
    if(j==1){
      lwr<-0
    }else{
      lwr<-maf[j-1]
    }
    upr<-maf[j]
    het[j]<-sum(plink[plink$V14>lwr&plink$V14<=upr,]$O.HET.)
  }
  thetas[,i]<-het
}
theta<-rowSums(thetas)/s
save.image("het.Rdata")


