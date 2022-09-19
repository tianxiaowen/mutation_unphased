## remove mutation count greater and equal to 50
gettrio<-function(x){
  trio7<-x*7
  trio1<-trio7-6
  trio1:trio7
}

allcount<-list()

for (i in 1:22){
  #nj<-10
  #nj<-11
  cmaxs<-list()
  trios<-read.table(paste0("chr",i,"/triosgt.maf0.25.2.5-6cM.txt"),fill=T,stringsAsFactors=F)
  len<-nrow(trios)/7
  cmax<-numeric(len)
  for(m in 1:len){
    count<-trios[(m*7-2):(m*7),2:3]
    cmax[m]<-max(count)
  }
  keeptrio<-which(cmax<50)
  if(length(keeptrio)<len){
    print(paste0('chr',i))
    keeprow<-sapply(keeptrio,gettrio)
    finaltrios<-trios[keeprow,]
    write.table(trios,paste0("chr",i,"/triosgt.maf0.25.2.5-6cM.old50.txt"),quote=F,col.names=F,row.names=F,na="")
    write.table(finaltrios,paste0("chr",i,"/triosgt.maf0.25.2.5-6cM.txt"),quote=F,col.names=F,row.names=F,na="")
  }
}

 
