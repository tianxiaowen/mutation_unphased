## usage
## Rscript search.eps2.R prefix
## prefix is the prefix of .like file
## for example, if .like file is saved in /projects/browning/tmp/fhs.chrX.like, prefix=/projects/browning/tmp/fhs

arg <- commandArgs(TRUE)
prefix <- arg[[1]]
n<-22
l.i<-read.table('chr1/mle_epsilon2.like',header=T)
L<-array(NA,c(nrow(l.i),n))
for (i in 1:n){
  L[,i]<-read.table(paste('chr',i,'/mle_epsilon2.like',sep=""),header=T)$log.likelihood
}
l<-rowSums(L)
mle<-l.i[l==max(l),1:4]

boot<-function(){
  index<-sample(1:n,n,replace=T)
  L.star<-L[,index]
  l.star<-rowSums(L.star)
  mle.index<-which(l.star==max(l.star))
  if (length(mle.index)>1){mle.index<-sample(mle.index,1)}
  mle.star<-l.i[mle.index,1:4]
  return(mle.star)
}
set.seed(99)
nboot=10000
do.boot<-replicate(nboot,boot())
ci.err1<-quantile(unlist(do.boot[1,1,]),c(0.025,0.975))
ci.err2<-quantile(unlist(do.boot[1,2,]),c(0.025,0.975))
ci.gc<-quantile(unlist(do.boot[1,3,]),c(0.025,0.975))
ci.mu<-quantile(unlist(do.boot[1,4,]),c(0.025,0.975))


print(mle)
print(ci.err1)
print(ci.err2)
print(ci.gc)
print(ci.mu)
