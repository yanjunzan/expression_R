#1. Function to screen the zero inflated expression trait.
# data is an expression data matrix,with each row as a individuals and each colom as a gene. the output is a list with two vector 
# num.ex stand for the number of genes expressing in each individual and gen.ex stands for the number of accession 
# in which a given gene is expressing. a candidate genelist is return if a lower and upper bundary is specified
zeroscan<-function(data,cutoff_lower=F,cutoff_upper=F){
  num.ex<-rep(NA,nrow(data))
  for( i in 1:nrow(data)){
    num.ex<-sum(which(data[i,]==0))
  }
  
  gen.ex<-rep(NA,nrow(data))
  for( j in 1:ncol(data)){
    gen.ex<-sum(which(data[,j]==0)) 
  }
  
  if(all(cutoff_lower,cutoff_upper)){
    if( cutoff_lower > cutoff_upper){
      cat(" cutoff_lower have to be greater than cutoff_upper")
    }else{
      genename<-colnames(data)[which(num.ex> cutoff_lower & num.ex<cutoff_upper)]
      return(list("num.ex"=num.ex,"gen.ex"=gen.ex,"trait"=genename))
    }
  }
  
  return(list("num.ex"=num.ex,"gen.ex"=gen.ex))
}

#2.Function to do genome wide scan and make mahattan plot for  loss of function allele
# input gwaadata is a GenABEL object, trait is the name of testing trait
scan_binary<-function(gwaadata,trait,plot=F){
  if(!require(GenABEL)){
    stop(" GenABEL has to be installed and loaded !") 
  }
  y<-gwaadata@phdata[,trait]
  y[which(y>0)]<-1
  qt<-qtscore(y,gwaadata,trait.type="binomial")
  
  if(plot==T){
    png(paste(trait,"_scan_loss.png",sep=""),width = 1024,height = 768)
    plot(qt,df="Pc1df") 
    thres<--log10(0.05/gwaadata@gtdata@nsnps)
    abline(h=thres,lty="dashed",col="red")
    dev.off()
  }
  return(qt)
}
#3.Function to do  permutation on the loss of function allele, 
#the return object is a matrix contain p-value adjusted for  genome wide test
permute_binary<-function(gwaadata,trait,times=1000,plot=T){
  if(!require(GenABEL)){
    stop(" GenABEL has to be installed and loaded !") 
  }
  y<-gwaadata@phdata[,trait]
  y[which(y>0)]<-1
  qt<-emp.qtscore(y,gwaadata,trait.type = "binomial",times = times)

  return(qt)
  if(plot==T){
    png(paste(trait,"permutation_loss.png",sep=""),width = 1024,height = 768)
    plot(qt) 
    abline(h=-log10(0.05),lty="dashed",col="red")
    dev.off()
  }
  return(qt)
}

##4. function to do genome wide scan and make mahattan plot for the second modifier allele

scan_normal<-function(gwaadata,trait,plot=F){
  y<-gwaadata@phdata[,trait]
  id<-which(y>0)
  idsubset<-gwaadata@gtdata@idnames[id]
  qt<-qtscore(zscore(y[id]),gwaadata,idsubset = idsubset)
  if(plot==T){
    png(paste(trait,"_scan_normal.png",sep=""),width = 1024,height = 768)
    plot(qt,df="Pc1df") 
    thres<--log10(0.05/gwaadata@gtdata@nsnps)
    abline(h=thres,lty="dashed",col="red")
    dev.off()
  }
  return(qt)
}
#5.Function to do  permutation on the second modifier allele, 
#the return object is a gwaa.scan data type contain p-value adjusted for genome wide test
permute_normal<-function(gwaadata,trait,times=1000,plot=F){
  y<-gwaadata@phdata[,trait]
  id<-which(y>0)
  idsubset<-gwaadata@gtdata@idnames[id]
  zscore<-function(x) qnorm((rank(x, na.last = "keep") - 0.5)/sum(!is.na(x)))
  qt<-emp.qtscore(zscore(y[id]),gwaadata,trait.type = "guess",idsubset = idsubset,times = times)  
  qt<-emp.qtscore(zscore(y[id]),gwaadata,trait.type = "guess",idsubset = idsubset,times = times)

  if(plot==T){
    png(paste(trait,"permutation_normal.png",sep=""),width = 1024,height = 768)
    plot(qt) 
    abline(h=-log10(0.05),lty="dashed",col="red")
    dev.off()
  }
  return(qt)
}
