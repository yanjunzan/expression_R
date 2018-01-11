c <- read.table("~/Dropbox/jy/C_matrix.txt",header=T,sep = "\t",stringsAsFactors = F,colClasses = c("character",rep("numeric",24)))
out <- array(NA,c(88,8))
c[,1] -> rownames(out)
colnames(out)<- c(0,1,2,4,5,7,10,14)
for( i in 1:nrow(c)){
  for (j in 0:7){
    out[i,j+1]<- mean(as.numeric(c[i,c(2,10,18)+j]),na.rm=T) 
  }
}
out.c <- out
id <- which(out.c[,1]< out.c[,2] & out.c[,1] < out.c[,3] & out.c[,4]< out.c[,3] & out.c[,2]>out.c[,4])
out.p <- out.c[id,]
p1 <- rownames(out.c)[id]

id2 <- which(out.c[,4]< out.c[,6] & out.c[,1]< out.c[,4] )
p2 <- rownames(out.c)[id2]

n <- read.table("~/Dropbox/jy/N_matrix.txt",header=T,sep = "\t",stringsAsFactors = F,colClasses = c("character",rep("numeric",24)))
out <- array(NA,c(23910,8))
n[,1] -> rownames(out)
colnames(out)<- c(0,1,2,4,5,7,10,14)
for( i in 1:nrow(n)){
  for (j in 0:7){
    out[i,j+1]<- mean(as.numeric(n[i,c(2,10,18)+j]),na.rm=T) 
  }
}
out.n <- out
# id
id <- read.table("~/Dropbox/jy/match.id.txt",stringsAsFactors = F,sep="\t",header = T)
#id <- gsub(pattern = "(AT.*)\\.*",replacement = "\\1",x = id)
out.n1 <- out.n[unique(id$ID),]
# 
ch <- cbind.data.frame(out.c,rep("ch",nrow(out.c)))
nc <- cbind.data.frame(out.n1,rep("nc",nrow(out.n1)))
lg <- function(x){
  return(log10(x))
}


get.log10 <- function(mat){
  return(apply(mat, 2, FUN= function(x) log(x = x,base = 2)))
}
ch.lg <- get.log10(mat = data.matrix(ch)+1)[,1:6]

nc.lg <- get.log10(mat = data.matrix(nc)+1)[,1:6]

gene.rm <- rownames(nc.lg)[ which(nc.lg[,3] <6)]

par(xpd=T,mar=c( 5.1,4.1,4.1,2.1))
x <- seq(0.5,5.5)

plot(1, type="n",xlim = c(0,6),ylim =c(-1,22),xlab = "",ylab = "",main = "",xaxt='n',frame.plot = F)

#plot(,y,frame.plot = F,type = "l",col=rgb(0,0,0,0.14),xaxs = "i",yaxs = "i")
legend(1.5,4,legend = c("Chloroplast encoded genes_NEP","Chloroplast encoded genes_PEP","Nuclear encoded genes"),fill = c("red","green","blue"),border="white",box.lwd = 0,box.col = "white",bg = "white")

for(i in 1:nrow(ch.lg)){
      if(rownames(ch.lg)[i] %in% p1){
        points(c(0:5),ch.lg[i,c(1:6)],type = "l",col=rgb(1,0,0,0.5),lwd=1)
      }else if(rownames(ch.lg)[i] %in% p2){
        points(c(0:5),ch.lg[i,c(1:6)],type = "l",col=rgb(0,1,0,0.5),lwd=1)
        
      }else{
        #points(c(0:7),ch.lg[i,c(1:8)],type = "l",col=rgb(0,0,0,0.14),lwd=1)
      }

}


for(i in 1:nrow(nc.lg)){
  if(!(rownames(nc.lg)[i] %in% gene.rm))
  points(c(0:5),nc.lg[i,c(1:6)],type = "l",col=rgb(0.2,0.2,1,0.3),lwd=1)
}

axis(1,at =0:5,labels = F,tck=-0.015,lwd = 1,pos = -0.5)
text(x=seq(0.5,5.5,1)-0.5,labels=as.numeric(colnames(out.c))[1:6],y=rep(-1.3,8),srt=0,cex=0.8)
par(xpd=T)
text(x=2.5,labels="Time (days)",y=-1.9,srt=0,cex=1)

#axis(2,at = 0:7,labels = F,tck=0.01,lwd=0.33)
#text(x=-0.3,y=0:7,srt=90,labels = 0:7,cex=0.8)
text(x=-0.8,12,srt=90,labels = "Log2(Normalised expression)",cex=1)

# mean 
ch.lg <- get.log10(mat = data.matrix(ch)+1)[,1:6]
ch.lg.p1 <- ch.lg[p1,]
ch.lg.p2 <- ch.lg[p2,]

nc.lg <- get.log10(mat = data.matrix(nc)+1)[,1:6]

gene.rm <- rownames(nc.lg)[ which(nc.lg[,3] <6)]

nc.lg.re <- nc.lg[rownames(nc.lg)[ which(nc.lg[,3] >=6)],]

ch.lg.p1.m <- apply(ch.lg.p1, 2, mean)
ch.lg.p2.m <- apply(ch.lg.p2, 2, mean)
nc.lg.re.m <- apply(nc.lg.re, 2, mean)

par(xpd=T,mar=c( 5.1,4.1,4.1,2.1))
x <- seq(0.5,5.5)
plot(1, type="n",xlim = c(0,6),ylim =c(-1,15),xlab = "",ylab = "",main = "",xaxt='n',frame.plot = F)

#plot(,y,frame.plot = F,type = "l",col=rgb(0,0,0,0.14),xaxs = "i",yaxs = "i")
legend(1.8,5.3,legend = c("Chloroplast encoded genes_NEP","Chloroplast encoded genes_PEP","Nuclear encoded genes"),fill = c("red","green","blue"),border="white",box.lwd = 0,box.col = "white",bg = "white")

points(c(0:5),ch.lg.p1.m,type = "l",col=rgb(1,0,0,1),lwd=2)
points(c(0:5),ch.lg.p2.m,type = "l",col=rgb(0,1,0,1),lwd=2)
points(c(0:5),nc.lg.re.m,type = "l",col=rgb(0.2,0.2,1,1),lwd=2)

axis(1,at =0:5,labels = F,tck=-0.015,lwd = 1,pos = -0.5)
text(x=seq(0.5,5.5,1)-0.5,labels=as.numeric(colnames(out.c))[1:6],y=rep(-1.3,8),srt=0,cex=0.8)
par(xpd=T)
text(x=2.5,labels="Time (days)",y=-1.9,srt=0,cex=1)

#axis(2,at = 0:7,labels = F,tck=0.01,lwd=0.33)
#text(x=-0.3,y=0:7,srt=90,labels = 0:7,cex=0.8)
text(x=-0.8,7.5,srt=90,labels = "Log2(Normalised expression)",cex=1)


## old 
# ch.lg <- get.log10(mat = data.matrix(ch)+1)  
# 
# nc.lg <- get.log10(mat = data.matrix(nc)+1)  
# 
# par(xpd=T,mar=c( 5.1,4.1,4.1,2.1))
# x <- seq(0.5,5.5)
# 
# plot(1, type="n",xlim = c(0,9),ylim =c(-1,22),xlab = "",ylab = "",main = "",xaxt='n',frame.plot = F)
# 
# #plot(,y,frame.plot = F,type = "l",col=rgb(0,0,0,0.14),xaxs = "i",yaxs = "i")
# 
# for(i in 1:nrow(ch.lg)){
#   if(rownames(ch.lg)[i] %in% p1){
#     points(c(0:7),ch.lg[i,c(1:8)],type = "l",col=rgb(1,0,0,0.5),lwd=1)
#   }else if(rownames(ch.lg)[i] %in% p2){
#     points(c(0:7),ch.lg[i,c(1:8)],type = "l",col=rgb(0,1,0,0.5),lwd=1)
#     
#   }else{
#     #points(c(0:7),ch.lg[i,c(1:8)],type = "l",col=rgb(0,0,0,0.14),lwd=1)
#   }
#   
# }
# 
# 
# for(i in 1:nrow(nc.lg)){
#   
#   points(c(0:7),nc.lg[i,c(1:8)],type = "l",col=rgb(0.2,0.2,1,0.3),lwd=1)
#   
# }
# 
# axis(1,at =0:7,labels = F,tck=-0.02,lwd = 1,pos = -0.5)
# text(x=seq(0.5,7.5,1)-0.5,labels=as.numeric(colnames(out.c))[1:8],y=rep(-1.3,8),srt=0,cex=0.8)
# par(xpd=T)
# text(x=3.5,labels="Time (days)",y=-1.9,srt=0,cex=1)
# 
# #axis(2,at = 0:7,labels = F,tck=0.01,lwd=0.33)
# #text(x=-0.3,y=0:7,srt=90,labels = 0:7,cex=0.8)
# text(x=-1,12,srt=90,labels = "Log2(absolute expression)",cex=1)
# legend(0,26.4,legend = c("Chloroplast encoded genes_NEP","Chloroplast encoded genes_PEP","Chloroplast encoded genes_other","Nuclear encoded genes"),fill = c("red","green","black","blue"),border="white",box.lwd = 0,box.col = "white",bg = "white")
# 
