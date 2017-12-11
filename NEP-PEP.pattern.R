c <- read.table("~/Desktop/C_matrix.txt",header=T,sep = "\t",stringsAsFactors = F,colClasses = c("character",rep("numeric",24)))
out <- array(NA,c(88,8))
c[,1] -> rownames(out)
colnames(out)<- c(0,1,2,4,5,7,10,14)
for( i in 1:nrow(c)){
  for (j in 0:7){
   out[i,j+1]<- mean(as.numeric(c[i,c(2,10,18)+j]),na.rm=T) 
  }
}
write.table(out,file="~/Desktop/mean.txt",sep = "\t",row.names = T,quote = F)
out <- out[,c(1:6)]
#id <- which(out[,1]< out[,2] & out[,1]< out[,3] & out[,3]>out[,4] & out[,2]>out[,4])
id <- which(out[,1]< out[,2] & out[,1] < out[,3] & out[,4]< out[,3] & out[,2]>out[,4])
out.p <- out[id,]
p1 <- rownames(out)[id]
id2 <- which(out[,4]< out[,6] & out[,1]< out[,4] )
p2 <- rownames(out)[id2]


# dir.create("~/Desktop/jy/")
# 
# for(i in 1:nrow(out.p)){
#   name<- paste("~/Desktop/jy/",rownames(out.p)[i],"pattern1.png",sep="")
#   png(name)
#   plot(out.p[i,],type="l",col="red",main=rownames(out.p)[i],xlab = "",ylab = "",frame.plot = F)
#   dev.off()
# }
# write.table(rownames(out.p),file="~/Desktop/list1.txt",sep = "\t",row.names = T,quote = F)
# write.table(p1,file="~/Desktop/jy/p1.txt",sep = "\t",row.names = T,quote = F)
# write.table(p2,file="~/Desktop/jy/p2.txt",sep = "\t",row.names = T,quote = F)
# 
# ##
#pdf(file = "~/Desktop/jy/test.pdf",width = 6,height = 6)
par(xpd=T,mar=c( 5.1,4.1,4.1,2.1))
x <- seq(0.5,5.5)
for(i in 1:nrow(out)){
  y<- jitter(rank(out[i,]))
  
  if(i ==1){
    #plot(0,frame.plot = F,type = "n",col=rgb(0,0,0,0.14),ylim =c(0,8),xlab = 1:7,ylab = 1:8,main = "")
    plot(x,y,frame.plot = F,type = "l",col=rgb(0,0,0,0.14),xlim = c(0,7),ylim =c(0,7),xlab = "",ylab = "",main = "",xaxt='n',yaxt='n',xaxs = "i",yaxs = "i")
  }else if(rownames(out)[i] %in% p1 ){
    points(x,y,type = "l",col=rgb(1,0,0,0.5),lwd=1)
  }else if(rownames(out)[i] %in% p2 ){
    points(x,y,type = "l",col=rgb(0,1,0,0.5),lwd=1)
  }else{
    points(x,y,type = "l",col=rgb(0,0,0,0.14),lwd=0.33)
    }
}

axis(1,at =0:6,labels = F,tck=-0.00,lwd = 0.33)
text(x=seq(0.5,5.5,1),labels=as.numeric(colnames(out)),y=rep(-0.5,8),srt=0,cex=0.8)
text(x=2.5,labels="Time (days)",y=-1,srt=0,cex=1)
par(xpd=T)
axis(2,at = 0:7,labels = F,tck=0.01,lwd=0.33)
text(x=-0.3,y=0:7,srt=90,labels = 0:7,cex=0.8)
text(x=-0.6,4,srt=90,labels = "Expression rank",cex=1)
legend(5.5,4.5,legend = c("NEP profile","PEP profile"),fill = c("red","green"),border="white",box.lwd = 0,box.col = "white",bg = "white")
#dev.off()

# 
# id <- which(out[,1]< out[,2] & out[,2]< out[,3] & out[,3]<out[,4] & out[,4]<out[,5])
# out.p <- out[id,]
# 
# for(i in 1:nrow(out.p)){
#   name<- paste("~/Desktop/jy/",rownames(out.p)[i],"pattern2.png",sep="")
#   png(name)
#   plot(out.p[i,],type="l",col="red",main=rownames(out.p)[i],xlab = "",ylab = "",frame.plot = F)
#   dev.off()
# }
# 
# ####all 
# dir.create("~/Desktop/jy/all/")
# for(i in 1:nrow(c)){
#   name<- paste("~/Desktop/jy/","all/",rownames(out)[i],"all.png",sep="")
#   png(name)
#   plot(out[i,],type="l",col="red",main=rownames(out)[i],xlab = "",ylab = "",frame.plot = F)
#   dev.off()
# }
