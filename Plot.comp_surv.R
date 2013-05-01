##############################
##############################
## EXPLORE SHAPE OF FUNCTION LINKING COMPETITIVE HIERARCHY AND MORTALITY ALONG CLIMATIC GRADIENT

## function describing competitive interaction and survival along climatic gradient
pdf("./figs/compet_hierarchy.pdf")
par(mfrow=c(1,2))
plot(seq(-1,1,length=100),1/(1+exp(-10*(seq(-1,1,length=100)))),type="l",xlab="Diff in competitive ability c",ylab="Probability of wins")
text(-0.75,1,labels="k=10")
c <-  1
plot(seq(0,1,length=100), 1-c*seq(0,1,length=100),type="l",lty=1,xlab="Climat severity gradient",ylab="Probabilty of survival")
text(0.9,1-c*0.9,labels=paste("c=",c,sep=""))
c <-  0.75
lines(seq(0,1,length=100), 1-c*seq(0,1,length=100),lty=2)
text(0.9,1-c*0.9,labels=paste("c=",c,sep=""))

c <-  0.25
lines(seq(0,1,length=100), 1-c*seq(0,1,length=100),lty=3)
text(0.9,1-c*0.9,labels=paste("c=",c,sep=""))
c <-  0.0
lines(seq(0,1,length=100), 1-c*seq(0,1,length=100),lty=4)
text(0.9,1-c*0.9,labels=paste("c=",c,sep=""))
dev.off()

## resulting trade off
pdf("./figs/trade_off.pdf")
c <-  seq(from=0,to=1,length=100)
plot(c,1-c,xlab="Competitive ability",ylab="Survival in extrem climate",type="l")
dev.off()
