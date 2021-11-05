#proportion.to.select
#s=proportion.to.select
CI<-function(x,y,s=s,top=T){
  
    x2<-as.data.frame(x)
    y2<-as.data.frame(y)
    size<-ceiling(nrow(x2)*s)
    x2<-x2[order(x2[,2], decreasing=top),]
    y2<-y2[order(y2[,2], decreasing=top),]
    both<-sum(x2[1:size,1]%in%y2[1:size,1])
    random<-both*s
    ci <-(both-random)/(size-random)
 
  return(ci)
}


#CI <- (B-R)/(T-R)
