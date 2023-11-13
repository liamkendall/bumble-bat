###functions

#basic summary for binary data
#can handle non-integers  (by rounding) 

binary.summary=function (x,n){
  
  x=round(n*x)
  out=c(rep(1,x),rep(0,n-x))
  mean=mean(out)
  sd=sd(out)
  se=sd(out)/sqrt(n)
  
  return(cbind(n,mean,sd,se))
  
}


