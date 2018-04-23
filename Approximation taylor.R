################################################################################
#                             Taylor approximation                             #
################################################################################

# Modulo:
modulo<-function(n,k)return(ifelse(k<0,n-ceiling(n/k)*k,n-floor(n/k)*k))
# Factorial:
factorial_relative<-function(n)return(factorial(abs(n))*ifelse(n>0,1,(-1)^n))

# Arrangement:
arrangement<-function(k,n){
  if(k>n)print("Error : k is larger than n")
  else return(factorial_relative(n)/factorial_relative(n-k))
}

# Polynom:
polynomial<-function(type,vector) return(function(x) ifelse("root"==type,prod(x-vector),sum(vector*x^(0:(length(vector)-1)))))
f<-polynomial("root",1:3)
g<-polynomial("coef",1:3)