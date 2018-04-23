################################################################################
#              NONPARAMETRIC ESTIMATION OF RISK_NEUTRAL DENSITY                #
################################################################################

######################### LOCAL POLYNOMIAL REGRESSION ##########################
################################################################################

# Data:
dir1<-"~/Devoirs/Mémoire/Data"
setwd(dir1)
data2<-read.table("Data3.txt",header=T,sep="\t",dec=",")
data1<-read.table("Data1.txt",header=T,sep="\t",dec=",")
data<-read.table("Data2.txt",header=T,sep="\t",dec=",")
S<-data[,1];Y<-data[,2];Z<-data[,3]
S1<-data1[,1];Y1<-data1[,4];Z1<-data1[,5];Strike<-data1$Strike;Coef<-data1[,3]
S2<-data2[,1];Y2<-data2[,3];Z2<-data2[,4];Strike<-data2$Strike
n<-length(Z);n1<-length(Z1);n2<-length(Z2)
z<-seq(min(Z),max(Z),length=n)
View(data)
n<-nrow(Z);d<-ncol(Z)
Z_min<-Z_max<-vector()
for(k in 1:d){
Z_max[k]<-max(Z[,k]);Z_min[k]<-min(Z[k])
}
z1<-seq(Z_min[1],Z_max[1],length=100)
z2<-seq(Z_min[2],Z_max[2],length=100)
z3<-seq(Z_min[3],Z_max[3],length=100)
z<-cbind(z1,z2,z3)

# Integrate:
int_function<-function(fun,lower,upper){
  if(lower>upper)stop("Bounded no-adapted")
  a<-ifelse(lower!=-Inf,lower,-10^3)
  b<-ifelse(upper!=Inf,upper,10^3)
  matx<-as.matrix(seq(a,b,length=10^6))
  vecy<-apply(matx,1,fun)
  stdx<-matx[2,1]-matx[1,1]
  return(sum(stdx*vecy[-length(vecy)]))
}

# Square integrated:
R<-function(fun,lower,upper){
  local_function<-function(x) (fun(x))^2
  return(int_function(local_function,lower,upper))
}

# Moments:
mu_n<-function(fun,order,lower,upper){
  local_function<-function(x) x^order*fun(x)
  return(int_function(local_function,lower,upper))
}

# Convolution:                                                                  # Produit de convolution
Conv<-function(fun_1,fun_2,x){
  fun_local<-function(y) fun_1(x-y)*fun_2(y)
  return(int_function(fun_local,-Inf,Inf))
}

# Combinatorial Coefficients:
kAn<-function(k,n) ifelse(k>n,0,factorial(n)/factorial(n-k))                    # Arrangement sans répétition
kCn<-function(k,n) kAn(k,n)/factorial(k)                                        # Combinaison sans répétition
kGn<-function(k,n) kCn(k,n+k-1)                                                 # Combinaison avec répétition

# Characteristic functions:
caract_left<-function(x)   ifelse((x[1]>=x[2])&&(x[1]<x[3]),1,0)                # Borne gauche appartient
caract_right<-function(x)  ifelse((x[1]>x[2])&&(x[1]<=x[3]),1,0)                # Borne droite appartient
caract_strict<-function(x) ifelse((x[1]>=x[2])&&(x[1]<x[3]),1,0)                # Aucune borne n'appartient
caract_large<-function(x)  ifelse((x[1]>=x[2])&&(x[1]<=x[3]),1,0)               # Les deux bornes appartiennent

# Univariate kernel functions:
Uniform<-function(x)       .5*caract_large(c(x,-1,1))
Triangle<-function(x)      (1-abs(x))*caract_large(c(x,-1,1))
Epanechnikov<-function(x)  .75*(1-x^2)*caract_large(c(x,-1,1))
Quartic<-function(x)       (15/16)*(1-x^2)^2*caract_large(c(x,-1,1))
Triweight<-function(x)     (35/32)*(1-x^2)^3*caract_large(c(x,-1,1))
Gaussian<-function(x)      (1/sqrt(2*pi))*exp(-.5*x^2)
Cosinus<-function(x)       (pi/4)*cos(pi*x/2)*caract_large(c(x,-1,1))

# List of univariate Kernel functions:
Kernel_functions<-c(Uniform,Triangle,Epanechnikov,Quartic,Triweight,Gaussian,
                    Cosinus)

# Kernel density estimator:
kernel_density<-function(x)return((1/(n*h))*sum(apply(as.matrix((x-sample)/h),2,
                                                      Gaussian)))

# Univariate Local Polynomial Regression:
LPM_smooth<-function(x,p,Kernel_function=Gaussian){
  if(p<0) stop("L'ordre de la régression est incorrecte")
  local<-matrix(NA,length(X),p+1)
  for(j in 1:(p+1))local[,j]<-(X-x[1])^(j-1)
  W<-diag(sapply((X-x[1])/h,Kernel_function))
  return(solve(t(local) %*% W %*% local, t(local) %*% W %*% Yy))
}

# Simulations:
set.seed(1)
n<-1000
X<-sort(sample(runif(n^2),n))
Zz<-sort(sample(rexp(n^2),n))
zz<-seq(min(Zz),max(Zz),length=n)
Obs<-cbind(X,Zz)
o<-cbind(x,zz)
sigma<-sd(X)
m<-function(x) sin(2*pi*(x^3)/3)
u<-function(x) sin(2*pi*(x[1]^3)/3)+sqrt(.5)*exp(-.5*x[2]^2)
Yy<-m(X)+sigma*rnorm(n)
U<-u(Obs)+sigma*rnorm(n)
h<-sigma/2
par(mfrow=c(4,2))
X<-Z;Yy<-Y;x<-z
h<-sd(X)/10

# Visualization:
p<-3
LPM_list<-list(NA)
for (k in 1:length(Kernel_functions)){
  LPM_list[[k]]<-matrix(NA,p+1,length(x))
  for (j in 1:length(x)){
    LPM_list[[k]][,j]<-LPM_smooth(x[j],p,Kernel_functions[[k]])
  }
}
plot(X,Yy,pch=20)
lines(x,m(x),col=2,lwd=2)
for (k in 1:length(Kernel_functions)){
  plot(x,LPM_list[[k]][3,],col=2,lwd=2,type='l')
}

# Multivariate kernel function:                                                 # Sans observation de X 
MK<-function(x,Kernel_function=Gaussian) prod(sapply(x,Kernel_function))

# Multivariate kernel density estimator:
h<-110
# Multivariate LPR estimator:
MLPM<-function(z,Kernel_function=Gaussian){
  f<-function(z) MK(z,Kernel_function=Gaussian)
  W<-diag(apply(t(z-t(Z))/h,1,f))
  X<-rep(1,nrow(Z))
  L<-t(t(Z)-z)
  X<-cbind(X,L)
  M<-list(L)                                                                    # Les sous-matrices de X
  M[[2]]<-L[,1]*M[[1]]
  for(j in 2:3) M[[2]]<-cbind(M[[2]],L[,j]*M[[1]][,j:3])
  X<-cbind(X,M[[2]])
  M[[3]]<-L[,1]*M[[2]]
  M[[3]]<-cbind(M[[3]],L[,2]*M[[2]][,4:kGn(2,3)])                               # k-ième sous-matrice de X
  M[[3]]<-cbind(M[[3]],L[,3]*M[[2]][,6:kGn(2,3)])
  X<-cbind(X,M[[3]])
  return(solve(t(X) %*% W %*% X, t(X) %*% W %*% Y))                             # Vecteur des betas estimés
}
det(t(X) %*% W %*% X)^(-1)
dim(t(X))
t(matrix(1:9,3,3))-(1:3)

# Visualization:
MLPM_vec<-vector()
for(i in 1:nrow(z)) MLPM_vec[i]<-MLPM(z[i,])[5]

hCV<-regCVBwSelC(Z,Y,deg=3,gaussK)*90
res<-locpoly(Z,Y,drv=2,degree=3,kernel="normal",bandwidth=hCV)
plot(res$x,res$y,type='l',ylab="q",xlab=expression(S[T]),col="blue",lwd=2)

hCV1<-regCVBwSelC(Z1,Y1,deg=3,gaussK)
res1<-locpoly(Z1,Y1,drv=1,degree=3,kernel="normal",bandwidth=hCV1)
res2<-locpoly(Z1,Y1,drv=2,degree=3,kernel="normal",bandwidth=hCV1)
x1<-res1$x;v1<-res1$y;x2<-res2$x;v2<-res2$y

hCVvi<-regCVBwSelC(Z2,Y2,deg=3,gaussK)
resvi0<-locpoly(Z2,Y2,drv=0,degree=3,kernel="normal",bandwidth=hCVvi)
resvi1<-locpoly(Z2,Y2,drv=1,degree=3,kernel="normal",bandwidth=hCVvi)
resvi2<-locpoly(Z2,Y2,drv=2,degree=3,kernel="normal",bandwidth=hCVvi)
xvi0<-resvi0$x;xvi1<-resvi1$x;xvi2<-resvi2$x;
yvi0<-resvi0$y;yvi1<-resvi1$y;yvi2<-resvi2$y

plot(xvi0,yvi0,type='l',ylab=expression(sigma),xlab="M",col="blue",lwd=2)
par(mfrow=c(1,2))
plot(xvi1,yvi1,type='l',ylab=expression(sigma),xlab="M",col="blue",lwd=2)
plot(xvi2,yvi2,type='l',ylab=expression(sigma),xlab="M",col="blue",lwd=2)

par(mfrow=c(1,1))
hCVi<-regCVBwSelC(Z,Y,deg=3,gaussK)*141.53
resvi<-locpoly(Z,Y,drv=2,degree=3,kernel="normal",bandwidth=hCVi)
plot(resvi$x,resvi$y,type='l',ylab="q",xlab=expression(S[T]),col="blue",lwd=2)

########################## Tikhonov regularization #############################
################################################################################
par(mfrow=c(1,1))
# Visualization of error's explosion into the PLMP:
n_ex<-100;a_ex<-0;b_ex<-6;step_ex<-(b_ex-a_ex)/(n_ex-1)
leg1.txt<-c("Données sans bruit","Observations")
leg2.txt<-c("Dérivée réelle","Dérivée approximative")
fun<-function(x) sin(x)
fun_d<-function(x) cos(x)
u_ex<-seq(a_ex,b_ex,length=n_ex)
v_ex<-fun(u_ex)+rnorm(length(u_ex),sd=.005)
dlt<-max(abs(v_ex-fun(u_ex)))
dv1_ex<-(v_ex[-c(1,2)]-v_ex[-c(n_ex-1,n_ex)])/(2*step_ex)
dv_ex<-c((v_ex[2]-v_ex[1])/step_ex,dv1_ex,(v_ex[n_ex]-v_ex[n_ex-1])/step_ex)
par(mfrow=c(1,2))
plot(u_ex,fun(u_ex),type="l",xlab="",ylab="")
lines(u_ex,v_ex,col=4)
legend("topright",leg1.txt,col = c(1,4),lty=c(1,1))
plot(u_ex,fun_d(u_ex),type="l",xlab="",ylab="")
lines(u_ex,dv_ex,col="red")
legend("topright",leg2.txt,col=1:2,lty=c(1,1))

sd_ex<-seq(0.01,3,by=.01)
total_error<-((sd_ex^2)/6)+(dlt/sd_ex)
sd_ex_m<-(3*dlt)^(1/3);n_ex_m<-floor(1+(b_ex-a_ex)/sd_ex_m)
u_ex_m<-seq(a_ex,b_ex,by=sd_ex_m)
v_ex_m<-fun(u_ex_m)+rnorm(n_ex_m,sd=.005)
dv1_ex_m<-(v_ex_m[-c(1,2)]-v_ex_m[-c(n_ex_m-1,n_ex_m)])/(2*sd_ex_m)
dv_ex_m<-c((v_ex_m[2]-v_ex_m[1])/sd_ex_m,dv1_ex_m,
           (v_ex_m[n_ex_m]-v_ex_m[n_ex_m-1])/sd_ex_m)
plot(sd_ex,total_error,type="l",ylab="erreur totale",
     xlab="pas de discrétisation")
points(sd_ex_m,((sd_ex_m^2)/6)+(dlt/sd_ex_m),col="blue",pch=20)
plot(u_ex,fun_d(u_ex),type="l",xlab="",ylab="")
lines(u_ex_m,dv_ex_m,col="red")
legend("topright",leg2.txt,col=1:2,lty=c(1,1))

al<-length(Z)-(2*7)
pc<-exp(-2*0.0054)*(Z[3]-Z[2])/6

bounds<-sort(c((1:7)*25,((1:7)*25)-1))
Cbut<-(Y[-(1:2)]-2*Y[-c(1,length(Y))]+Y[-c(length(Y)-1,length(Y))])/50
Cbut<-Cbut[-bounds[-c(length(bounds)-1,length(bounds))]]
H<-matrix(0,al,al)
for(i in 1:nrow(H)){
  for(j in 1:ncol(H)) H[i,j]<-ifelse(j==i+1,1,ifelse(j==i-1,1,ifelse(i==j,4,0)))
}
H<-pc*H
alp<-2
lamb<--solve(H + diag(alp,al,al),Cbut)

Kb<-Z[-bounds]
qalpha<-function(x){
  vecteurH<-c()
  for(i in 1:length(Z)){
    vecteurH[i]<-ifelse(max(i==bounds)==1,NA,caract_right(c(x,Z[i],Z[i+1]))*
                      (x-Z[i])+caract_right(c(x,Z[i+1],Z[i+2]))*(Z[i+2]-x))/50
  }
  vecteurh<-exp(-0.0054)*na.omit(vecteurH)
  return(-sum(lamb*vecteurh))
}

vecteur.test<-seq(min(Z)/2,max(Z)*2,by=95)
qreg1<-sapply(vecteur.test,qalpha)
plot(vecteur.test,qreg1,type='l',ylab="q",xlab=expression(S[T]),col="blue",lwd=2)
