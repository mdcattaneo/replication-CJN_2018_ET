################################################################################
## R CODE FOR CATTANEO-JANSSON-NEWEY (2014): Many Series Terms -- ET paper
## DATE: 15-Dec-2014
################################################################################
################################################################################
# TO RUN: source("main.R");
################################################################################
rm(list=ls(all=TRUE))
library(MASS);library(Hmisc);

################################################################################
## QR-based (X'X)^(-1)
################################################################################
#XXinv = function(x) tcrossprod(solve(qr.R(qr(x))))
XXinv = function(x) chol2inv(chol(crossprod(x)))

################################################################################
## DGP
################################################################################
## Function g(z)
g = function(Z) apply(Z,1,function(x)(exp(sqrt(crossprod(x)))))

## Function h(z)
h = function(Z) apply(Z,1,function(x)(exp(sqrt(crossprod(x)))))

## Function sigma(z)
sigma.v.fn = function(Z,gamma) return(abs(1+rowSums(Z))^gamma);

## Distributions
gen.Z = function(n,dz) return(matrix(runif(n*dz,-1,1),nrow=n,ncol=dz))

rmixnorm = function(n,m1,m2,s1,s2,alpha) {
    I = runif(n)<alpha; rnorm(n,mean=ifelse(I,m1,m2),sd=ifelse(I,s1,s2));}

dmixnorm = function(x,m1,m2,s1,s2,alpha) return(alpha * dnorm(x,m1,s1) + (1-alpha)*dnorm(x,m2,s2))

dgp.fn = function(m,n){
    ## MODELS 1-2: GAUSSIAN
    if (m==1){u = rnorm(n);
              v = rnorm(n); gamma.v = 0;}
    if (m==2){u = rnorm(n);
              v = rnorm(n); gamma.v = 1;}
    ## MODELS 3-4: ASYMMETRIC
    if (m==3){u = rmixnorm(n,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2);
              v = rmixnorm(n,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2); gamma.v = 0;}
    if (m==4){u = rmixnorm(n,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2);
              v = rmixnorm(n,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2); gamma.v = 1;}
    ## MODELS 5-6: BIMODAL
    if (m==5){u = rmixnorm(n,m1=-3/2,m2=3/2,s1=1/2,s2=1,alpha=1/2);
              v = rmixnorm(n,m1=-3/2,m2=3/2,s1=1/2,s2=1,alpha=1/2); gamma.v = 0;}
    if (m==6){u = rmixnorm(n,m1=-3/2,m2=3/2,s1=1/2,s2=1,alpha=1/2);
              v = rmixnorm(n,m1=-3/2,m2=3/2,s1=1/2,s2=1,alpha=1/2); gamma.v = 1;}

    return(list(v=v,u=u,gamma.v=gamma.v))
}

## Figure Distributions for Paper
postscript(file=paste0("Figure1.ps"), horizontal=F)
    par(mfrow=c(3,1))
    curve(dnorm, from=-3, to=3, main="Normal Distribution")
    fun = function(x) dmixnorm(x,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2); curve(fun, from=-3, to=3, main="Asymmetric Distribution")
    fun = function(x) dmixnorm(x,m1=-3/2,m2=3/2,s1=1/2,s2=1/2,alpha=1/2); curve(fun, from=-3, to=3, main="Bimodal Distribution")
dev.off()

################################################################################
## BASIS OF APPROXIMATION
################################################################################
gen.P = function(Z,K) {
    if (K==0)   out = NULL;
    if (K==1)   out = poly(Z,degree=1,raw=TRUE);
    if (K==2)  {out = poly(Z,degree=1,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^2);}
    if (K==2.5) out = poly(Z,degree=2,raw=TRUE);
    if (K==3)  {out = poly(Z,degree=2,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^3);}
    if (K==3.5) out = poly(Z,degree=3,raw=TRUE);
    if (K==4)  {out = poly(Z,degree=3,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^4);}
    if (K==4.5) out = poly(Z,degree=4,raw=TRUE);
    if (K>=5)  {out = poly(Z,degree=4,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^5);}
    if (K==5.5) out = poly(Z,degree=5,raw=TRUE);
    if (K>=6)  {out = poly(Z,degree=5,raw=TRUE); for (k in 6:K) for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^k);}
    ## RETURN POLYNOMIAL BASIS
    return(out)
}

################################################################################
## MONTECARLO SETUP
################################################################################
models = 6
## Number of simulations (S)
S = 5000
## Sample size (n)
n = 500
## dimensions & parameters
d = 1; beta = 1; dz = 5
## K grid; K = dim(gamma)
K.grid = c(1,2,2.5,3,3.5,4,4.5,5,5.5,6:10)

## Output Tables
## 
## SE0 = HOM/n | SE1 = HOM/(n-K-d)
col.names=c("m","n","d","beta","s","K.id","K","beta.hat","SE0","SE1")
out = list(NULL)
for (m in 1:models) out[[m]] = matrix(NA, nrow=S*length(K.grid), ncol=length(col.names), dimnames=list(NULL,col.names))

################################################################################
# Run Monte Carlo Experiment
################################################################################
## Heteroskedasticity Scales
if (TRUE){
    message("\nComputing Scaling Constants.\n")
    dimnames = list(NULL,c("m","scale.sigma.v"))
    popval = matrix(NA,nrow=models,ncol=length(dimnames[[2]]),dimnames=dimnames)
    I = 1000000; Z = gen.Z(n=I,dz=dz);

    for (m in 1:models){
        dgp = dgp.fn(m=m,n=I);
        popval[m,] = c(m,sqrt(as.numeric(crossprod(sigma.v.fn(Z,dgp$gamma.v)*dgp$v)/I)));
    }
    write.csv(popval, file=paste0("output/popval.csv"))
}
popval = read.csv(paste0("output/popval.csv"), row.names=1)

message("Simulation began. Time: ", Sys.time(),"\n")
showevery=1000
row=1;
# m=1; s=1; K.i=5;
for (s in 1:S) {
    if (max(s==seq.int(1000,S,showevery))==1) {message("Simulations Completed: ",s," of ",S," - ", Sys.time())}

    Z = gen.Z(n,dz); H = h(Z); G = g(Z);

    for (K.i in 1:length(K.grid)) {
        P = cbind(rep.int(1,n),gen.P(Z,K.grid[K.i])); qr.P = qr(P)
        K = ncol(P);

        if (qr.P$rank==K){

        Q = tcrossprod(P%*%tcrossprod(solve(qr.R(qr.P))),P); M = diag(n) - Q;

        for (m in 1:models){
            dgp = dgp.fn(m,n);

            X = H + sigma.v.fn(Z,dgp$gamma.v)*dgp$v/popval[m,"scale.sigma.v"];
            Y = X*beta + G + dgp$u;

            MX = M%*%X; XMX_inv = XXinv(MX); XMY = crossprod(MX,Y);

            beta.hat = XMX_inv%*%XMY; u2.hat = (M%*%Y - MX%*%beta.hat)^2

            SE0 = sum(u2.hat)/n * XMX_inv
            SE1 = SE0 * n/(n-d-K)

            out[[m]][row,] = c(m,n,d,beta,s,K,ncol(P),beta.hat,SE0,SE1);
        }
        row=row+1;
        }
    }
}

## Save final table
for (m in 1:models) write.csv(out[[m]], file=paste0("output/table_m",m,".csv"));
message("\nSimulation completed. Time: ", Sys.time())


################################################################################
## LOAD DATA & CONSTRUCT DATA & TABLES LATEX
################################################################################
#m=1
for (m in 1:6) {
    data = read.csv(paste0("output/table_m",m,".csv"), row.names=1)

    data1 = cbind(data[,"s"], data[,"K"]/data[,"n"], data[,"beta.hat"]-data[,"beta"])
    colnames(data1) = c("s","K/n","Bias")
    out1 = apply(data1,2,function(x) by(x,data1[,c("K/n")],mean,na.rm=TRUE))

    data2 = cbind(data[,"s"], data[,"K"]/data[,"n"], data[,"beta.hat"])
    colnames(data2) = c("s","K/n","SD")
    out2 = apply(data2,2,function(x) by(x,data2[,c("K/n")],sd,na.rm=TRUE))

    data3 = cbind(data[,"s"], data[,"K"]/data[,"n"], (data[,"beta.hat"]-data[,"beta"])^2,
                  abs(data[,"beta.hat"]-data[,"beta"])/sqrt(data[,"SE0"])<=1.96,
                  abs(data[,"beta.hat"]-data[,"beta"])/sqrt(data[,"SE1"])<=1.96,
                  sqrt(data[,"SE0"]),sqrt(data[,"SE1"]))
    colnames(data3) = c("s","K/n","MSE","CI$_0$","CI$_1$","$\\hat\\sigma$","$s$")
    out3 = apply(data3,2,function(x) by(x,data3[,c("K/n")],mean,na.rm=TRUE))

    results = cbind(out1[,-1],out2[,c(-1,-2)],sqrt(out3[,3]),out1[,3]/out2[,3],out3[,c(-1,-2,-3)])

    table.latex = round(results[,],3);
    colnames(table.latex) = c("$K/n$","$\\mathsf{Bias}$","$\\mathsf{SD}$",
                              "$\\mathsf{RMSE}$","$\\frac{\\mathsf{Bias}}{\\mathsf{SD}}$",
                              "CI$_0$","CI$_1$","$\\hat\\sigma$","$s$")

    latex(table.latex, file=paste0("table_m",m,".txt"),
          append=FALSE, table.env=FALSE, center="none",
          title="", insert.bottom="", rowname=NULL)
}





