
##
#
# Logistic regression
# 
# Y_{i} | \beta \sim \textrm{Bin}\left(n_{i},e^{x_{i}^{T}\beta}/(1+e^{x_{i}^{T}\beta})\right)
# \beta \sim N\left(\beta_{0},\Sigma_{0}\right)
#
##

library(MASS)
########################################################################################
########################################################################################
## Handle batch job arguments:

# 1-indexed version is used now.

####
# sim_start ==> Lowest simulation number to be analyzed by this particular batch job
###

#######################
#######################

# Simulation datasets numbered 1001-1200

########################################################################################
########################################################################################
bayes.logreg = function (m,y,X,beta.0,Sigma.0.inv,niter,burnin,print.every,retune,ini,upd){
cat('this is a update for variables: ',upd,"\n")
re=floor(burnin/retune);
p=length(upd)
x.tmp=X[,upd];
p=length(upd)
sigma.beta=solve(Sigma.0.inv)
d=dim(X)
d2=d[2]
tmp.data=X[,2:d2]
f=glm(y~tmp.data,family=binomial)
tmp.f=vcov(f)
var.ini=tmp.f[upd,upd]#solve(t(x.tmp)%*%x.tmp)*min(diag(sigma.beta))
#print(var.ini)
#cat("\n")
#diag(var.ini)=rep(1,p)
v.square=c(seq(16,0.8,-0.8),-seq(0,16,0.8)); #variance for proposal dist

p.beta=function(b){
	inn=X%*%b
	tmp=sum(inn*y)-sum(m*log(1+exp(inn)))-0.5*t(b-beta.0)%*%Sigma.0.inv%*%(b-beta.0);
	return(tmp);	
}

#burnin stage

beta.update=matrix(ini[upd],p,1)
beta.curr.cond=ini
beta.prop.cond=ini
t.stage=0;
up=1;
up.ind=21;
var.proposal=var.ini
var.curr=var.ini
s=dim(var.ini)
s1=s[1]
s2=s[2]
I=matrix(0,s1,s2)
rate.accept=0
for (i in 1:re){
	u.threshold=runif(retune,0,1);
	if(up==1){
	#var.proposal=rate.accept*var.curr+(1-rate.accept)*var.proposal
	#print(var.proposal)
	#cat("\n")
	var.proposal=var.proposal#var.ini*exp(v.square[(up.ind)])
	} else if(up==0) {
	#var.proposal=rate.accept*var.curr+(1-rate.accept)*var.proposal
	#var.proposal=var.proposal*exp(0.1)
	up.ind=up.ind-1
	#ad=I
	#diag(ad)=diag(var.proposal)*(exp(0.5)-1)	
	var.proposal=var.proposal*exp(0.85)#var.ini*exp(v.square[(up.ind)])
	} else if(up==2) {
	#var.proposal=rate.accept*var.curr+(1-rate.accept)*var.proposal
	#var.proposal=var.proposal*exp(-0.1)
	up.ind=up.ind+1
	#ad=I
        #diag(ad)=diag(var.proposal)*(exp(-0.8)-1)
	var.proposal=var.proposal*exp(-0.85)#var.ini*exp(v.square[(up.ind)])
	}

	rate.accept=0;
	for (j in 1:retune){
		beta.current=matrix(beta.update[,(t.stage+1)])
		beta.proposal=(mvrnorm(1,beta.current,var.proposal));
		#tmp2=which.max(apply(tmp1,1,function(x) p.beta(matrix(x))))
		#beta.proposal=tmp1[tmp2,];
		beta.curr.cond[upd]=beta.current
		beta.prop.cond[upd]=beta.proposal
		tmp=min(0,(p.beta(beta.prop.cond)-p.beta(beta.curr.cond)));
		#cat(p.beta(beta.prop.cond),p.beta(beta.curr.cond),"\n")
		#cat(tmp,log((u.threshold[j])),"\n")
		#criterion=min(1,tmp);
		if(tmp > log(u.threshold[j])){
		
		t.stage=t.stage+1;
		beta.update=cbind(beta.update,beta.proposal);
		rate.accept=rate.accept+1;
		#cat(beta.proposal,rate.accept,"\n")
		} else {
		beta.update=cbind(beta.update,beta.current);
		t.stage=t.stage+1;
		}
	}
	#tmp1=t(beta.update)
	#tmp2=tail(tmp1,n=retune)
	#var.curr=var(tmp2)
	#print(var.curr)
	#cat("\n")
	rate.accept=rate.accept/retune;
	#cat(rate.accept)
	cat("Acceptance Rate:",rate.accept,"\n");
	if(rate.accept>0.35 && rate.accept<0.55) {
	#cat("Tuning for proposal distributions is successful, v.square: ",v.square[i],"\n");
	up=1
	} else if(rate.accept>0.55) {
	up=0
	} else{
	up=2
	}
	#cat(up,"\n")
}


#MCMC
step.upt=0;
rate.accept=0;
repeat{
	step.upt=step.upt+1;
	beta.current=beta.update[,(t.stage+1)];
	beta.proposal=(mvrnorm(1,beta.current,var.proposal));
	beta.curr.cond[upd]=beta.current
	beta.prop.cond[upd]=beta.proposal
	tmp=p.beta(beta.prop.cond)-p.beta(beta.curr.cond);
	u.threshold=runif(1,0,1)
	#criterion=min(1,tmp);
	if(tmp>log(u.threshold)){
	t.stage=t.stage+1;
	beta.update=cbind(beta.update,beta.proposal);
	rate.accept=rate.accept+1;
	} else {
	beta.update=cbind(beta.update,beta.current);
	t.stage=t.stage+1;
	}
	if((step.upt%%print.every)==0){
	cat("The acceptance rate is: ",(rate.accept/step.upt),"\n")
	}
	
	if(t.stage==(niter+burnin)){
	break;
	}
}
	return(t(beta.update))
}

#################################################
# Set up the specifications:
beta.0 <- matrix(rep(0,11));
Sigma.0.inv <- diag(rep(0.001,11));
niter <- 20000;
burnin=2000;
print.every=1000;
retune=100;
verbose=TRUE;

upd=1:11
# etc... (more needed here)
#################################################
# Read data corresponding to appropriate sim_num:
cancer=as.matrix(read.table(file="breast_cancer.txt",header=TRUE));
# Extract X and y:
size=dim(cancer)
tmp=matrix(as.double(cancer[,1:10]),ncol=10)
X=cbind(rep(1.0,size[1]),tmp)
y=matrix(0,size[1],1)
y[which(cancer[,11]=='M')]=1
m=matrix(rep(1,size[1]))
# Fit the Bayesian model:
d=dim(X)
d2=d[2]
tmp.data=X[,2:d2]
f=glm(y~tmp.data,family=binomial)
tmp=vcov(f)
ini=matrix(f$coefficients)+0.005
#matrix(c(0.61,15,0.2,2.29,0.33,-0.3,0.14,-9.73,1.07,0.469,1.675))#matrix(f$coefficients)+0.005
# Fit the Bayesian model:
# beta.post=bayes.logreg(m,y,X,beta.0,Sigma.0.inv,niter,burnin,
                           # print.every,retune,verbose)
						   

beta.post=bayes.logreg(m,y,X,beta.0,Sigma.0.inv,niter,burnin,print.every,retune,ini,upd)

#standardized
#col=colMeans(X)
#col=matrix(col[2:11],1,10)
#colvar=apply(X,2,function(x) var(x))
#colvar=matrix(colvar[2:11],1,10)
#s.X=apply(tmp.data,1,function(x) (x-col)/colvar)
#s.X=cbind(rep(1.0,size[1]),t(s.X))
#tmp.data=s.X[,2:d2]
#f.s=glm(y~tmp.data,family=binomial)
#ini.s=matrix(f$coefficients)

#beta.post.s=bayes.logreg(m,y,s.X,beta.0,Sigma.0.inv,niter,burnin,print.every,retune,ini.s,upd)
# Extract posterior quantiles...
# Write results to a (99 x p) csv file...
# Go celebrate.
cat("done. :)\n")
save.image(file="mh_breast_test11.RData")






