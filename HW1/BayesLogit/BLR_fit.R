
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
args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)

####
# sim_start ==> Lowest simulation number to be analyzed by this particular batch job
###

#######################
sim_start <- 1000
length.datasets <- 200
#######################

if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}

# Simulation datasets numbered 1001-1200

########################################################################################
########################################################################################
bayes.logreg = function (m,y,X,beta.0,Sigma.0.inv,niter,burnin,print.every,retune,ini,upd){
#cat('this is a update for variables: ',upd,"\n")
re=floor(burnin/retune);
p=length(upd)
x.tmp=X[,upd];
p=length(upd)
sigma.beta=solve(Sigma.0.inv)
d=dim(X)
d2=d[2]
tmp.data=X[,2:d2]
f=glm(cbind(y,(m-y))~tmp.data,family="binomial")
tmp.f=vcov(f)
var.ini=tmp.f[upd,upd]#solve(t(x.tmp)%*%x.tmp)*min(diag(sigma.beta))
#diag(var.ini)=rep(1,p)
v.square=c(seq(16,0.8,-0.8),-seq(0,16,0.8)); #variance for proposal dist

p.beta=function(b){
	inn=X%*%b
	if(max(inn)<11){
	tmp=sum(inn*y)-sum(m*log(1+exp(inn)))-1/2*t(b-beta.0)%*%Sigma.0.inv%*%(b-beta.0);
	#cat(b,"\n")
	} else{
	tmp=sum(inn*y)-sum(m*inn)-1/2*t(b-beta.0)%*%Sigma.0.inv%*%(b-beta.0);
	}
	#cat(sum(X%*%b*y),-sum(m*log(1+exp(X%*%b))),-1/2*t(b-beta.0)%*%Sigma.0.inv%*%(b-beta.0),"\n")
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
for (i in 1:re){
	u.threshold=runif(retune,0,1);
	if(up==1){
	var.proposal=var.ini*exp(v.square[(up.ind)])
	} else if(up==0) {
	up.ind=up.ind-1
	var.proposal=var.ini*exp(v.square[(up.ind)])
	} else if(up==2) {
	up.ind=up.ind+1
	var.proposal=var.ini*exp(v.square[(up.ind)])
	}

	rate.accept=0;
	for (j in 1:retune){
		beta.current=matrix(beta.update[,(t.stage+1)])
		beta.proposal=(mvrnorm(1,beta.current,var.proposal));
		beta.curr.cond[upd]=beta.current
		beta.prop.cond[upd]=beta.proposal
		tmp=min((p.beta(beta.prop.cond)-p.beta(beta.curr.cond)),0);
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
	rate.accept=rate.accept/retune;
	cat("Acceptance Rate:",rate.accept,"\n");
	if(rate.accept>0.3 && rate.accept<0.6) {
	#cat("Tuning for proposal distributions is successful, v.square: ",v.square[i],"\n");
	up=1
	} else if(rate.accept>0.6) {
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
	tmp=min((p.beta(beta.prop.cond)-p.beta(beta.curr.cond)),0);
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
beta.0 <- matrix(c(0,0));
Sigma.0.inv <- diag(rep(1.0,2));
niter <- 10000;
burnin=1800;
print.every=1000;
retune=200;
verbose=TRUE;

upd=c(1,2)
# etc... (more needed here)
#################################################
# Read data corresponding to appropriate sim_num:
file.name=paste('./data/blr_data_',sim_num,'.csv',sep='');
data.tmp=as.matrix(read.table(file.name,header=TRUE,sep=','));
# Extract X and y:
y=data.tmp[,1];
X=data.tmp[,c(3,4)];
m=data.tmp[,2];
x.tmp=X[,2]
ini.tmp=glm(y~x.tmp)
ini=c(0,0)#matrix(ini.tmp$coeff)
# Fit the Bayesian model:
# beta.post=bayes.logreg(m,y,X,beta.0,Sigma.0.inv,niter,burnin,
                           # print.every,retune,verbose)
						   

beta.post=bayes.logreg(m,y,X,beta.0,Sigma.0.inv,niter,burnin,print.every,retune,ini,upd)
# Extract posterior quantiles...
l1=burnin+1
tmp=dim(beta.post)
l2=tmp[1]
b0.percent=matrix(quantile(beta.post[l1:l2,1],probs=seq(0.01,0.99,0.01)))
b1.percent=matrix(quantile(beta.post[l1:l2,2],probs=seq(0.01,0.99,0.01)))
# Write results to a (99 x p) csv file...
file.name=paste('./results/blr_res_',sim_num,".csv",sep='')
write.table(cbind(b0.percent,b1.percent),file=file.name,sep=',',
row.names = FALSE,col.names = FALSE)
# Go celebrate.
cat("done. :)\n")







