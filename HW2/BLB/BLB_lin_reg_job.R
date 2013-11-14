
mini <- FALSE

#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 1000
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  sim_num <- sim_start + as.numeric(args[1])
  sim_seed <- (762*(sim_num-1) + 121231)
}

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))

# Find r and s indices:

#============================== Run the simulation study ==============================#

# Load packages:
library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)
library(MASS)
# I/O specifications:
datapath <- "/home/pdbaines/data"
outpath <- "output/"

# mini or full?
if (mini){
	rootfilename <- "blb_lin_reg_mini"
} else {
	rootfilename <- "blb_lin_reg_data"
}

# Filenames:

# Set up I/O stuff:
tmp=sprintf("%s.desc",rootfilename)
#desc=dget(tmp)
# Attach big.matrix :
foo=attach.big.matrix(tmp,path=datapath)
# Remaining BLB specs:
n=1000000
d=1000
r=50
s=5
gam=0.7
b=ceiling(n^gam)
# Extract the subset:
s_index=as.integer((sim_num-1)/r)-(sim_start/r)+1

tmp1=(sim_num %% r)
if(tmp1!=0){
	r_index=tmp1
} else {
	r_index=50
}
set.seed(s_index*300)
b.index=sample(1:n,b , replace = FALSE)
b.x=foo[b.index,1:d]
b.y=foo[b.index,(d+1)]
# Reset simulation seed:
set.seed(sim_seed)
# Bootstrap dataset:
bs.sample=rmultinom(1,n,prob=rep(1/b,b))
# Fit lm:
fit=lm(b.y~b.x-1,weights=bs.sample)
beta.k=fit$coeff
# Output file:
#file.name=sprintf('./%s/blb_lin_reg_%04d',outpath,sim_num)
#outfile = paste("./output/","coef_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt",sep='')
outfile = sprintf("./output/coef_%02d_%02d.txt",s_index,r_index)
# Save estimates to file:
write.table(beta.k,file=outfile,row.names=FALSE)

