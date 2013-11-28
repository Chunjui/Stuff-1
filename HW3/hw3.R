##Question 1
#function of bisection algorithm
#f: target function; l_ini=lower bound of interval;tol: tolerance
#u_ini: upper bound of interval; max_iter: maximum number of iteration
#debugging: print out information
bisection=function(f,l_ini,u_ini,tol,max_iter,debugging=TRUE){
	a=l_ini
	b=u_ini
	converge=0
	root=0
	for (i in 1:max_iter){
		p=a+(b-a)/2
		if(f(p)==0 | (b-a)/2< tol){ #if p is root or the interval 
			output=p                #smaller than tolerance, stop
			converge=1
			break
		}
		
		fa=f(a);fp=f(p)
		if(fa*fp>0){
			a=p
		} else {
			b=p
			root=1
		}
		if(debugging){
			cat(i,"th updated [",a,",",b,"]",'\t',"sign of f(a)*f(b):",sign(f(a)*f(b)),"\n");
		}
	
	}
	

	if(converge==1 && root==1){ #returning output or p depends on converge or not
		print("the algorithm achieved convergence")
		cat('the iteration number is:',i,'\n')
		result=list(root=output,funvalue=f(output),Step=i,converge="TRUE",root="TRUE")
		return(result)
	} else if(converge==1 && root==0 ){
		print("the algorithm achieved convergence, but there is no root in this interval")
		cat('the iteration number is:',i,'\n')
		result=list(root=output,funvalue=f(output),Step=i,converge="TRUE",root="FALSE")
		return(result)
	} else if(converge==0 && root==1 ){
		print("the maximum iteration number is attained, but haven't converge to root")
		result=list(root=p,funvalue=f(p),Step=i,converge="FALSE",root="TRUE")
		return(result)
	} else {
		print("the maximum iteration number is attained, there is no root in this interval")
		result=list(root=p,funvalue=f(p),Step=i,converge="FALSE",root="FALSE")
		return(result)
	}
	

}


#function of Newton-Raphson
#f: target function; d_f: derivative of f
#st: starting point; tol: tolerance
#max_iter: maximum number of iteration
#debugging: print out information
NR=function(f,d_f,st,tol,max_iter,debugging=TRUE){
	p0=st
	converge=0
	for (i in 1:max_iter){
		p=p0-f(p0)/d_f(p0)
		
		if(debugging){
			cat("Updated value:",p,"\t","f(p):",f(p),"\n");
		}

		if(abs(f(p))<tol ){
			output=p
			converge=1
			break
		}
		p0=p

	}
	if(converge==1){
		print("the algorithm achieved convergence")	
		cat('the iteration number is:',i,'\n')
		result=list(root=output,funvalue=f(output),Step=i,converge="TRUE")
		return(result)
	} else {
		print("the maximum iteration number is attained")
		result=list(root=p0,funvalue=f(p0),Step=i,converge="FALSE")
		return(result)
	}
}

l_1=function(x){
	tmp=125/(2+x)-38/(1-x)+34/x
	return(tmp)
}
l_2=function(x){
	tmp=-(125/(2+x)^2)-(38/(1-x)^2)-(34/x^2)
	return(tmp)
}

curve(l_1(x),seq(-2,3,by=0.01),xlab='lambda',ylab='log likelihood',
      main="Log likelihood",lwd=4,lty=2,col='blue',
      xlim=c(-2,3)) # theta in [0,1], thus lambda in [0,1] too.
#savePlot('loglikelihood',type='png')
x1=bisection(l_1,0,1,0.000001,100) 
x2=NR(l_1,l_2,0.1,0.000001,100)


#para=function(x){ return((x-1)^2-4)}
#dpara=function(x){return(2*(x-1))}

##Question 2
sa=function(x){return(1/(1+x))}
aa=function(x){return(x/(1+x))}
plot(1:5,xlab='y^bar',ylab='converge rates',main='Converge rate for both SA and AA',
     type='n',xlim=c(0,5),ylim=c(0,1))
curve(sa(x),seq(0,5,by=0.01),col='red',lwd=4,add=TRUE)
curve(aa(x),seq(0,5,by=0.01),col='blue',lwd=4,lty=2,add=TRUE)
