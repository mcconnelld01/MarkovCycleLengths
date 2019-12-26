require(MCMCpack)
require(rlist)
require(expm)
require(pracma)
require(tidyverse)
require(matrixcalc)
require(complexplus)
require(expm)
require(pracma)
require(tidyverse)


## Function to generate a random 3*3 stochasti matrix
## By default, values with be drawn from dirichlet distributions 

rand33<-function(){
  r1<-rdirichlet(1,c(3,1,1))
  r2<-rdirichlet(1,c(1,3,1))
  r3<-rdirichlet(1,c(1,1,3))
  return(matrix(c(r1,r2,r3),nrow=3,byrow=T))
}



# Take the p-th root of the matrix A
# Only works for non-repeated eigenvalues, need to generalise this
# If root taking fails, returns the 1*1 zero matrix

matroot<-function(A,p)
{
  
  if(is.integer(1/p))
  {
    return(A%^%(1/p))
  }
  else if (round(Det(eigen(A)$vectors),15)!=0)
  {
    return(eigen(A)$vectors %*% diag(eigen(A)$values^{1/p}) %*% solve(eigen(A)$vectors))
  }
  else {return(matrix(0,1,1))}
}


# Check if a matrix is stochastix (up to 15 decimal places)
# Probably should allow this to take complex entries and return false,
# or maybe true if they are sufficiently close to zero.


is.stochastic<-function(A=diag(1,nrow=2,ncol=2))
{
  if (any(is.na(A)))
  {
    return(FALSE)
  }
  
  # maybe if any imaginary parts rounded to 15 places are nonzero return false here?
  
  rowsums<-round(rowSums(A),15)
  rowmins<-numeric()
  if(any(Im(A)!=0))
  {return(FALSE)}
  else{
    
    if (all(Re(rowsums)==1) & min(Re(A)) >= 0 )
    {return(TRUE)}
    else {return(FALSE)}
  }}



# Generate a list of stochastic matrices with non-stochastic 12-th roots
# Could probably change this to a while loop which runs until the list has been filled.


makelist<-function(n=100,r=12)
{
  j=1
  matlist<-list()
  for (i in 1:n)
  {
  M<-rand33()
  
   if (!is.stochastic(matroot(M,r)) & !any(is.na(matroot(M,r))))
      {
      matlist[[j]]<-M
      j<-j+1
      }
  
  
  }
  return(matlist)
}

## Do the wrong conversion, i.e. the traditional 'constant-rate' formula

pconvert<-function(A,p)
{
  n=nrow(A)
  m=ncol(A)
  B=matrix(rep(0,n*m),nrow=n,ncol=m)
  for (i  in 1:n){
    for (j in 1:m){
      if (i != j)
      {
        B[i,j] <- 1-(1-A[i,j])^{p}
      }
    }
    B[i,i]<-1-sum(B[i,-i])
  }
  return(B)
}




# The function stochroot takes a (stochastic) matrix A and returns either the p-th root
# of that matrix, if it is stochastic, or the closest stochastic matrix to its p-th root
# in the Frobenius norm.


# Can probably insert my own matroot function here

stochroot<-function(A,p)
{
  # Calculate the p-th root of A
  
  if (is.numeric(eigen(A)$values))
  {
    Ap = expm(logm(A)*(1/p))
  }
  else
  {
    Ap = eigen(A)$vectors %*% diag(sapply(eigen(A)$values, function (x) x^p)) %*% solve(eigen(A)$vectors)
  }
  n=nrow(A)
  b=matrix(0,nrow=n,ncol=n)
  for (j in 1:n)
  {
    a=Re(Ap[j,])
    cond<-TRUE
    
    while (cond==TRUE)
      
    { if (round(sum(a),15)==1 & min(a)>=0)
    { b[j,]<-a
    
    cond<-FALSE
    }
      
      lambda = (sum(a)-1)/n
      x=a-lambda*rep(1,n)
      
      if (min(x)>=0)
      {    b[j,]<-x 
      
      cond<-FALSE }
      
      for (k in 1:n)
      {
        x[k]<-max(0,x[k])
      }
      a<-x
      
      
    }
  }
  return(b)
  
}


## Return matrix, bad approximation 1, bad approximation 2

# Rename this, tidy up output

# Integrate it with Markov Comparison App?

wronglist1<-function(n=100,p=12)
{
  inputlist<-makelist(n,p)
 outlist<- lapply(inputlist,function(x){
        return(list("Matrix"=x,
                    "Att1"=x-pconvert(x,1/p)%^%p,
                    "Att2"=x-stochroot(x,p)%^%p,
                    "NormAtt1"=norm(x-pconvert(x,1/p)%^%p,type="F"),
                    "NormAtt2"=norm(x-stochroot(x,p)%^%p,type="F"),
                    "RealEigenvalues"=is.numeric(eigen(x)$values)
        ))
    }
    )
 return(list.sort(outlist,NormAtt2-NormAtt1))
}

## Perhaps change the rand33 function to take as input a matrix of counts,
# then use random draws from the corresponding dirichlet distribution.

## Maybe we would also like to obtain an estimate of how often these errors arise
