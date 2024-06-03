#################################
# 3. Když výběr existuje.
#################################

make.X <- function(){
  Xpocet <- matrix(
    c(1,1,1, 6, 1,1,2, 6, 1,1,3,26, 1,1,4,10, 1,2,1, 2,
      1,2,2,25, 1,2,3,17, 1,2,4, 9, 1,3,2,12, 1,3,3,16,
      1,3,4, 3, 1,4,1, 2, 1,4,2, 1, 1,4,3,12, 1,4,4, 9,
      2,1,1,16, 2,1,2,46, 2,1,3,67, 2,1,4,11, 2,2,1, 4,
      2,2,2,44, 2,2,3,50, 2,2,4,26, 2,3,2,20, 2,3,3,25,
      2,3,4, 9, 2,4,1, 2, 2,4,2,12, 2,4,3, 7, 2,4,4, 5),
    ncol=4, byrow=T)
  Xdf <- as.data.frame(Xpocet[rep(1:nrow(Xpocet),Xpocet[,4]),1:3])
  Xdffac <- as.data.frame(lapply(Xdf,as.factor))
  names(Xdffac)<-c('sex','agecat','educat')
  X <- with(Xdffac,
            {cbind(
              model.matrix(~sex-1),
              model.matrix(~agecat-1),
              model.matrix(~educat-1),
              total = 1)
  })
  }
  X <- make.X()
  N <- dim(X)[1]
  p <- dim(X)[2]
  s <- c(44,56,8,38,37,17,2,26,43,29,100)
  
  
  #################################
  
  library('Rsymphony')
  result <- Rsymphony_solve_LP(obj=rep(1,N),mat=t(X),
                               dir=rep('==',p),rhs=s, types=rep('B',N), max=F)
  w <- result$solution[1:N]
  tabvysledek <- cbind(t(X)%*%w, s)
  colnames(tabvysledek) <- c('vyber','s')
  tabvysledek
  
  #################################
  # 4. Když výběr neexistuje.
  #################################
  
  
  s <- c(111,139,21,95,92,42,5,66,107,72,250)
  obj <- c(rep(1,N),rep(1,p))
  A <- matrix(0, nrow=2*p, ncol=p+N)
  A[1:p,1:N] <- t(X)
  A[1:p,(N+1):(N+p)] <- -diag(rep(1,p))
  A[(p+1):(2*p),1:N] <- -t(X)
  A[(p+1):(2*p),(N+1):(N+p)] <- -diag(rep(1,p))
  op <- c(rep('<=',2*p))
  rhs <- c(s,-s)
  bounds <- list(
    lower = list(ind=c((N+1):(N+p)), val=c(rep(0,p))),
    upper = list(ind=c((N+1):(N+p)), val=c(s[1:p-1],0))
  )
  result <- Rsymphony_solve_LP(obj, A, op, rhs,
                               types = c(rep('B',N),rep('I',p)),
                               max = F, bounds=bounds)
  w <- result$solution[1:N]
  tabvysledek <- cbind(t(X)%*%w, s)
  colnames(tabvysledek) <- c('vyber','s')
  tabvysledek   
    

  #################################
  # 4. Kvadratická minimalizace.
  #################################
  
  library('Rcplex')
  s <- c(111,139,21,95,92,42,5,66,107,72,250)
  A <- matrix(0,nrow=2*p,ncol=p+N)
  A[1:p,1:N] <- t(X)
  A[1:p,(N+1):(N+p)] <- -diag(rep(1,p))
  A[(p+1):(2*p),1:N] <- -t(X)
  A[(p+1):(2*p),(N+1):(N+p)] <- -diag(rep(1,p))
  rhs <- c(s,-s)
  result <- Rcplex(cvec = c(rep(0,N+p)), Amat=A, bvec=rhs,
                   Qmat = diag(c(rep(1,N),rep(1,p))),lb = 0,
                   ub = c(rep(1,N),s[1:(p-1)],0),
                   objsense = c('min'), sense = 'L',
                   vtype = c(rep('B',N),rep('I',p)))
  w <- result$xopt[1:N]
  tabvysledek <- cbind(t(X)%*%w,s)
  colnames(tabvysledek) <- c('vyber','struktura')
  tabvysledek
  
  
  ##############################
  
  library('gurobi')
  model<-list(Q = diag(c(rep(1,N), rep(1,p))),
              obj = 0, A = A, rhs = rhs,
              sense = c(rep('<=',2*p)),
              lb=c(rep(0,N+p)), ub=c(rep(1,N),s[1:(p-1)],0),
              vtype= c(rep("B",N),rep("I",p)) )
  params<-list(OutputFlag=0)
  result<-gurobi(model,params)
  w<-result$x[1:N]
  tabvysledek <- cbind(t(X)%*%w,s)
  colnames(tabvysledek) <- c('vyber','struktura')
  tabvysledek
  
  #################################
  # 5. Největší možný výběr.
  #################################
  
  zaklad <- s[11]
  r <- s/zaklad
  obj <- c(rep(1,N)) 
  A <- matrix(0, nrow=2*p, ncol=N)
  J <- matrix(1, nrow=p, ncol=N)
  dr <- diag(r) %*% J
  A[1:p,1:N] <- t(X) - dr
  A[(p+1):(2*p),1:N] <- t(X) - dr
  op <- c(rep('>=',p),rep('<=',p))
  b <- c(rep(0.9,p))
  rhs <- c(rep(0,p)-b,rep(0,p)+b)
  result <- Rsymphony_solve_LP(obj, A, op, rhs,
                               types = c(rep('B',N)),
                               max = T)
  w <- result$solution
  tabvysledek <- as.data.frame(cbind(c(t(X)%*%w), c(round(s/zaklad*sum(w),1))))
  colnames(tabvysledek) <- c('vyber','zadani')
  tabvysledek
  

  ####################################################
  # 6. Výběr s použitím knihovny ROI a jejích doplňků.
  ####################################################
  
  library(ROI)
  library(ROI.plugin.neos)
  
  s <- c(111,139,21,95,92,42,5,66,107,72,250)
  
  obj <- c(rep(1,N),rep(1,p))
  A <- matrix(0, nrow=2*p, ncol=p+N)
  A[1:p,1:N] <- t(X)
  A[1:p,(N+1):(N+p)] <- -diag(rep(1,p))
  A[(p+1):(2*p),1:N] <- -t(X)
  A[(p+1):(2*p),(N+1):(N+p)] <- -diag(rep(1,p))
  op <- c(rep('<=',2*p))
  rhs <- c(s,-s)
  lp_bound <- V_bound(li=c((N+1):(N+p)),ui=c((N+1):(N+p)), lb=c(rep(0,p)), ub=round(c(s[1:p-1],0),0)) 
  types <- c(rep('B',N),rep('I',p))
  
  lp <- OP(objective = obj,
           L_constraint( L = A,
                         dir = op,
                         rhs = round(rhs,0)),
           bounds=lp_bound,
           types=types,
           maximum = FALSE)
  sol <- ROI_solve(lp, solver = "neos", method="cplex",email="emailova adresa")
  
  w <- sol$solution[1:N]
  tabvysledek <- cbind(t(X)%*%w,s)
  colnames(tabvysledek) <- c('vyber','struktura')
  tabvysledek
  
  
  
  
  
  