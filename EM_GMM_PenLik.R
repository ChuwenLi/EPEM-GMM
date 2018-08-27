###### NOTE: Functions are based on a very specific set up of the initial data matrix
###### The main function is "EMfunc" which calls all other functions defined here

##### Function to get sigmas and betas for each gene

### Not averaging over replicates (TR observations per gene)
init1<-function(Y,p){
  m         <- lm(Y~poly(timepoints,p,raw=T)) 
  beta0T    <- m$coef
  sig20T    <- summary(m)$sigma^2
  p.one     <- sum(summary(m)$coef[2:(p+1),4] <= 0.05)
  
  return(c(beta0T,s2=sig20T,p.one))
}

### Averaging over replicates (T observations per gene)
init2<-function(Y,p){
  time2     <- rep(timepoints,length(Y)/length(timepoints))
  m         <- lm(Y~poly(time2,p,raw=T)) 
  beta0T    <- m$coef
  sig20T    <- summary(m)$sigma^2
  p.one     <- sum(summary(m)$coef[2:(p+1),4] <= 0.05)
  
  return(c(beta0T,s2=sig20T,p.one))
}

### FUNCTIONS FOR PPFUN
covfcn.FE    <- function(s2,X){s2*diag(dim(X)[1])}
covfcn.ME1NR <- function(g,s2,X,U){U%*%g%*%t(U)+s2*diag(dim(X)[1])}
covfcn.ME    <- function(g,c2,s2,X,U,V){U%*%g%*%t(U)+c2*V%*%t(V)+s2*diag(dim(X)[1])}

### Function to obtain posterior probabilities (wik)
ppfun<-function(X,U,V,dat.sh,beta,S2,G,C2,alpha,model){
  
  X.L<-lapply(1:length(alpha), function(j) X)
  U.L<-lapply(1:length(alpha), function(j) U)
  V.L<-lapply(1:length(alpha), function(j) V)
  
  if(model == 'NOM'){
    print(model)
    mu = beta
    cov.y=lapply(S2,function(s){s*diag(dim(dat.sh)[2])})
  }
  else{
    mu        <- apply(beta,2,function(b){X%*%b})

    if(model == 'FE'){
      print(model)
      cov.y     <- mapply(covfcn.FE,s2=as.list(S2),X=X.L,SIMPLIFY=FALSE)
    }
    else if(model%in% c('ME1NR','ME2NR')){
      print("No Rep")
      cov.y     <- mapply(covfcn.ME1NR,g=G,s2=as.list(S2),X=X.L,U=U.L,SIMPLIFY=FALSE)
    }
    else{
      print("Rep")
      cov.y     <- mapply(covfcn.ME,g=G,c2=as.list(C2),s2=as.list(S2),X=X.L,U=U.L,V=V.L,SIMPLIFY=FALSE)
    }
  }
  
  wk1.T      <- matrix(0,ncol=length(alpha),nrow=dim(dat.sh)[1])
  for(k in 1:length(alpha)){
    wk1.T[,k]<- alpha[k]*dmvn(X=as.matrix(dat.sh), mu=mu[,k], sigma=cov.y[[k]], log=FALSE) 
  }
  
  wik<-wk1.T/apply(wk1.T,1,sum)
  
  whichfcn<-function(w){which(w==max(w))}
  cluster<-apply(wik,1,whichfcn)
  
  l<-list("post.prob" = wik,
          "cluster"   = cluster)
  
  return(l)
}

### Function to obtain Mixing probabilities of each component of the mixture model
alphafun<-function(pp,lambda,a.last){
  (colSums(pp)/dim(pp)[1])+lambda*a.last*(log(a.last)-sum(a.last*log(a.last)))
}

### Function to obtain lambda
lambda.fun<-function(eta, alpha.cur, alpha.last, pp){
  l1T<-exp(-1*eta*dim(pp)[1]*abs(alpha.cur-alpha.last))
  l2T<-colSums(pp)/dim(pp)[1]
  l1 <-sum(l1T)/length(alpha.cur)
  l2 <-(1-max(l2T))/(-1*max(alpha.last)*sum(alpha.last*log(alpha.last)))
  return(min(l1,l2))
}

### Function to obtain Maximum likelihood estimates
mlefun<-function(X,U,V,pp,dat.sh,beta,S2,G,C2,alpha,model,R,ntimepts,lambda){
  
  ppk <- colSums(pp)
  X.L<-lapply(1:length(alpha), function(j) X)
  U.L<-lapply(1:length(alpha), function(j) U)
  V.L<-lapply(1:length(alpha), function(j) V)
  
  
  if(model == 'NOM'){
    #print('Estimated Parameters, mu and sigma- NO MODEL')

    covy=lapply(S2,function(s){s*diag(dim(dat.sh)[2])})
    
    betaT  <- matrix(NA,ncol=dim(pp)[2],nrow=dim(dat.sh)[2])  ##Not really beta- this is a mu
    sig2T  <- matrix(NA,ncol=dim(pp)[2],nrow=1)
    loglik <- matrix(NA,ncol=length(alpha),nrow=dim(dat.sh)[1])
    
    for(k in 1:dim(pp)[2]){
      betaT[,k]   <- (pp[,k]%*%as.matrix(dat.sh))/ppk[k]
      Ym          <- t(apply(dat.sh,1,function(Y){Y-betaT[,k]}))
      sig2T[k]    <- (pp[,k]%*% apply(Ym,1,function(y){matrix(y,nrow=1)%*%matrix(y,ncol=1)}))/ (dim(dat.sh)[2]*ppk[k])
      loglik[,k]  <- alpha[k]*dmvn(X=as.matrix(dat.sh), mu=betaT[,k], sigma=covy[[k]], log=FALSE)
      
    }
    wik<-loglik/apply(loglik,1,sum)
    
    whichfcn<-function(w){which(w==max(w))}
    cluster<-apply(wik,1,whichfcn)
    
    ll<-sum(log(apply(loglik,1,sum)))+ lambda*dim(pp)[1]*sum(alpha*log(alpha))
    
    return(list(beta=betaT,
                sig2=sig2T,
                loglik=ll,
                post.prob = wik,
                cluster   = cluster))
  }
  else{
    
    mu        <- apply(beta,2,function(b){X%*%b})
    
  if(model=='FE'){
    dat.sh.L <-lapply(seq_len(nrow(dat.sh)), function(i) matrix(as.numeric(dat.sh[i,]),ncol=1))
    XY <-lapply(seq_len(nrow(dat.sh)), function(i) t(X)%*%matrix(as.numeric(dat.sh[i,]),ncol=1))
    XX <-lapply(1:dim(pp)[1], function(j) t(X)%*%X)

    betaT  <- matrix(NA,ncol=dim(pp)[2],nrow=dim(beta)[1])
    sig2T  <- matrix(NA,ncol=dim(pp)[2],nrow=1)
    loglik <- matrix(NA,ncol=length(alpha),nrow=dim(dat.sh)[1])

    for(k in 1:dim(pp)[2]){
      num  <-Reduce('+',mapply(function(a,b){a*b},a=XX,b=as.list(pp[,k]),SIMPLIFY=FALSE))
      denom<-Reduce('+',mapply(function(a,b){a*b},a=XY,b=as.list(pp[,k]),SIMPLIFY=FALSE))

      betaT[,k]   <- solve(num)%*%denom
      Ym          <- mapply(function(Y){crossprod(Y-X%*%betaT[,k])},Y=dat.sh.L,SIMPLIFY=FALSE)
      sig2T[k]    <- Reduce('+',mapply(function(a,b){a*b},a=Ym,b=as.list(pp[,k]),SIMPLIFY=FALSE))/ (ppk[k]*ntimepts)
      
      covy        <- sig2T[k]*diag(ntimepts)
      loglik[,k]  <- alpha[k]*dmvn(X=as.matrix(dat.sh), mu=X%*%betaT[,k], sigma=covy, log=FALSE)
      
    }
    wik<-loglik/apply(loglik,1,sum)
    
    whichfcn<-function(w){which(w==max(w))}
    cluster<-apply(wik,1,whichfcn)
    
    ll<-sum(log(apply(loglik,1,sum)))+ lambda*dim(pp)[1]*sum(alpha*log(alpha))
    
    return(list(beta=betaT,
                sig2=sig2T,
                loglik=ll,
                post.prob = wik,
                cluster   = cluster))
    
  }
  else if(model %in% c('ME1NR','ME2NR')){
    #print('MLEs for models without random-effect for replicates')
    
    covfcn     <- function(g,c2,s2){U%*%g%*%t(U)+s2*diag(dim(X)[1])}
    cov.y      <- mapply(covfcn,g=G,s2=as.list(S2),SIMPLIFY=FALSE) 
    inv.cov.y  <- lapply(cov.y,function(c){solve(c)})
    covDfcn    <- function(g,inv.err,s2){U%*%(g-g%*%t(U)%*%inv.err%*%U%*%g)%*%t(U)}
    cov.d      <- mapply(covDfcn,g=G,inv.err=inv.cov.y,s2=as.list(S2),SIMPLIFY=FALSE)
    
    B1 <-lapply(inv.cov.y,function(c) {solve(t(X)%*%X)%*%t(X)%*%c})
    S1 <-lapply(inv.cov.y,function(c){c%*%c})
    GF1<-matrix(mapply(function(g,c){g%*%t(U)%*%c},g=G,c=inv.cov.y,SIMPLIFY=F))
    GF2<-mapply(function(g,c){g - g%*%t(U)%*%c%*%U%*%g},g=G,c=inv.cov.y,SIMPLIFY=F)
    
    b.fcn1     <- function(m){beta[,k] + S2[k]*B1[[k]]%*%m}
    G.fcn1     <- function(m){quad.tform(m%*%t(m),GF1[[k]]) + GF2[[k]]}
    s2.fcn1    <- function(m){S2[k]^2*quad.form(S1[[k]],m) + sum(diag(cov.d[[k]]))}

    betaT  <- matrix(NA,ncol=dim(pp)[2],nrow=dim(beta)[1])
    sig2T  <- matrix(NA,ncol=dim(pp)[2],nrow=1)
    G.lstT <- list()
    loglik <- matrix(NA,ncol=length(alpha),nrow=dim(dat.sh)[1])
    
    for(k in 1:dim(pp)[2]){
      Ym          <- t(apply(dat.sh,1,function(Y){Y-mu[,k]}))
      betaT[,k]   <- (pp[,k]%*%t(apply(Ym,1,b.fcn1))) / ppk[k] 
      G.lstT[[k]] <- matrix((pp[,k]%*%matrix(t(apply(Ym,1,G.fcn1)),ncol=dim(G[[k]])[1]^2))/ ppk[k],ncol=(dim(G[[k]])[1])) 
      sig2T[k]    <- (pp[,k]%*%matrix(t(apply(Ym,1,s2.fcn1)),ncol=1))/ (ppk[k]*R*ntimepts) 

      covy        <- U%*%G.lstT[[k]]%*%t(U)+sig2T[k]*diag(ntimepts*R)
      loglik[,k]  <- alpha[k]*dmvn(X=as.matrix(dat.sh), mu=X%*%betaT[,k], sigma=covy, log=FALSE)
    }
    wik<-loglik/apply(loglik,1,sum)
    
    whichfcn<-function(w){which(w==max(w))}
    cluster<-apply(wik,1,whichfcn)
    ll<-sum(log(apply(loglik,1,sum))) + lambda*dim(pp)[1]*sum(alpha*log(alpha))
    
    return(list(beta      = betaT,
                G         = G.lstT,
                sig2      = sig2T,
                loglik    = ll,
                post.prob = wik,
                cluster   = cluster))
  }
  
  else{
    
    #pp=pp.cur;dat.sh=dat.sh.all;beta=beta;S2=sig2; G=G; C2=sigc2
    #X=X;U=U;V=V;pp=pp.cur;dat.sh=dat.sh.all;
    #beta=beta.last;S2=sig2.last; G=G.last; C2=sigc2.last; alpha=alpha.cur;
    #model=model;R=R;ntimepts=ntimepts;lambda=lambdaN[e]
    #print("MLEs for FULLY specified model")
    covfcn     <- function(g,c2,s2){U%*%g%*%t(U)+c2*V%*%t(V)+s2*diag(dim(X)[1])}
    cov.y      <- mapply(covfcn,g=G,c2=as.list(C2),s2=as.list(S2),SIMPLIFY=FALSE) 
    inv.cov.y  <- lapply(cov.y,function(c){solve(c)})
    covDfcn    <- function(g,inv.err,c2,s2){U%*%(g-g%*%t(U)%*%inv.err%*%U%*%g)%*%t(U) + V%*%(c2*diag(R)-c2^2*t(V)%*%inv.err%*%V)%*%t(V)}
    cov.d      <- mapply(covDfcn,g=G,inv.err=inv.cov.y,c2=as.list(C2),s2=as.list(S2),SIMPLIFY=FALSE)
    
    B1 <-lapply(inv.cov.y,function(c) {solve(t(X)%*%X)%*%t(X)%*%c})
    T1 <-lapply(inv.cov.y,function(c){c%*%V%*%t(V)%*%c})
    T2 <-lapply(inv.cov.y,function(c){sum(diag(c%*%V%*%t(V)))})
    S1 <-lapply(inv.cov.y,function(c){c%*%c})
    GF1<-matrix(mapply(function(g,c){g%*%t(U)%*%c},g=G,c=inv.cov.y,SIMPLIFY=F))
    GF2<-mapply(function(g,c){g - g%*%t(U)%*%c%*%U%*%g},g=G,c=inv.cov.y,SIMPLIFY=F)
    
    b.fcn1     <- function(m){beta[,k] + S2[k]*B1[[k]]%*%m}
    G.fcn1     <- function(m){quad.tform(m%*%t(m),GF1[[k]]) + GF2[[k]]}
    s2.fcn1    <- function(m){S2[k]^2*quad.form(S1[[k]],m) + sum(diag(cov.d[[k]]))}
    c2.fcn1    <- function(m){C2[k]^2*quad.form(T1[[k]],m) + R*C2[k] - C2[k]^2*T2[[k]]}
    
    betaT  <- matrix(NA,ncol=dim(pp)[2],nrow=dim(beta)[1])
    sig2T  <- matrix(NA,ncol=dim(pp)[2],nrow=1)
    sigc2T <- matrix(NA,ncol=dim(pp)[2],nrow=1)
    G.lstT <- list()
    loglik <- matrix(NA,ncol=length(alpha),nrow=dim(dat.sh)[1])
    
    for(k in 1:dim(pp)[2]){
      Ym          <- t(apply(dat.sh,1,function(Y){Y-mu[,k]}))
      betaT[,k]   <- (pp[,k]%*%t(apply(Ym,1,b.fcn1))) / ppk[k] 
      G.lstT[[k]] <- matrix((pp[,k]%*%matrix(t(apply(Ym,1,G.fcn1)),ncol=dim(G[[k]])[1]^2))/ ppk[k],ncol=dim(G[[k]])[1]) 
      sig2T[k]    <- (pp[,k]%*%matrix(t(apply(Ym,1,s2.fcn1)),ncol=1))/ (ppk[k]*R*ntimepts) 
      sigc2T[k]   <- (pp[,k]%*%matrix(t(apply(Ym,1,c2.fcn1)),ncol=1))/ (ppk[k]*R)
      
      covy        <- U%*%G.lstT[[k]]%*%t(U)+sigc2T[k]*V%*%t(V)+sig2T[k]*diag(ntimepts*R)
      loglik[,k]  <- alpha[k]*dmvn(X=as.matrix(dat.sh), mu=X%*%betaT[,k], sigma=covy, log=FALSE)
    }
    wik<-loglik/apply(loglik,1,sum)
    
    whichfcn<-function(w){which(w==max(w))}
    cluster<-apply(wik,1,whichfcn)
    ll<-sum(log(apply(loglik,1,sum)))+ lambda*dim(pp)[1]*sum(alpha*log(alpha))
   
    return(list(beta=betaT,
                G=G.lstT,
                sig2=sig2T,
                sigc2=sigc2T,
                loglik=ll,
                post.prob = wik,
                cluster   = cluster))
  }
  }
}


###############################  FINAL EM Algorithm function #######################################
EPEM_function<-function(dat.use,    ## Input dataset
                        timepoints, ## vector of unique time-points used in the dataset
                        Ruse,       ## Nubmer of replicates
                        iterview)        
{
  
  dat.sh.RepLvl   <- dat.use[,c(paste("V1.",timepoints,sep=""),"GeneID","RepNum")]
  dat.sh.GeneLvl  <- aggregate(dat.use[,c(paste("V1.",timepoints,sep=""))], by=list(dat.use$GeneID), FUN=mean)
  names(dat.sh.GeneLvl)<-c("GeneID",paste("V1.",timepoints,sep=""))
  head(dat.sh.GeneLvl)
  
  gathercols      <-c(paste("V1.",timepoints,sep=""))
  
  ## Create long form of the dataset (either Genelvl or Replvl)
  if(Ruse==1){
    ## Long form of the Gene-Level data and short form where each row corresponds to the T observations for each gene
    dat.l.T        <-gather_(dat.sh.GeneLvl, "Day", "Y", gathercols)
    dat.l.T$Day    <-as.numeric(gsub('\\V1.', '', dat.l.T$Day))
    dat.l          <-dat.l.T[order(dat.l.T$GeneID,dat.l.T$Day),]
    dat.l.T$plotid <-dat.l.T$GeneID
    dat.Temp       <- reshape(cbind(dat.l[,c("Y","GeneID")],Day2=paste(dat.l$Day,sep='.')), 
                              idvar = c("GeneID"), timevar = "Day2", direction = "wide")
    
  } else if(Ruse > 1) {
    ## Long form of the Rep-Level data and short form where each row corresponds to the T*R observations for each gene
    dat.l.T        <-gather_(dat.sh.RepLvl, "Day", "Y", gathercols)
    dat.l.T$Day    <-as.numeric(gsub('\\V1.', '', dat.l.T$Day))
    dat.l.T$plotid <-paste(dat.l.T$GeneID,dat.l.T$RepNum,sep='_')
    dat.l        <-dat.l.T[order(dat.l.T$GeneID,dat.l.T$RepNum,dat.l.T$Day),]
    
    dat.Temp     <- reshape(cbind(dat.l[,c("Y","GeneID")],Day2=paste(dat.l$Day,dat.l$RepNum,sep='.')), 
                            idvar = c("GeneID"), timevar = "Day2", direction = "wide")
  }
  
  ## Data matrix with only expression values 
  dat.sh.all     <- dat.Temp[, !(names(dat.Temp) %in% c("GeneID"))]
  
  
  dat.sh.dim  <- dim(dat.sh.all)
  dat.sh.nrow <- dat.sh.dim[1]
  dat.sh.ncol <- dat.sh.dim[2]
  
  
  ngenes     = dat.sh.nrow   #Number of genes in dataset
  R          = Ruse          #Number of replicates (strains)
  time       = c(unique(dat.l$Day))                              #Unique time-points- we are assuming they are the same for each gene
  eta        = min(1,0.5^floor((length(time)/2)-1))              #Input to calculate lambda
  ntimepts   = length(time)
  
  ##------------------------------------------------------------------
  
  drops             <- c("RepNum","GeneID")
  PolyReg.RepLvl    <- data.frame(GeneID=dat.sh.RepLvl$GeneID,
                                  t(apply(dat.sh.RepLvl[,!(names(dat.sh.RepLvl) %in% drops)],1,init1,p=p)))
  
  PolyReg.GeneLvl   <- data.frame(GeneID=dat.Temp$GeneID,
                                  t(apply(dat.Temp[,!(names(dat.Temp) %in% drops)],1,init2,p=p)))
  
  names(PolyReg.RepLvl)   <-c('GeneID',paste('b',seq(0,p,1),sep=''),'s2','p.one')  
  names(PolyReg.GeneLvl)  <-c('GeneID',paste('b',seq(0,p,1),sep=''),'s2','p.one')
  
  #################################################
  ## Initializations to EM algorithm
  s2.temp<-PolyReg.GeneLvl$s2
  K0     <- ngenes
  true.beta<-NA
  
  beta0  <- t(PolyReg.GeneLvl[,paste('b',0:p,sep='')])
  cov0   <- rep(median(s2.temp),K0)
  alpha0 <- rep(1/K0,K0)
  
  X       <- vandermonde.matrix(rep(time,R), p+1)
  U       <- vandermonde.matrix(rep(time,R), q+1) 
  G.T     <- solve(crossprod(vandermonde.matrix(time, q+1)))
  G0.lst  <- lapply(1:K0, function(j) cov0[1]*G.T)#cov0[1]*GMat)#*G.T)
  c20     <- rep(median(PolyReg.RepLvl$s2),K0)
  T       <- length(time) 
  V       <- rep(rep(c(1,0), c(T, T*R)), R)  
  V       <- V[1:(T*R*R)]
  dim(V)  <- c(T*R, R) 
  tcrossprodV <- tcrossprod(V) 
  
  # define some commom calc. terms: 
  crossprodX  <- crossprod(X)
  solveX      <- solve(crossprodX)%*%t(X)
  XY          <- as.matrix(dat.sh.all)%*%X
  tcrossprodV <- tcrossprod(V) 
  
  id        <- dat.sh.GeneLvl$GeneID
  group     <- NA
  
  if(R == 1){model = "ME2NR"}
  if(R > 1) {model = "ME2"}
  
  ########## INITIALIZE posterior probability matrix
  pp0       <- ppfun(X=X,U=U,V=V,dat.sh=dat.sh.all,beta=beta0,S2=cov0,G=G0.lst,C2=c20,alpha=alpha0,model=model)
  
  post.prob0<-pp0$post.prob  
  summary(c(post.prob0))  
  
  diff.vector<- rep(NA,B)
  K          <- rep(NA,B)
  lambdaN    <- rep(NA,B)
  loglik.iter<- rep(NA,B)
  
  ptm1       <- proc.time()
  
  for(e in 1:B){
    if(e==1){
      alpha         <- alphafun(pp=post.prob0,lambda=1,a.last=alpha0)
      lambdaN[e]    <- lambda.fun(eta           = eta,
                                  alpha.cur     = alpha,
                                  alpha.last    = alpha0,
                                  pp            = post.prob0)
      i.lt          <- alpha< (1/K0)
      #which(i.lt)
      if(sum(i.lt,na.rm=TRUE) > 0) {
        K[e]         <- K0-sum(i.lt)
        pp.curT      <- post.prob0[,-c(which(i.lt))]
        alpha.curT   <- alpha[-c(which(i.lt))]
        pp.cur       <- pp.curT/apply(pp.curT,1,sum)
        alpha.cur    <- alpha.curT/sum(alpha.curT)
        
        G.last         <-G0.lst[-c(which(i.lt))]
        sig2.last      <-  cov0[-c(which(i.lt))]
        sigc2.last     <-   c20[-c(which(i.lt))]
        beta.last      <- beta0[,-c(which(i.lt))]
        
      } else {
        K[e]       <- length(alpha)
        pp.cur     <- post.prob0
        alpha.cur  <- alpha
        G.last     <- G0.lst
        sig2.last  <- cov0
        sigc2.last <- c20
        beta.last  <- beta0
      }
      
      
      mle            <- mlefun(X=X,U=U,V=V,pp=pp.cur,dat.sh=dat.sh.all,
                               beta=beta.last,S2=sig2.last, G=G.last, C2=sigc2.last, alpha=alpha.cur,
                               model=model,R=R,ntimepts=ntimepts,lambda=lambdaN[e])
      loglik.iter[e] <- mle$loglik
      post.prob     <- mle$post.prob
      summary(c(post.prob))
      cur.beta      <- mle$beta
      
      diff          <- mle$beta-beta0[,-c(which(i.lt))]
      norm.beta<-NULL
      for(k in 1:dim(mle$beta)[2]){
        norm.beta[k] <- sqrt(sum(diff[,k]^2))
      }
      diff.vector[e] <- max(norm.beta)
      loglikdiff<-NA
      
    }  else {
      alpha        <- alphafun(pp=post.prob,lambda=lambdaN[e-1],a.last=alpha.cur) #
      if(e > 61 && (K[e-60]-K[e-1]) == 0) { 
        lambdaN[e]   <- 0 
      } else {
        lambdaN[e]   <- lambda.fun(eta           = eta,
                                   alpha.cur     = alpha,
                                   alpha.last    = alpha.cur,
                                   pp            = post.prob)
      }
      i.lt         <- alpha < (1/K0)
      which(i.lt)
      
      
      if(sum(i.lt,na.rm=TRUE) > 0) {
        K[e]         <- length(alpha)-sum(i.lt)
        pp.curT      <- post.prob[,-c(which(i.lt))]
        alpha.curT   <- alpha[-c(which(i.lt))]
        pp.cur       <- pp.curT/apply(pp.curT,1,sum)
        alpha.cur    <- alpha.curT/sum(alpha.curT)
        G.last       <- mle$G[-c(which(i.lt))]
        sig2.last    <- mle$sig2[-c(which(i.lt))]
        sigc2.last   <- mle$sigc2[-c(which(i.lt))]
        beta.last    <- mle$beta[,-c(which(i.lt))]
        
      } else {
        K[e]       <- length(alpha)
        pp.cur     <- post.prob
        alpha.cur  <- alpha
        G.last     <- mle$G
        sig2.last  <- mle$sig2
        sigc2.last <- mle$sigc2
        beta.last  <- mle$beta
      }
      
      ## If we get rid of all non-zero values for one row/gene, then need to reassign a small number
      ## to each post.prob to prevent errors
      notok<-is.na(apply(pp.cur,1,sum))==TRUE
      pp.cur[notok,]
      pp.cur[notok,]<-rep(1E-6,dim(pp.cur)[2]) 
      
      mle            <- mlefun(X=X,U=U,V=V,pp=pp.cur,dat.sh=dat.sh.all,
                               beta=beta.last,S2=sig2.last, G=G.last, C2=sigc2.last, alpha=alpha.cur, 
                               model=model,R=R,ntimepts=ntimepts,lambda=lambdaN[e])
      loglik.iter[e] <- mle$loglik
      post.prob      <- mle$post.prob
      cluster        <- unlist(mle$cluster)
      
      if(sum(i.lt,na.rm=TRUE) > 0) {cur.beta <- cur.beta[,-c(which(i.lt))]}
      
      diff           <- mle$beta-cur.beta
      norm.beta<- NULL
      for(k in 1:dim(mle$beta)[2]){
        norm.beta[k] <- sqrt(sum(diff[,k]^2))
      }
      diff.vector[e] <- max(norm.beta)
      
      alpha.F<-NULL
      for(k in 1:K[e]){ alpha.F[k]<-sum(post.prob[,k])/K0}
      
      loglikdiff<-abs(loglik.iter[e]-loglik.iter[e-1])
      
      if(max(norm.beta)<epsilon && lambdaN[e]==0){  # loglikdiff<epsilon
        lambdaN
        converge=1
        print(paste("converged! Num Iter=",e,sep=""))
        break
      } else {
        cur.beta <- mle$beta
        converge=0
        
      }
    }
    
    if(iterview==TRUE){
      print(paste('K=',K[e]))
      print(paste('betaDiff=',diff.vector[e]))
      print(paste('lambda=',lambdaN[e]))
      print(paste('LogLik Diff=',loglikdiff))
      print(paste('IterNum=',e))
      if(K[e] < 100) table(cluster)
    }
  } 
  ptm2<-proc.time()
  
  
  ## output result after convergence
  list(beta         =mle$beta,
       sig2         =mle$sig2,
       tau2         =mle$sigc2,
       G            =mle$G,
       alpha        =alpha.F,
       K            =K[1:e],
       loglik.iter  =loglik.iter[1:e],
       IterNum      =seq(1,e,1),
       diff.vector  =diff.vector[1:e],
       e.last       =e,
       lambda       =lambdaN[1:e],
       convtime.sec =(ptm2-ptm1)[3],
       converge     = converge
  )
}
