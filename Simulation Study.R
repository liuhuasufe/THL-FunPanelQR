###################################################
#### Simulation Functional QR Panel Data  ###
###################################################

rm(list=ls())

library(conquer)
library(phtt)
library(progress)
library(xtable)
library(NbClust)
library(mclust)

######################################
## Setup: 
######################################

Sc<-2
mcreps=500
n<-100
T<-50

######################################
## Setup (End)
######################################

######################################
### Paramters & Scenarios:
######################################

xgrid<-seq(0,1,length.out=100)
tau<-c(0.25,0.5,0.75)
alpha0=5*cos(seq(1,T)/pi)
beta=5*sin(seq(1,T)/pi) # needs to be time invariant!!
alpha2<-function(s) {-s*2+8*s^2+5*s^3+2*sin(s*8)}
alpha3<-function(s) {-s*2+cos(6*s)}
if (Sc==1) {alpha1<-function(s) {sqrt(2)*sin(pi*s/2)-s^3/2+sqrt(18)*sin(3*pi*s/2)}}
if (Sc==2) {alpha1<-function(s) {s*8-4*s^2-5*s^3+2*sin(s*8)} }

######################################
### Paramters & Scenarios (End)
######################################


#######################################
## Data Generating Propcess:
#######################################
phi<-function(j,s) {sqrt(2)*sin((j-1/2)*pi*s)}
product1<-function(j,s) {phi(j=j,s)*alpha1(s)}
product2<-function(j,s) {phi(j=j,s)*alpha2(s)}
product3<-function(j,s) {phi(j=j,s)*alpha3(s)}
inner1<-function(j) {integrate(function(s) product1(j=j,s),0,1)$value}
inner2<-function(j) {integrate(function(s) product2(j=j,s),0,1)$value}
inner3<-function(j) {integrate(function(s) product3(j=j,s),0,1)$value}
M=20 # that many summands for X 
I1<-mapply(inner1,j=1:M)
I2<-mapply(inner2,j=1:M)
I3<-mapply(inner3,j=1:M)
sqrtlambda<-function(j) {((j-1/2)*pi)^(-1)} # Eigenvalues
sqL<-mapply(sqrtlambda,j=1:M)
Th<- function(it) {sqL*rnorm(M)} # M-vector of scores / function of it creates ability to repeat draws
#######################################
## Data Generating Propcess (End)
#######################################

###################################################
######## Monte Carlo 
###################################################

pb <- progress_bar$new(total = mcreps)

for (k in 1:3){
  #######################################
  ### Storage of Results:
  #######################################
  k=1
  CError<-rep(NA,mcreps)
  CError_TH<-rep(NA,mcreps)
  ErKm<-rep(NA,mcreps)
  K_hat<-rep(NA,mcreps)
  InterceptMC<-matrix(NA,mcreps,T)
  BetaMC<-matrix(NA,mcreps,T)
  AlphaMC<-array(NA,dim=c(mcreps,3,length(xgrid)))
  #CI_asy_up<-array(NA,dim=c(mcreps,3,length(xgrid)))
  #CI_asy_low<-array(NA,dim=c(mcreps,3,length(xgrid)))
  CI_boot_up<-array(NA,dim=c(mcreps,3,length(xgrid)))
  CI_boot_low<-array(NA,dim=c(mcreps,3,length(xgrid)))
  CI_up<-array(NA,dim=c(mcreps,T,length(xgrid)))
  CI_low<-array(NA,dim=c(mcreps,T,length(xgrid)))
  
  A1 <- matrix(NA,nrow = mcreps,ncol = length(xgrid))
  A2 <- matrix(NA,nrow = mcreps,ncol = length(xgrid))
  A3 <- matrix(NA,nrow = mcreps,ncol = length(xgrid))
  
  dd=0
  #mcreps
  for (r in 1:mcreps){
    print(r)
    ###############################################################################
    ####### Random Number Generation (Beginning)
    ###############################################################################
    ThM<-array(NA,dim=c(T,n,M))
    for (t in 1:T) {ThM[t,,]<-t(mapply(Th,it=1:n))}
    
    Xobs<-array(NA,dim=c(T,n,length(xgrid))) # 'Observed' X-values
    PhiEval<-matrix(NA,length(xgrid),M)
    fctl<-matrix(NA,n,T)  # Functional (including regimes already)
    
    for (j in 1:M){ PhiEval[,j]<-phi(j=j,s=xgrid) } # True Basis Functions evaluated at xgrid
    for (t in 1:T){
      Xobs[t,,]<-ThM[t,,]%*%t(PhiEval)    # observed X-values on xgrid
      if (t<=(T/3)) {fctl[,t]<-t(I1)%*%t(ThM[t,,])} # functionals in 1st regime
      if (t>(T/3) & t<=(2*T/3) ) {fctl[,t]<-t(I2)%*%t(ThM[t,,])} # functionals in 2nd regime
      if (t>(2*T/3)) {fctl[,t]<-t(I3)%*%t(ThM[t,,])} # functionals in 3rd regime
    }
    
    #case1  eps follows a standard normal distribution
    eps<-matrix(rnorm(n*T),n,T)
    epsilon = eps - qnorm(tau[k])
    z<-matrix(rnorm(n*T),n,T)
    y<-matrix(rep(alpha0,n),n,T,byrow=T)+fctl+z%*%diag(beta,T)+epsilon
    ###############################################################################
    ####### Random Number Generation (End)
    ###############################################################################
    
    
    ###############################################################################
    ####### 1st Stage Estimation (Beginning)
    ###############################################################################  
    alpha0_hat<-colMeans(y)  
    InterceptMC[r,]<-(alpha0_hat-alpha0)/alpha0
    #ybreve<-y-matrix(rep(alpha0_hat,n),n,T,byrow=T) # subtract timespecific means (note: rows=i / cols= t)
    ybreve<-y
    z<-z-matrix(colSums(z)/nrow(z),n,T,byrow=T) ## Demean z
    #espilon<-epsilon-matrix(colSums(epsilon)/nrow(epsilon),n,T,byrow=T) ## Demean eps
    espilon<-epsilon
    for (s in 1:T){Xobs[s,,]<- Xobs[s,,]-matrix(colSums(Xobs[s,,])/nrow(Xobs[s,,]),n,length(xgrid),byrow=T)} ## Demean X
    
    Sig<-array(NA,dim=c(T,length(xgrid),length(xgrid)))
    ## Implementation of Information Criteria (t-wise):
    m_data<-rep(NA,T)
    for (t in 1:T) {
      Sig[t,,]<-t(Xobs[t,,])%*%Xobs[t,,]/(n*length(xgrid))
      suppressWarnings(m_data[t]<-c(as.list(OptDim(t(Xobs[t,1:n,]),criteria = "ER"))[[1]]$`Optimal Dimension`))
    } 
    
    #plot(mapply(w_t,p=1:15))
    #summary(prcomp(Sig[t,,]))
    
    m_hat<-max(round(mean(m_data)),2)
    #m_hat=3
    Phi_hat<-array(NA,dim=c(T,length(xgrid),m_hat))
    Lmbd<-matrix(NA,T,m_hat)
    for (t in 1:T) {  
      Phi_hat[t,,]<-eigen(Sig[t,,])$vectors[,1:m_hat]*(length(xgrid)^(1/2))
      Lmbd[t,]<-eigen(Sig[t,,])$values[1:m_hat]
    }
    Theta_hat<-array(NA,dim=c(T,n,m_hat))
    for (t in 1:T) {
      for (i in 1:n){
        for (j in 1:m_hat){
          Theta_hat[t,i,j]<- t(Xobs[t,i,])%*%Phi_hat[t,,j]/length(xgrid) # Riemann Approx to integral
        }
      }
    }
    
    Mod<-function(t){conquer(cbind(z[,t],Theta_hat[t,,]),t(t(ybreve[,t])),
                             tau = tau[k],kernel = "parabolic",ci=TRUE)}
    SigEps<-mapply(function(t){sum(Mod(t=t)$residual^2)/(n-m_hat-1)},t=1:T)
    coef_hat<-function(t) { Mod(t=t)$coeff} 
    Coef_Hat<-mapply(coef_hat,1:T)
    a_hat<-Coef_Hat[3:(m_hat+2),]
    beta_hat<-Coef_Hat[2,]
    
    a_hat_d<-matrix(NA,m_hat,T)
    for (t in 1:T) {  a_hat_d[,t]<-a_hat[,t]*(Lmbd[t,]^(1/2))/(SigEps[t]^(1/2))}
    
    alph_d<-matrix(NA,length(xgrid),T)
    for(t in 1:T) {  alph_d[,t]<-Phi_hat[t,,]%*%a_hat_d[,t]}
    
    alpha_hat<-matrix(NA,length(xgrid),T)
    sgm <- array(data=NA,T)
    for(t in 1:T) {  
      alpha_hat[,t]<-Phi_hat[t,,]%*%a_hat[,t]
      sgm <- sqrt(diag(Phi_hat[t,,] %*% diag((a_hat[,t]-Mod(t=t)$perCI[3:(m_hat+2),1])/1.96) %*% t(Phi_hat[t,,])))
      CI_up[r,t,] <- alpha_hat[,t] + 1.96 * sgm
      CI_low[r,t,] <- alpha_hat[,t] - 1.96 * sgm
    }
    
    ################################################################
    ## 2nd stage  (Beginning)
    ################################################################
    ###########    BIC    #############
    Ds_bs2<-function(s2) {sum((alpha_hat[,s2+1]-alpha_hat[,s2])^2)/30}
    dst_bs2<-mapply(Ds_bs2,s2=1:(T-1))
    ord = order(dst_bs2)
    G1<-c(1:min(ord[T-1],ord[T-2]));n1<-length(G1)
    if(length(G1)==1){print("G1==1");dd = dd+1;next()}
    G2<-c((min(ord[T-1],ord[T-2])+1):max(ord[T-1],ord[T-2]));n2<-length(G2)
    if(length(G2)==1){print("G2==1");dd = dd+1;next()}
    G3<-c((max(ord[T-1],ord[T-2])+1):T);n3<-length(G3)
    if(length(G3)==1){print("G3==1");dd = dd+1;next()}
    
    #######################  END  ###################################
    
    g1star<-1:floor((T/3))
    g2star<-(tail(g1star,1)+1):floor(2*T/3)  #(T/2+1):(T)
    g3star<-(tail(g2star,1)+1):T
    
    CError[r]<-length(which(!G1%in%g1star))+length(which(!G2%in%g2star))+length(which(!G3%in%g3star))
    K_hat[r]<-NbClust(t(alpha_hat[1:length(xgrid),]), distance = "euclidean", min.nc=2, max.nc=8,method = "complete", index = "ch")$Best.nc["Number_clusters"]
    kmclst<-Mclust(t(a_hat),3,verbose=FALSE)$classification
    g1km<-which(kmclst==1)
    g2km<-which(kmclst==2)
    g3km<-which(kmclst==3)
    e1km<-min(c(length(which(!g1km%in%g1star)),length(which(!g1km%in%g2star)),length(which(!g1km%in%g3star))))
    e2km<-min(c(length(which(!g2km%in%g1star)),length(which(!g2km%in%g2star)),length(which(!g2km%in%g3star))))
    e3km<-min(c(length(which(!g3km%in%g1star)),length(which(!g3km%in%g2star)),length(which(!g3km%in%g3star))))
    ErKm[r]<-(e1km+e2km+e3km)
    
    Ds_TH<-function(s2) {sum((alph_d[,1]-alph_d[,s2])^2)/length(xgrid)} # Riemann Approx
    dst_TH<-mapply(Ds_TH,s2=1:T)
    tau_TH<-qchisq(p=.99,df=m_hat)/(n/2)
    #plot(dst)
    G1_TH<-which(dst_TH<tau_TH) ; n1_TH<-length(G1_TH)
    G2or3_TH<-which(dst_TH>=tau_TH) 
    Ds2or3_TH<-function(s2) {sum((alph_d[,g2star[1]]-alph_d[,s2])^2)/length(xgrid)} # Riemann Approx
    dst2or3_TH<-mapply(Ds2or3_TH,s2=G2or3_TH)
    #plot(dst2or3) ; abline(h=2*tau)
    G2_TH<-G2or3_TH[which(dst2or3_TH<tau_TH)] ; n2_TH<-length(G2_TH)
    G3_TH<-G2or3_TH[which(dst2or3_TH>=tau_TH)] ; n3_TH<-length(G3_TH)
    
    CError_TH[r]<-length(which(!G1_TH%in%g1star))+length(which(!G2_TH%in%g2star))+length(which(!G3_TH%in%g3star))
    
    
    ################################################################
    ## 2nd stage  (End)
    ################################################################
    
    ################################################################
    ## 3rd stage  (Beginning)
    ################################################################
    
    m_tilde<-c(NA,NA,NA); #c(NA,NA)
    mean_t = array(data = NA,c(T,n,length(xgrid)))
    for(s in 1:T){mean_t[s,,] = matrix(colSums(Xobs[s,,])/nrow(Xobs[s,,]),n,length(xgrid),byrow=T)}
    mean1 = mean_t[G1,,]
    mean2 = mean_t[G2,,]
    mean3 = mean_t[G3,,]
    meanG1 = array(data = NA,c(1,n,length(xgrid)))
    meanG2 = array(data = NA,c(1,n,length(xgrid)))
    meanG3 = array(data = NA,c(1,n,length(xgrid)))
    for(i in 1:n){meanG1[,i,] = matrix((colSums(mean1[,i,]))/nrow(mean1[,i,]),1,length(xgrid),byrow=T)}
    for(i in 1:n){meanG2[,i,] = matrix((colSums(mean2[,i,]))/nrow(mean2[,i,]),1,length(xgrid),byrow=T)}
    for(i in 1:n){meanG3[,i,] = matrix((colSums(mean3[,i,]))/nrow(mean3[,i,]),1,length(xgrid),byrow=T)}
   
    Xc = array(data = NA,c(T,n,length(xgrid)))
    for (s in 1:T){
      if(s %in% G1){Xc[s,,]<- Xobs[s,,]-meanG1[1,,]}
      if(s %in% G2){Xc[s,,]<- Xobs[s,,]-meanG2[1,,]}
      if(s %in% G3){Xc[s,,]<- Xobs[s,,]-meanG3[1,,]}
    }
    
    if(length(G1)>1) {Xobs1<-apply(aperm(Xc[G1,,],c(2,1,3)), 3, rbind)}
    if(length(G2)>1) {Xobs2<-apply(aperm(Xc[G2,,],c(2,1,3)), 3, rbind)}
    if(length(G3)>1) {Xobs3<-apply(aperm(Xc[G3,,],c(2,1,3)), 3, rbind)}
   
    if(length(G1)==1) {Xobs1<-Xobs[G1,,]}
    if(length(G2)==1) {Xobs2<-Xobs[G2,,]}
    if(length(G3)==1) {Xobs3<-Xobs[G3,,]}
    
    #T*T
    Sig1tilde<-t(Xobs1)%*%Xobs1/(dim(Xobs1)[1]*length(xgrid)) ; 
    Sig2tilde<-t(Xobs2)%*%Xobs2/(dim(Xobs2)[1]*length(xgrid)) ; 
    Sig3tilde<-t(Xobs3)%*%Xobs3/(dim(Xobs3)[1]*length(xgrid)) ; 
 

    m_tilde[1]<-max(2,suppressWarnings(c(as.list(OptDim(t(Xobs1),criteria = "ER"))[[1]]$`Optimal Dimension`)))
    m_tilde[2]<-max(2,suppressWarnings(c(as.list(OptDim(t(Xobs2),criteria = "ER"))[[1]]$`Optimal Dimension`)))
    m_tilde[3]<-max(2,suppressWarnings(c(as.list(OptDim(t(Xobs3),criteria = "ER"))[[1]]$`Optimal Dimension`)))

    
    Phi1tilde<-eigen(Sig1tilde)$vectors[,1:(m_tilde[1])]*(length(xgrid)^(1/2)) #  
    Phi2tilde<-eigen(Sig2tilde)$vectors[,1:(m_tilde[2])]*(length(xgrid)^(1/2)) # 
    Phi3tilde<-eigen(Sig3tilde)$vectors[,1:(m_tilde[3])]*(length(xgrid)^(1/2)) #

    #Theta
    Th1tilde<-Xobs1%*%Phi1tilde/length(xgrid)
    Th2tilde<-Xobs2%*%Phi2tilde/length(xgrid)
    Th3tilde<-Xobs3%*%Phi3tilde/length(xgrid)
    
    ### 1st Regime
    y1tilde<-c((ybreve[,G1]-z[,G1]%*%diag(beta_hat[G1],length(G1))))
    Mod1<-conquer(Th1tilde,y1tilde,tau = tau[k],kernel = "parabolic",ci=TRUE)
    SigEps<-mapply(function(t){sum(Mod1$residual^2)/(n-m_hat-1)},t=1:T)
    a1_hat<-Mod1$coeff[2:(m_tilde[1]+1)]
    alpha1tilde<-Phi1tilde%*%a1_hat
    sgm1 <- sqrt(diag(Phi1tilde %*% diag((a1_hat-Mod1$perCI[2:(m_tilde[1]+1),1])/1.96) %*% t(Phi1tilde)))
    CI_boot_up[r,1,] <- alpha1tilde + 1.96 * sgm1
    CI_boot_low[r,1,] <- alpha1tilde - 1.96 * sgm1

    ### 2nd Regime
    y2tilde<-c((ybreve[,G2]-z[,G2]%*%diag(beta_hat[G2],length(G2))))
    Mod2<-conquer(Th2tilde,y2tilde,tau = tau[k],kernel = "parabolic",ci=TRUE)
    SigEps<-mapply(function(t){sum(Mod2$residual^2)/(n-m_hat-1)},t=1:T)
    a2_hat<-Mod2$coeff[2:(m_tilde[2]+1)]
    alpha2tilde<-Phi2tilde%*%a2_hat
    sgm2 <- sqrt(diag(Phi2tilde %*% diag((a2_hat-Mod2$perCI[2:(m_tilde[2]+1),1])/1.96) %*% t(Phi2tilde)))
    CI_boot_up[r,2,] <- alpha2tilde + 1.96 * sgm2
    CI_boot_low[r,2,] <- alpha2tilde - 1.96 * sgm2
    
    ### 3rd Regime
    y3tilde<-c((ybreve[,G3]-z[,G3]%*%diag(beta_hat[G3],length(G3))))
    Mod3<-conquer(Th3tilde,y3tilde,tau = tau[k],kernel = "parabolic",ci=TRUE)
    SigEps<-mapply(function(t){sum(Mod3$residual^2)/(n-m_hat-1)},t=1:T)
    a3_hat<-Mod3$coeff[2:(m_tilde[3]+1)]
    alpha3tilde<-Phi3tilde%*%a3_hat
    sgm3 <- sqrt(diag(Phi3tilde %*% diag((a3_hat-Mod3$perCI[2:(m_tilde[3]+1),1])/1.96) %*% t(Phi3tilde)))
    CI_boot_up[r,3,] <- alpha3tilde + 1.96 * sgm3
    CI_boot_low[r,3,] <- alpha3tilde - 1.96 * sgm3
    ################################################################
    ## 3rd stage  (End)
    ################################################################
    
    
    ################################################################
    ## 4th stage  (Beginning)
    ################################################################
    
    #BetaMC[r,]<-(c(b1tilde,b2tilde)-beta[c(G1,G2)])/beta[c(G1,G2)]
    BetaMC[r,]<-(beta_hat-beta)
    AlphaMC[r,1,]<-alpha1tilde
    AlphaMC[r,2,]<-alpha2tilde
    AlphaMC[r,3,]<-alpha3tilde
    ################################################################
    ## 4th stage  (End)
    ################################################################
    #print(r)
    
  }
  
  A1nrm<-integrate(function(u) {(alpha1(u))^2},lower=0,upper=1)$value
  A2nrm<-integrate(function(u) {(alpha2(u))^2},lower=0,upper=1)$value
  A3nrm<-integrate(function(u) {(alpha3(u))^2},lower=0,upper=1)$value
  
  ############### Save Results:
  
  mCerror<-(CError)/T
  mCerror_TH<-(CError_TH)/T
  mfc<-function(x)  {sum(x^2)/length(x)}
  Bbias<-(mapply(function(s){sum(BetaMC[s,])/length(BetaMC[s,])}, s=1:nrow(BetaMC)))
  BError<-(mapply(function(s) {mfc(BetaMC[s,])},s=1:nrow(BetaMC)))
  A0Error<-(mapply(function(s) {mfc(InterceptMC[[s]])},s=1:length(InterceptMC)))^2
  A1 = na.omit(AlphaMC[,1,])
  A2 = na.omit(AlphaMC[,2,])
  A3 = na.omit(AlphaMC[,3,])
  Alph1Err<-A1-matrix(rep(c(alpha1(xgrid)),mcreps-dd),mcreps-dd,length(xgrid),byrow=T)
  Alph2Err<-A2-matrix(rep(c(alpha2(xgrid)),mcreps-dd),mcreps-dd,length(xgrid),byrow=T)
  Alph3Err<-A3-matrix(rep(c(alpha3(xgrid)),mcreps-dd),mcreps-dd,length(xgrid),byrow=T)
  
  Alph1ErrV<-(rowSums(Alph1Err^2)/length(xgrid))/A1nrm
  Alph2ErrV<-(rowSums(Alph2Err^2)/length(xgrid))/A2nrm
  Alph3ErrV<-(rowSums(Alph3Err^2)/length(xgrid))/A3nrm
  
  #KmK<-K_hat-2
  ErrMat<-matrix(NA,8,5)
  ErrMat[1,]<-round(c(quantile(Bbias,p=c(.25,.5),na.rm = T),mean(Bbias,na.rm = T),quantile(Bbias,p=c(.75),na.rm = T),sqrt(var(Bbias,na.rm = T))),3)
  ErrMat[2,]<-round(c(quantile(BError,p=c(.25,.5),na.rm = T),mean(BError,na.rm = T),quantile(BError,p=c(.75),na.rm = T),sqrt(var(BError,na.rm = T))),3)
  #ErrMat[2,]<-round(c(quantile(A0Error,p=c(.25,.5)),mean(A0Error),quantile(A0Error,p=c(.75)),sqrt(var(A0Error))),2)
  ErrMat[3,]<-round(c(quantile(mCerror,p=c(.25,.5),na.rm = T),mean(mCerror,na.rm = T),quantile(mCerror,p=c(.75),na.rm = T),sqrt(var(mCerror,na.rm = T))),3)
  ErrMat[4,]<-round(c(quantile(mCerror_TH,p=c(.25,.5),na.rm = T),mean(mCerror_TH,na.rm = T),quantile(mCerror_TH,p=c(.75),na.rm = T),sqrt(var(mCerror_TH,na.rm = T))),3)
  #ErrMat[3,]<-round(c(quantile(KmK,p=c(.25,.5)),mean(KmK),quantile(KmK,p=c(.75)),sqrt(var(KmK))),2)
  ErrMat[5,]<-round(c(quantile(ErKm/T,p=c(.25,.5),na.rm = T),mean(ErKm/T,na.rm = T),quantile(ErKm/T,p=c(.75),na.rm = T),sqrt(var(ErKm/T,na.rm = T))),3)
  ErrMat[6,]<-round(c(quantile(Alph1ErrV,p=c(.25,.5)),mean(Alph1ErrV),quantile(Alph1ErrV,p=c(.75)),sqrt(var(Alph1ErrV))),3)
  ErrMat[7,]<-round(c(quantile(Alph2ErrV,p=c(.25,.5)),mean(Alph2ErrV),quantile(Alph2ErrV,p=c(.75)),sqrt(var(Alph2ErrV))),3)
  ErrMat[8,]<-round(c(quantile(Alph3ErrV,p=c(.25,.5)),mean(Alph3ErrV),quantile(Alph3ErrV,p=c(.75)),sqrt(var(Alph3ErrV))),3)
  colnames(ErrMat)<-c("Q(0.25)","Q(0.5)","Avg","Q(0.75)","Std")
  print(k)
  print(Sc)
  print(ErrMat)
  
  A1_mean = colSums(A1)/(mcreps-dd)
  A2_mean = colSums(A2)/(mcreps-dd)
  A3_mean = colSums(A3)/(mcreps-dd)
  
  A1low = na.omit(CI_asy_low[,1,])
  A2low =  na.omit(CI_asy_low[,2,])
  A3low =  na.omit(CI_asy_low[,3,])
  A1up = na.omit(CI_asy_up[,1,])
  A2up =  na.omit(CI_asy_up[,2,])
  A3up =  na.omit(CI_asy_up[,3,])
  
  A1low = na.omit(CI_boot_low[,1,])
  A2low =  na.omit(CI_boot_low[,2,])
  A3low =  na.omit(CI_boot_low[,3,])
  A1up = na.omit(CI_boot_up[,1,])
  A2up =  na.omit(CI_boot_up[,2,])
  A3up =  na.omit(CI_boot_up[,3,])
  
  A1l_mean = colSums(A1low)/(mcreps-dd)
  A2l_mean = colSums(A2low)/(mcreps-dd)
  A3l_mean = colSums(A3low)/(mcreps-dd)
  
  A1u_mean = colSums(A1up)/(mcreps-dd)
  A2u_mean = colSums(A2up)/(mcreps-dd)
  A3u_mean = colSums(A3up)/(mcreps-dd)
  
  
  
  par(mfrow=c(1,1))
  par(mgp=c(2.2,1,0))
  plot(xgrid[5:100],A1_mean[5:100],lwd=2,type = "l",col="black", xlab = "u",ylab=expression(A[1]),
       ylim=c(-5,15)
       ,cex.axis=1.7,cex.lab=1.7)
  lines(xgrid[5:100],smooth(A1l_mean[5:100]),lty="dashed",col="blue")
  lines(xgrid[5:100],smooth(A1u_mean[5:100]),lty="dashed",col="blue")
  polygon(c(xgrid[5:100],rev(xgrid[5:100])),c(CI_low[1,8,5:100],rev(CI_up[1,12,5:100])),col="grey95",lty=0,lwd=3)
  polygon(c(xgrid[5:100],rev(xgrid[5:100])),c(A1l_mean[5:100],rev(A1u_mean[5:100])),col="grey80",lty=0,lwd=3)
  lines(xgrid[5:100],A1_mean[5:100],lwd=2,type = "l",col="black")
  lines(xgrid[5:100],alpha1(xgrid[5:100]),lwd=2,cex.axis=1.2,cex.lab=1.2,col="red",lty="dashed")  
  legend(x=0.0005,y=14,cex=1.2,lwd=3,c(expression(tilde(A)[1]),expression(A[1])),lty=c(1,2),bty="n",col = c("black","red"))
  
    
  plot(xgrid[5:100],A2_mean[5:100],lwd=2,type = "l",col="black", xlab = "u",ylab=expression(A[2]),
       ylim=c(-5,15)
       ,cex.axis=1.7,cex.lab=1.7)
  lines(xgrid[5:100],smooth(A2l_mean[5:100]),lty="dashed",col="blue")
  lines(xgrid[5:100],smooth(A2u_mean[5:100]),lty="dashed",col="blue")
  lines(xgrid[5:100],smooth(CI_low[1,30,5:100]),lty="dashed",col="red")
  lines(xgrid[5:100],smooth(CI_up[1,30,5:100]),lty="dashed",col="red")
  polygon(c(xgrid[5:100],rev(xgrid[5:100])),c(CI_low[1,29,5:100],rev(CI_up[1,29,5:100])),col="grey95",lty=0,lwd=3)
  polygon(c(xgrid[5:100],rev(xgrid[5:100])),c(A2l_mean[5:100],rev(A2u_mean[5:100])),col="grey80",lty=0,lwd=3)
  lines(xgrid[5:100],A2_mean[5:100],lwd=2,type = "l",col="black")
  lines(xgrid[5:100],alpha2(xgrid[5:100]),lwd=2,cex.axis=1.2,cex.lab=1.2,col="red",lty="dashed")  
  legend(x=0.0005,y=14,cex=1.2,lwd=3,c(expression(tilde(A)[2]),expression(A[2])),lty=c(1,2),bty="n",col = c("black","red"))
  
  plot(xgrid[5:100],A3_mean[5:100],lwd=2,type = "l",col="black", xlab = "u",ylab=expression(A[3]),
       ylim=c(-6.5,15)
       ,cex.axis=1.7,cex.lab=1.7)
  lines(xgrid[5:100],smooth(A3l_mean[5:100]),lty="dashed",col="blue")
  lines(xgrid[5:100],smooth(A3u_mean[5:100]),lty="dashed",col="blue")
  polygon(c(xgrid[5:100],rev(xgrid[5:100])),c(CI_low[1,39,5:100],rev(CI_up[1,39,5:100])),col="grey95",lty=0,lwd=3)
  polygon(c(xgrid[5:100],rev(xgrid[5:100])),c(A3l_mean[5:100],rev(A3u_mean[5:100])),col="grey80",lty=0,lwd=3)
  lines(xgrid[5:100],A3_mean[5:100],lwd=2,type = "l",col="black")
  lines(xgrid[5:100],alpha3(xgrid[5:100]),lwd=2,cex.axis=1.2,cex.lab=1.2,col="red",lty="dashed")  
  legend(x=0.0005,y=14,cex=1.2,lwd=3,c(expression(tilde(A)[3]),expression(A[3])),lty=c(1,2),bty="n",col = c("black","red"))

}
