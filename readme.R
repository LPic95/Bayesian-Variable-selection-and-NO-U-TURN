#Project code extract referring to the fifth simulation scenario.
#Simulation design with 10 out of 50 non-zero predictors and orthogonal regressors by varing sample size from n>p to n<<p (high dimensional contexts)

install.packages("MASS")
install.packages("BoomSpikeSlab")
install.packages("statmod")
install.packages("glmnet")
install.packages("tictoc")
install.packages("ggpubr")
install.packages("spikeslab")
install.packages("ggrepel")
install.packages("SSLASSO")

library(MASS)
library(statmod)
library(glmnet)
library(tictoc)
library(BoomSpikeSlab)
library(ggpubr)
library(spikeslab)
library(stringr)
library(ggrepel)
library(SSLASSO)

getmode= function(v) {
  v=round(v,4)
  uniqv=unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

T=6000 #Gibbs Sampling interactions
Q=200 #lasso chains
burnin=800  
K=50 #number of predictors
set.seed(1234)
lasso_frequentist_matrix=matrix(0,Q,K)
beta_bayesian_matrix_mode=matrix(0,Q,K)
beta_bayesian_matrix_median=matrix(0,Q,K)
beta_bayesian_matrix_mean=matrix(0,Q,K)

truevslasso=data.frame(specificity=rep(0,Q),sensitivity=rep(0,Q),accuracy=rep(0,Q))
truevsbayesian_mean=data.frame(specificity=rep(0,Q),sensitivity=rep(0,Q),accuracy=rep(0,Q))
truevsbayesian_mode=data.frame(specificity=rep(0,Q),sensitivity=rep(0,Q),accuracy=rep(0,Q))
truevsbayesian_median=data.frame(specificity=rep(0,Q),sensitivity=rep(0,Q),accuracy=rep(0,Q))

obs=c(1000,100,60,50,30,20)

#True model:

Sigma=diag(1,K,K)
set.seed(1)
true_betas=matrix(rep(0,K),K,1)
index=sample(1:K,10,replace=F)
true_betas[index,]=matrix(rnorm(length(index),0,2),length(index),1)
X=list()
Ys=list()

for(nn in obs){
  passo=which(obs == nn)
  X[[passo]]=matrix(mvrnorm(nn,rep(0,K),Sigma),nn,K)
  Ys[[passo]]=X[[passo]] %*% true_betas+cbind(rnorm(nn,0,2))
}

expected_value_bayesian_mean=matrix(0,length(obs),K)
expected_value_bayesian_median=matrix(0,length(obs),K)
expected_value_bayesian_mode=matrix(0,length(obs),K)
expected_value_lasso=matrix(0,length(obs),K)


expected_truevsbayesian_median=data.frame(specificity=rep(0,length(obs)),sensitivity=rep(0,length(obs)),accuracy=rep(0,length(obs)))
expected_truevsbayesian_mean=data.frame(specificity=rep(0,length(obs)),sensitivity=rep(0,length(obs)),accuracy=rep(0,length(obs)))
expected_truevsbayesian_mode=data.frame(specificity=rep(0,length(obs)),sensitivity=rep(0,length(obs)),accuracy=rep(0,length(obs)))
expected_truevslasso=data.frame(specificity=rep(0,length(obs)),sensitivity=rep(0,length(obs)),accuracy=rep(0,length(obs)))

#Hamming distance
true_model_h=matrix(0,length(nns),K)
true_model_h[,abs(true_betas)>0]=1

bayesian_mean_h=c(rep(0,Q))
bayesian_mode_h=c(rep(0,Q))
bayesian_median_h=c(rep(0,Q))
frequentist_h=c(rep(0,Q))

expected_hamming_lasso_bayesian_mean=matrix(rep(0,length(nns)),length(nns),1)
expected_hamming_lasso_bayesian_mode=matrix(rep(0,length(nns)),length(nns),1)
expected_hamming_lasso_bayesian_median=matrix(rep(0,length(nns)),length(nns),1)
expected_hamming_lasso_frequentist=matrix(rep(0,length(nns)),length(nns),1)

lista_bayesian_mean=list()
lista_bayesian_mode=list()
lista_bayesian_median=list()
tic("DESIGN SIMULATION")

for(nn in obs){
  passo=which(obs == nn)
  Y_model=Ys[[passo]]
  betas_model=true_betas
  x=X[[passo]]
  for(r in 1:Q){
    
    ######glmnet lasso frequentist####
    
    fit.glmnet=glmnet(x,Y_model, 
                      lambda=cv.glmnet(x,Y_model,nfolds=5)$lambda.min,intercept=FALSE)
    lasso_frequentist_matrix[r,]=t(as.matrix(fit.glmnet$beta[,1]))
    
    ######GIBBS SAMPLING#############
    
    beta=matrix(0,T+1,K)
    lambda=rep(sqrt(1/2),T+1.75) #frazione dei due parametri della gamma prior.
    start=glmnet(x,Y_model,lambda=c(lambda[1]))
    beta[1,]=matrix(start$beta)[,1]
    sigma=rep(0.1,T+1)
    sigma[1]=0.1
    tau=matrix(1,T+1,K)
    for(i in 1:T){
      lambda[i+1]=sqrt(rgamma(1,shape=K+1,rate=sum(tau[i,])/2+1.75))
      for(j in 1:K){
        tau[i+1,j]=1/(rinvgauss(1,sqrt((lambda[i]^2)*sigma[i]/beta[i+1,j]),lambda[i]^2))
      }
      sigma[i+1]=1/rgamma(1,(nn+K-1)/2,1/(t(Y_model-x%*%beta[i,])%*%(Y_model-x%*%beta[i,])/2+t(beta[i,])%*%solve(diag(tau[i,]),tol = 1e-22)%*%(beta[i,]/2)))
      beta[i+1,]=mvrnorm(1,(solve(t(x)%*%x+solve(diag(tau[i,]),tol = 1e-22),tol = 1e-22))%*%t(x)%*%Y_model,sigma[i]*solve(t(x)%*%x+solve(diag(tau[i,]),tol = 1e-22),tol = 1e-22))
    }
    
    beta_bayesian=apply(beta[-(1:(burnin+1)),],2,getmode)
    beta_bayesian_median=apply(beta[-(1:(burnin+1)),],2,median)
    beta_bayesian_mean=apply(beta[-(1:(burnin+1)),],2,mean)
    beta_bayesian_matrix_mode[r,]=beta_bayesian
    beta_bayesian_matrix_median[r,]=beta_bayesian_median
    beta_bayesian_matrix_mean[r,]=beta_bayesian_mean
    
    #Sensitivity and Specificity
    
    coeff=data.frame(vero=betas_model,lasso=fit.glmnet$beta[,1],bayesian_mode=beta_bayesian,bayesian_median=beta_bayesian_median,bayesian_mean=beta_bayesian_mean)
    coeff=round(coeff,4)
    coeff$count.lasso=ifelse(coeff$lasso==0,0,1)
    coeff$count.true_model=ifelse(coeff$vero==0,0,1)
    coeff$count.bayesian_mode_lasso=ifelse(coeff$bayesian_mode==0,0,1)
    coeff$count.bayesian_mean_lasso=ifelse(coeff$bayesian_mean==0,0,1)
    coeff$count.bayesian_median_lasso=ifelse(coeff$bayesian_median==0,0,1)
    
    ####### true model  vs frequentist lasso #####
    o=table(factor(coeff$count.true_model,level=0:1),factor(coeff$count.lasso,level=0:1))
    names(dimnames(o)) <- c("true model","lasso")
    
    #Specificity
    
    truevslasso$specificity[r]=(o[2,2])/(o[2,2]+o[2,1])
    #Sensitivity
    truevslasso$sensitivity[r]=(o[1,1])/(o[1,1]+o[1,2])
    #Accuracy
    truevslasso$accuracy[r]=((o[1,1]+o[2,2])/sum(o))
    
    ###############################################
    ###### true model vs lasso Bayesian mode ######
    ###############################################
    
    lbmode=table(factor(coeff$count.true_model,level=0:1),factor(coeff$count.bayesian_mode_lasso,level=0:1))
    names(dimnames(lbmode)) <- c("true model","bayesian lasso")
    #Specificità
    truevsbayesian_mode$specificity[r]=(lbmode[2,2])/(lbmode[2,2]+lbmode[2,1])
    #Sensitività
    truevsbayesian_mode$sensitivity[r]=(lbmode[1,1])/(lbmode[1,1]+lbmode[1,2])
    #Accuratezza
    truevsbayesian_mode$accuracy[r]=((lbmode[1,1]+lbmode[2,2])/sum(lbmode))
    
    ##############################################
    ###### true model vs lasso Bayesian mean #####
    ##############################################
    lbmean=table(factor(coeff$count.true_model,level=0:1),factor(coeff$count.bayesian_mean_lasso,level=0:1))
    names(dimnames(lbmean)) <- c("true model","bayesian lasso")
    #Specificità
    truevsbayesian_mean$specificity[r]=(lbmean[2,2])/(lbmean[2,2]+lbmean[2,1])
    #Sensitività
    truevsbayesian_mean$sensitivity[r]=(lbmean[1,1])/(lbmean[1,1]+lbmean[1,2])
    #Accuratezza
    truevsbayesian_mean$accuracy[r]=((lbmean[1,1]+lbmean[2,2])/sum(lbmean))
    
    ############################################
    ###### true model vs lasso Bayesian median #
    ############################################
    lbmedian=table(factor(coeff$count.true_model,level=0:1),factor(coeff$count.bayesian_median_lasso,level=0:1))
    names(dimnames(lbmedian)) <- c("true model","bayesian lasso")
    #Specificità
    truevsbayesian_median$specificity[r]=(lbmedian[2,2])/(lbmedian[2,2]+lbmedian[2,1])
    #Sensitività
    truevsbayesian_median$sensitivity[r]=(lbmedian[1,1])/(lbmedian[1,1]+lbmedian[1,2])
    #Accuratezza
    truevsbayesian_median$accuracy[r]=((lbmedian[1,1]+lbmedian[2,2])/sum(lbmedian))
    
    ##########################################
    ###### Hamming distance ##################
    ##########################################
    
    bayesian_mean_h[r]=sum(abs(true_model_h[passo,]-coeff$count.bayesian_mean_lasso)) 
    bayesian_mode_h[r]=sum(abs(true_model_h[passo,]-coeff$count.bayesian_mode_lasso)) 
    bayesian_median_h[r]=sum(abs(true_model_h[passo,]-coeff$count.bayesian_median_lasso)) 
    frequentist_h[r]=sum(abs(true_model_h[passo,]-coeff$count.lasso)) 
  }
  
  lista_bayesian_mean[[passo]]=beta_bayesian_matrix_mean
  lista_bayesian_mode[[passo]]=beta_bayesian_matrix_mode
  lista_bayesian_median[[passo]]=beta_bayesian_matrix_median
  
  expected_hamming_lasso_bayesian_mean[passo,]=mean(bayesian_mean_h)
  expected_hamming_lasso_bayesian_mode[passo,]=mean(bayesian_mode_h)
  expected_hamming_lasso_bayesian_median[passo,]=mean(bayesian_median_h)
  expected_hamming_lasso_frequentist[passo,]=mean(frequentist_h)
  
  #Expected specificity/sensitivity/accuracy
  
  expected_truevsbayesian_median[passo,]=apply(truevsbayesian_median,2,mean)
  expected_truevsbayesian_mean[passo,]=apply(truevsbayesian_mean,2,mean)
  expected_truevsbayesian_mode[passo,]=apply(truevsbayesian_mode,2,mean)
  expected_truevslasso[passo,]=apply(truevslasso,2,mean)
  
  #Expected values
  expected_value_bayesian_mean[passo,]=round(apply(beta_bayesian_matrix_mean,2,mean),4)
  expected_value_bayesian_median[passo,]=round(apply(beta_bayesian_matrix_median,2,mean),4)
  expected_value_bayesian_mode[passo,]=round(apply(beta_bayesian_matrix_mode,2,mean),4)
  expected_value_lasso[passo,]=round(apply(lasso_frequentist_matrix,2,mean),4)
}
toc() 


##################################################OTHER METHODOLOGIES#######################################################

##########SPIKE AND SLAB  #########################

########George and McCulloch

ss_George_McCulloch_oracle=list()#Oracle
ss_George_McCulloch=list()
niter=10000
truevsSS_George_McCulloch_oracle=data.frame(specificity=rep(0,length(obs)),sensitivity=rep(0,length(obs)),accuracy=rep(0,length(obs)))
truevsSS_George_McCulloch=data.frame(specificity=rep(0,length(obs)),sensitivity=rep(0,length(obs)),accuracy=rep(0,length(obs)))
SS_hamming=matrix(rep(0,length(obs)),length(obs),1)
set.seed(1)
R=30
seed=sample(1:10000000,R,rep=F)
expected_hamming_SS=matrix(0,length(obs),R)
expected_sensitivity_SS=matrix(0,length(obs),R)
expected_specificity_SS=matrix(0,length(obs),R)
expected_accuracy_SS=matrix(0,length(obs),R)
SS_mean_coeff=matrix(0,length(obs),K+1)
expected_mean_coeff=matrix(0,R*length(obs),K+1)
dim(Ys[[1]])

for(r in 1:R){
  for(nn in obs){
    passo=which(obs == nn)
    prior=SpikeSlabPrior(X[[passo]], Ys[[passo]])
    ss_George_McCulloch[[passo]]=lm.spike(Ys[[passo]] ~ X[[passo]] -1, niter=niter,prior=prior, error.distribution = c("gaussian"), ping=1000,seed=seed)
    coefs=data.frame(vero=true_betas,ss_GM=apply(ss_George_McCulloch[[passo]]$beta,2,mean))
    coefs=round(coefs,4)
    coefs$count.true_model=ifelse(coefs$vero==0,0,1)
    coefs$count.SS=ifelse(coefs$ss_GM==0,0,1)
    
    ss=table(factor(coefs$count.true_model,level=0:1),factor(coefs$count.SS,level=0:1))
    names(dimnames(ss)) <- c("true model","SS George McCulloch")
    #Specificity
    truevsSS_George_McCulloch$specificity[passo]=(ss[2,2])/(ss[2,2]+ss[2,1])
    #Sensitivity
    truevsSS_George_McCulloch$sensitivity[passo]=(ss[1,1])/(ss[1,1]+ss[1,2])
    #Accuracy
    truevsSS_George_McCulloch$accuracy[passo]=((ss[1,1]+ss[2,2])/sum(ss))
    
    SS_hamming[passo,]=sum(abs(true_model_h[passo,]-coefs$count.SS)) 
    SS_mean_coeff[passo,]=c(apply(ss_George_McCulloch[[passo]]$beta[-c(1:500),],2,mean),passo=passo)
  }
  expected_hamming_SS[,r]=SS_hamming
  expected_sensitivity_SS[,r]=truevsSS_George_McCulloch$sensitivity
  expected_specificity_SS[,r]=truevsSS_George_McCulloch$specificity
  expected_accuracy_SS[,r]=truevsSS_George_McCulloch$accuracy
  expected_mean_coeff[c(r+length(obs)*(r-1)-(r-1)):c(r*length(obs)),]=SS_mean_coeff
}

#################################


expected_truevsSS_George_McCulloch=data.frame(specificity=apply(expected_specificity_SS,1,mean),sensitivity=apply(expected_sensitivity_SS,1,mean),accuracy=apply(expected_accuracy_SS,1,mean))
expected_truevsSS_George_McCulloch

expected_hamming_SS=apply(expected_hamming_SS,1,mean)
expected_hamming_SS

expected_SS_George_McCulloch=matrix(0,length(obs),K)

expected_SS_George_McCulloch[1,]=apply(expected_mean_coeff[which(expected_mean_coeff[,51]==1),c(1:K)],2,mean)
expected_SS_George_McCulloch[2,]=apply(expected_mean_coeff[which(expected_mean_coeff[,51]==2),c(1:K)],2,mean)
expected_SS_George_McCulloch[3,]=apply(expected_mean_coeff[which(expected_mean_coeff[,51]==3),c(1:K)],2,mean)
expected_SS_George_McCulloch[4,]=apply(expected_mean_coeff[which(expected_mean_coeff[,51]==4),c(1:K)],2,mean)
expected_SS_George_McCulloch[5,]=apply(expected_mean_coeff[which(expected_mean_coeff[,51]==5),c(1:K)],2,mean)
expected_SS_George_McCulloch[6,]=apply(expected_mean_coeff[which(expected_mean_coeff[,51]==6),c(1:K)],2,mean)
expected_SS_George_McCulloch=round(expected_SS_George_McCulloch,5)
expected_SS_George_McCulloch



##1000 obs

beta_model=data.frame(ss_George_McCulloch[[1]]$beta)
beta_model$count=c(1:dim(beta_model)[1])
true_betas[which(abs(true_betas)>0)][1]

b1=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[1]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[3]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][1],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[1]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][1]), ymax = max(beta_model[,which(abs(true_betas)>0)][1]), fill = "#FFCC99", alpha = 0.4)

b2=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[2]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[10]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][2],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[2]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][2]), ymax = max(beta_model[,which(abs(true_betas)>0)][2]), fill = "#FFCC99", alpha = 0.4)

b3=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[3]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[14]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][3],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[3]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][3]), ymax = max(beta_model[,which(abs(true_betas)>0)][3]), fill = "#FFCC99", alpha = 0.4)

b4=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[4]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[19]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][4],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[4]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][4]), ymax = max(beta_model[,which(abs(true_betas)>0)][4]), fill = "#FFCC99", alpha = 0.4)

b5=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[5]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[27]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][5],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[5]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][5]), ymax = max(beta_model[,which(abs(true_betas)>0)][5]), fill = "#FFCC99", alpha = 0.4)


b6=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[6]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[28]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][6],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[6]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][6]), ymax = max(beta_model[,which(abs(true_betas)>0)][6]), fill = "#FFCC99", alpha = 0.4)

b7=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[7]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[29]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][7],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[7]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][7]), ymax = max(beta_model[,which(abs(true_betas)>0)][7]), fill = "#FFCC99", alpha = 0.4)

b8=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[8]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[41]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][8],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[8]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][8]), ymax = max(beta_model[,which(abs(true_betas)>0)][8]), fill = "#FFCC99", alpha = 0.4)

b9=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[9]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[42]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][9],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[9]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][9]), ymax = max(beta_model[,which(abs(true_betas)>0)][9]), fill = "#FFCC99", alpha = 0.4)

b10=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[10]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[43]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][10],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[10]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][10]), ymax = max(beta_model[,which(abs(true_betas)>0)][10]), fill = "#FFCC99", alpha = 0.4)

library(ggpubr)
figure=ggarrange(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,
                 ncol = 2, nrow = 5)
annotate_figure(figure,top = text_grob("Trace plots of Gibbs Sampling draws", color = "Black", face = "bold", size = 14))



#100 obs

beta_model=data.frame(ss_George_McCulloch[[2]]$beta)
beta_model$count=c(1:dim(beta_model)[1])
true_betas[which(abs(true_betas)>0)][1]

b1=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[1]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[3]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][1],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[1]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][1]), ymax = max(beta_model[,which(abs(true_betas)>0)][1]), fill = "#FFCC99", alpha = 0.4)

b2=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[2]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[10]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][2],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[2]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][2]), ymax = max(beta_model[,which(abs(true_betas)>0)][2]), fill = "#FFCC99", alpha = 0.4)

b3=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[3]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[14]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][3],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[3]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][3]), ymax = max(beta_model[,which(abs(true_betas)>0)][3]), fill = "#FFCC99", alpha = 0.4)

b4=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[4]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[19]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][4],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[4]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][4]), ymax = max(beta_model[,which(abs(true_betas)>0)][4]), fill = "#FFCC99", alpha = 0.4)

b5=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[5]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[27]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][5],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[5]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][5]), ymax = max(beta_model[,which(abs(true_betas)>0)][5]), fill = "#FFCC99", alpha = 0.4)

b6=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[6]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[28]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][6],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[6]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][6]), ymax = max(beta_model[,which(abs(true_betas)>0)][6]), fill = "#FFCC99", alpha = 0.4)

b7=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[7]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[29]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][7],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[7]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][7]), ymax = max(beta_model[,which(abs(true_betas)>0)][7]), fill = "#FFCC99", alpha = 0.4)

b8=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[8]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[41]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][8],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[8]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][8]), ymax = max(beta_model[,which(abs(true_betas)>0)][8]), fill = "#FFCC99", alpha = 0.4)

b9=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[9]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[42]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][9],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[9]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][9]), ymax = max(beta_model[,which(abs(true_betas)>0)][9]), fill = "#FFCC99", alpha = 0.4)

b10=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[10]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[43]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][10],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[10]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][10]), ymax = max(beta_model[,which(abs(true_betas)>0)][10]), fill = "#FFCC99", alpha = 0.4)

figure=ggarrange(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,
                 ncol = 2, nrow = 5)
annotate_figure(figure,top = text_grob("Trace plots of Gibbs Sampling draws", color = "Black", face = "bold", size = 14))


#60 obs

beta_model=data.frame(ss_George_McCulloch[[3]]$beta)
beta_model$count=c(1:dim(beta_model)[1])
true_betas[which(abs(true_betas)>0)][1]

b1=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[1]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[3]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][1],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[1]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][1]), ymax = max(beta_model[,which(abs(true_betas)>0)][1]), fill = "#FFCC99", alpha = 0.4)

b2=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[2]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[10]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][2],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[2]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][2]), ymax = max(beta_model[,which(abs(true_betas)>0)][2]), fill = "#FFCC99", alpha = 0.4)

b3=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[3]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[14]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][3],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[3]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][3]), ymax = max(beta_model[,which(abs(true_betas)>0)][3]), fill = "#FFCC99", alpha = 0.4)

b4=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[4]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[19]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][4],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[4]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][4]), ymax = max(beta_model[,which(abs(true_betas)>0)][4]), fill = "#FFCC99", alpha = 0.4)

b5=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[5]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[27]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][5],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[5]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][5]), ymax = max(beta_model[,which(abs(true_betas)>0)][5]), fill = "#FFCC99", alpha = 0.4)

b6=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[6]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[28]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][6],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[6]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][6]), ymax = max(beta_model[,which(abs(true_betas)>0)][6]), fill = "#FFCC99", alpha = 0.4)

b7=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[7]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[29]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][7],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[7]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][7]), ymax = max(beta_model[,which(abs(true_betas)>0)][7]), fill = "#FFCC99", alpha = 0.4)

b8=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[8]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[41]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][8],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[8]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][8]), ymax = max(beta_model[,which(abs(true_betas)>0)][8]), fill = "#FFCC99", alpha = 0.4)

b9=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[9]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[42]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][9],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[9]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][9]), ymax = max(beta_model[,which(abs(true_betas)>0)][9]), fill = "#FFCC99", alpha = 0.4)

b10=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[10]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[43]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][10],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[10]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][10]), ymax = max(beta_model[,which(abs(true_betas)>0)][10]), fill = "#FFCC99", alpha = 0.4)

library(ggpubr)
figure=ggarrange(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,
                 ncol = 2, nrow = 5)
annotate_figure(figure,top = text_grob("Trace plots of Gibbs Sampling draws", color = "Black", face = "bold", size = 14))



#50 obs

beta_model=data.frame(ss_George_McCulloch[[4]]$beta)
beta_model$count=c(1:dim(beta_model)[1])
true_betas[which(abs(true_betas)>0)][1]

b1=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[1]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[3]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][1],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[1]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][1]), ymax = max(beta_model[,which(abs(true_betas)>0)][1]), fill = "#FFCC99", alpha = 0.4)

b2=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[2]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[10]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][2],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[2]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][2]), ymax = max(beta_model[,which(abs(true_betas)>0)][2]), fill = "#FFCC99", alpha = 0.4)

b3=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[3]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[14]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][3],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[3]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][3]), ymax = max(beta_model[,which(abs(true_betas)>0)][3]), fill = "#FFCC99", alpha = 0.4)

b4=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[4]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[19]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][4],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[4]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][4]), ymax = max(beta_model[,which(abs(true_betas)>0)][4]), fill = "#FFCC99", alpha = 0.4)

b5=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[5]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[27]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][5],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[5]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][5]), ymax = max(beta_model[,which(abs(true_betas)>0)][5]), fill = "#FFCC99", alpha = 0.4)

b6=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[6]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[28]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][6],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[6]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][6]), ymax = max(beta_model[,which(abs(true_betas)>0)][6]), fill = "#FFCC99", alpha = 0.4)

b7=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[7]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[29]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][7],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[7]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][7]), ymax = max(beta_model[,which(abs(true_betas)>0)][7]), fill = "#FFCC99", alpha = 0.4)

b8=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[8]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[41]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][8],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[8]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][8]), ymax = max(beta_model[,which(abs(true_betas)>0)][8]), fill = "#FFCC99", alpha = 0.4)

b9=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[9]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[42]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][9],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[9]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][9]), ymax = max(beta_model[,which(abs(true_betas)>0)][9]), fill = "#FFCC99", alpha = 0.4)

b10=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[10]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[43]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][10],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[10]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][10]), ymax = max(beta_model[,which(abs(true_betas)>0)][10]), fill = "#FFCC99", alpha = 0.4)

library(ggpubr)
figure=ggarrange(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,
                 ncol = 2, nrow = 5)
annotate_figure(figure,top = text_grob("Trace plots of Gibbs Sampling draws", color = "Black", face = "bold", size = 14))


#30 obs

beta_model=data.frame(ss_George_McCulloch[[5]]$beta)
beta_model$count=c(1:dim(beta_model)[1])
true_betas[which(abs(true_betas)>0)][1]

b1=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[1]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[3]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][1],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[1]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][1]), ymax = max(beta_model[,which(abs(true_betas)>0)][1]), fill = "#FFCC99", alpha = 0.4)

b2=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[2]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[10]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][2],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[2]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][2]), ymax = max(beta_model[,which(abs(true_betas)>0)][2]), fill = "#FFCC99", alpha = 0.4)

b3=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[3]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[14]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][3],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[3]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][3]), ymax = max(beta_model[,which(abs(true_betas)>0)][3]), fill = "#FFCC99", alpha = 0.4)

b4=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[4]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[19]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][4],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[4]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][4]), ymax = max(beta_model[,which(abs(true_betas)>0)][4]), fill = "#FFCC99", alpha = 0.4)

b5=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[5]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[27]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][5],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[5]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][5]), ymax = max(beta_model[,which(abs(true_betas)>0)][5]), fill = "#FFCC99", alpha = 0.4)

b6=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[6]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[28]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][6],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[6]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][6]), ymax = max(beta_model[,which(abs(true_betas)>0)][6]), fill = "#FFCC99", alpha = 0.4)

b7=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[7]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[29]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][7],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[7]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][7]), ymax = max(beta_model[,which(abs(true_betas)>0)][7]), fill = "#FFCC99", alpha = 0.4)

b8=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[8]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[41]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][8],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[8]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][8]), ymax = max(beta_model[,which(abs(true_betas)>0)][8]), fill = "#FFCC99", alpha = 0.4)

b9=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[9]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[42]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][9],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[9]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][9]), ymax = max(beta_model[,which(abs(true_betas)>0)][9]), fill = "#FFCC99", alpha = 0.4)

b10=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[10]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[43]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][10],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[10]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][10]), ymax = max(beta_model[,which(abs(true_betas)>0)][10]), fill = "#FFCC99", alpha = 0.4)

library(ggpubr)
figure=ggarrange(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,
                 ncol = 2, nrow = 5)
annotate_figure(figure,top = text_grob("Trace plots of Gibbs Sampling draws", color = "Black", face = "bold", size = 14))


#20 obs

beta_model=data.frame(ss_George_McCulloch[[6]]$beta)
beta_model$count=c(1:dim(beta_model)[1])
true_betas[which(abs(true_betas)>0)][1]

b1=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[1]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[3]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][1],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[1]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][1]), ymax = max(beta_model[,which(abs(true_betas)>0)][1]), fill = "#FFCC99", alpha = 0.4)

b2=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[2]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[10]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][2],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[2]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][2]), ymax = max(beta_model[,which(abs(true_betas)>0)][2]), fill = "#FFCC99", alpha = 0.4)

b3=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[3]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[14]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][3],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[3]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][3]), ymax = max(beta_model[,which(abs(true_betas)>0)][3]), fill = "#FFCC99", alpha = 0.4)

b4=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[4]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[19]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][4],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[4]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][4]), ymax = max(beta_model[,which(abs(true_betas)>0)][4]), fill = "#FFCC99", alpha = 0.4)

b5=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[5]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[27]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][5],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[5]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][5]), ymax = max(beta_model[,which(abs(true_betas)>0)][5]), fill = "#FFCC99", alpha = 0.4)

b6=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[6]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[28]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][6],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[6]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][6]), ymax = max(beta_model[,which(abs(true_betas)>0)][6]), fill = "#FFCC99", alpha = 0.4)

b7=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[7]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[29]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][7],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[7]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][7]), ymax = max(beta_model[,which(abs(true_betas)>0)][7]), fill = "#FFCC99", alpha = 0.4)

b8=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[8]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[41]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][8],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[8]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][8]), ymax = max(beta_model[,which(abs(true_betas)>0)][8]), fill = "#FFCC99", alpha = 0.4)

b9=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[9]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[42]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][9],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[9]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][9]), ymax = max(beta_model[,which(abs(true_betas)>0)][9]), fill = "#FFCC99", alpha = 0.4)

b10=ggplot(beta_model,aes(x=beta_model$count,y=beta_model[,which(abs(true_betas)>0)[10]]))+geom_line(size=0.2)+xlab("Interaction")+ylab("Beta[43]")+theme_bw()+scale_x_continuous(limits=c(-1, 10000))+geom_hline(yintercept=true_betas[which(abs(true_betas)>0)][10],linetype="dashed", color = "red")+
  geom_hline(yintercept=mean(beta_model[-c(1:500),which(abs(true_betas)>0)[10]]),linetype="dashed", color = "orange")+annotate("rect", xmin = 0, xmax = 500, ymin = min(beta_model[,which(abs(true_betas)>0)][10]), ymax = max(beta_model[,which(abs(true_betas)>0)][10]), fill = "#FFCC99", alpha = 0.4)

library(ggpubr)
figure=ggarrange(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,
                 ncol = 2, nrow = 5)
annotate_figure(figure,top = text_grob("Trace plots of Gibbs Sampling draws", color = "Black", face = "bold", size = 14))



####SPIKE AND SLAB: Ishwaran and Rao  ###################


ss_Ishwaran_Rao=list()
stability=list()
stb=list()
regexp="[[:digit:]]+"#serve per estrarre dal vettore stringa solo i numeri
truevsSS_Ishwaran_Rao=data.frame(specificity=rep(0,length(obs)),sensitivity=rep(0,length(obs)),accuracy=rep(0,length(obs)))
SS_hamming=matrix(0,length(obs),1)
SS_mean_coeff_Ishwaran_Rao=matrix(0,length(obs),K+1)
expected_hamming_SS_Ishwaran_Rao=matrix(0,length(obs),R)
expected_sensitivity_SS_Ishwaran_Rao=matrix(0,length(obs),R)
expected_specificity_SS_Ishwaran_Rao=matrix(0,length(obs),R)
expected_accuracy_SS_Ishwaran_Rao=matrix(0,length(obs),R)
expected_mean_coeff_SS_Ishwaran_Rao=matrix(0,R*length(obs),K+1)


for(r in 1:R){
  for(nn in obs){
    passo=which(obs == nn)
    if(nn<50){
      ss_Ishwaran_Rao[[passo]]=spikeslab(y=Ys[[passo]],x= X[[passo]], n.iter1=500,n.iter2=1000,bigp.smalln=T,screen=F,verbose=TRUE,intercept=F)
      cv.obj=cv.spikeslab(x = X[[passo]], y = Ys[[passo]], K = 10,bigp.smalln=T,plot.it=F)
    } else{
      ss_Ishwaran_Rao[[passo]]=spikeslab(y=Ys[[passo]],x= X[[passo]], n.iter1=500,n.iter2=1000,bigp.smalln=F,screen=F,verbose=TRUE,intercept=F)
      cv.obj=cv.spikeslab(x = X[[passo]], y = Ys[[passo]], K = 10,plot.it=F)
    }
    cv.stb=as.data.frame(cv.obj$stability) 
    cv.stb$coeff=as.numeric(str_extract(rownames(cv.stb), regexp))
    cv.stb=cv.stb[order(cv.stb$coeff),]
    stb[[passo]]=cv.stb
    ss_Ishwaran_Rao[[passo]]$bma.scale[which(cv.stb$stability<90)]=0
    stability[[passo]]=cv.stb$stability
    coefs=data.frame(vero=true_betas,SS=ss_Ishwaran_Rao[[passo]]$bma.scale)
    coefs=round(coefs,4)
    coefs$count.true_model=ifelse(coefs$vero==0,0,1)
    coefs$count.SS=ifelse(coefs$SS==0,0,1)
    ss=table(factor(coefs$count.true_model,level=0:1),factor(coefs$count.SS,level=0:1))
    
    names(dimnames(ss)) <- c("true model","ss Ishwaran Rao")
    #Specificità
    truevsSS_Ishwaran_Rao$specificity[passo]=(ss[2,2])/(ss[2,2]+ss[2,1])
    #Sensitività
    truevsSS_Ishwaran_Rao$sensitivity[passo]=(ss[1,1])/(ss[1,1]+ss[1,2])
    #Accuratezza
    truevsSS_Ishwaran_Rao$accuracy[passo]=((ss[1,1]+ss[2,2])/sum(ss))
    SS_hamming[passo,]=sum(abs(coefs$count.true_model-coefs$count.SS))
    SS_mean_coeff_Ishwaran_Rao[passo,]=c(ss_Ishwaran_Rao[[passo]]$bma.scale,passo=passo)
  }
  expected_hamming_SS_Ishwaran_Rao[,r]=SS_hamming
  expected_sensitivity_SS_Ishwaran_Rao[,r]=truevsSS_Ishwaran_Rao$sensitivity
  expected_specificity_SS_Ishwaran_Rao[,r]=truevsSS_Ishwaran_Rao$specificity
  expected_accuracy_SS_Ishwaran_Rao[,r]=truevsSS_Ishwaran_Rao$accuracy
  expected_mean_coeff_SS_Ishwaran_Rao[c(r+length(obs)*(r-1)-(r-1)):c(r*length(obs)),]=SS_mean_coeff_Ishwaran_Rao
}

expected_truevsSS_Ishwaran_Rao=data.frame(specificity=apply(expected_specificity_SS_Ishwaran_Rao,1,mean),sensitivity=apply(expected_sensitivity_SS_Ishwaran_Rao,1,mean),accuracy=apply(expected_accuracy_SS_Ishwaran_Rao,1,mean))
expected_truevsSS_Ishwaran_Rao

expected_hamming_SS_Ishwaran=apply(expected_hamming_SS_Ishwaran_Rao,1,mean)
expected_hamming_SS_Ishwaran

expected_SS_Ishwaran_Rao=matrix(0,length(obs),K)

expected_SS_Ishwaran_Rao[1,]=apply(expected_mean_coeff_SS_Ishwaran_Rao[which(expected_mean_coeff_SS_Ishwaran_Rao[,51]==1),c(1:K)],2,mean)
expected_SS_Ishwaran_Rao[2,]=apply(expected_mean_coeff_SS_Ishwaran_Rao[which(expected_mean_coeff_SS_Ishwaran_Rao[,51]==2),c(1:K)],2,mean)
expected_SS_Ishwaran_Rao[3,]=apply(expected_mean_coeff_SS_Ishwaran_Rao[which(expected_mean_coeff_SS_Ishwaran_Rao[,51]==3),c(1:K)],2,mean)
expected_SS_Ishwaran_Rao[4,]=apply(expected_mean_coeff_SS_Ishwaran_Rao[which(expected_mean_coeff_SS_Ishwaran_Rao[,51]==4),c(1:K)],2,mean)
expected_SS_Ishwaran_Rao[5,]=apply(expected_mean_coeff_SS_Ishwaran_Rao[which(expected_mean_coeff_SS_Ishwaran_Rao[,51]==5),c(1:K)],2,mean)
expected_SS_Ishwaran_Rao[6,]=apply(expected_mean_coeff_SS_Ishwaran_Rao[which(expected_mean_coeff_SS_Ishwaran_Rao[,51]==6),c(1:K)],2,mean)
expected_SS_Ishwaran_Rao=round(expected_SS_Ishwaran_Rao,5)
expected_SS_Ishwaran_Rao


############### STABILITY ASSESSMENT ##############

library(tidyverse)


#1000 obs
s=ggplot(data=stb[[1]],aes(x=stb[[1]]$bma,y=stb[[1]]$stability))
s+geom_point(aes(color = stb[[1]]$stability))+xlab("Bma")+ylab("Stability (%)")+ 
  labs(color="Stability (%)") +theme_bw()+geom_text_repel(aes(label=rownames(stb[[1]]) %>% str_replace("x.","X")), size = 2,segment.color="tan2",segment.size = 0.25,segment.alpha =0.4)+
  ggtitle("Stability analysis.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 1000 observations.")+
  theme( panel.background = element_blank(),plot.title = element_text(size = 14, face = "bold", hjust = 0.5),axis.text.x=element_text(colour="black", size = 9), axis.text.y=element_text(colour="black", size = 9),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))



#100 obs
s=ggplot(data=stb[[2]],aes(x=stb[[2]]$bma,y=stb[[2]]$stability))
s+geom_point(aes(color = stb[[2]]$stability))+xlab("Bma")+ylab("Stability (%)")+ 
  labs(color="Stability (%)") +theme_bw()+geom_text_repel(aes(label=rownames(stb[[2]]) %>% str_replace("x.","X")), size = 2,segment.color="tan2",segment.size = 0.25,segment.alpha =0.4)+
  ggtitle("Stability analysis.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 100 observations.")+theme( panel.background = element_blank(),plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
                                                                                                                                                              axis.text.x=element_text(colour="black", size = 9), axis.text.y=element_text(colour="black", size = 9),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))


#60 obs
s=ggplot(data=stb[[3]],aes(x=stb[[3]]$bma,y=stb[[3]]$stability))
s+geom_point(aes(color = stb[[3]]$stability))+xlab("Bma")+ylab("Stability (%)")+ 
  labs(color="Stability (%)") +theme_bw()+geom_text_repel(aes(label=rownames(stb[[3]]) %>% str_replace("x.","X")), size = 2,segment.color="tan2",segment.size = 0.25,segment.alpha =0.4)+
  ggtitle("Stability analysis.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 60 observations.")+theme( panel.background = element_blank(),plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
                                                                                                                                                              axis.text.x=element_text(colour="black", size = 9), axis.text.y=element_text(colour="black", size = 9),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))

#50 obs
s=ggplot(data=stb[[4]],aes(x=stb[[4]]$bma,y=stb[[4]]$stability))
s+geom_point(aes(color = stb[[4]]$stability))+xlab("Bma")+ylab("Stability (%)")+ 
  labs(color="Stability (%)") +theme_bw()+geom_text_repel(aes(label=rownames(stb[[4]]) %>% str_replace("x.","X")), size = 2,segment.color="tan2",segment.size = 0.25,segment.alpha =0.4)+
  ggtitle("Stability analysis.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 50 observations.")+theme( panel.background = element_blank(),plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
                                                                                                                                                              axis.text.x=element_text(colour="black", size = 9), axis.text.y=element_text(colour="black", size = 9),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))

#30 obs
s=ggplot(data=stb[[5]],aes(x=stb[[5]]$bma,y=stb[[5]]$stability))
s+geom_point(aes(color = stb[[5]]$stability))+xlab("Bma")+ylab("Stability (%)")+ 
  labs(color="Stability (%)") +theme_bw()+geom_text_repel(aes(label=rownames(stb[[5]]) %>% str_replace("x.","X")), size = 2,segment.color="tan2",segment.size = 0.25,segment.alpha =0.4)+
  ggtitle("Stability analysis.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 30 observations.")+theme( panel.background = element_blank(),plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
                                                                                                                                                              axis.text.x=element_text(colour="black", size = 9), axis.text.y=element_text(colour="black", size = 9),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))

#20 obs
s=ggplot(data=stb[[6]],aes(x=stb[[6]]$bma,y=stb[[6]]$stability))
s+geom_point(aes(color = stb[[6]]$stability))+xlab("Bma")+ylab("Stability (%)")+ 
  labs(color="Stability (%)") +theme_bw()+geom_text_repel(aes(label=rownames(stb[[6]]) %>% str_replace("x.","X")), size = 2,segment.color="tan2",segment.size = 0.25,segment.alpha =0.4)+
  ggtitle("Stability analysis.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 20 observations.")+theme( panel.background = element_blank(),plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
                                                                                                                                                              axis.text.x=element_text(colour="black", size = 9), axis.text.y=element_text(colour="black", size = 9),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))




##########SPIKE AND SLAB LASSO #########################

ssl_lasso_known_variance_separable=list()#Oracle
ssl_lasso_known_variance_notseparable=list()
ssl_lasso_unknown_variance=list()

truevsSSL_oracle=data.frame(specificity=rep(0,length(obs)),sensitivity=rep(0,length(obs)),accuracy=rep(0,length(obs)))
truevsSSL_unknown_variance=data.frame(specificity=rep(0,length(obs)),sensitivity=rep(0,length(obs)),accuracy=rep(0,length(obs)))
SSL_oracle_hamming=matrix(rep(0,length(obs)),length(obs),1)
SSL_hamming=matrix(rep(0,length(obs)),length(obs),1)

obs
d=seq(1,50,1)
lambda0=c(1 + d)
lambda0[1]=1.5
lambda0[2]=2.3
lambda0[3]=5
lambda0[4]=5.9
lambda0[5]=6.5
lambda0[6]=7.8

for(nn in obs){
  passo=which(obs == nn)
  ssl_lasso_known_variance_separable[[passo]]=SSLASSO(X[[passo]], Ys[[passo]], penalty = "separable", variance = "fixed",lambda1=1,lambda0 = lambda0, theta = 10/50,sigma=2,eps = 0.001, warn = T)
  ssl_lasso_known_variance_notseparable[[passo]]=SSLASSO(X[[passo]], Ys[[passo]], penalty = "separable", variance = "fixed",lambda1=1,lambda0 = lambda0,a=1,b=K,sigma=2,eps = 0.001, warn = T) 
  ssl_lasso_unknown_variance[[passo]]=SSLASSO(X[[passo]], Ys[[passo]], penalty = "adaptive", variance = "unknown",lambda1=1,lambda0 = lambda0,eps = 0.001,a=1,b=K, warn = T)
  
  coefs=data.frame(vero=true_betas,oracle=apply(ssl_lasso_known_variance_separable[[passo]]$beta,1,getmode),SSL=apply(ssl_lasso_unknown_variance[[passo]]$beta,1,getmode))
  coefs=round(coefs,4)
  coefs$count.true_model=ifelse(coefs$vero==0,0,1)
  coefs$count.oracle=ifelse(coefs$oracle==0,0,1)
  coefs$count.SSL=ifelse(coefs$SSL==0,0,1)
  
  ssl_or=table(factor(coefs$count.true_model,level=0:1),factor(coefs$count.oracle,level=0:1))
  names(dimnames(ssl_or)) <- c("true model","ssl oracle")
  #Specificità
  truevsSSL_oracle$specificity[passo]=(ssl_or[2,2])/(ssl_or[2,2]+ssl_or[2,1])
  #Sensitività
  truevsSSL_oracle$sensitivity[passo]=(ssl_or[1,1])/(ssl_or[1,1]+ssl_or[1,2])
  #Accuratezza
  truevsSSL_oracle$accuracy[passo]=((ssl_or[1,1]+ssl_or[2,2])/sum(ssl_or))
  
  ssl_uvariance=table(factor(coefs$count.true_model,level=0:1),factor(coefs$count.SSL,level=0:1))
  names(dimnames(ssl_uvariance)) <- c("true model","ssl")
  #Specificità
  truevsSSL_unknown_variance$specificity[passo]=(ssl_uvariance[2,2])/(ssl_uvariance[2,2]+ssl_uvariance[2,1])
  #Sensitività
  truevsSSL_unknown_variance$sensitivity[passo]=(ssl_uvariance[1,1])/(ssl_uvariance[1,1]+ssl_uvariance[1,2])
  #Accuratezza
  truevsSSL_unknown_variance$accuracy[passo]=((ssl_uvariance[1,1]+ssl_uvariance[2,2])/sum(ssl_uvariance))
  SSL_oracle_hamming[passo,]=sum(abs(true_model_h[passo,]-coefs$count.oracle.SS)) 
  SSL_hamming[passo,]=sum(abs(true_model_h[passo,]-coefs$count.SSL)) 
}


#Unknown variance SSL, 10 not null coefficients. 1000 obs

lambda0=rep(ssl_lasso_unknown_variance[[1]]$lambda0, each=50)
coeff=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
        "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37","X38","X39",
        "X40","X41","X42","X43","X44","X45","X46","X47","X48","X49","X50"
)
coefficient=rep(coeff,length(ssl_lasso_unknown_variance[[1]]$lambda0))
length(lambda0)==length(coefficient)
value=c(ssl_lasso_unknown_variance[[1]]$beta)
length(value)
true_value=c(rep(0,50))#utile per forma linea se 0 o differente da zero nel modello originale.
true_value[which(true_betas!=0)]=1


coeff_sslasso=data.frame(coefficient,value,lambda0,true=rep(true_value,length(ssl_lasso_unknown_variance[[1]]$lambda0)))
Coefficients=factor(coeff_sslasso$coefficient,levels=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
                                                       "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37","X38","X39",
                                                       "X40","X41","X42","X43","X44","X45","X46","X47","X48","X49","X50"))
position=rep(90,50)
position[ssl_lasso_unknown_variance[[1]]$beta[,100]==0]=seq(1,90,2)
ggplot(coeff_sslasso, aes(x=lambda0,y=value)) +
  geom_line(aes(x=lambda0,y=value,linetype =as.factor(true), color=Coefficients)) +
  scale_linetype_manual("True value",values =c("dashed" ,"solid"))+
  xlab(TeX("$\\zeta_0$")) +
  ylab("Estimated coefficients")+
  geom_text(data = subset(coeff_sslasso, lambda0 == 101), aes(label = coefficient,color=coefficient, x =position, y = value-0.09), hjust = -.1, size=2) +
  ggtitle("Unknown variance SSL estimated coefficients",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 1000 observations.")+ theme_bw()+
  theme(legend.key.width = unit(2,"lines"),
        plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))




#Unknown variance SSL, 10 not null coefficients. 100 obs

lambda0=rep(ssl_lasso_unknown_variance[[2]]$lambda0, each=50)
coeff=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
        "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37","X38","X39",
        "X40","X41","X42","X43","X44","X45","X46","X47","X48","X49","X50"
)
coefficient=rep(coeff,length(ssl_lasso_unknown_variance[[2]]$lambda0))
length(lambda0)==length(coefficient)
value=c(ssl_lasso_unknown_variance[[2]]$beta)
length(value)
true_value=c(rep(0,50))#utile per forma linea se 0 o differente da zero nel modello originale.
true_value[which(true_betas!=0)]=1


coeff_sslasso=data.frame(coefficient,value,lambda0,true=rep(true_value,length(ssl_lasso_unknown_variance[[2]]$lambda0)))
Coefficients=factor(coeff_sslasso$coefficient,levels=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
                                                       "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37","X38","X39",
                                                       "X40","X41","X42","X43","X44","X45","X46","X47","X48","X49","X50"))
position=rep(90,50)
position[ssl_lasso_unknown_variance[[2]]$beta[,100]==0]=seq(20,99,2)

ggplot(coeff_sslasso, aes(x=lambda0,y=value)) +
  geom_line(aes(x=lambda0,y=value,linetype =as.factor(true), color=Coefficients)) +
  scale_linetype_manual("True value",values =c("dashed" ,"solid"))+
  xlab(TeX("$\\zeta_0$")) +
  ylab("Estimated coefficients")+
  geom_text(data = subset(coeff_sslasso, lambda0 == 101), aes(label = coefficient,color=coefficient, x =position, y = value-0.09), hjust = -.1, size=2) +
  ggtitle("Unknown variance SSL estimated coefficients",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 100 observations.")+ theme_bw()+
  theme(legend.key.width = unit(2,"lines"),
        plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))



#Unknown variance SSL, 10 not null coefficients. 60 obs

lambda0=rep(ssl_lasso_unknown_variance[[3]]$lambda0, each=50)
coeff=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
        "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37","X38","X39",
        "X40","X41","X42","X43","X44","X45","X46","X47","X48","X49","X50"
)
coefficient=rep(coeff,length(ssl_lasso_unknown_variance[[3]]$lambda0))
length(lambda0)==length(coefficient)
value=c(ssl_lasso_unknown_variance[[3]]$beta)
length(value)
true_value=c(rep(0,50))#utile per forma linea se 0 o differente da zero nel modello originale.
true_value[which(true_betas!=0)]=1


coeff_sslasso=data.frame(coefficient,value,lambda0,true=rep(true_value,length(ssl_lasso_unknown_variance[[3]]$lambda0)))
Coefficients=factor(coeff_sslasso$coefficient,levels=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
                                                       "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37","X38","X39",
                                                       "X40","X41","X42","X43","X44","X45","X46","X47","X48","X49","X50"))
position=rep(90,50)
position[ssl_lasso_unknown_variance[[3]]$beta[,100]==0]=seq(20,99,2)
position_y=c(subset(coeff_sslasso, lambda0 == 101)$value-0.09)
position[which(ssl_lasso_unknown_variance[[3]]$beta[,20]!=0 & ssl_lasso_unknown_variance[[3]]$beta[,100]==0)]=rep(20,length(which(ssl_lasso_unknown_variance[[3]]$beta[,20]!=0 & ssl_lasso_unknown_variance[[3]]$beta[,100]==0)))
position_y[which(ssl_lasso_unknown_variance[[3]]$beta[,20]!=0 & ssl_lasso_unknown_variance[[3]]$beta[,100]==0)]=c(subset(coeff_sslasso, lambda0 == 20)$value[which(ssl_lasso_unknown_variance[[3]]$beta[,20]!=0 & ssl_lasso_unknown_variance[[3]]$beta[,100]==0)])-0.09
ggplot(coeff_sslasso, aes(x=lambda0,y=value)) +
  geom_line(aes(x=lambda0,y=value,linetype =as.factor(true), color=Coefficients)) +
  scale_linetype_manual("True value",values =c("dashed" ,"solid"))+
  xlab(TeX("$\\zeta_0$")) +
  ylab("Estimated coefficients")+
  geom_text(data = subset(coeff_sslasso, lambda0 == 101), aes(label = coefficient,color=coefficient, x =position, y = position_y), hjust = -.1, size=2) +
  ggtitle("Unknown variance SSL estimated coefficients",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 60 observations.")+ theme_bw()+
  theme(legend.key.width = unit(2,"lines"),
        plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))



#Unknown variance SSL, 10  not null coefficients. 20 obs

lambda0=rep(ssl_lasso_unknown_variance[[6]]$lambda0, each=50)
coeff=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
        "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37","X38","X39",
        "X40","X41","X42","X43","X44","X45","X46","X47","X48","X49","X50"
)
coefficient=rep(coeff,length(ssl_lasso_unknown_variance[[6]]$lambda0))
length(lambda0)==length(coefficient)
value=c(ssl_lasso_unknown_variance[[6]]$beta)
length(value)
true_value=c(rep(0,50))#utile per forma linea se 0 o differente da zero nel modello originale.
true_value[which(true_betas!=0)]=1


coeff_sslasso=data.frame(coefficient,value,lambda0,true=rep(true_value,length(ssl_lasso_unknown_variance[[6]]$lambda0)))
Coefficients=factor(coeff_sslasso$coefficient,levels=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
                                                       "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37","X38","X39",
                                                       "X40","X41","X42","X43","X44","X45","X46","X47","X48","X49","X50"))
position_x=rep(90,50)
position_x[ssl_lasso_unknown_variance[[6]]$beta[,5]==0]=seq(5,99,2)
length(position_x[ssl_lasso_unknown_variance[[6]]$beta[,5]==0])
position_x[which(ssl_lasso_unknown_variance[[6]]$beta[,5]!=0)]=rep(5,2)
position_y=subset(coeff_sslasso, lambda0 == 101)$value -0.1
position_y[which(ssl_lasso_unknown_variance[[6]]$beta[,5]!=0)]=subset(coeff_sslasso, lambda0 == 5)$value[which(ssl_lasso_unknown_variance[[6]]$beta[,5]!=0)] +0.18

ggplot(coeff_sslasso, aes(x=lambda0,y=value)) +
  geom_line(aes(x=lambda0,y=value,linetype =as.factor(true), color=Coefficients)) +
  scale_linetype_manual("True value",values =c("dashed" ,"solid"))+
  xlab(TeX("$\\zeta_0$")) +
  ylab("Estimated coefficients")+
  geom_text(data = subset(coeff_sslasso, lambda0 == 101), aes(label = coefficient,color=coefficient, x =position_x, y = position_y), hjust = -.1, size=2) +
  ggtitle("Unknown variance SSL estimated coefficients",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 20 observations.")+ theme_bw()+
  theme(legend.key.width = unit(2,"lines"),
        plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))



##################################################

#Expected Specificity/sensitivity/accurancy con 1000 obs

expected_truevsbayesian_median[1,]
expected_truevsbayesian_mean[1,]
expected_truevsbayesian_mode[1,]
expected_truevslasso[1,]
expected_truevsSS_George_McCulloch[1,]
expected_truevsSS_Ishwaran_Rao[1,]
truevsSSL_unknown_variance[1,]

#Expected Specificity/sensitivity con 100 obs

expected_truevsbayesian_median[2,]
expected_truevsbayesian_mean[2,]
expected_truevsbayesian_mode[2,]
expected_truevslasso[2,]
expected_truevsSS_George_McCulloch[2,]
expected_truevsSS_Ishwaran_Rao[2,]
truevsSSL_unknown_variance[2,]

#Expected Specificity/sensitivity con 60 obs

expected_truevsbayesian_median[3,]
expected_truevsbayesian_mean[3,]
expected_truevsbayesian_mode[3,]
expected_truevslasso[3,]
expected_truevsSS_George_McCulloch[3,]
expected_truevsSS_Ishwaran_Rao[3,]
truevsSSL_unknown_variance[3,]

#Expected Specificity/sensitivity con 50 obs

expected_truevsbayesian_median[4,]
expected_truevsbayesian_mean[4,]
expected_truevsbayesian_mode[4,]
expected_truevslasso[4,]
expected_truevsSS_George_McCulloch[4,]
expected_truevsSS_Ishwaran_Rao[4,]
truevsSSL_unknown_variance[4,]

#Expected Specificity/sensitivity con 30 obs

expected_truevsbayesian_median[5,]
expected_truevsbayesian_mean[5,]
expected_truevsbayesian_mode[5,]
expected_truevslasso[5,]
expected_truevsSS_George_McCulloch[5,]
expected_truevsSS_Ishwaran_Rao[5,]
truevsSSL_unknown_variance[5,]

#Expected Specificity/sensitivity con 20 obs

expected_truevsbayesian_median[6,]
expected_truevsbayesian_mean[6,]
expected_truevsbayesian_mode[6,]
expected_truevslasso[6,]
expected_truevsSS_George_McCulloch[6,]
expected_truevsSS_Ishwaran_Rao[6,]
truevsSSL_unknown_variance[6,]

########################################################################################################################
########################## HAMMING DISTANCE COMPARISON  ################################################################
########################################################################################################################

Hamming_distance=cbind(expected_hamming_lasso_bayesian_mode,expected_hamming_lasso_bayesian_mean,expected_hamming_lasso_bayesian_median,expected_hamming_lasso_frequentist,expected_hamming_SS,expected_hamming_SS_Ishwaran,SSL_hamming)


colnames(Hamming_distance)=c("Bayesian Lasso Mode GS","Bayesian Lasso Mean GS","Bayesian Lasso Median GS","Frequentist Lasso","George and McCulloch SS","Ishwaran and Rao SS","SSL")
rownames(Hamming_distance)=c("1000 obs.","100 obs.","60 obs.","50 obs.","30 obs.","20 obs.")
Hamming_distance


############# PLOTS ###############

#1000 obs
set.seed(1)
casual=sample(1:50,6,rep=F)
casual=casual[order(casual)]
casual
beta_graph=matrix(0,48,1)
beta_graph[1:8]=rep("beta[10]",8)
beta_graph[9:16]=rep("beta[14]",8)
beta_graph[17:24]=rep("beta[19]",8)
beta_graph[25:32]=rep("beta[28]",8)
beta_graph[33:40]=rep("beta[41]",8)
beta_graph[41:48]=rep("beta[43]",8)

value=c(true_betas[casual[1]],expected_value_bayesian_mode[1,casual[1]],expected_value_bayesian_median[1,casual[1]],expected_value_bayesian_mean[1,casual[1]],expected_value_lasso[1,casual[1]],expected_SS_George_McCulloch[1,casual[1]],expected_SS_Ishwaran_Rao[1,casual[1]],apply(ssl_lasso_unknown_variance[[1]]$beta,1,getmode)[casual[1]],
        true_betas[casual[2]],expected_value_bayesian_mode[1,casual[2]],expected_value_bayesian_median[1,casual[2]],expected_value_bayesian_mean[1,casual[2]],expected_value_lasso[1,casual[2]],expected_SS_George_McCulloch[1,casual[2]],expected_SS_Ishwaran_Rao[1,casual[2]],apply(ssl_lasso_unknown_variance[[1]]$beta,1,getmode)[casual[2]],
        true_betas[casual[3]],expected_value_bayesian_mode[1,casual[3]],expected_value_bayesian_median[1,casual[3]],expected_value_bayesian_mean[1,casual[3]],expected_value_lasso[1,casual[3]],expected_SS_George_McCulloch[1,casual[3]],expected_SS_Ishwaran_Rao[1,casual[3]],apply(ssl_lasso_unknown_variance[[1]]$beta,1,getmode)[casual[3]],
        true_betas[casual[4]],expected_value_bayesian_mode[1,casual[4]],expected_value_bayesian_median[1,casual[4]],expected_value_bayesian_mean[1,casual[4]],expected_value_lasso[1,casual[4]],expected_SS_George_McCulloch[1,casual[4]],expected_SS_Ishwaran_Rao[1,casual[4]],apply(ssl_lasso_unknown_variance[[1]]$beta,1,getmode)[casual[4]],
        true_betas[casual[5]],expected_value_bayesian_mode[1,casual[5]],expected_value_bayesian_median[1,casual[5]],expected_value_bayesian_mean[1,casual[5]],expected_value_lasso[1,casual[5]],expected_SS_George_McCulloch[1,casual[5]],expected_SS_Ishwaran_Rao[1,casual[5]],apply(ssl_lasso_unknown_variance[[1]]$beta,1,getmode)[casual[5]],
        true_betas[casual[6]],expected_value_bayesian_mode[1,casual[6]],expected_value_bayesian_median[1,casual[6]],expected_value_bayesian_mean[1,casual[6]],expected_value_lasso[1,casual[6]],expected_SS_George_McCulloch[1,casual[6]],expected_SS_Ishwaran_Rao[1,casual[6]],apply(ssl_lasso_unknown_variance[[1]]$beta,1,getmode)[casual[6]])

statistics=matrix(0,48,1)
statistics[c(1,9,17,25,33,41),]="True_betas"
statistics[c(2,10,18,26,34,42),]="Bayesian Lasso mode"
statistics[c(3,11,19,27,35,43),]="Bayesian Lasso median"
statistics[c(4,12,20,28,36,44),]="Bayesian Lasso mean"
statistics[c(5,13,21,29,37,45),]="Frequentist Lasso"
statistics[c(6,14,22,30,38,46),]="George and McCulloch SS"
statistics[c(7,15,23,31,39,47),]="Ishwaran and Rao SS"
statistics[c(8,16,24,32,40,48),]="SSL"

data_graph=data.frame(beta_graph=beta_graph,values=value,"estimates and true value"=statistics)
data_graph$estimates.and.true.value   

library(bayestestR)

c_low=c(true_betas[casual[1]],ci(lista_bayesian_mode[[1]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[1]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[1]][,casual[1]],.95,method="HDI")$CI_low,expected_value_lasso[1,casual[1]],expected_SS_George_McCulloch[1,casual[1]],expected_SS_Ishwaran_Rao[1,casual[1]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[1],],.95,method="HDI")$CI_low,
        true_betas[casual[2]],ci(lista_bayesian_mode[[1]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[1]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[1]][,casual[2]],.95,method="HDI")$CI_low,expected_value_lasso[1,casual[2]],expected_SS_George_McCulloch[1,casual[2]],expected_SS_Ishwaran_Rao[1,casual[2]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[2],],.95,method="HDI")$CI_low,
        true_betas[casual[3]],ci(lista_bayesian_mode[[1]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[1]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[1]][,casual[3]],.95,method="HDI")$CI_low,expected_value_lasso[1,casual[3]],expected_SS_George_McCulloch[1,casual[3]],expected_SS_Ishwaran_Rao[1,casual[3]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[3],],.95,method="HDI")$CI_low,
        true_betas[casual[4]],ci(lista_bayesian_mode[[1]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[1]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[1]][,casual[4]],.95,method="HDI")$CI_low,expected_value_lasso[1,casual[4]],expected_SS_George_McCulloch[1,casual[4]],expected_SS_Ishwaran_Rao[1,casual[4]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[4],],.95,method="HDI")$CI_low,
        true_betas[casual[5]],ci(lista_bayesian_mode[[1]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[1]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[1]][,casual[5]],.95,method="HDI")$CI_low,expected_value_lasso[1,casual[5]],expected_SS_George_McCulloch[1,casual[5]],expected_SS_Ishwaran_Rao[1,casual[5]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[5],],.95,method="HDI")$CI_low,
        true_betas[casual[6]],ci(lista_bayesian_mode[[1]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[1]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[1]][,casual[6]],.95,method="HDI")$CI_low,expected_value_lasso[1,casual[6]],expected_SS_George_McCulloch[1,casual[6]],expected_SS_Ishwaran_Rao[1,casual[6]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[6],],.95,method="HDI")$CI_low)

c_high=c(true_betas[casual[1]],ci(lista_bayesian_mode[[1]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[1]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[1]][,casual[1]],.95,method="HDI")$CI_high,expected_value_lasso[1,casual[1]],expected_SS_George_McCulloch[1,casual[1]],expected_SS_Ishwaran_Rao[1,casual[1]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[1],],.95,method="HDI")$CI_high,
         true_betas[casual[2]],ci(lista_bayesian_mode[[1]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[1]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[1]][,casual[2]],.95,method="HDI")$CI_high,expected_value_lasso[1,casual[2]],expected_SS_George_McCulloch[1,casual[2]],expected_SS_Ishwaran_Rao[1,casual[2]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[2],],.95,method="HDI")$CI_high,
         true_betas[casual[3]],ci(lista_bayesian_mode[[1]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[1]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[1]][,casual[3]],.95,method="HDI")$CI_high,expected_value_lasso[1,casual[3]],expected_SS_George_McCulloch[1,casual[3]],expected_SS_Ishwaran_Rao[1,casual[3]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[3],],.95,method="HDI")$CI_high,
         true_betas[casual[4]],ci(lista_bayesian_mode[[1]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[1]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[1]][,casual[4]],.95,method="HDI")$CI_high,expected_value_lasso[1,casual[4]],expected_SS_George_McCulloch[1,casual[4]],expected_SS_Ishwaran_Rao[1,casual[4]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[4],],.95,method="HDI")$CI_high,
         true_betas[casual[5]],ci(lista_bayesian_mode[[1]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[1]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[1]][,casual[5]],.95,method="HDI")$CI_high,expected_value_lasso[1,casual[5]],expected_SS_George_McCulloch[1,casual[5]],expected_SS_Ishwaran_Rao[1,casual[5]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[5],],.95,method="HDI")$CI_high,
         true_betas[casual[6]],ci(lista_bayesian_mode[[1]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[1]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[1]][,casual[6]],.95,method="HDI")$CI_high,expected_value_lasso[1,casual[6]],expected_SS_George_McCulloch[1,casual[6]],expected_SS_Ishwaran_Rao[1,casual[6]],ci(ssl_lasso_unknown_variance[[1]]$beta[casual[6],],.95,method="HDI")$CI_high)

data_graph$c_low=c_low   
data_graph$c_high=c_high 
levels(data_graph$beta_graph)
data_graph$beta_graph=factor(data_graph$beta_graph, levels = c("beta[10]", "beta[14]", "beta[19]", "beta[28]", "beta[41]", "beta[43]"))
ggplot(data_graph, aes(colour = estimates.and.true.value))+ geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+ 
  geom_linerange(aes(x = beta_graph, ymin = c_low,
                     ymax =c_high),
                 lwd = 1, position = position_dodge(width = 1/2))+ 
  geom_pointrange(aes(x = beta_graph, y = values, ymin = c_low ,
                      ymax =c_high ),
                  lwd = 1/2, position = position_dodge(width = 1/2),
                  shape = 21, fill = "WHITE") + coord_flip() + labs(title="Comparison of Bayesian Lasso estimates, Frequentist Lasso estimates, George and McCulloch Spike and Slab, \n Ishwaran and Rao Spike and Slab, Spike-and- Slab Lasso and true coefficients.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 1000 observations.") +xlab("Coefficients")+ylab("Values")+labs(colour="Estimates and true values")+
  theme_bw()+theme(plot.title=element_text(size=10, hjust=0.5,face="bold"),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))

#100 obs
set.seed(1)
casual=sample(1:50,6,rep=F)
casual=casual[order(casual)]
casual
beta_graph=matrix(0,48,1)
beta_graph[1:8]=rep("beta[10]",8)
beta_graph[9:16]=rep("beta[14]",8)
beta_graph[17:24]=rep("beta[19]",8)
beta_graph[25:32]=rep("beta[28]",8)
beta_graph[33:40]=rep("beta[41]",8)
beta_graph[41:48]=rep("beta[43]",8)

value=c(true_betas[casual[1]],expected_value_bayesian_mode[2,casual[1]],expected_value_bayesian_median[2,casual[1]],expected_value_bayesian_mean[2,casual[1]],expected_value_lasso[2,casual[1]],expected_SS_George_McCulloch[2,casual[1]],expected_SS_Ishwaran_Rao[2,casual[1]],apply(ssl_lasso_unknown_variance[[2]]$beta,1,getmode)[casual[1]],
        true_betas[casual[2]],expected_value_bayesian_mode[2,casual[2]],expected_value_bayesian_median[2,casual[2]],expected_value_bayesian_mean[2,casual[2]],expected_value_lasso[2,casual[2]],expected_SS_George_McCulloch[2,casual[2]],expected_SS_Ishwaran_Rao[2,casual[2]],apply(ssl_lasso_unknown_variance[[2]]$beta,1,getmode)[casual[2]],
        true_betas[casual[3]],expected_value_bayesian_mode[2,casual[3]],expected_value_bayesian_median[2,casual[3]],expected_value_bayesian_mean[2,casual[3]],expected_value_lasso[2,casual[3]],expected_SS_George_McCulloch[2,casual[3]],expected_SS_Ishwaran_Rao[2,casual[3]],apply(ssl_lasso_unknown_variance[[2]]$beta,1,getmode)[casual[3]],
        true_betas[casual[4]],expected_value_bayesian_mode[2,casual[4]],expected_value_bayesian_median[2,casual[4]],expected_value_bayesian_mean[2,casual[4]],expected_value_lasso[2,casual[4]],expected_SS_George_McCulloch[2,casual[4]],expected_SS_Ishwaran_Rao[2,casual[4]],apply(ssl_lasso_unknown_variance[[2]]$beta,1,getmode)[casual[4]],
        true_betas[casual[5]],expected_value_bayesian_mode[2,casual[5]],expected_value_bayesian_median[2,casual[5]],expected_value_bayesian_mean[2,casual[5]],expected_value_lasso[2,casual[5]],expected_SS_George_McCulloch[2,casual[5]],expected_SS_Ishwaran_Rao[2,casual[5]],apply(ssl_lasso_unknown_variance[[2]]$beta,1,getmode)[casual[5]],
        true_betas[casual[6]],expected_value_bayesian_mode[2,casual[6]],expected_value_bayesian_median[2,casual[6]],expected_value_bayesian_mean[2,casual[6]],expected_value_lasso[2,casual[6]],expected_SS_George_McCulloch[2,casual[6]],expected_SS_Ishwaran_Rao[2,casual[6]],apply(ssl_lasso_unknown_variance[[2]]$beta,1,getmode)[casual[6]])

statistics=matrix(0,48,1)
statistics[c(1,9,17,25,33,41),]="True_betas"
statistics[c(2,10,18,26,34,42),]="Bayesian Lasso mode"
statistics[c(3,11,19,27,35,43),]="Bayesian Lasso median"
statistics[c(4,12,20,28,36,44),]="Bayesian Lasso mean"
statistics[c(5,13,21,29,37,45),]="Frequentist Lasso"
statistics[c(6,14,22,30,38,46),]="George and McCulloch SS"
statistics[c(7,15,23,31,39,47),]="Ishwaran and Rao SS"
statistics[c(8,16,24,32,40,48),]="SSL"

data_graph=data.frame(beta_graph=beta_graph,values=value,"estimates and true value"=statistics)
data_graph$estimates.and.true.value   

library(bayestestR)

c_low=c(true_betas[casual[1]],ci(lista_bayesian_mode[[2]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[2]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[2]][,casual[1]],.95,method="HDI")$CI_low,expected_value_lasso[2,casual[1]],expected_SS_George_McCulloch[2,casual[1]],expected_SS_Ishwaran_Rao[2,casual[1]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[1],],.95,method="HDI")$CI_low,
        true_betas[casual[2]],ci(lista_bayesian_mode[[2]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[2]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[2]][,casual[2]],.95,method="HDI")$CI_low,expected_value_lasso[2,casual[2]],expected_SS_George_McCulloch[2,casual[2]],expected_SS_Ishwaran_Rao[2,casual[2]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[2],],.95,method="HDI")$CI_low,
        true_betas[casual[3]],ci(lista_bayesian_mode[[2]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[2]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[2]][,casual[3]],.95,method="HDI")$CI_low,expected_value_lasso[2,casual[3]],expected_SS_George_McCulloch[2,casual[3]],expected_SS_Ishwaran_Rao[2,casual[3]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[3],],.95,method="HDI")$CI_low,
        true_betas[casual[4]],ci(lista_bayesian_mode[[2]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[2]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[2]][,casual[4]],.95,method="HDI")$CI_low,expected_value_lasso[2,casual[4]],expected_SS_George_McCulloch[2,casual[4]],expected_SS_Ishwaran_Rao[2,casual[4]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[4],],.95,method="HDI")$CI_low,
        true_betas[casual[5]],ci(lista_bayesian_mode[[2]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[2]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[2]][,casual[5]],.95,method="HDI")$CI_low,expected_value_lasso[2,casual[5]],expected_SS_George_McCulloch[2,casual[5]],expected_SS_Ishwaran_Rao[2,casual[5]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[5],],.95,method="HDI")$CI_low,
        true_betas[casual[6]],ci(lista_bayesian_mode[[2]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[2]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[2]][,casual[6]],.95,method="HDI")$CI_low,expected_value_lasso[2,casual[6]],expected_SS_George_McCulloch[2,casual[6]],expected_SS_Ishwaran_Rao[2,casual[6]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[6],],.95,method="HDI")$CI_low)

c_high=c(true_betas[casual[1]],ci(lista_bayesian_mode[[2]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[2]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[2]][,casual[1]],.95,method="HDI")$CI_high,expected_value_lasso[2,casual[1]],expected_SS_George_McCulloch[2,casual[1]],expected_SS_Ishwaran_Rao[2,casual[1]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[1],],.95,method="HDI")$CI_high,
         true_betas[casual[2]],ci(lista_bayesian_mode[[2]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[2]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[2]][,casual[2]],.95,method="HDI")$CI_high,expected_value_lasso[2,casual[2]],expected_SS_George_McCulloch[2,casual[2]],expected_SS_Ishwaran_Rao[2,casual[2]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[2],],.95,method="HDI")$CI_high,
         true_betas[casual[3]],ci(lista_bayesian_mode[[2]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[2]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[2]][,casual[3]],.95,method="HDI")$CI_high,expected_value_lasso[2,casual[3]],expected_SS_George_McCulloch[2,casual[3]],expected_SS_Ishwaran_Rao[2,casual[3]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[3],],.95,method="HDI")$CI_high,
         true_betas[casual[4]],ci(lista_bayesian_mode[[2]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[2]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[2]][,casual[4]],.95,method="HDI")$CI_high,expected_value_lasso[2,casual[4]],expected_SS_George_McCulloch[2,casual[4]],expected_SS_Ishwaran_Rao[2,casual[4]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[4],],.95,method="HDI")$CI_high,
         true_betas[casual[5]],ci(lista_bayesian_mode[[2]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[2]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[2]][,casual[5]],.95,method="HDI")$CI_high,expected_value_lasso[2,casual[5]],expected_SS_George_McCulloch[2,casual[5]],expected_SS_Ishwaran_Rao[2,casual[5]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[5],],.95,method="HDI")$CI_high,
         true_betas[casual[6]],ci(lista_bayesian_mode[[2]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[2]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[2]][,casual[6]],.95,method="HDI")$CI_high,expected_value_lasso[2,casual[6]],expected_SS_George_McCulloch[2,casual[6]],expected_SS_Ishwaran_Rao[2,casual[6]],ci(ssl_lasso_unknown_variance[[2]]$beta[casual[6],],.95,method="HDI")$CI_high)

data_graph$c_low=c_low   
data_graph$c_high=c_high 
levels(data_graph$beta_graph)
data_graph$beta_graph=factor(data_graph$beta_graph, levels = c("beta[10]", "beta[14]", "beta[19]", "beta[28]", "beta[41]", "beta[43]"))
ggplot(data_graph, aes(colour = estimates.and.true.value))+ geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+ 
  geom_linerange(aes(x = beta_graph, ymin = c_low,
                     ymax =c_high),
                 lwd = 1, position = position_dodge(width = 1/2))+ 
  geom_pointrange(aes(x = beta_graph, y = values, ymin = c_low ,
                      ymax =c_high ),
                  lwd = 1/2, position = position_dodge(width = 1/2),
                  shape = 21, fill = "WHITE") + coord_flip() + labs(title="Comparison of Bayesian Lasso estimates, Frequentist Lasso estimates, George and McCulloch Spike and Slab, \n Ishwaran and Rao Spike and Slab, Spike-and- Slab Lasso and true coefficients.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 100 observations.") +xlab("Coefficients")+ylab("Values")+labs(colour="Estimates and true values")+
  theme_bw()+theme(plot.title=element_text(size=10, hjust=0.5,face="bold"),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))

####60 obs
set.seed(1)
casual=sample(1:50,6,rep=F)
casual=casual[order(casual)]
casual
beta_graph=matrix(0,48,1)
beta_graph[1:8]=rep("beta[10]",8)
beta_graph[9:16]=rep("beta[14]",8)
beta_graph[17:24]=rep("beta[19]",8)
beta_graph[25:32]=rep("beta[28]",8)
beta_graph[33:40]=rep("beta[41]",8)
beta_graph[41:48]=rep("beta[43]",8)

value=c(true_betas[casual[1]],expected_value_bayesian_mode[3,casual[1]],expected_value_bayesian_median[3,casual[1]],expected_value_bayesian_mean[3,casual[1]],expected_value_lasso[3,casual[1]],expected_SS_George_McCulloch[3,casual[1]],expected_SS_Ishwaran_Rao[3,casual[1]],apply(ssl_lasso_unknown_variance[[3]]$beta,1,getmode)[casual[1]],
        true_betas[casual[2]],expected_value_bayesian_mode[3,casual[2]],expected_value_bayesian_median[3,casual[2]],expected_value_bayesian_mean[3,casual[2]],expected_value_lasso[3,casual[2]],expected_SS_George_McCulloch[3,casual[2]],expected_SS_Ishwaran_Rao[3,casual[2]],apply(ssl_lasso_unknown_variance[[3]]$beta,1,getmode)[casual[2]],
        true_betas[casual[3]],expected_value_bayesian_mode[3,casual[3]],expected_value_bayesian_median[3,casual[3]],expected_value_bayesian_mean[3,casual[3]],expected_value_lasso[3,casual[3]],expected_SS_George_McCulloch[3,casual[3]],expected_SS_Ishwaran_Rao[3,casual[3]],apply(ssl_lasso_unknown_variance[[3]]$beta,1,getmode)[casual[3]],
        true_betas[casual[4]],expected_value_bayesian_mode[3,casual[4]],expected_value_bayesian_median[3,casual[4]],expected_value_bayesian_mean[3,casual[4]],expected_value_lasso[3,casual[4]],expected_SS_George_McCulloch[3,casual[4]],expected_SS_Ishwaran_Rao[3,casual[4]],apply(ssl_lasso_unknown_variance[[3]]$beta,1,getmode)[casual[4]],
        true_betas[casual[5]],expected_value_bayesian_mode[3,casual[5]],expected_value_bayesian_median[3,casual[5]],expected_value_bayesian_mean[3,casual[5]],expected_value_lasso[3,casual[5]],expected_SS_George_McCulloch[3,casual[5]],expected_SS_Ishwaran_Rao[3,casual[5]],apply(ssl_lasso_unknown_variance[[3]]$beta,1,getmode)[casual[5]],
        true_betas[casual[6]],expected_value_bayesian_mode[3,casual[6]],expected_value_bayesian_median[3,casual[6]],expected_value_bayesian_mean[3,casual[6]],expected_value_lasso[3,casual[6]],expected_SS_George_McCulloch[3,casual[6]],expected_SS_Ishwaran_Rao[3,casual[6]],apply(ssl_lasso_unknown_variance[[3]]$beta,1,getmode)[casual[6]])

statistics=matrix(0,48,1)
statistics[c(1,9,17,25,33,41),]="True_betas"
statistics[c(2,10,18,26,34,42),]="Bayesian Lasso mode"
statistics[c(3,11,19,27,35,43),]="Bayesian Lasso median"
statistics[c(4,12,20,28,36,44),]="Bayesian Lasso mean"
statistics[c(5,13,21,29,37,45),]="Frequentist Lasso"
statistics[c(6,14,22,30,38,46),]="George and McCulloch SS"
statistics[c(7,15,23,31,39,47),]="Ishwaran and Rao SS"
statistics[c(8,16,24,32,40,48),]="SSL"

data_graph=data.frame(beta_graph=beta_graph,values=value,"estimates and true value"=statistics)
data_graph$estimates.and.true.value   

library(bayestestR)

c_low=c(true_betas[casual[1]],ci(lista_bayesian_mode[[3]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[3]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[3]][,casual[1]],.95,method="HDI")$CI_low,expected_value_lasso[3,casual[1]],expected_SS_George_McCulloch[3,casual[1]],expected_SS_Ishwaran_Rao[3,casual[1]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[1],],.95,method="HDI")$CI_low,
        true_betas[casual[2]],ci(lista_bayesian_mode[[3]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[3]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[3]][,casual[2]],.95,method="HDI")$CI_low,expected_value_lasso[3,casual[2]],expected_SS_George_McCulloch[3,casual[2]],expected_SS_Ishwaran_Rao[3,casual[2]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[2],],.95,method="HDI")$CI_low,
        true_betas[casual[3]],ci(lista_bayesian_mode[[3]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[3]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[3]][,casual[3]],.95,method="HDI")$CI_low,expected_value_lasso[3,casual[3]],expected_SS_George_McCulloch[3,casual[3]],expected_SS_Ishwaran_Rao[3,casual[3]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[3],],.95,method="HDI")$CI_low,
        true_betas[casual[4]],ci(lista_bayesian_mode[[3]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[3]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[3]][,casual[4]],.95,method="HDI")$CI_low,expected_value_lasso[3,casual[4]],expected_SS_George_McCulloch[3,casual[4]],expected_SS_Ishwaran_Rao[3,casual[4]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[4],],.95,method="HDI")$CI_low,
        true_betas[casual[5]],ci(lista_bayesian_mode[[3]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[3]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[3]][,casual[5]],.95,method="HDI")$CI_low,expected_value_lasso[3,casual[5]],expected_SS_George_McCulloch[3,casual[5]],expected_SS_Ishwaran_Rao[3,casual[5]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[5],],.95,method="HDI")$CI_low,
        true_betas[casual[6]],ci(lista_bayesian_mode[[3]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[3]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[3]][,casual[6]],.95,method="HDI")$CI_low,expected_value_lasso[3,casual[6]],expected_SS_George_McCulloch[3,casual[6]],expected_SS_Ishwaran_Rao[3,casual[6]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[6],],.95,method="HDI")$CI_low)

c_high=c(true_betas[casual[1]],ci(lista_bayesian_mode[[3]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[3]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[3]][,casual[1]],.95,method="HDI")$CI_high,expected_value_lasso[3,casual[1]],expected_SS_George_McCulloch[3,casual[1]],expected_SS_Ishwaran_Rao[3,casual[1]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[1],],.95,method="HDI")$CI_high,
         true_betas[casual[2]],ci(lista_bayesian_mode[[3]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[3]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[3]][,casual[2]],.95,method="HDI")$CI_high,expected_value_lasso[3,casual[2]],expected_SS_George_McCulloch[3,casual[2]],expected_SS_Ishwaran_Rao[3,casual[2]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[2],],.95,method="HDI")$CI_high,
         true_betas[casual[3]],ci(lista_bayesian_mode[[3]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[3]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[3]][,casual[3]],.95,method="HDI")$CI_high,expected_value_lasso[3,casual[3]],expected_SS_George_McCulloch[3,casual[3]],expected_SS_Ishwaran_Rao[3,casual[3]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[3],],.95,method="HDI")$CI_high,
         true_betas[casual[4]],ci(lista_bayesian_mode[[3]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[3]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[3]][,casual[4]],.95,method="HDI")$CI_high,expected_value_lasso[3,casual[4]],expected_SS_George_McCulloch[3,casual[4]],expected_SS_Ishwaran_Rao[3,casual[4]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[4],],.95,method="HDI")$CI_high,
         true_betas[casual[5]],ci(lista_bayesian_mode[[3]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[3]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[3]][,casual[5]],.95,method="HDI")$CI_high,expected_value_lasso[3,casual[5]],expected_SS_George_McCulloch[3,casual[5]],expected_SS_Ishwaran_Rao[3,casual[5]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[5],],.95,method="HDI")$CI_high,
         true_betas[casual[6]],ci(lista_bayesian_mode[[3]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[3]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[3]][,casual[6]],.95,method="HDI")$CI_high,expected_value_lasso[3,casual[6]],expected_SS_George_McCulloch[3,casual[6]],expected_SS_Ishwaran_Rao[3,casual[6]],ci(ssl_lasso_unknown_variance[[3]]$beta[casual[6],],.95,method="HDI")$CI_high)

data_graph$c_low=c_low  
data_graph$c_high=c_high 
levels(data_graph$beta_graph)
data_graph$beta_graph=factor(data_graph$beta_graph, levels = c("beta[10]", "beta[14]", "beta[19]", "beta[28]", "beta[41]", "beta[43]"))
ggplot(data_graph, aes(colour = estimates.and.true.value))+ geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+ 
  geom_linerange(aes(x = beta_graph, ymin = c_low,
                     ymax =c_high),
                 lwd = 1, position = position_dodge(width = 1/2))+ 
  geom_pointrange(aes(x = beta_graph, y = values, ymin = c_low ,
                      ymax =c_high ),
                  lwd = 1/2, position = position_dodge(width = 1/2),
                  shape = 21, fill = "WHITE") + coord_flip() + labs(title="Comparison of Bayesian Lasso estimates, Frequentist Lasso estimates, George and McCulloch Spike and Slab, \n Ishwaran and Rao Spike and Slab, Spike-and- Slab Lasso and true coefficients.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 60 observations.") +xlab("Coefficients")+ylab("Values")+labs(colour="Estimates and true values")+
  theme_bw()+theme(plot.title=element_text(size=10, hjust=0.5,face="bold"),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))


####50 obs
set.seed(1)
casual=sample(1:50,6,rep=F)
casual=casual[order(casual)]
casual
beta_graph=matrix(0,48,1)
beta_graph[1:8]=rep("beta[10]",8)
beta_graph[9:16]=rep("beta[14]",8)
beta_graph[17:24]=rep("beta[19]",8)
beta_graph[25:32]=rep("beta[28]",8)
beta_graph[33:40]=rep("beta[41]",8)
beta_graph[41:48]=rep("beta[43]",8)

value=c(true_betas[casual[1]],expected_value_bayesian_mode[4,casual[1]],expected_value_bayesian_median[4,casual[1]],expected_value_bayesian_mean[4,casual[1]],expected_value_lasso[4,casual[1]],expected_SS_George_McCulloch[4,casual[1]],expected_SS_Ishwaran_Rao[4,casual[1]],apply(ssl_lasso_unknown_variance[[4]]$beta,1,getmode)[casual[1]],
        true_betas[casual[2]],expected_value_bayesian_mode[4,casual[2]],expected_value_bayesian_median[4,casual[2]],expected_value_bayesian_mean[4,casual[2]],expected_value_lasso[4,casual[2]],expected_SS_George_McCulloch[4,casual[2]],expected_SS_Ishwaran_Rao[4,casual[2]],apply(ssl_lasso_unknown_variance[[4]]$beta,1,getmode)[casual[2]],
        true_betas[casual[3]],expected_value_bayesian_mode[4,casual[3]],expected_value_bayesian_median[4,casual[3]],expected_value_bayesian_mean[4,casual[3]],expected_value_lasso[4,casual[3]],expected_SS_George_McCulloch[4,casual[3]],expected_SS_Ishwaran_Rao[4,casual[3]],apply(ssl_lasso_unknown_variance[[4]]$beta,1,getmode)[casual[3]],
        true_betas[casual[4]],expected_value_bayesian_mode[4,casual[4]],expected_value_bayesian_median[4,casual[4]],expected_value_bayesian_mean[4,casual[4]],expected_value_lasso[4,casual[4]],expected_SS_George_McCulloch[4,casual[4]],expected_SS_Ishwaran_Rao[4,casual[4]],apply(ssl_lasso_unknown_variance[[4]]$beta,1,getmode)[casual[4]],
        true_betas[casual[5]],expected_value_bayesian_mode[4,casual[5]],expected_value_bayesian_median[4,casual[5]],expected_value_bayesian_mean[4,casual[5]],expected_value_lasso[4,casual[5]],expected_SS_George_McCulloch[4,casual[5]],expected_SS_Ishwaran_Rao[4,casual[5]],apply(ssl_lasso_unknown_variance[[4]]$beta,1,getmode)[casual[5]],
        true_betas[casual[6]],expected_value_bayesian_mode[4,casual[6]],expected_value_bayesian_median[4,casual[6]],expected_value_bayesian_mean[4,casual[6]],expected_value_lasso[4,casual[6]],expected_SS_George_McCulloch[4,casual[6]],expected_SS_Ishwaran_Rao[4,casual[6]],apply(ssl_lasso_unknown_variance[[4]]$beta,1,getmode)[casual[6]])



statistics=matrix(0,48,1)
statistics[c(1,9,17,25,33,41),]="True_betas"
statistics[c(2,10,18,26,34,42),]="Bayesian Lasso mode"
statistics[c(3,11,19,27,35,43),]="Bayesian Lasso median"
statistics[c(4,12,20,28,36,44),]="Bayesian Lasso mean"
statistics[c(5,13,21,29,37,45),]="Frequentist Lasso"
statistics[c(6,14,22,30,38,46),]="George and McCulloch SS"
statistics[c(7,15,23,31,39,47),]="Ishwaran and Rao SS"
statistics[c(8,16,24,32,40,48),]="SSL"

data_graph=data.frame(beta_graph=beta_graph,values=value,"estimates and true value"=statistics)
data_graph$estimates.and.true.value   

library(bayestestR)

c_low=c(true_betas[casual[1]],ci(lista_bayesian_mode[[4]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[4]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[4]][,casual[1]],.95,method="HDI")$CI_low,expected_value_lasso[4,casual[1]],expected_SS_George_McCulloch[4,casual[1]],expected_SS_Ishwaran_Rao[4,casual[1]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[1],],.95,method="HDI")$CI_low,
        true_betas[casual[2]],ci(lista_bayesian_mode[[4]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[4]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[4]][,casual[2]],.95,method="HDI")$CI_low,expected_value_lasso[4,casual[2]],expected_SS_George_McCulloch[4,casual[2]],expected_SS_Ishwaran_Rao[4,casual[2]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[2],],.95,method="HDI")$CI_low,
        true_betas[casual[3]],ci(lista_bayesian_mode[[4]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[4]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[4]][,casual[3]],.95,method="HDI")$CI_low,expected_value_lasso[4,casual[3]],expected_SS_George_McCulloch[4,casual[3]],expected_SS_Ishwaran_Rao[4,casual[3]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[3],],.95,method="HDI")$CI_low,
        true_betas[casual[4]],ci(lista_bayesian_mode[[4]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[4]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[4]][,casual[4]],.95,method="HDI")$CI_low,expected_value_lasso[4,casual[4]],expected_SS_George_McCulloch[4,casual[4]],expected_SS_Ishwaran_Rao[4,casual[4]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[4],],.95,method="HDI")$CI_low,
        true_betas[casual[5]],ci(lista_bayesian_mode[[4]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[4]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[4]][,casual[5]],.95,method="HDI")$CI_low,expected_value_lasso[4,casual[5]],expected_SS_George_McCulloch[4,casual[5]],expected_SS_Ishwaran_Rao[4,casual[5]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[5],],.95,method="HDI")$CI_low,
        true_betas[casual[6]],ci(lista_bayesian_mode[[4]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[4]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[4]][,casual[6]],.95,method="HDI")$CI_low,expected_value_lasso[4,casual[6]],expected_SS_George_McCulloch[4,casual[6]],expected_SS_Ishwaran_Rao[4,casual[6]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[6],],.95,method="HDI")$CI_low)

c_high=c(true_betas[casual[1]],ci(lista_bayesian_mode[[4]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[4]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[4]][,casual[1]],.95,method="HDI")$CI_high,expected_value_lasso[4,casual[1]],expected_SS_George_McCulloch[4,casual[1]],expected_SS_Ishwaran_Rao[4,casual[1]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[1],],.95,method="HDI")$CI_high,
         true_betas[casual[2]],ci(lista_bayesian_mode[[4]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[4]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[4]][,casual[2]],.95,method="HDI")$CI_high,expected_value_lasso[4,casual[2]],expected_SS_George_McCulloch[4,casual[2]],expected_SS_Ishwaran_Rao[4,casual[2]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[2],],.95,method="HDI")$CI_high,
         true_betas[casual[3]],ci(lista_bayesian_mode[[4]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[4]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[4]][,casual[3]],.95,method="HDI")$CI_high,expected_value_lasso[4,casual[3]],expected_SS_George_McCulloch[4,casual[3]],expected_SS_Ishwaran_Rao[4,casual[3]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[3],],.95,method="HDI")$CI_high,
         true_betas[casual[4]],ci(lista_bayesian_mode[[4]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[4]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[4]][,casual[4]],.95,method="HDI")$CI_high,expected_value_lasso[4,casual[4]],expected_SS_George_McCulloch[4,casual[4]],expected_SS_Ishwaran_Rao[4,casual[4]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[4],],.95,method="HDI")$CI_high,
         true_betas[casual[5]],ci(lista_bayesian_mode[[4]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[4]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[4]][,casual[5]],.95,method="HDI")$CI_high,expected_value_lasso[4,casual[5]],expected_SS_George_McCulloch[4,casual[5]],expected_SS_Ishwaran_Rao[4,casual[5]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[5],],.95,method="HDI")$CI_high,
         true_betas[casual[6]],ci(lista_bayesian_mode[[4]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[4]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[4]][,casual[6]],.95,method="HDI")$CI_high,expected_value_lasso[4,casual[6]],expected_SS_George_McCulloch[4,casual[6]],expected_SS_Ishwaran_Rao[4,casual[6]],ci(ssl_lasso_unknown_variance[[4]]$beta[casual[6],],.95,method="HDI")$CI_high)

data_graph$c_low=c_low   
data_graph$c_high=c_high 
levels(data_graph$beta_graph)
data_graph$beta_graph=factor(data_graph$beta_graph, levels = c("beta[10]", "beta[14]", "beta[19]", "beta[28]", "beta[41]", "beta[43]"))
ggplot(data_graph, aes(colour = estimates.and.true.value))+ geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+ 
  geom_linerange(aes(x = beta_graph, ymin = c_low,
                     ymax =c_high),
                 lwd = 1, position = position_dodge(width = 1/2))+ 
  geom_pointrange(aes(x = beta_graph, y = values, ymin = c_low ,
                      ymax =c_high ),
                  lwd = 1/2, position = position_dodge(width = 1/2),
                  shape = 21, fill = "WHITE") + coord_flip() + labs(title="Comparison of Bayesian Lasso estimates, Frequentist Lasso estimates, George and McCulloch Spike and Slab, \n Ishwaran and Rao Spike and Slab, Spike-and- Slab Lasso and true coefficients.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 50 observations.") +xlab("Coefficients")+ylab("Values")+labs(colour="Estimates and true values")+
  theme_bw()+theme(plot.title=element_text(size=10, hjust=0.5,face="bold"),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))
####30 obs
set.seed(1)
casual=sample(1:50,6,rep=F)
casual=casual[order(casual)]
casual
beta_graph=matrix(0,48,1)
beta_graph[1:8]=rep("beta[10]",8)
beta_graph[9:16]=rep("beta[14]",8)
beta_graph[17:24]=rep("beta[19]",8)
beta_graph[25:32]=rep("beta[28]",8)
beta_graph[33:40]=rep("beta[41]",8)
beta_graph[41:48]=rep("beta[43]",8)

value=c(true_betas[casual[1]],expected_value_bayesian_mode[5,casual[1]],expected_value_bayesian_median[5,casual[1]],expected_value_bayesian_mean[5,casual[1]],expected_value_lasso[5,casual[1]],expected_SS_George_McCulloch[5,casual[1]],expected_SS_Ishwaran_Rao[5,casual[1]],apply(ssl_lasso_unknown_variance[[5]]$beta,1,getmode)[casual[1]],
        true_betas[casual[2]],expected_value_bayesian_mode[5,casual[2]],expected_value_bayesian_median[5,casual[2]],expected_value_bayesian_mean[5,casual[2]],expected_value_lasso[5,casual[2]],expected_SS_George_McCulloch[5,casual[2]],expected_SS_Ishwaran_Rao[5,casual[2]],apply(ssl_lasso_unknown_variance[[5]]$beta,1,getmode)[casual[2]],
        true_betas[casual[3]],expected_value_bayesian_mode[5,casual[3]],expected_value_bayesian_median[5,casual[3]],expected_value_bayesian_mean[5,casual[3]],expected_value_lasso[5,casual[3]],expected_SS_George_McCulloch[5,casual[3]],expected_SS_Ishwaran_Rao[5,casual[3]],apply(ssl_lasso_unknown_variance[[5]]$beta,1,getmode)[casual[3]],
        true_betas[casual[4]],expected_value_bayesian_mode[5,casual[4]],expected_value_bayesian_median[5,casual[4]],expected_value_bayesian_mean[5,casual[4]],expected_value_lasso[5,casual[4]],expected_SS_George_McCulloch[5,casual[4]],expected_SS_Ishwaran_Rao[5,casual[4]],apply(ssl_lasso_unknown_variance[[5]]$beta,1,getmode)[casual[4]],
        true_betas[casual[5]],expected_value_bayesian_mode[5,casual[5]],expected_value_bayesian_median[5,casual[5]],expected_value_bayesian_mean[5,casual[5]],expected_value_lasso[5,casual[5]],expected_SS_George_McCulloch[5,casual[5]],expected_SS_Ishwaran_Rao[5,casual[5]],apply(ssl_lasso_unknown_variance[[5]]$beta,1,getmode)[casual[5]],
        true_betas[casual[6]],expected_value_bayesian_mode[5,casual[6]],expected_value_bayesian_median[5,casual[6]],expected_value_bayesian_mean[5,casual[6]],expected_value_lasso[5,casual[6]],expected_SS_George_McCulloch[5,casual[6]],expected_SS_Ishwaran_Rao[5,casual[6]],apply(ssl_lasso_unknown_variance[[5]]$beta,1,getmode)[casual[6]])


statistics=matrix(0,48,1)
statistics[c(1,9,17,25,33,41),]="True_betas"
statistics[c(2,10,18,26,34,42),]="Bayesian Lasso mode"
statistics[c(3,11,19,27,35,43),]="Bayesian Lasso median"
statistics[c(4,12,20,28,36,44),]="Bayesian Lasso mean"
statistics[c(5,13,21,29,37,45),]="Frequentist Lasso"
statistics[c(6,14,22,30,38,46),]="George and McCulloch SS"
statistics[c(7,15,23,31,39,47),]="Ishwaran and Rao SS"
statistics[c(8,16,24,32,40,48),]="SSL"

data_graph=data.frame(beta_graph=beta_graph,values=value,"estimates and true value"=statistics)
data_graph$estimates.and.true.value   


c_low=c(true_betas[casual[1]],ci(lista_bayesian_mode[[5]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[5]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[5]][,casual[1]],.95,method="HDI")$CI_low,expected_value_lasso[5,casual[1]],expected_SS_George_McCulloch[5,casual[1]],expected_SS_Ishwaran_Rao[5,casual[1]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[1],],.95,method="HDI")$CI_low,
        true_betas[casual[2]],ci(lista_bayesian_mode[[5]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[5]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[5]][,casual[2]],.95,method="HDI")$CI_low,expected_value_lasso[5,casual[2]],expected_SS_George_McCulloch[5,casual[2]],expected_SS_Ishwaran_Rao[5,casual[2]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[2],],.95,method="HDI")$CI_low,
        true_betas[casual[3]],ci(lista_bayesian_mode[[5]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[5]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[5]][,casual[3]],.95,method="HDI")$CI_low,expected_value_lasso[5,casual[3]],expected_SS_George_McCulloch[5,casual[3]],expected_SS_Ishwaran_Rao[5,casual[3]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[3],],.95,method="HDI")$CI_low,
        true_betas[casual[4]],ci(lista_bayesian_mode[[5]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[5]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[5]][,casual[4]],.95,method="HDI")$CI_low,expected_value_lasso[5,casual[4]],expected_SS_George_McCulloch[5,casual[4]],expected_SS_Ishwaran_Rao[5,casual[4]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[4],],.95,method="HDI")$CI_low,
        true_betas[casual[5]],ci(lista_bayesian_mode[[5]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[5]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[5]][,casual[5]],.95,method="HDI")$CI_low,expected_value_lasso[5,casual[5]],expected_SS_George_McCulloch[5,casual[5]],expected_SS_Ishwaran_Rao[5,casual[5]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[5],],.95,method="HDI")$CI_low,
        true_betas[casual[6]],ci(lista_bayesian_mode[[5]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[5]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[5]][,casual[6]],.95,method="HDI")$CI_low,expected_value_lasso[5,casual[6]],expected_SS_George_McCulloch[5,casual[6]],expected_SS_Ishwaran_Rao[5,casual[6]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[6],],.95,method="HDI")$CI_low)

c_high=c(true_betas[casual[1]],ci(lista_bayesian_mode[[5]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[5]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[5]][,casual[1]],.95,method="HDI")$CI_high,expected_value_lasso[5,casual[1]],expected_SS_George_McCulloch[5,casual[1]],expected_SS_Ishwaran_Rao[5,casual[1]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[1],],.95,method="HDI")$CI_high,
         true_betas[casual[2]],ci(lista_bayesian_mode[[5]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[5]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[5]][,casual[2]],.95,method="HDI")$CI_high,expected_value_lasso[5,casual[2]],expected_SS_George_McCulloch[5,casual[2]],expected_SS_Ishwaran_Rao[5,casual[2]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[2],],.95,method="HDI")$CI_high,
         true_betas[casual[3]],ci(lista_bayesian_mode[[5]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[5]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[5]][,casual[3]],.95,method="HDI")$CI_high,expected_value_lasso[5,casual[3]],expected_SS_George_McCulloch[5,casual[3]],expected_SS_Ishwaran_Rao[5,casual[3]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[3],],.95,method="HDI")$CI_high,
         true_betas[casual[4]],ci(lista_bayesian_mode[[5]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[5]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[5]][,casual[4]],.95,method="HDI")$CI_high,expected_value_lasso[5,casual[4]],expected_SS_George_McCulloch[5,casual[4]],expected_SS_Ishwaran_Rao[5,casual[4]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[4],],.95,method="HDI")$CI_high,
         true_betas[casual[5]],ci(lista_bayesian_mode[[5]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[5]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[5]][,casual[5]],.95,method="HDI")$CI_high,expected_value_lasso[5,casual[5]],expected_SS_George_McCulloch[5,casual[5]],expected_SS_Ishwaran_Rao[5,casual[5]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[5],],.95,method="HDI")$CI_high,
         true_betas[casual[6]],ci(lista_bayesian_mode[[5]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[5]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[5]][,casual[6]],.95,method="HDI")$CI_high,expected_value_lasso[5,casual[6]],expected_SS_George_McCulloch[5,casual[6]],expected_SS_Ishwaran_Rao[5,casual[6]],ci(ssl_lasso_unknown_variance[[5]]$beta[casual[6],],.95,method="HDI")$CI_high)

data_graph$c_low=c_low   
data_graph$c_high=c_high  

levels(data_graph$beta_graph)
data_graph$beta_graph=factor(data_graph$beta_graph, levels = c("beta[10]", "beta[14]", "beta[19]", "beta[28]", "beta[41]", "beta[43]"))
ggplot(data_graph, aes(colour = estimates.and.true.value))+ geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+ 
  geom_linerange(aes(x = beta_graph, ymin = c_low,
                     ymax =c_high),
                 lwd = 1, position = position_dodge(width = 1/2))+ 
  geom_pointrange(aes(x = beta_graph, y = values, ymin = c_low ,
                      ymax =c_high ),
                  lwd = 1/2, position = position_dodge(width = 1/2),
                  shape = 21, fill = "WHITE") + coord_flip() + labs(title="Comparison of Bayesian Lasso estimates, Frequentist Lasso estimates, George and McCulloch Spike and Slab, \n Ishwaran and Rao Spike and Slab, Spike-and- Slab Lasso and true coefficients.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 30 observations.") +xlab("Coefficients")+ylab("Values")+labs(colour="Estimates and true values")+
  theme_bw()+theme(plot.title=element_text(size=10, hjust=0.5,face="bold"),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))

####20 obs.

set.seed(1)
casual=sample(1:50,6,rep=F)
casual=casual[order(casual)]
casual
beta_graph=matrix(0,48,1)
beta_graph[1:8]=rep("beta[10]",8)
beta_graph[9:16]=rep("beta[14]",8)
beta_graph[17:24]=rep("beta[19]",8)
beta_graph[25:32]=rep("beta[28]",8)
beta_graph[33:40]=rep("beta[41]",8)
beta_graph[41:48]=rep("beta[43]",8)

value=c(true_betas[casual[1]],expected_value_bayesian_mode[6,casual[1]],expected_value_bayesian_median[6,casual[1]],expected_value_bayesian_mean[6,casual[1]],expected_value_lasso[6,casual[1]],expected_SS_George_McCulloch[6,casual[1]],expected_SS_Ishwaran_Rao[6,casual[1]],apply(ssl_lasso_unknown_variance[[6]]$beta,1,getmode)[casual[1]],
        true_betas[casual[2]],expected_value_bayesian_mode[6,casual[2]],expected_value_bayesian_median[6,casual[2]],expected_value_bayesian_mean[6,casual[2]],expected_value_lasso[6,casual[2]],expected_SS_George_McCulloch[6,casual[2]],expected_SS_Ishwaran_Rao[6,casual[2]],apply(ssl_lasso_unknown_variance[[6]]$beta,1,getmode)[casual[2]],
        true_betas[casual[3]],expected_value_bayesian_mode[6,casual[3]],expected_value_bayesian_median[6,casual[3]],expected_value_bayesian_mean[6,casual[3]],expected_value_lasso[6,casual[3]],expected_SS_George_McCulloch[6,casual[3]],expected_SS_Ishwaran_Rao[6,casual[3]],apply(ssl_lasso_unknown_variance[[6]]$beta,1,getmode)[casual[3]],
        true_betas[casual[4]],expected_value_bayesian_mode[6,casual[4]],expected_value_bayesian_median[6,casual[4]],expected_value_bayesian_mean[6,casual[4]],expected_value_lasso[6,casual[4]],expected_SS_George_McCulloch[6,casual[4]],expected_SS_Ishwaran_Rao[6,casual[4]],apply(ssl_lasso_unknown_variance[[6]]$beta,1,getmode)[casual[4]],
        true_betas[casual[5]],expected_value_bayesian_mode[6,casual[5]],expected_value_bayesian_median[6,casual[5]],expected_value_bayesian_mean[6,casual[5]],expected_value_lasso[6,casual[5]],expected_SS_George_McCulloch[6,casual[5]],expected_SS_Ishwaran_Rao[6,casual[5]],apply(ssl_lasso_unknown_variance[[6]]$beta,1,getmode)[casual[5]],
        true_betas[casual[6]],expected_value_bayesian_mode[6,casual[6]],expected_value_bayesian_median[6,casual[6]],expected_value_bayesian_mean[6,casual[6]],expected_value_lasso[6,casual[6]],expected_SS_George_McCulloch[6,casual[6]],expected_SS_Ishwaran_Rao[6,casual[6]],apply(ssl_lasso_unknown_variance[[6]]$beta,1,getmode)[casual[6]])


statistics=matrix(0,48,1)
statistics[c(1,9,17,25,33,41),]="True_betas"
statistics[c(2,10,18,26,34,42),]="Bayesian Lasso mode"
statistics[c(3,11,19,27,35,43),]="Bayesian Lasso median"
statistics[c(4,12,20,28,36,44),]="Bayesian Lasso mean"
statistics[c(5,13,21,29,37,45),]="Frequentist Lasso"
statistics[c(6,14,22,30,38,46),]="George and McCulloch SS"
statistics[c(7,15,23,31,39,47),]="Ishwaran and Rao SS"
statistics[c(8,16,24,32,40,48),]="SSL"

data_graph=data.frame(beta_graph=beta_graph,values=value,"estimates and true value"=statistics)
data_graph$estimates.and.true.value   



c_low=c(true_betas[casual[1]],ci(lista_bayesian_mode[[6]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[6]][,casual[1]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[6]][,casual[1]],.95,method="HDI")$CI_low,expected_value_lasso[6,casual[1]],expected_SS_George_McCulloch[6,casual[1]],expected_SS_Ishwaran_Rao[6,casual[1]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[1],],.95,method="HDI")$CI_low,
        true_betas[casual[2]],ci(lista_bayesian_mode[[6]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[6]][,casual[2]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[6]][,casual[2]],.95,method="HDI")$CI_low,expected_value_lasso[6,casual[2]],expected_SS_George_McCulloch[6,casual[2]],expected_SS_Ishwaran_Rao[6,casual[2]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[2],],.95,method="HDI")$CI_low,
        true_betas[casual[3]],ci(lista_bayesian_mode[[6]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[6]][,casual[3]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[6]][,casual[3]],.95,method="HDI")$CI_low,expected_value_lasso[6,casual[3]],expected_SS_George_McCulloch[6,casual[3]],expected_SS_Ishwaran_Rao[6,casual[3]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[3],],.95,method="HDI")$CI_low,
        true_betas[casual[4]],ci(lista_bayesian_mode[[6]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[6]][,casual[4]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[6]][,casual[4]],.95,method="HDI")$CI_low,expected_value_lasso[6,casual[4]],expected_SS_George_McCulloch[6,casual[4]],expected_SS_Ishwaran_Rao[6,casual[4]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[4],],.95,method="HDI")$CI_low,
        true_betas[casual[5]],ci(lista_bayesian_mode[[6]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[6]][,casual[5]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[6]][,casual[5]],.95,method="HDI")$CI_low,expected_value_lasso[6,casual[5]],expected_SS_George_McCulloch[6,casual[5]],expected_SS_Ishwaran_Rao[6,casual[5]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[5],],.95,method="HDI")$CI_low,
        true_betas[casual[6]],ci(lista_bayesian_mode[[6]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_median[[6]][,casual[6]],.95,method="HDI")$CI_low,ci(lista_bayesian_mean[[6]][,casual[6]],.95,method="HDI")$CI_low,expected_value_lasso[6,casual[6]],expected_SS_George_McCulloch[6,casual[6]],expected_SS_Ishwaran_Rao[6,casual[6]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[6],],.95,method="HDI")$CI_low)

c_high=c(true_betas[casual[1]],ci(lista_bayesian_mode[[6]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[6]][,casual[1]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[6]][,casual[1]],.95,method="HDI")$CI_high,expected_value_lasso[6,casual[1]],expected_SS_George_McCulloch[6,casual[1]],expected_SS_Ishwaran_Rao[6,casual[1]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[1],],.95,method="HDI")$CI_high,
         true_betas[casual[2]],ci(lista_bayesian_mode[[6]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[6]][,casual[2]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[6]][,casual[2]],.95,method="HDI")$CI_high,expected_value_lasso[6,casual[2]],expected_SS_George_McCulloch[6,casual[2]],expected_SS_Ishwaran_Rao[6,casual[2]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[2],],.95,method="HDI")$CI_high,
         true_betas[casual[3]],ci(lista_bayesian_mode[[6]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[6]][,casual[3]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[6]][,casual[3]],.95,method="HDI")$CI_high,expected_value_lasso[6,casual[3]],expected_SS_George_McCulloch[6,casual[3]],expected_SS_Ishwaran_Rao[6,casual[3]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[3],],.95,method="HDI")$CI_high,
         true_betas[casual[4]],ci(lista_bayesian_mode[[6]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[6]][,casual[4]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[6]][,casual[4]],.95,method="HDI")$CI_high,expected_value_lasso[6,casual[4]],expected_SS_George_McCulloch[6,casual[4]],expected_SS_Ishwaran_Rao[6,casual[4]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[4],],.95,method="HDI")$CI_high,
         true_betas[casual[5]],ci(lista_bayesian_mode[[6]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[6]][,casual[5]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[6]][,casual[5]],.95,method="HDI")$CI_high,expected_value_lasso[6,casual[5]],expected_SS_George_McCulloch[6,casual[5]],expected_SS_Ishwaran_Rao[6,casual[5]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[5],],.95,method="HDI")$CI_high,
         true_betas[casual[6]],ci(lista_bayesian_mode[[6]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_median[[6]][,casual[6]],.95,method="HDI")$CI_high,ci(lista_bayesian_mean[[6]][,casual[6]],.95,method="HDI")$CI_high,expected_value_lasso[6,casual[6]],expected_SS_George_McCulloch[6,casual[6]],expected_SS_Ishwaran_Rao[6,casual[6]],ci(ssl_lasso_unknown_variance[[6]]$beta[casual[6],],.95,method="HDI")$CI_high)

data_graph$c_low=c_low   
data_graph$c_high=c_high   

levels(data_graph$beta_graph)
data_graph$beta_graph=factor(data_graph$beta_graph, levels = c("beta[10]", "beta[14]", "beta[19]", "beta[28]", "beta[41]", "beta[43]"))
ggplot(data_graph, aes(colour = estimates.and.true.value))+ geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+ 
  geom_linerange(aes(x = beta_graph, ymin = c_low,
                     ymax =c_high),
                 lwd = 1, position = position_dodge(width = 1/2))+ 
  geom_pointrange(aes(x = beta_graph, y = values, ymin = c_low ,
                      ymax =c_high ),
                  lwd = 1/2, position = position_dodge(width = 1/2),
                  shape = 21, fill = "WHITE") + coord_flip() + labs(title="Comparison of Bayesian Lasso estimates, Frequentist Lasso estimates, George and McCulloch Spike and Slab, \n Ishwaran and Rao Spike and Slab, Spike-and- Slab Lasso and true coefficients.",subtitle="Simulation design based on a sparse context with ten non-zero coefficients out of 50 and 20 observations.") +xlab("Coefficients")+ylab("Values")+labs(colour="Estimates and true values")+
  theme_bw()+theme(plot.title=element_text(size=10, hjust=0.5,face="bold"),plot.subtitle=element_text(size=7, hjust=0.5, face="italic"))




############################################## MSE ##############################################

set.seed(2)
obs

R=200 #number of resamples

mse_lasso_mode=matrix(0,R,5)
mse_lasso_mean=matrix(0,R,5)
mse_lasso_median=matrix(0,R,5)
mse_lasso=matrix(0,R,5)
mse_GMSS=matrix(0,R,5)
mse_IRSS=matrix(0,R,5)
mse_SSL=matrix(0,R,5)
mse_mean=matrix(0,7,5)

rownames(mse_mean)=c("Bayesian Lasso Mode GS","Bayesian Lasso Mean GS","Bayesian Lasso Median GS","Frequentist Lasso","George and McCulloch SS","Ishwaran and Rao SS","SSL")

#n_test=100 #sample size
n_test=100

for(r in 1:R){

  x_test=matrix(mvrnorm(n_test,rep(0,K),Sigma),n_test,K)
  y_test=matrix(0,n_test,1)
  y_test=x_test%*% true_betas+cbind(rnorm(n_test,0,2))

#1000 obbs

y_pred_lasso_bayesian_gibbs_sampling_mode=x_test%*%expected_value_bayesian_mode[1,]
y_pred_lasso_bayesian_gibbs_sampling_mean=x_test%*%expected_value_bayesian_mean[1,]
y_pred_lasso_bayesian_gibbs_sampling_median=x_test%*%expected_value_bayesian_median[1,]
y_pred_lasso=x_test%*%expected_value_lasso[1,]
y_pred_SS_George_McCulloch=x_test%*%expected_SS_George_McCulloch[1,]
y_pred_SS_Ishwaran_Rao=x_test%*%expected_SS_Ishwaran_Rao[1,]
y_pred_ssl=x_test%*%apply(ssl_lasso_unknown_variance[[1]]$beta,1,getmode)

mse_lasso_mode[r,1]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_mode)^2)
mse_lasso_mean[r,1]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_mean)^2)
mse_lasso_median[r,1]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_median)^2)
mse_lasso[r,1]=mean((y_test-y_pred_lasso)^2)
mse_GMSS[r,1]=mean((y_test-y_pred_SS_George_McCulloch)^2)
mse_IRSS[r,1]=mean((y_test-y_pred_SS_Ishwaran_Rao)^2)
mse_SSL[r,1]=mean((y_test-y_pred_ssl)^2)

#100 obs
y_pred_lasso_bayesian_gibbs_sampling_mode=x_test%*%expected_value_bayesian_mode[2,]
y_pred_lasso_bayesian_gibbs_sampling_mean=x_test%*%expected_value_bayesian_mean[2,]
y_pred_lasso_bayesian_gibbs_sampling_median=x_test%*%expected_value_bayesian_median[2,]
y_pred_lasso=x_test%*%expected_value_lasso[2,]
y_pred_SS_George_McCulloch=x_test%*%expected_SS_George_McCulloch[2,]
y_pred_SS_Ishwaran_Rao=x_test%*%expected_SS_Ishwaran_Rao[2,]
y_pred_ssl=x_test%*%apply(ssl_lasso_unknown_variance[[2]]$beta,1,getmode)

mse_lasso_mode[r,2]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_mode)^2)
mse_lasso_mean[r,2]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_mean)^2)
mse_lasso_median[r,2]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_median)^2)
mse_lasso[r,2]=mean((y_test-y_pred_lasso)^2)
mse_GMSS[r,2]=mean((y_test-y_pred_SS_George_McCulloch)^2)
mse_IRSS[r,2]=mean((y_test-y_pred_SS_Ishwaran_Rao)^2)
mse_SSL[r,2]=mean((y_test-y_pred_ssl)^2)

#50 obs
y_pred_lasso_bayesian_gibbs_sampling_mode=x_test%*%expected_value_bayesian_mode[4,]
y_pred_lasso_bayesian_gibbs_sampling_mean=x_test%*%expected_value_bayesian_mean[4,]
y_pred_lasso_bayesian_gibbs_sampling_median=x_test%*%expected_value_bayesian_median[4,]
y_pred_lasso=x_test%*%expected_value_lasso[4,]
y_pred_SS_George_McCulloch=x_test%*%expected_SS_George_McCulloch[4,]
y_pred_SS_Ishwaran_Rao=x_test%*%expected_SS_Ishwaran_Rao[4,]
y_pred_ssl=x_test%*%apply(ssl_lasso_unknown_variance[[4]]$beta,1,getmode)

mse_lasso_mode[r,3]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_mode)^2)
mse_lasso_mean[r,3]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_mean)^2)
mse_lasso_median[r,3]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_median)^2)
mse_lasso[r,3]=mean((y_test-y_pred_lasso)^2)
mse_GMSS[r,3]=mean((y_test-y_pred_SS_George_McCulloch)^2)
mse_IRSS[r,3]=mean((y_test-y_pred_SS_Ishwaran_Rao)^2)
mse_SSL[r,3]=mean((y_test-y_pred_ssl)^2)

#30 obs
y_pred_lasso_bayesian_gibbs_sampling_mode=x_test%*%expected_value_bayesian_mode[5,]
y_pred_lasso_bayesian_gibbs_sampling_mean=x_test%*%expected_value_bayesian_mean[5,]
y_pred_lasso_bayesian_gibbs_sampling_median=x_test%*%expected_value_bayesian_median[5,]
y_pred_lasso=x_test%*%expected_value_lasso[5,]
y_pred_SS_George_McCulloch=x_test%*%expected_SS_George_McCulloch[5,]
y_pred_SS_Ishwaran_Rao=x_test%*%expected_SS_Ishwaran_Rao[5,]
y_pred_ssl=x_test%*%apply(ssl_lasso_unknown_variance[[5]]$beta,1,getmode)

mse_lasso_mode[r,4]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_mode)^2)
mse_lasso_mean[r,4]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_mean)^2)
mse_lasso_median[r,4]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_median)^2)
mse_lasso[r,4]=mean((y_test-y_pred_lasso)^2)
mse_GMSS[r,4]=mean((y_test-y_pred_SS_George_McCulloch)^2)
mse_IRSS[r,4]=mean((y_test-y_pred_SS_Ishwaran_Rao)^2)
mse_SSL[r,4]=mean((y_test-y_pred_ssl)^2)

#20 obs
obs
y_pred_lasso_bayesian_gibbs_sampling_mode=x_test%*%expected_value_bayesian_mode[6,]
y_pred_lasso_bayesian_gibbs_sampling_mean=x_test%*%expected_value_bayesian_mean[6,]
y_pred_lasso_bayesian_gibbs_sampling_median=x_test%*%expected_value_bayesian_median[6,]
y_pred_lasso=x_test%*%expected_value_lasso[6,]
y_pred_SS_George_McCulloch=x_test%*%expected_SS_George_McCulloch[6,]
y_pred_SS_Ishwaran_Rao=x_test%*%expected_SS_Ishwaran_Rao[6,]
y_pred_ssl=x_test%*%apply(ssl_lasso_unknown_variance[[6]]$beta,1,getmode)

mse_lasso_mode[r,5]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_mode)^2)
mse_lasso_mean[r,5]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_mean)^2)
mse_lasso_median[r,5]=mean((y_test-y_pred_lasso_bayesian_gibbs_sampling_median)^2)
mse_lasso[r,5]=mean((y_test-y_pred_lasso)^2)
mse_GMSS[r,5]=mean((y_test-y_pred_SS_George_McCulloch)^2)
mse_IRSS[r,5]=mean((y_test-y_pred_SS_Ishwaran_Rao)^2)
mse_SSL[r,5]=mean((y_test-y_pred_ssl)^2)

}

mse_mean=rbind(apply(mse_lasso_mode,2,mean),apply(mse_lasso_mean,2,mean),apply(mse_lasso_median,2,mean),apply(mse_lasso,2,mean),apply(mse_GMSS,2,mean),apply(mse_IRSS,2,mean),apply(mse_SSL,2,mean))
rownames(mse_mean)=c("Bayesian Lasso Mode GS","Bayesian Lasso Mean GS","Bayesian Lasso Median GS","Frequentist Lasso","George and McCulloch SS","Ishwaran and Rao SS","SSL")
colnames(mse_mean)=c("1000 obs.","100 obs.","50 obs.","30 obs.","20 obs.")
mse_mean
