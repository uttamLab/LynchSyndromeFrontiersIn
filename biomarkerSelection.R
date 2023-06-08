#Biomarker Selection 
#Shikhar Uttam
#Dept of Computational and Systems Biology
#UPMC Hillman Cancer Center
#University of Pittsburgh
rm(list=ls());
library(ggplot2);
library(ggpubr)
library(reshape2)
library(Matrix)
library(glmnet)
library(splitstackshape)
library(zoo)
load('./patientData.RData')
#===============
select.func <- function(perc.df,p.select,bs.out){
  perc.df[perc.df <= (1- p.select)] <- 1000
  perc.df[perc.df != 1000] <- 0
  perc.df[perc.df == 1000] <- 1
  select.df <- apply(perc.df,1, sum)/bs.out;
  return(select.df)
}

#trapezoidal quadrature rule based AUC;
trap.quad.auc <- function(x,y){
  id <- order(x);
  return(sum(diff(x[id])*rollmean(y[id],2)))
}
#==========================
#CHOOSE 1,2 OR 3 for sections 3.4, 3.6, and 3.5 respectively in the paper.
cntg <- 3; 
#========================

#penalized logistic regression=======

#outer bootstrap
bs.out <- 100;# .
perc <- 0.70;
#inner bootstrap
bs.in <- 200;#

cv.fitvec.stim.visit1 <- vector("list",length(data.mats.visit1));
beta.coefvec.stim.visit1 <- vector("list",length(data.mats.visit1));
names.coefvec.stim.visit1 <-  vector("list",length(data.mats.visit1));

cv.fitvec.stim.visit2 <- vector("list",length(data.mats.visit2));
beta.coefvec.stim.visit2 <- vector("list",length(data.mats.visit2));
names.coefvec.stim.visit2 <-  vector("list",length(data.mats.visit2)); 


thresh.reject <- 0;
rangevec <- 1:length(data.mats.visit1);
for(k1 in rangevec[cntg]){
  set.seed(2);
  print(k1)
  #visit 1
  cv.fitvec.stim.visit1[[k1]] <- vector("list",bs.out)
  beta.coefvec.stim.visit1[[k1]] <- vector("list",bs.out);
  names.coefvec.stim.visit1[[k1]] <-  vector("list",bs.out);
  data.mat.cv.visit1 <- cbind(data.mats.visit1[[k1]],group=factor.condition.visit1[[k1]]);
  #visit 2
  cv.fitvec.stim.visit2[[k1]] <- vector("list",bs.out)
  beta.coefvec.stim.visit2[[k1]] <- vector("list",bs.out);
  names.coefvec.stim.visit2[[k1]] <-  vector("list",bs.out);
  data.mat.cv.visit2 <- cbind(data.mats.visit2[[k1]],group=factor.condition.visit2[[k1]]);

  #visit 1
  if(length(which(ratio.cv.df[k1,]<0))>thresh.reject){ 
    idx.rem.visit1 <- which(ratio.cv.df[k1,]<thresh.reject)
  }else{
    idx.rem.visit1 <-32
  }
  #visit 2
  if(length(which(ratio.cv.df[k1,]<0))>thresh.reject){ 
    idx.rem.visit2 <- which(ratio.cv.df[k1,]<thresh.reject)
  }else{
    idx.rem.visit2 <-32
  }

  #===========
  data.mat.cv.visit2 <- data.mat.cv.visit2[,-idx.rem.visit2];
  #==========
  
  cnt.visit1 <-0;
  cnt.visit2 <-0;
  for(t1 in 1:bs.out){#loop through outer bootstrap
    time1 <- Sys.time();
    #
    data.mat.cv2.visit1 <- stratified(data.mat.cv.visit1,group='group',size=perc,replace = TRUE)
    print(paste(c('...',t1), collapse=' '));
    cv.fitvec.stim.visit1[[k1]][[t1]] <- vector("list",bs.in)
    beta.coefvec.stim.visit1[[k1]][[t1]] <- vector("list",bs.in)
    names.coefvec.stim.visit1[[k1]][[t1]] <- vector("list",bs.in)
    
    #
    data.mat.cv2.visit2 <- stratified(data.mat.cv.visit2,group='group',size=perc,replace = TRUE)
    print(paste(c('...',t1), collapse=' '));
    cv.fitvec.stim.visit2[[k1]][[t1]] <- vector("list",bs.in)
    beta.coefvec.stim.visit2[[k1]][[t1]] <- vector("list",bs.in)
    names.coefvec.stim.visit2[[k1]][[t1]] <- vector("list",bs.in)
    

    for(i1 in 1:bs.in){#loop through inner bootstrap
      print(paste(c('....',i1),collapse=' '))
      data.mat.cv3.visit1 <- stratified(data.mat.cv2.visit1,group='group',size=0.999,replace = TRUE)
      data.mat.cv3.visit2 <- stratified(data.mat.cv2.visit2,group='group',size=0.999,replace = TRUE)
      #initialize: simulation
      cv.fitvec.stim.visit1[[k1]][[t1]][[i1]] <- cv.glmnet(x=as.matrix(data.mat.cv3.visit1[,1:(dim(data.mat.cv3.visit1)[2]-1)]),y=data.mat.cv3.visit1$group,family="binomial",
                                                    maxit=100000,alpha=0.5,nlambda = 1000,type.measure = "class",nfolds = dim(data.mat.cv3.visit1)[1],keep=TRUE,
                                                    grouped= FALSE,standardize=TRUE);#dim(data.mat.cv)[1]
      
      devpos.visit1 <- which(cv.fitvec.stim.visit1[[k1]][[t1]][[i1]]$lambda.min == cv.fitvec.stim.visit1[[k1]][[t1]][[i1]]$lambda);#
      devpos2.visit1 <- which.min(cv.fitvec.stim.visit1[[k1]][[t1]][[i1]]$cvm);
      
      if(devpos.visit1 != devpos2.visit1){
        cnt.visit1 <- 1;
        break
      }
      
      betacoeftmp.visit1 <- cv.fitvec.stim.visit1[[k1]][[t1]][[i1]]$glmnet.fit$beta[,devpos.visit1]
      beta.coefvec.stim.visit1[[k1]][[t1]][[i1]] <- betacoeftmp.visit1;
      names.coefvec.stim.visit1[[k1]][[t1]][[i1]] <- names(which(abs(betacoeftmp.visit1) > 0));
      

        cv.fitvec.stim.visit2[[k1]][[t1]][[i1]] <- cv.glmnet(x=as.matrix(data.mat.cv3.visit2[,1:(dim(data.mat.cv3.visit2)[2]-1)]),y=data.mat.cv3.visit2$group,family="binomial",
                                                    maxit=100000,alpha=0.5,nlambda = 1000,type.measure = "class",nfolds = dim(data.mat.cv3.visit2)[1],keep=TRUE,
                                                    grouped= FALSE,standardize=TRUE);#dim(data.mat.cv)[1]
        
        devpos.visit2 <- which(cv.fitvec.stim.visit2[[k1]][[t1]][[i1]]$lambda.min == cv.fitvec.stim.visit2[[k1]][[t1]][[i1]]$lambda);#
        devpos2.visit2 <- which.min(cv.fitvec.stim.visit2[[k1]][[t1]][[i1]]$cvm);
        
        
        if(devpos.visit2 != devpos2.visit2){
          cnt.visit2 <- 1;
          break
        }
        
        betacoeftmp.visit2 <- cv.fitvec.stim.visit2[[k1]][[t1]][[i1]]$glmnet.fit$beta[,devpos.visit2]
        beta.coefvec.stim.visit2[[k1]][[t1]][[i1]] <- betacoeftmp.visit2;
        names.coefvec.stim.visit2[[k1]][[t1]][[i1]] <- names(which(abs(betacoeftmp.visit2) > 0));
    }
    if(cnt.visit1){
      break;
    }
    if(cnt.visit2){
      break;
    }
    time2<-Sys.time();
    print(time2-time1)
  }
  if(cnt.visit1){
    break
  }
  if(cnt.visit2){
    break
  }
}

##biomarker selection probability
#visit 1
alpha.conf.vec.visit1 <- seq(0.5,1,0.01); #initialize
betacoef.bs.visit1 <- vector("list",length(data.mats.visit1));
conf.type.visit1 <- vector("list",length(data.mats.visit1));
bm.names.visit1 <- vector("list",length(data.mats.visit1));
for(i1 in rangevec[cntg]){
  betacoef.bs.visit1[[i1]] <- vector("list",bs.out);
  conf.type.visit1[[i1]] <- vector("list",bs.out);
  bm.names.visit1[[i1]] <- vector("list",bs.out);
  for(k1 in 1:bs.out){#loop through outer bootstraps
    betacoef.bs.visit1[[i1]][[k1]] <- sapply(beta.coefvec.stim.visit1[[i1]][[k1]], function(x) x);
    conf.type.visit1[[i1]][[k1]] <- vector("numeric",dim(betacoef.bs.visit1[[i1]][[k1]])[1]) #initialize
    brks <- seq(-.Machine$double.eps,.Machine$double.eps,by=.Machine$double.eps);
    for(t1 in 1:length(conf.type.visit1[[i1]][[k1]])){#looping through chemokines/cytokines
      cdf.type.visit1 <- ecdf(betacoef.bs.visit1[[i1]][[k1]][t1,]);
      cdf.specific.visit1 <- cdf.type.visit1(brks);
      conf.type.visit1[[i1]][[k1]][t1] <- cdf.specific.visit1[length(cdf.specific.visit1)]-cdf.specific.visit1[1] 
    }
    names(conf.type.visit1[[i1]][[k1]]) <- rownames(betacoef.bs.visit1[[i1]][[k1]])
    bm.names.visit1[[i1]][[k1]] <- vector("list",length(alpha.conf.vec.visit1));#initialize
    for(j1 in 1:length(alpha.conf.vec.visit1)){
      bm.names.visit1[[i1]][[k1]][[j1]] <- names(conf.type.visit1[[i1]][[k1]])[conf.type.visit1[[i1]][[k1]] <= (1 - alpha.conf.vec.visit1[j1])];
    }
  }
}
names(betacoef.bs.visit1)<-names(data.mats.visit1);
names(bm.names.visit1)<-names(data.mats.visit1);
names(conf.type.visit1)<-names(data.mats.visit1);

#visit 2
alpha.conf.vec.visit2 <- seq(0.5,1,0.01); #initialize
betacoef.bs.visit2 <- vector("list",length(data.mats.visit2));
conf.type.visit2 <- vector("list",length(data.mats.visit2));
bm.names.visit2 <- vector("list",length(data.mats.visit2));
for(i1 in rangevec[cntg]){
  betacoef.bs.visit2[[i1]] <- vector("list",bs.out);
  conf.type.visit2[[i1]] <- vector("list",bs.out);
  bm.names.visit2[[i1]] <- vector("list",bs.out);
  for(k1 in 1:bs.out){#loop through outer bootstraps
    betacoef.bs.visit2[[i1]][[k1]] <- sapply(beta.coefvec.stim.visit2[[i1]][[k1]], function(x) x);
    conf.type.visit2[[i1]][[k1]] <- vector("numeric",dim(betacoef.bs.visit2[[i1]][[k1]])[1]) #initialize
    brks <- seq(-.Machine$double.eps,.Machine$double.eps,by=.Machine$double.eps);
    for(t1 in 1:length(conf.type.visit2[[i1]][[k1]])){#looping through chemokines/cytokines
      cdf.type.visit2 <- ecdf(betacoef.bs.visit2[[i1]][[k1]][t1,]);
      cdf.specific.visit2 <- cdf.type.visit2(brks);
      conf.type.visit2[[i1]][[k1]][t1] <- cdf.specific.visit2[length(cdf.specific.visit2)]-cdf.specific.visit2[1] 
    }
    names(conf.type.visit2[[i1]][[k1]]) <- rownames(betacoef.bs.visit2[[i1]][[k1]])
    bm.names.visit2[[i1]][[k1]] <- vector("list",length(alpha.conf.vec.visit2));#initialize
    for(j1 in 1:length(alpha.conf.vec.visit2)){
      bm.names.visit2[[i1]][[k1]][[j1]] <- names(conf.type.visit2[[i1]][[k1]])[conf.type.visit2[[i1]][[k1]] <= (1 - alpha.conf.vec.visit2[j1])];
    }
  }
}
names(betacoef.bs.visit2)<-names(data.mats.visit2);
names(bm.names.visit2)<-names(data.mats.visit2);
names(conf.type.visit2)<-names(data.mats.visit2);


#====Overall selection probability
#visit 1
conf.df.visit1 <- sapply(conf.type.visit1[[cntg]], function(x) x)
select.prob.df.visit1 <- matrix(0,nrow=dim(conf.df.visit1)[1],ncol=length(alpha.conf.vec.visit1))
for(i1 in 1:length(alpha.conf.vec.visit1)){
  select.prob.df.visit1[,i1] <- select.func(conf.df.visit1, alpha.conf.vec.visit1[i1],bs.out);
}
rownames(select.prob.df.visit1) <- rownames(conf.df.visit1)

par(mar=c(5, 4, 4, 8), xpd=TRUE)
colorvec <- sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],size=length(alpha.conf.vec.visit1),replace = FALSE)
plot(alpha.conf.vec.visit1,select.prob.df.visit1[1,],col=colorvec[1],type='l',lwd=3,ylim=c(0,1),ylab='Selection probability',xlab='conf threshold',main='Visit 1');
for(i1 in 2:dim(select.prob.df.visit1)[1]){
  lines(alpha.conf.vec.visit1,select.prob.df.visit1[i1,],col=colorvec[i1],type='l',lwd=3);
}
lines(alpha.conf.vec.visit1,1-(alpha.conf.vec.visit1-0.5)*2,col='black',lwd=3,lty=2)
legend("topright", inset=c(-0.3, 0), legend=paste(rownames(select.prob.df.visit1), sep = ','), lwd=2,col=colorvec)

auc.vec.visit1 <- apply(select.prob.df.visit1,1,function(x) trap.quad.auc(alpha.conf.vec.visit1,x))/0.5
names(auc.vec.visit1) <- rownames(select.prob.df.visit1)

#visit 2
conf.df.visit2 <- sapply(conf.type.visit2[[cntg]], function(x) x)
select.prob.df.visit2 <- matrix(0,nrow=dim(conf.df.visit2)[1],ncol=length(alpha.conf.vec.visit2))
for(i1 in 1:length(alpha.conf.vec.visit2)){
  select.prob.df.visit2[,i1] <- select.func(conf.df.visit2, alpha.conf.vec.visit2[i1],bs.out);
}
rownames(select.prob.df.visit2) <- rownames(conf.df.visit2)


par(mar=c(5, 4, 4, 8), xpd=TRUE)

plot(alpha.conf.vec.visit2,select.prob.df.visit2[1,],col=colorvec[1],type='l',lwd=3,xlim=c(0.5,1), ylim=c(0,1),ylab='Selection probability',xlab='conf threshold', main='Visit 2');
for(i1 in 2:dim(select.prob.df.visit2)[1]){
  lines(alpha.conf.vec.visit2,select.prob.df.visit2[i1,],col=colorvec[i1],type='l',lwd=3);
}
lines(alpha.conf.vec.visit2,1-(alpha.conf.vec.visit2-0.5)*2,col='black',lwd=3,lty=2)
legend("topright", inset=c(-0.3, 0), legend=paste(rownames(select.prob.df.visit2), sep = ','), lwd=2,col=colorvec)

auc.vec.visit2 <- apply(select.prob.df.visit2,1,function(x) trap.quad.auc(alpha.conf.vec.visit2,x))/0.5
names(auc.vec.visit2) <- rownames(select.prob.df.visit2)


#barplot
ttl <- '';
df.auc.visit1.visit2 <- data.frame(visit1 = auc.vec.visit1, visit2 = auc.vec.visit2, analyte=names(auc.vec.visit1))
df.auc.visit1.visit2.m <- melt(df.auc.visit1.visit2)
colorvec <- c('orange', 'red');
df.auc.visit1.visit2.m$variable <- as.factor(df.auc.visit1.visit2.m$variable)
colnames(df.auc.visit1.visit2.m) <- c('Analyte','Visit','Selection.Probability')
#df.auc.visit1.visit2.m <- cbind(df.auc.visit1.visit2.m, color=colorvec[df.auc.visit1.visit2.m$Visit])
p <- ggplot(df.auc.visit1.visit2.m, aes(x=Analyte,y=Selection.Probability,fill=Visit))+
  geom_bar(position = 'dodge',stat='identity')+ 
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='bold'), 
        axis.text.y = element_text(face='bold'))+
  scale_fill_manual(values=colorvec)+
  geom_hline(yintercept=0.5, linetype="dashed", color = "black",lwd=1)+
  labs(title=ttl)+
  ylab('Overall Selection Probability')+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
print(p)



