# Network analysis and outcome prediction 

# setwd("C:/Users/SAMSUNG/Desktop/Research/Computers_In_Biology_And_Medicine/ST000355_COHPLASMA")
#install.packages('glmnet')
library(ggfortify)
library(glmnet)
library(glasso)
library(igraph)
library(factoextra)
set.seed(100)
a= read.csv('OS2.csv', check.names = F)

a[a[,2]=='Hypoxic',2]='Hypoxic'; a[grepl('Normoxic',a[,2]),2]='Normoxic';
a[a=="ND"]=0


forPCA= a[,c(2, (4:dim(a)[2]))]
forPCA_Met= forPCA[,2:(dim(forPCA)[2])]
forPCA_Met  = matrix(as.numeric(as.matrix(  forPCA_Met  )),    # Convert to numeric matrix
       ncol = ncol(forPCA_Met))
pca_res <- prcomp(forPCA_Met, scale. = T)
autoplot(pca_res, data=forPCA, colour = 'Treatment.group')

# Network Analysis 

#1 glasso + community detection 

X= a[,4:dim(a)[2]]
X=matrix(as.numeric(as.matrix(X)),    # Convert to numeric matrix
       ncol = ncol(X))

metnames= colnames(a)[4:dim(a)[2]]

tableX= as.data.frame(X)
s<- var(X)
ai<-glasso(s, rho=0.01)
aa<-glasso(s,rho=0.01, w.init=ai$w, wi.init=ai$wi)

g=graph_from_adjacency_matrix(aa$wi, mode='undirected', diag=F)

V(g)$name=metnames
plot(simplify(g),vertex.label.cex = 0.5, vertex.size=8)

# greedy method (hiearchical, fast method)
c1 = cluster_fast_greedy(simplify(g))

# modularity measure


#2 K-means clustering 

fviz_nbclust(t(X), kmeans, method = "silhouette")
k2 <- kmeans(t(X), centers = 8, nstart = 25)
kmgroup = k2$cluster



# Dissimilarity matrix
d <- dist(X, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
hcgroup = cutree(hc1, k=2)



# Load the packages
library(gtsummary)
library(dplyr)
library(tableone)
# The script

# Baseline Characteristics 
adultX=as.data.frame( X[1:38,] )
fetusX=as.data.frame( X[39:66,])

hypoxicX = as.data.frame(X[a[,2]==1,])
normoxicX = as.data.frame(X[a[,2]==0,])

hi=CreateTableOne(data = adultX); write.csv(print(hi$ContTable),'BC_adultX.csv')
hi2=CreateTableOne(data = fetusX); write.csv(print(hi2$ContTable),'BC_fetusX.csv')
hi3=CreateTableOne(data = hypoxicX); write.csv(print(hi3$ContTable),'BC_hypoxicX.csv')
hi4=CreateTableOne(data = normoxicX); write.csv(print(hi4$ContTable),'BC_normoxicX.csv')


tableX %>% 
  CreateTableOne(data =tableX) %>% 
  kableone()



set.seed(100)
coefEN75 = rep(0,71)
errEN75 = 0 
fEN75=0

for(i in 1:100){

  setindex= sample(1:dim(X)[1], replace = F)
  tempX= X[setindex,]
trainX=(tempX[c(1:53),])
y= a[,2] 
y= y[setindex]
trainy=y[1:53]

#p.fac <- rep(1, dim(X)[2])
#p.fac[c(dim(X)[2]-1, dim(X)[2])] 
pfit <- cv.glmnet(trainX, trainy, family='binomial',
                  type.measure = 'class', alpha=0.75)
plot(pfit, label = TRUE)

pfit$lambda.min
print(i)
print(which(pfit$glmnet.fit$lambda==pfit$lambda.min))

predy=predict(pfit, newx = X[54:66, ], s = "lambda.min")
predy=(predy>0)+0
obsvy=(y[54:66]=='Hypoxic')

print(as.numeric(coef(pfit)))
coefEN75 = coefEN75+as.numeric(coef(pfit))
errEN75 = errEN75+ sum(predy==obsvy)/13

tempF= predy+obsvy
fEN75 = fEN75 + 2*sum(tempF==2)/(2*sum(tempF==2)+ sum(tempF==1))

}


## Lasso 


set.seed(100)
coefLA = rep(0,71)
errLA = 0 
fLA=0

for(i in 1:100){
  
  setindex= sample(1:dim(X)[1], replace = F)
  tempX= X[setindex,]
  trainX=(tempX[c(1:53),])
  y= a[,2] 
  y= y[setindex]
  trainy=y[1:53]
  
  #p.fac <- rep(1, dim(X)[2])
  #p.fac[c(dim(X)[2]-1, dim(X)[2])] 
  pfit <- cv.glmnet(trainX, trainy, family='binomial',
                    type.measure = 'class', alpha=1)
  plot(pfit, label = TRUE)
  
  pfit$lambda.min
  print(i)
  print(which(pfit$glmnet.fit$lambda==pfit$lambda.min))
  
  predy=predict(pfit, newx = X[54:66, ], s = "lambda.min")
  predy=(predy>0)+0
  obsvy=(y[54:66]=='Hypoxic')
  
  print(as.numeric(coef(pfit)))
  coefLA = coefLA+as.numeric(coef(pfit))
  errLA = errLA+ sum(predy==obsvy)/13
  
  tempF = predy+obsvy 
  fLA= fLA + 2*sum(tempF==2)/(2*sum(tempF==2)+ sum(tempF==1))
  
}


## Elastic Net 0.5

set.seed(100)
coefEN50 = rep(0,71)
errEN50 = 0 
fEN50 = 0 


for(i in 1:100){
  
  setindex= sample(1:dim(X)[1], replace = F)
  tempX= X[setindex,]
  trainX=(tempX[c(1:53),])
  y= a[,2] 
  y= y[setindex]
  trainy=y[1:53]
  
  #p.fac <- rep(1, dim(X)[2])
  #p.fac[c(dim(X)[2]-1, dim(X)[2])] 
  pfit <- cv.glmnet(trainX, trainy, family='binomial',
                    type.measure = 'class', alpha=0.5)
  plot(pfit, label = TRUE)
  
  pfit$lambda.min
  print(i)
  print(which(pfit$glmnet.fit$lambda==pfit$lambda.min))
  
  predy=predict(pfit, newx = X[54:66, ], s = "lambda.min")
  predy=(predy>0)+0
  obsvy=(y[54:66]=='Hypoxic')
  
  print(as.numeric(coef(pfit)))
  coefEN50 = coefEN50+as.numeric(coef(pfit))
  errEN50 = errEN50+ sum(predy==obsvy)/13
  
  tempF = predy+obsvy 
  fEN50= fEN50 + 2*sum(tempF==2)/(2*sum(tempF==2)+ sum(tempF==1))
  
}

#Elastic Net 0.25 


set.seed(100)
coefEN25 = rep(0,71)
errEN25 = 0 
fEN25=0 

for(i in 1:100){
  
  setindex= sample(1:dim(X)[1], replace = F)
  tempX= X[setindex,]
  trainX=(tempX[c(1:53),])
  y= a[,2] 
  y= y[setindex]
  trainy=y[1:53]
  
  #p.fac <- rep(1, dim(X)[2])
  #p.fac[c(dim(X)[2]-1, dim(X)[2])] 
  pfit <- cv.glmnet(trainX, trainy, family='binomial',
                    type.measure = 'class', alpha=0.25)
  plot(pfit, label = TRUE)
  
  pfit$lambda.min
  print(i)
  print(which(pfit$glmnet.fit$lambda==pfit$lambda.min))
  
  predy=predict(pfit, newx = X[54:66, ], s = "lambda.min")
  predy=(predy>0)+0
  obsvy=(y[54:66]=='Hypoxic')
  
  print(as.numeric(coef(pfit)))
  coefEN25 = coefEN25+as.numeric(coef(pfit))
  errEN25 = errEN25+ sum(predy==obsvy)/13
  
  tempF = predy+obsvy 
  fEN25= fEN25 + 2*sum(tempF==2)/(2*sum(tempF==2)+ sum(tempF==1))
  
  
}


# Ridge 



set.seed(100)
coefRI = rep(0,71)
errRI = 0 
fRI =0 

for(i in 1:100){
  
  setindex= sample(1:dim(X)[1], replace = F)
  tempX= X[setindex,]
  trainX=(tempX[c(1:53),])
  y= a[,2] 
  y= y[setindex]
  trainy=y[1:53]
  
  #p.fac <- rep(1, dim(X)[2])
  #p.fac[c(dim(X)[2]-1, dim(X)[2])] 
  pfit <- cv.glmnet(trainX, trainy, family='binomial',
                    type.measure = 'class', alpha=0)
  plot(pfit, label = TRUE)
  
  pfit$lambda.min
  print(i)
  print(which(pfit$glmnet.fit$lambda==pfit$lambda.min))
  
  predy=predict(pfit, newx = X[54:66, ], s = "lambda.min")
  predy=(predy>0)+0
  obsvy=(y[54:66]=='Hypoxic')
  
  print(as.numeric(coef(pfit)))
  coefRI = coefRI+as.numeric(coef(pfit))
  errRI = errRI+ sum(predy==obsvy)/13
  
  tempF = predy+obsvy 
fRI= fRI + 2*sum(tempF==2)/(2*sum(tempF==2)+ sum(tempF==1))
}


metnames[which(coefLA>0) -1 ]

metnames[which(coefLA<0) -1]

coefLAMET= coefLA[-1]
POS=cbind(metnames[which(coefLAMET>0) ], round(coefLAMET[coefLAMET>0],3) )
write.csv(POS, "POS.csv")

NEG=cbind(metnames[which(coefLAMET<0)  ], round(coefLAMET[coefLAMET<0],3))
write.csv(NEG, 'NEG.csv')


