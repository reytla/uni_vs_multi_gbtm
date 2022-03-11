### Simulate according to this perhaps:
### https://blogs.sas.com/content/iml/2017/09/13/simulate-clusters-multivariate-sas.html

### Simulate according to this perhaps:
### https://blogs.sas.com/content/iml/2017/09/13/simulate-clusters-multivariate-sas.html


### Load required packages
#.libPaths("D:/R_library/4.0")
#.libPaths( c( "~/userLibrary" , .libPaths() ) )
# Package names
# install.packages("MplusAutomation",dependencies = TRUE,Ncpus=4)

#install.packages("tidyverse", Ncpus = 6)

packages <- c("MplusAutomation","MASS","mvrnorm")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages],Ncpus=4)
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

### Check if required fonts are loaded and load if necessary
#font_import()
loadfonts()
extrafont::loadfonts(device="win")

truedgp <- "Multivariate 3 class"
## Set working directory

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

### Defining population parameters

model = "GBTM"                      # Population model
n=200                               # Number of replications (i.e. independent data sets)
# Change
N= 1000                             # Number of subjects per data set
m = 5                               # Number of repeated measures
start = 0                           # Starting time of observations i.e. t=0
Nt = N*m                            # Dimension of data frame - i.e. X matrix
rp = 2                              # Number of random effects
k = 3                               # Number of classes in outcome Y1
l = 3                               # Number of classes in outcome Y2
cond ="3 class univariate"
# Change
pr <- rep(1/k,k) #c(.5,.5)   # Class sizes desired: 30%,20%,30%,20%
#
# Creating quadratic fixed trend and linear RE

ltime = cbind(seq(0,7,length=m))             # Time points of repeated measurements
ltime2 = (ltime)**2       # Second order polynomial time points
int <- cbind(rep(1,length.out=m))               # Generate intercept

fdm <- cbind(int,ltime,ltime2)                  # Fixed design matrix - quadratic
rdm <- cbind(int,ltime)                         # Random design matrix - linear

# Parameters of individual classes
# Informed by class separation measure - keep total variation constant between simulations
# Ratio between RS and RQ set to 0.8

# Fixed effects mean matrix 
# Fixed effects mean matrix 

# Change
bmeanmat =matrix(4*c(0,.438,-0.035,
                     0,0,0,
                     0,-.438,0.035),nrow=3,ncol=3)

# class separation from GBTM2

# Change
bmeanmat2 =matrix(4*c(.214,.125,.01,
                      0,-.027,0.002,
                      1.688,.125,-.01),nrow=3,ncol=3)

# Giving column names
colnames(bmeanmat)<-c("k=1","k=2","k=3")
rownames(bmeanmat)<-c("B0","B1","B2")

colnames(bmeanmat2)<-c("k=1","k=2","k=3")
rownames(bmeanmat2)<-c("B0","B1","B2")


### Mean structure
# Fixed effects per latent class i.e. XB
fe1 <- fdm%*%bmeanmat[,1]
fe2 <- fdm%*%bmeanmat[,2]
fe3 <- fdm%*%bmeanmat[,3]


fe1_2 <- fdm%*%bmeanmat2[,1]
fe2_2 <- fdm%*%bmeanmat2[,2]
fe3_2 <- fdm%*%bmeanmat2[,3]

# Residual error matrix
residualmean <- rep(0,length.out=m)

### Covariance structure
# Here assumed diagonal and fixed over time - can modify

### Outcome 1 

vary1 <-3.125 
vary2 <-3.16

sig2I1_1 <-rep(vary1,m)*diag(m)
sig2I1_2 <-rep(vary1,m)*diag(m)
sig2I1_3 <-rep(vary1,m)*diag(m)

sig2I1<-rep(vary1,m)*diag(m)

### Outcome 2

sig2I2_1 <-rep(vary2,m)*diag(m)
sig2I2_2 <-rep(vary2,m)*diag(m)
sig2I2_3 <-rep(vary2,m)*diag(m)

sig2I2<- rep(vary2,m)*diag(m)


################################################################
################################################################
#### /CLASS SEPARATION CALCULATIONS 
################################################################
################################################################

#### Individual point overlap
## NB takes standard deviation
int_f <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  pmin(f1, f2)
}

### Outcome Y1 - CLASSES
sdev1=sqrt(diag(sig2I1_1))
sdev2=sqrt(diag(sig2I1_2))
sdev3 =sqrt(diag(sig2I1_3))

### Class 1 and 2
a0 <- integrate(int_f, -Inf, Inf, mu1=fe1[1], mu2=fe2[1], sd1=sdev1[1], sd2=sdev2[1])$value
a1 <- integrate(int_f, -Inf, Inf, mu1=fe1[2], mu2=fe2[2], sd1=sdev1[2], sd2=sdev2[2])$value
a2 <- integrate(int_f, -Inf, Inf, mu1=fe1[3], mu2=fe2[3], sd1=sdev1[3], sd2=sdev2[3])$value
a3 <- integrate(int_f, -Inf, Inf, mu1=fe1[4], mu2=fe2[4], sd1=sdev1[4], sd2=sdev2[4])$value
a4 <- integrate(int_f, -Inf, Inf, mu1=fe1[5], mu2=fe2[5], sd1=sdev1[5], sd2=sdev2[5])$value

a<-mean(c(a0,a1,a2,a3,a4))

### Class 1 and 3

b0 <- integrate(int_f, -Inf, Inf, mu1=fe1[1], mu2=fe3[1], sd1=sdev1[1], sd2=sdev3[1])$value
b1 <- integrate(int_f, -Inf, Inf, mu1=fe1[2], mu2=fe3[2], sd1=sdev1[2], sd2=sdev3[2])$value
b2 <- integrate(int_f, -Inf, Inf, mu1=fe1[3], mu2=fe3[3], sd1=sdev1[3], sd2=sdev3[3])$value
b3 <- integrate(int_f, -Inf, Inf, mu1=fe1[4], mu2=fe3[4], sd1=sdev1[4], sd2=sdev3[4])$value
b4 <- integrate(int_f, -Inf, Inf, mu1=fe1[5], mu2=fe3[5], sd1=sdev1[5], sd2=sdev3[5])$value

b<- mean(c(b0,b1,b2,b3,b4))

### Class 2 and 3

c0 <- integrate(int_f, -Inf, Inf, mu1=fe2[1], mu2=fe3[1], sd1=sdev2[1], sd2=sdev3[1])$value
c1 <- integrate(int_f, -Inf, Inf, mu1=fe2[2], mu2=fe3[2], sd1=sdev2[2], sd2=sdev3[2])$value
c2 <- integrate(int_f, -Inf, Inf, mu1=fe2[3], mu2=fe3[3], sd1=sdev2[3], sd2=sdev3[3])$value
c3 <- integrate(int_f, -Inf, Inf, mu1=fe2[4], mu2=fe3[4], sd1=sdev2[4], sd2=sdev3[4])$value
c4 <- integrate(int_f, -Inf, Inf, mu1=fe2[5], mu2=fe3[5], sd1=sdev2[5], sd2=sdev3[5])$value

c<- mean(c(c0,c1,c2,c3,c4))



OVL<-matrix(round(c(a,1,b,c),3),nrow=2,byrow=TRUE)
colnames(OVL)<-c("k=1","k=2")
rownames(OVL)<-c("k=2","k=3")
OVL

### Outcome Y2
sdev1a=sqrt(diag(sig2I2_1))
sdev2a=sqrt(diag(sig2I2_2))
sdev3a=sqrt(diag(sig2I2_3))

### Class 1 and 2
aa0 <- integrate(int_f, -Inf, Inf, mu1=fe1_2[1], mu2=fe2_2[1], sd1=sdev1a[1], sd2=sdev2a[1])$value
aa1 <- integrate(int_f, -Inf, Inf, mu1=fe1_2[2], mu2=fe2_2[2], sd1=sdev1a[2], sd2=sdev2a[2])$value
aa2 <- integrate(int_f, -Inf, Inf, mu1=fe1_2[3], mu2=fe2_2[3], sd1=sdev1a[3], sd2=sdev2a[3])$value
aa3 <- integrate(int_f, -Inf, Inf, mu1=fe1_2[4], mu2=fe2_2[4], sd1=sdev1a[4], sd2=sdev2a[4])$value
aa4 <- integrate(int_f, -Inf, Inf, mu1=fe1_2[5], mu2=fe2_2[5], sd1=sdev1a[5], sd2=sdev2a[5])$value

aa<-mean(c(aa0,aa1,aa2,aa3,aa4))

### Class 1 and 3

bb0 <- integrate(int_f, -Inf, Inf, mu1=fe1_2[1], mu2=fe3_2[1], sd1=sdev1a[1], sd2=sdev3a[1])$value
bb1 <- integrate(int_f, -Inf, Inf, mu1=fe1_2[2], mu2=fe3_2[2], sd1=sdev1a[2], sd2=sdev3a[2])$value
bb2 <- integrate(int_f, -Inf, Inf, mu1=fe1_2[3], mu2=fe3_2[3], sd1=sdev1a[3], sd2=sdev3a[3])$value
bb3 <- integrate(int_f, -Inf, Inf, mu1=fe1_2[4], mu2=fe3_2[4], sd1=sdev1a[4], sd2=sdev3a[4])$value
bb4 <- integrate(int_f, -Inf, Inf, mu1=fe1_2[5], mu2=fe3_2[5], sd1=sdev1a[5], sd2=sdev3a[5])$value

bb<- mean(c(bb0,bb1,bb2,bb3,bb4))

### Class 2 and 3

cc0 <- integrate(int_f, -Inf, Inf, mu1=fe2_2[1], mu2=fe3_2[1], sd1=sdev2a[1], sd2=sdev3a[1])$value
cc1 <- integrate(int_f, -Inf, Inf, mu1=fe2_2[2], mu2=fe3_2[2], sd1=sdev2a[2], sd2=sdev3a[2])$value
cc2 <- integrate(int_f, -Inf, Inf, mu1=fe2_2[3], mu2=fe3_2[3], sd1=sdev2a[3], sd2=sdev3a[3])$value
cc3 <- integrate(int_f, -Inf, Inf, mu1=fe2_2[4], mu2=fe3_2[4], sd1=sdev2a[4], sd2=sdev3a[4])$value
cc4 <- integrate(int_f, -Inf, Inf, mu1=fe2_2[5], mu2=fe3_2[5], sd1=sdev2a[5], sd2=sdev3a[5])$value

cc<- mean(c(cc0,cc1,cc2,cc3,cc4))

OVLa<-matrix(round(c(aa,1,bb,cc),3),nrow=2,byrow=TRUE)
colnames(OVLa)<-c("k=1","k=2")
rownames(OVLa)<-c("k=2","k=3")
OVLa

### Multivariate OVL

#### Create Multivariate Sigma
# Create correlation of 0.5 between contemporaneous outcomes cov(Y1,Y2)=p*sigmaY1*sigmaY2

### Outcome 1



q <- diag((sig2I1_1))     # C1 Y var
s <- diag((sig2I2_1))     # C1 Z var
# Change
r <- rep(0.5,m)         # corr between Y and Z

q2 <- diag((sig2I1_2))     # C1 Y var
s2 <- diag((sig2I2_2))     # C1 Z var
# Change
r2 <- rep(0.5,m)         # corr between Y and Z

q3 <- diag((sig2I1_3))     # C1 Y var
s3 <- diag((sig2I2_3))     # C1 Z var
# Change
r3 <- rep(0.5,m)         # corr between Y and Z

covc1 <- r*matrix(c(diag(c(sqrt(q[1])*sqrt(s[1]),sqrt(q[2])*sqrt(s[2]),sqrt(q[3])*sqrt(s[3]),sqrt(q[4])*sqrt(s[4]),sqrt(q[5])*sqrt(s[5])))),nrow=5)
covc2 <- matrix(c(r2*diag(c(sqrt(q2[1])*sqrt(s2[1]),sqrt(q2[2])*sqrt(s2[2]),sqrt(q2[3])*sqrt(s2[3]),sqrt(q2[4])*sqrt(s2[4]),sqrt(q2[5])*sqrt(s2[5])))),nrow=5)
covc3 <- matrix(c(r3*diag(c(sqrt(q3[1])*sqrt(s3[1]),sqrt(q3[2])*sqrt(s3[2]),sqrt(q3[3])*sqrt(s3[3]),sqrt(q3[4])*sqrt(s3[4]),sqrt(q3[5])*sqrt(s3[5])))),nrow=5)


# Calculate class covariance matrices

sig2I1 <-rbind((cbind(sig2I1_1,covc1)),cbind(covc1,sig2I2_1)) ### Since sqrt(sigmax)*sqrt(sigmax)=sigmax since both equal
sig2I2 <-rbind((cbind(sig2I1_2,covc2)),cbind(covc2,sig2I2_2))
sig2I3 <-rbind((cbind(sig2I1_3,covc3)),cbind(covc3,sig2I2_3))

eigen(sig2I1)                # Check if positive semi-definite (all eigenvalues !=0) to ensure invertibility
eigen(sig2I2)                # Check if positive semi-definite (all eigenvalues !=0) to ensure invertibility
eigen(sig2I3)                # Check if positive semi-definite (all eigenvalues !=0) to ensure invertibility

# ### Create Multivariate Mu
# 
fe11 <- rbind(fe1,fe1_2)
fe22 <- rbind(fe2,fe2_2)
fe33 <- rbind(fe3,fe3_2)

####
id <- 1:N  # subject ids

### Fucntion to Generate Multivariate data

mult_mixture <- function(x){
  set.seed(x)
  samp_size_comp1 <- rmultinom(1,N,prob = pr)                 # Vector of sample sizes for classes in outcome Y1
  y1 <- cbind(mvrnorm(n=samp_size_comp1[1],mu = fe11, Sigma = sig2I1),class_m=c(rep(1,samp_size_comp1[1])))
  y2 <- cbind(mvrnorm(n=samp_size_comp1[2],mu = fe22, Sigma = sig2I2),class_m=c(rep(2,samp_size_comp1[2])))
  y3 <- cbind(mvrnorm(n=samp_size_comp1[3],mu = fe33, Sigma = sig2I3),class_m=c(rep(3,samp_size_comp1[3])))  
  y <- rbind(y1,y2,y3)
  colnames(y)<-c("Y1","Y2","Y3","Y4","Y5","Z1","Z2","Z3","Z4","Z5","class_m")
  mydata <-data.frame(repl=x,y)
  rows <-(sample(nrow(mydata))) # Random sorting
  mydata<- mydata[rows,]
  mydata<- data.frame(mydata,id)
}

rep <-(1:n)

res<-lapply(rep, mult_mixture)  # Generate n independent data sets

#### To Mplus
### Generate n different data sets for Mplus
dir.create("Data")
library(MplusAutomation)
for (i in 1:n){
  prepareMplusData(as.data.frame(res[[i]]), paste0(wd,"/Data/Mplus","_",i,".dat"),
                   overwrite = TRUE)}


### Create models according to Mplus .txt script
### Create models according to Mplus .txt script
### First set working directory again
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
getwd()

filenames <-noquote(list.files(paste0(wd,"/Data")))

write.table(noquote(filenames),file=paste0(wd,"/Data/Data.txt"),col.names = FALSE,row.names=FALSE,sep="\t", quote = FALSE)



filenames <-list.files(paste0(wd,"/Templates"))

for (j in filenames){
  createModels(paste0(wd,"/Templates/",j))
}

#Change number of starts etc
### Run all models within a directory and log output
### Skipping models with existing output files
### replaceOutfile ="modifiedDate"
### Checks whether there is an existing output file for a given input file
### If there is, it checks to see whether the date the input file was modified is newer than output
runModels(wd,recursive=TRUE,showOutput = TRUE,replaceOutfile = "never",
          logFile = paste0(wd,"/ComparisonLog.txt"))


### Extract summaries - this is where fit statistics reside
library(plyr)
allModelParams <- readModels(wd, recursive = TRUE, what="summaries")
justSummaries <- do.call("rbind.fill",sapply(allModelParams,"[","summaries"))

## Remove incomplete runs - relic of Mplus Automation
nrow(justSummaries)
justSummaries <- justSummaries[complete.cases(justSummaries$NGroups),]
nrow(justSummaries)


### custom functions for string extraction

left = function(text, num_char) {
  substr(text, 1, num_char)
}

mid = function(text, start_num, num_char) {
  substr(text, start_num, start_num + num_char - 1)
}

right = function(text, num_char) {
  substr(text, nchar(text) - (num_char-1), nchar(text))
}

library(stringr)
### Extract replication number
justSummaries$rep <- right(justSummaries$Title,3)
table(justSummaries$rep)
justSummaries$rep <-stringr::str_remove_all(justSummaries$rep," ")
table(justSummaries$rep)
View(justSummaries)

#### Keep specific columns

df <- justSummaries[c(2,6,7,8,10,11,12,13,15,16,17,18,19,20,21)]

### Extract model type
df$model <- mid(df$Title,9,10)
### Modify

df$model <- str_remove_all(df$model,"rep")
df$model <- str_remove_all(df$model,"re")
df$model <- str_remove_all(df$model," ")

table(df$model)

### Set model and replication as factors
df$model <-as.factor(df$model)
df$rep <- as.factor(df$rep)

### Want to create estimated class number (K) variable
df$class <-left(df$Title,2)
df$class <- str_remove_all(df$class,"-")

df$class <-as.numeric(df$class)
table(df$class)

### Save Data frame

save(df,file="df.Rda")