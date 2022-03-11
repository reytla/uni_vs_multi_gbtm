packages <- c("MplusAutomation","MASS","stringr")
# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages],Ncpus=4)
# }

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Change
wd <- "/data/volume_8/U_GBTM_1"

# Change

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

pr <- rep(1/k,k)

y2_y1_1 <- c(.6,.2,.2)  #Given Y1 K=1
y2_y1_2 <- c(.2,.6,.2)  # Given Y1 K=2
y2_y1_3 <- c(.2,.2,.6)  # Given Y1 K=2


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
bmeanmat =matrix(4*c(0,.438,-0.035,
                     0,0,0,
                     0,-.438,0.035),nrow=3,ncol=3)

# class separation from GBTM2
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

### Outcome 1 - based on paper

vary1 <-3.125 #.01
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
#### CLASS SEPARATION CALCULATIONS 
################################################################
################################################################


### Conditional probabilities
### Ensure joint probabilities sum to 1
### Conditional probabilities
### Ensure joint probabilities sum to 1


####
id <- 1:N


### alternative

mult_mixture <- function(x){
  set.seed(x)
  samp_size_comp1 <- rmultinom(1,N,prob = pr)                   # Vector of sample sizes for classes in outcome Y1
  y2_y11 <-rmultinom(1,samp_size_comp1[1],prob = y2_y1_1)       # Vector of sample sizes for classes in outcome Y2 given class=1 in outcome Y1
  y2_y12 <-rmultinom(1,samp_size_comp1[2],prob = y2_y1_2)       # Vector of sample sizes for classes in outcome Y2 given class=2 in outcome Y1
  y2_y13 <-rmultinom(1,samp_size_comp1[3],prob = y2_y1_3)       # Vector of sample sizes for classes in outcome Y2 given class=2 in outcome Y1
  
  
  y1 <- rbind(mvrnorm(n = samp_size_comp1[1], mu = fe1, Sigma = sig2I1), 
              mvrnorm(n = samp_size_comp1[2], mu = fe2, Sigma = sig2I1),
              mvrnorm(n = samp_size_comp1[3], mu = fe3, Sigma = sig2I1))
  colnames(y1)<-c("Y1","Y2","Y3","Y4","Y5")
  class_y = c(c(rep(1,samp_size_comp1[1])),
              c(rep(2,samp_size_comp1[2])),
              c(rep(3,samp_size_comp1[3])))
  y2 <-rbind(mvrnorm(n=y2_y11[1],mu =fe1_2,Sigma= sig2I2),mvrnorm(n=y2_y11[2],mu = fe2_2, Sigma = sig2I2),mvrnorm(n=y2_y11[3],mu = fe3_2, Sigma = sig2I2),
             mvrnorm(n=y2_y12[1], mu = fe1_2, Sigma = sig2I2),mvrnorm(n=y2_y12[2], mu = fe2_2, Sigma = sig2I2),mvrnorm(n=y2_y12[3], mu = fe3_2, Sigma = sig2I2),
             mvrnorm(n=y2_y13[1], mu = fe1_2, Sigma = sig2I2),mvrnorm(n=y2_y13[2], mu = fe2_2, Sigma = sig2I2),mvrnorm(n=y2_y13[3], mu = fe3_2, Sigma = sig2I2))
  class_z<- c(c(rep(1,y2_y11[1])),c(rep(2,y2_y11[2])),c(rep(3,y2_y11[3])),
              c(rep(1,y2_y12[1])),c(rep(2,y2_y12[2])),c(rep(3,y2_y12[3])),
              c(rep(1,y2_y13[1])),c(rep(2,y2_y13[2])),c(rep(3,y2_y13[3])))
  colnames(y2)<-c("Z1","Z2","Z3","Z4","Z5") 
  mydata <-data.frame(repl=x,y1,class_y,y2,class_z)
  # Randomly sort the df
  rows <-sample(nrow(mydata))
  mydata<- mydata[rows,]
  mydata<- data.frame(mydata,id)
}


rep <-(1:n)

res<-lapply(rep,mult_mixture)

#### To Mplus
### Generate n different data sets for Mplus
dir.create(paste0(wd,"/Data"))
library(MplusAutomation)
for (i in 1:n){
  prepareMplusData(as.data.frame(res[[i]]), paste0(wd,"/Data/Mplus","_",i,".dat"),
                   overwrite = TRUE)}


### Create models according to Mplus .txt script
### Create models according to Mplus .txt script
### First set working directory again


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
### If input is newer then the model is run, otherwise it is skipped: replaceOutfile = "modifiedDate"
runModels(wd,recursive=TRUE,showOutput = FALSE,replaceOutfile = "never",
          logFile = paste0(wd,"/ComparisonLog.txt"))


### Extract summaries - this is where fit statistics reside
library(plyr)
allModelParams <- readModels(wd, recursive = TRUE, what="summaries")
justSummaries <- do.call("rbind.fill",sapply(allModelParams,"[","summaries"))

## Remove incomplete runs - relic of Mplus Automation
#nrow(justSummaries)
justSummaries <- justSummaries[complete.cases(justSummaries$NGroups),]
#nrow(justSummaries)


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
#table(justSummaries$rep)
justSummaries$rep <-stringr::str_remove_all(justSummaries$rep," ")
#table(justSummaries$rep)
#View(justSummaries)


#### Keep specific columns

#### Keep specific columns

df <- justSummaries[c(2,6,7,8,10,11,12,13,15,16,17,18,19,20,21)]

#View(df)

### Extract model type
df$model <- mid(df$Title,9,10)
### Modify
#df$model <-recode(df$model,"Y1 UNI"="Y1_UNI")
#df$model <-recode(df$model,"Y2 UNI"="Y2_UNI")
#table(df$model)
df$model <- str_remove_all(df$model,"rep")
df$model <- str_remove_all(df$model,"re")
df$model <- str_remove_all(df$model," ")
#df$model <- str_remove_all(df$model,"r")
#table(df$model)

### Set model and replication as factors
df$model <-as.factor(df$model)
df$rep <- as.factor(df$rep)

### Want to create estimated class number (K) variable
df$class <-left(df$Title,2)
df$class <- str_remove_all(df$class,"-")

df$class <-as.numeric(df$class)
#table(df$class)

### Save Data frame

save(df,file=paste0(wd,"/df.Rda"))



###########################################
###########################################
###########################################
