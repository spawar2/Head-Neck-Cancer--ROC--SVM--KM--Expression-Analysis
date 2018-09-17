true_Y = c(1,1,1,1,0,1,0,1,0,0)
probs = c(2,3.999,6.999,10.973,78.568,10.421,1.382,0.377,99.146,0.11)

getROC_AUC = function(probs, true_Y){
    probsSort = sort(probs, decreasing = TRUE, index.return = TRUE)
    val = unlist(probsSort$x)
    idx = unlist(probsSort$ix)  

    roc_y = true_Y[idx];
    stack_x = cumsum(roc_y == 0)/sum(roc_y == 0)
    stack_y = cumsum(roc_y == 1)/sum(roc_y == 1)    

    auc = sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
    return(list(stack_x=stack_x, stack_y=stack_y, auc=auc))
}

aList = getROC_AUC(probs, true_Y) 

stack_x = unlist(aList$stack_x)
stack_y = unlist(aList$stack_y)
auc = unlist(aList$auc)

plot(stack_x, stack_y, type = "l", col = "blue", xlab = "False Positive Rate", ylab = "True Positive Rate", main = "ROC")
axis(1, seq(0.0,1.0,0.1))
axis(2, seq(0.0,1.0,0.1))
abline(h=seq(0.0,1.0,0.1), v=seq(0.0,1.0,0.1), col="gray", lty=3)
legend(0.7, 0.3, sprintf("%3.3f",auc), lty=c(1,1), lwd=c(2.5,2.5), col="blue", title = "AUC")

-------------------Not used below------------------------------------------

url = "https://github.com/zhandong/TCGA2STAT/blob/master/TCGA2STAT_1.2.tar.gz"

destfile = "C:/Users/Bio-user/Desktop/NEWTCGA2STAT_1.2.tar.gz"
download.file(url, destfile, quiet = FALSE, mode = "w")

install.packages("NEWTCGA2STAT_1.2.tar.gz", repos = NULL, type = "source")

----------------------------MAIN-------------------------------------------

# A KM COX analysis was performed on survial and days, with significant P value
1477 genes were choosen. Further, patients were classifed as with tumor and
without tumor/control, finally ROC (AUC >0.65) was then calcualted to check
how effectievly can they prdict death or survival for final marker genes, there  same was validated with
SVM algorithm with 7, 5, 3 fold training and testing set. Finally plots were
made for grade, stages, etc to show high or low expression of these genes and
their functional annotation.

setwd("C:/Users/Bio-user/Desktop/TCGA Paper")
All_data <- read.csv("Final.csv", header=TRUE, sep=",")
library("survival")
library("Biobase")
library(pbapply)

samp2 <- All_data[,-1]
rownames(samp2) <- All_data[,1]

m <- as.numeric(as.matrix(samp2[8,]))
n <- as.numeric(as.matrix(samp2[5,]))
surv_object <- with(samp2, Surv(m, n))

j <- as.numeric(as.matrix(samp2[15,]))
mod <- coxph(surv_object~j)

sink('analysis-output.csv')

for (i in 15:20545){
  j <- as.numeric(as.matrix(samp2[i,]))
  print(coxph(surv_object~j))
  pbsapply(i, sqrt)
}

sink()

nums <- 15:20545
sqrt <- sapply(nums, print)

Final <- read.csv("Final_Final.csv", header=TRUE, sep=",")

FinalT <- t(Final)

with_tumor <- FinalT[which(FinalT[,1485]=='with tumor'),]

tumor_free <- FinalT[which(FinalT[,1485]=='tumor free'),]


for (i in 1:1476){
true_Y = as.numeric(as.matrix(with_tumor[,1481]))
probs = as.numeric(as.matrix(with_tumor[,i]))

getROC_AUC = function(probs, true_Y){
    probsSort = sort(probs, decreasing = TRUE, index.return = TRUE)
    val = unlist(probsSort$x)
    idx = unlist(probsSort$ix)  

    roc_y = true_Y[idx];
    stack_x = cumsum(roc_y == 0)/sum(roc_y == 0)
    stack_y = cumsum(roc_y == 1)/sum(roc_y == 1)    

    auc = sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
    return(list(stack_x=stack_x, stack_y=stack_y, auc=auc))
}

aList = getROC_AUC(probs, true_Y) 

stack_x = unlist(aList$stack_x)
stack_y = unlist(aList$stack_y)
auc = unlist(aList$auc)
print(c(i, auc))
}

#Plot for ROC
marker1 <- with_tumor[,c(1298, 826, 425, 1082, 846, 177)]

true_Y = as.numeric(as.matrix(with_tumor[,1481]))
probs = as.numeric(as.matrix(marker1[,1]))

getROC_AUC = function(probs, true_Y){
    probsSort = sort(probs, decreasing = TRUE, index.return = TRUE)
    val = unlist(probsSort$x)
    idx = unlist(probsSort$ix)  

    roc_y = true_Y[idx];
    stack_x = cumsum(roc_y == 0)/sum(roc_y == 0)
    stack_y = cumsum(roc_y == 1)/sum(roc_y == 1)    

    auc = sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
    return(list(stack_x=stack_x, stack_y=stack_y, auc=auc))
}

aList = getROC_AUC(probs, true_Y) 

stack_x = unlist(aList$stack_x)
stack_y = unlist(aList$stack_y)
auc = unlist(aList$auc)

par(mfrow=c(3,2))
plot(stack_x, stack_y, type = "l", col = "blue", xlab = "False Positive Rate", ylab = "True Positive Rate", main = "ROC")
axis(1, seq(0.0,1.0,0.1))
axis(2, seq(0.0,1.0,0.1))
abline(h=seq(0.0,1.0,0.1), v=seq(0.0,1.0,0.1), col="gray", lty=3)
legend(0.7, 0.3, sprintf("%3.3f",auc), lty=c(1,1), lwd=c(2.5,2.5), col="blue", title = "AUC")

Final gene signature:
RGMA|56963   col 1298
KLHL14|57565  col 826
DLG2|1740  col 425
NOVA1|4857  col 1082
KRTAP5-8|57830      col 846
C1orf190|541468   col 177

SVM for ROC validation
# Summation of all genes for marker
marker <- tumor_free[,c(1298, 826, 425, 1082, 846, 177)]
with_tumor[,1481]

# Making a numeric matrix for row sums
dim(marker)
[1] 110   6

N <- mapply(marker, FUN=as.numeric)
N <- matrix(data=N, ncol=6, nrow=110)

Index <- rowSums(N)
Status <- as.matrix(tumor_free[,1480])
DataSVM <- as.data.frame(cbind(Index, Status))
library(caret)
library(e1071)
set.seed(3033)

intrain <- createDataPartition(y = Index, p= 0.3, list = FALSE)
training <- DataSVM[intrain,]
testing <- DataSVM[-intrain,]

svm_model <- svm(V2~ ., data=DataSVM , method="C-classification", kernel="linear")

pred_train <-predict(svm_model,training)
pred_test <-predict(svm_model,testing)

#7 Fold
table(training$V2,pred_train)
  pred_train
        alive dead
  alive    68    0
  dead      5    5
((68+5)/(nrow(training)))*100
[1] 93.58974

table(testing$V2,pred_test)
pred_test
        alive dead
  alive   180    5
  dead     27    7
((180+7)/(nrow(testing)))*100
[1] 85.38813

#5 Fold
pred_train
        alive dead
  alive    51    1
  dead      2    2
((51+2)/(nrow(training)))*100
[1] 94.64286

pred_test
        alive dead
  alive   197    4
  dead     30   10
((197+10)/(nrow(testing)))*100
[1] 85.89212

#10 Fold
pred_train
        alive dead
  alive    34    0
  dead      1    1
((34+1)/(nrow(training)))*100
[1] 97.22222

table(testing$V2,pred_test)
pred_test
        alive dead
  alive   214    5
  dead     31   11
((214+11)/(nrow(testing)))*100
[1] 86.2069


marker <- tumor_free[,c(1298, 826, 425, 1082, 846, 177)]
N <- mapply(marker, FUN=as.numeric)
N <- matrix(data=N, ncol=6, nrow=297)
Index <- rowSums(N)

stage <- as.matrix(tumor_free[,1478])
gender <- as.matrix(tumor_free[,1477])
grade <- as.matrix(tumor_free[,1479])
vital <- as.matrix(tumor_free[,1480])
status <- as.matrix(tumor_free[,1481])
race <- as.matrix(tumor_free[,1482])
hpv <- as.matrix(tumor_free[,1483])
days_to_death <- as.matrix(tumor_free[,1484])
neoplasm <- as.matrix(tumor_free[,1485])
margin <- as.matrix(tumor_free[,1486])
lymphnodes <- as.matrix(tumor_free[,1487])
perineural <- as.matrix(tumor_free[,1488])
radiation.course <- as.matrix(tumor_free[,1489])
radiation.dosage <- as.matrix(tumor_free[,1490])

marker1 <- with_tumor[,c(1298, 826, 425, 1082, 846, 177)]
N <- mapply(marker1, FUN=as.numeric)
N <- matrix(data=N, ncol=6, nrow=110)
Index1 <- rowSums(N)

stage1 <- as.matrix(with_tumor[,1478])
gender1 <- as.matrix(with_tumor[,1477])
grade1 <- as.matrix(with_tumor[,1479])
vital1 <- as.matrix(with_tumor[,1480])   /////NO USE
status1 <- as.matrix(with_tumor[,1481])   //////Make KM for individual gene
race1 <- as.matrix(with_tumor[,1482])
hpv1 <- as.matrix(with_tumor[,1483])
days_to_death1 <- as.matrix(with_tumor[,1484]) //////Make KM for individual gene
neoplasm1 <- as.matrix(with_tumor[,1485])  /////NO USE
margin1 <- as.matrix(with_tumor[,1486])
lymphnodes1 <- as.matrix(with_tumor[,1487])
perineural1 <- as.matrix(with_tumor[,1488])
radiation.course1 <- as.matrix(with_tumor[,1489])
radiation.dosage1 <- as.matrix(with_tumor[,1490])

#Stage
counts <- table(Index, stage)
counts1 <- table(Index1, stage1)

par(mfrow=c(5,2))
barplot(counts, xlab="Number of Patients (N=297)", ylab="Marker Index", col=c("Red"))
barplot(counts1, xlab="Number of Patients (N=110)", ylab="Marker Index", col=c("Green"))

#Grade
counts <- table(Index, grade)
counts1 <- table(Index1, grade1)

#par(mfrow=c(2,1))
barplot(counts, xlab="Number of Patients (N=297)", ylab="Marker Index", col=c("Red"))
barplot(counts1, xlab="Number of Patients (N=110)", ylab="Marker Index", col=c("Green"))

#Gender
counts <- table(Index, gender)
counts1 <- table(Index1, gender1)

#par(mfrow=c(2,1))
barplot(counts, xlab="Number of Patients (N=297)", ylab="Marker Index", col=c("Red"))
barplot(counts1, xlab="Number of Patients (N=110)", ylab="Marker Index", col=c("Green"))

#Gender
counts <- table(Index, race)
counts1 <- table(Index1, race1)

#par(mfrow=c(2,1))
barplot(counts, xlab="Number of Patients (N=297)", ylab="Marker Index", col=c("Red"))
barplot(counts1, xlab="Number of Patients (N=110)", ylab="Marker Index", col=c("Green"))

#HPV Status
counts <- table(Index, hpv)
counts1 <- table(Index1, hpv1)

#par(mfrow=c(2,1))
barplot(counts, xlab="Number of Patients (N=297)", ylab="Marker Index", col=c("Red"))
barplot(counts1, xlab="Number of Patients (N=110)", ylab="Marker Index", col=c("Green"))

#Marginal Status
counts <- table(Index, margin)
counts1 <- table(Index1, margin1)

par(mfrow=c(5,2))
barplot(counts, xlab="Number of Patients (N=297)", ylab="Marker Index", col=c("Red"))
barplot(counts1, xlab="Number of Patients (N=110)", ylab="Marker Index", col=c("Green"))

#Lymphnode Status
counts <- table(Index, lymphnodes)
counts1 <- table(Index1, lymphnodes1)

#par(mfrow=c(2,1))
barplot(counts, xlab="Number of Patients (N=297)", ylab="Marker Index", col=c("Red"))
barplot(counts1, xlab="Number of Patients (N=110)", ylab="Marker Index", col=c("Green"))

#Perineural
counts <- table(Index, perineural)
counts1 <- table(Index1, perineural1)

#par(mfrow=c(2,1))
barplot(counts, xlab="Number of Patients (N=297)", ylab="Marker Index", col=c("Red"))
barplot(counts1, xlab="Number of Patients (N=110)", ylab="Marker Index", col=c("Green"))

#radiation.course
counts <- table(Index, radiation.course)
counts1 <- table(Index1, radiation.course1)

#par(mfrow=c(2,1))
barplot(counts, xlab="Number of Patients (N=297)", ylab="Marker Index", col=c("Red"))
barplot(counts1, xlab="Number of Patients (N=110)", ylab="Marker Index", col=c("Green"))


#radiation.dosage
counts <- table(Index, radiation.dosage)
counts1 <- table(Index1, radiation.dosage1)

#par(mfrow=c(2,1))
barplot(counts, xlab="Number of Patients (N=297)", ylab="Marker Index", col=c("Red"))
barplot(counts1, xlab="Number of Patients (N=110)", ylab="Marker Index", col=c("Green"))

#KM for Index genes 

a <- as.data.frame(with_tumor[,1:1476])
b <- as.data.frame(tumor_free[,1:1476])

m <- as.numeric(as.matrix(days_to_death))
n <- as.numeric(as.matrix(status))
surv_object <- Surv(m, n)
fit <- survfit(surv_object ~ n, data=b)
plot(fit, col = "Green", xlab="Survival time (days)", ylab="Survival Probability")
m <- as.numeric(as.matrix(days_to_death1))
n <- as.numeric(as.matrix(status1))
surv_object <- Surv(m, n)
fit2 <- survfit(surv_object ~ n, data=a)
lines(fit2, col="Red")
legend("topright", legend=c("Tumor", "Normal"),
       col=c("red", "Green"), lty=1:2, cex=0.8)

library(survminer)
library(survival)
fit.list <- list(Normal = fit, Tumor = fit2)
surv_pvalue(fit.list, combine = TRUE)

      id variable        pval   method   pval.txt
1 Normal        n 0.001752763 Log-rank p = 0.0018
2  Tumor        n 0.001752763 Log-rank p = 0.0018






































