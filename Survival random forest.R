library(readxl)
library(survival)
library(SurvMetrics)
options(warn = -1)
#seed=23
sim.data<- read_excel("~/Documents/Research/Factor PLANN new1/Cirrihosis/cirrhosis.xlsx")
sim.data<-sim.data[1:312,]
sim.data$Cholesterol<-as.numeric(sim.data$Cholesterol)
sim.data$Copper<-as.numeric(sim.data$Copper)
sim.data$Alk_Phos<-as.numeric(sim.data$Alk_Phos)
sim.data$SGOT<-as.numeric(sim.data$SGOT)
sim.data$Tryglicerides<-as.numeric(sim.data$Tryglicerides)
sim.data$Platelets<-as.numeric(sim.data$Platelets)
sim.data$Prothrombin<-as.numeric(sim.data$Prothrombin)
sim.data$Stage<-as.numeric(sim.data$Stage)
sim.data2<-sim.data[,-1]
sim.data1<-na.omit(sim.data2)
sim.data5<-sim.data1
colnames(sim.data5)[1]<-"time"
colnames(sim.data5)[2]<-"status"

edema_new=vector()
for (i in 1:nrow(sim.data1)) {
  if(sim.data1$Edema[i]=="S"){edema_new[i]=1}
  else if(sim.data1$Edema[i]=="Y"){edema_new[i]=2}
  else{edema_new[i]=0}
}
sim.data5$Edema<-edema_new
sim.data5$status<-ifelse(sim.data1$Status=="D",1,0)
sim.data5$Sex<-ifelse(sim.data1$Sex=="M",1,0)
sim.data5$Drug<-ifelse(sim.data1$Drug=="Placebo",1,0)
sim.data5$Ascites<-ifelse(sim.data1$Ascites=="Y",1,0)
sim.data5$Hepatomegaly<-ifelse(sim.data1$Hepatomegaly=="Y",1,0)
sim.data5$Spiders<-ifelse(sim.data1$Spiders=="Y",1,0)

data_scaled<-sim.data5


data_scaled<-as.data.frame(data_scaled)
id<-seq(1, nrow(data_scaled),1)
data_scaled<-cbind(data_scaled,id)
data_scaled<-data_scaled[order(data_scaled$id),]

# defining test set and training set
set.seed(seed)
N <- nrow (data_scaled)
index <- sample (1:N, round (N/3), replace = FALSE )
training_data<- data_scaled[-index,]
test_data<- data_scaled[index,]
cox.train<-training_data
cox.test<-test_data

cox.train.s<-cox.train
for (i in 9:18) {
  cox.train.s[,i]<-scale(cox.train[,i])
}

cox.test.s<-cox.test
for (i in 9:18) {
  cox.test.s[,i]<-scale(cox.test[,i],center = mean(cox.test[,i]),
                        scale = sd(cox.test[,i]))
}


srf=matrix(NA,10,3)
tree=c(50,100,150,200,250,300,400,500,700,1000)

for (k in 1:10) {
  set.seed(seed)
  library(randomForestSRC)
  v.obj <- rfsrc(Surv(time, status) ~ ., data = cox.train.s[,-20], 
                 ntree = tree[k], block.size = 1, ntime=10)
  
  # Predicting survival probabilities for new data
  predictions <- predict(v.obj, newdata = cox.test.s, estimate = "survival")
  
  ff1<-matrix(NA,10,2)
  for (i in 1:10) {
    data <- data.frame(time = cox.test$time, status = cox.test$status, 
                       predicted_scores = predictions$survival[,i])
    
    # Calculate C-index
    c_index <- concordance(Surv(time, status) ~ predicted_scores, data)
    c_index<-c_index$concordance
    brier<-Brier(Surv(data$time, data$status), data$predicted_scores)
    brier<-as.numeric(brier)
    ff1[i,]=cbind(c_index, brier)
  }
  object=Surv(data$time, data$status)
  time_interest=v.obj$time.interest
  ibs.sc=IBS(object, predictions$survival, time_interest)
  srf[k,]=c(apply(ff1, 2, mean), ibs.sc)
  colnames(srf)<- c("c-index","brier","IBS")
  rownames(srf)<-c(50,100,150,200,250,300,400,500,700,1000)
}
srf


