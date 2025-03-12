library(keras)
seed=23
source("~/Documents/Research/Factor PLANN new1 Revised 1 after sub 1/GitHub codes/Survival random forest.R")

# defining test set and training set

set.seed(seed)
N <- nrow (data_scaled)
index <- sample (1:N, round (N/3), replace = FALSE )
training_data1<- data_scaled[-index,]
test_data1<- data_scaled[index,]
cox.train1<-training_data1
cox.test1<-test_data1

#factor analysis
fa_result<-factanal(training_data1[,3:19], factors = 11, scores = "regression")

#Weight
W=solve(fa_result$correlation) %*% fa_result$loadings

#Standardize
mus<-colMeans(cox.train1[,3:19])
sigmas<-apply(cox.train1[,3:19],2,sd)

#Scores
cal_scores=scale(cox.test1[,3:19],mus,sigmas) %*% W

#New train and test set
training_data<-cox.train<-cbind(fa_result$scores,cox.train1[,1:2])
test_data<-cox.test<-cbind(cal_scores,cox.test1[,1:2])


#transform the training data into long format
data_train<-function(data){
  set.seed (seed)
  N<-nrow(data)
  data$interval<-as.numeric(cut(data$time, 
          breaks = c(40,v.obj$time.interest[-1],4556)))
  data$survival<-as.numeric(cut(data$time,breaks = c(40,v.obj$time.interest[-1],4556)))
  data$id<-1:N
  n.times<-data$interval
  data_long<-data[rep(seq_len(N), times = n.times),]
  
  #create the correct intervals
  for (i in unique(data_long$id)) {
    n_length<-length(data_long$interval[data_long$id == i])
    data_long$interval[data_long$id == i] <- 1:n_length
  }
  data_long$stat<-vector(mode = "numeric", length = nrow(data_long))
  #put indication 1 on status at the interval that patient dies
  for (i in 1:nrow(data_long)) {
    if (data_long$status[i] == 1 &&
        data_long$survival[i] <= data_long$interval[i])
      data_long$stat[i]<-1
  }
  return(data_long)
}

# function that creates data in the right long format for
# test set for each patient the interval 


data_test<-function(data){
  set.seed (seed)
  N<-nrow(data)
  #assign survival times to 10 intervals
  data$interval<-max(as.numeric(cut(data$time, breaks = c(40,v.obj$time.interest[-1],4556))))
  # the true interval survival
  data$survival<-as.numeric(cut(data$time,breaks = c(40,v.obj$time.interest[-1],4556)))
  data$id<-70001:(70000+N) # define the patient ids abstractly
  n.times<-data$interval
  data_long<-data[rep(seq_len(N), times = n.times),]
  # create the correct intervals
  for (i in unique(data_long$id)) {
    n_length<-length(data_long$interval[data_long$id == i])  
    data_long$interval[data_long$id == i] <- 1: n_length
  }
  data_long$stat<-vector(mode = "numeric", length = nrow(data_long))
  #put indication 1 on status at the intervals on which a patient has died
  for (i in 1:nrow(data_long)) {
    if (data_long$status[i] == 1 &&
        data_long$survival[i] <= data_long$interval[i])
      data_long$stat[i]<-1
  }
  return(data_long)
}


#Kaplian Meier plot and event time distribution plot

par(mfrow = (c(1,2)))
plot(survfit(Surv(data_scaled$time, data_scaled$status)~1), xlab = "Time", ylab = "Survival probability")
hist(data_scaled$time[which(data_scaled$status == 1)], xlab = "Time", main = "", col = "yellow",xlim = c(0,5000))


FEDeepSurv<- function(training_data, test_data, node.size, node.size1, node.size2){
  set.seed (seed)
  tensorflow::set_random_seed(seed)
  training<-data_train(training_data)
  test<-data_test(test_data)
  interval.test<- test$interval
  surv.test<-test$survival
  
  X.test<-test[,c(1:11,14)]
  
  status<- training$status
  time <- training$time
  interval<- scale(training$interval)
  stat<- training$stat
  
  X<-as.data.frame(training[,c(1:11,14)])
  
  # Create a sequential model
  model <- keras_model_sequential()
  
  # Add the first hidden layer with sigmoid activation
  model %>%
    layer_dense(units = node.size, activation = 'relu', input_shape = c(dim(X)[2])) %>%
    
    # Add the second hidden layer with ReLU activation
    layer_dense(units = node.size1, activation = 'relu') %>%
    
    # Add the third hidden layer with ReLU activation
    layer_dense(units = node.size2, activation = 'relu') %>%
    
    # Add the output layer (adjust units according to your problem)
    layer_dense(units = 1, activation = 'sigmoid')
  
  # Compile the model
  model %>% compile(
    #loss = 'mean_squared_error',
    loss = 'binary_crossentropy',
    optimizer = 'adam'
    #optimizer = optimizer_sgd(learning_rate = 0.01)
  )
  
  # Train the model with your data
  model %>% fit(
    x = as.matrix(X),
    y = as.vector(stat),
    epochs = 50,
    batch_size = 32
  )
  
  
  #data frame with prediction, prediction times, id, status
  # and original survival time in long format
  # Make predictions
  predictions <- predict(model, as.matrix(X.test))
  
  
  
  coll <- as.data.frame(cbind(predictions, test$survival, test$id, test$stat))
  
  ITs<- sort(unique(training_data$time), decreasing = FALSE)
  
  colnames(coll) <- c("prob", "surv", "id", "stat")
  
  stat <- test_data$status
  idd <- unique(coll$id)
  
  
  # obtaining survival estimates
  coll.s<- coll[, c(1,3)]
  coll.s<- split(coll.s, coll$id)
  coll.t<- lapply(coll.s, function(x) {
    x <-  cumprod(1-x$prob)
  })
  
  
  #obtaining prediction matrix
  pred.matrix <- do.call("rbind", coll.t)
  
  ff<-matrix(NA,10,2)
  for (i in 1:10) {
    data <- data.frame(time = test_data$time, status = test_data$status, 
                       predicted_scores = pred.matrix[,i])
    
    # Calculate C-index
    c_index <- concordance(Surv(time, status) ~ predicted_scores, data)
    c_index<-c_index$concordance
    brier<-Brier(Surv(data$time, data$status), data$predicted_scores)
    brier<-as.numeric(brier)
    
    ff[i,]=cbind(c_index, brier)
    colnames(ff)<-c("c-index","brier")
  }
  object=Surv(data$time, data$status)
  time_interest=v.obj$time.interest
  ibs.sc=IBS(object, pred.matrix, time_interest)
  lst <- list()
  lst[[1]]<- c(apply(ff, 2, mean), ibs.sc)
  lst[[2]]<-ff
  obj<- lst
  return(obj)
}

FEDeepSurv(training_data, test_data, 4, 17, 5)
