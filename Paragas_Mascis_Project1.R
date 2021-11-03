# Names: Spencer Paragas and Christopher Mascis




##################################### LOESS Regression Function #####################################

  


library(dplyr)
library(ggplot2)

# This function has the following inputs.
# 
# * x - a numeric input vector
# * y - a numeric response
#
# Note span and degree are shown with their default values.
# * degree should be 1 or 2 only
# * span can be any value in interval (0, 1) non-inclusive.
#
# If show.plot = TRUE then a plot of the final fit is shown

myloess <- function(x, y, span = 0.5, degree = 1, show.plot = TRUE) {
  
  # Getting range of values
  xrange <- diff(range(x))
  
  # Checking span meets requirements and setting width
  if (between(span, 0, 1)) {
    width <- span*xrange
  }
  else {
    stop("Span must be between 0 and 1 non-inclusive")
  }
  
  # Getting total number of points and windows (they'll be the same)
  N_total <- length(x)
  Win_total <- length(x)
  
  # Allocating space for vector of each window's population
  n_points <- vector(mode = "integer", length = length(x))
  
  # Allocating space for vector of fitted values
  yhat <- vector(mode = "numeric", length = length(x))
  
  # Combining our variables in data frame
  mydata <- cbind.data.frame(x, y)
  
  # Fitting each point
  for(x0 in x) {
    
    # Setting population of window
    sample <- subset(mydata, between(x, x0-width/2, x0+width/2))
    n <- length(sample[,1])
    
    # Getting weights into diagonal matrix
    weights <- (1-(abs(sample[,1] - x0)/width*2)^3)^3
    W_mat <- diag(weights, n, n)
    
    # Checking degree and completing our regression accordingly
    if (degree == 1) {
      
      X_mat <- cbind(rep(1, n), sample[,1])
      betahat <- solve( t(X_mat)%*%W_mat%*%X_mat ) %*% t(X_mat) %*% W_mat %*% sample[,2]
      
      # Getting fitted value
      yhat[which(x0 == x)] <-  betahat[1] + betahat[2]*x0
    }
    else if (degree == 2) {
      
      X_mat <- cbind(rep(1, n), sample[,1], sample[,1]^2)
      betahat <- solve( t(X_mat)%*%W_mat%*%X_mat ) %*% t(X_mat) %*% W_mat %*% sample[,2]
      
      # Getting fitted value
      yhat[which(x0 == x)] <-  betahat[1] + betahat[2]*x0 + betahat[3]*x0^2
    }
    else {
      stop("Degree must be set to either 1 or 2")
    }
    
    # Getting population of window
    n_points[which(x0 == x)] <- n
  }
  
  # Calculating SSE, MSE, and residual standard error
  SSE <- sum((y-yhat)^{2})
  MSE <- SSE/(N_total-2)
  RSE <- sqrt(MSE)
  
  # Creating plot of final fit
  loessplot <- ggplot(mydata, aes(x, y)) +
    geom_point(size = 3, alpha = 0.5, color = "grey") +
    geom_line(aes(x, yhat), color = "red", lty = 1) +
    xlab(deparse(substitute(x))) + ylab(deparse(substitute(y))) +
    ggtitle(paste("LOESS with degree =", degree,"and span =", span, sep = " "))
  
  # Checking whether to show plot or not
  if (show.plot == T) {
    print(loessplot)
  }
  
  # Returning named list
  return(invisible(list(span = span, degree = degree, N_total = N_total, Win_total = Win_total,
                        n_points = n_points, SSE = SSE, RSE = RSE, loessplot = loessplot)))
}




############################ kNN Classification and Regression Function ############################




library(caret)

# This function has the following inputs
#
# * train - matrix or data frame of training set cases
# * test - matrix or data frame of test set cases.  
#     (A vector will be interpreted as a row vector for a single case.)
# * y_train - Either a numeric vector, or factor vector for the responses in the training set
# * y_test - Either a numeric vector, or factor vector for the responses in the testing set
# * k - number of neighbors considered, the default value is 3
#
# If weighted = TRUE, then function uses the distance weighted kNN,
#  otherwise it does the default kNN method.

mykNN <- function(train, test, y_train, y_test, k = 3, weighted = TRUE) {
  
  # Allocating space for vector of fitted values
  yhat <- vector(mode = "numeric", length = length(y_test))
  
  # Fitting each point
  for (row in 1:nrow(test)) {
    
    # Getting features
    x0 <- test[row,]
    
    # Getting Euclidean distance and weights
    distmat <- as.matrix(dist(rbind(x0, train), method = "euclidean", diag = T, upper = T))
    dist <- distmat[1,]
    weights <- 1/dist
    
    # Getting kth nearest neighbor
    a.order <- order(weights, decreasing = T)[-1]
    last <- weights[a.order[k]]
    
    # Setting weights to 0 except for kNN
    y_train_sub <- y_train[-which(!is.finite(weights))]
    weights <- weights[-which(!is.finite(weights))]
    weights[weights < last] <- 0
    
    # Handling unweighted case
    if(weighted == F) {
      weights[weights != 0] <- 1
    }
    
    tot_weight <- sum(weights)
    
    # Handling classification
    if (is.factor(y_train)) {
      
      # Allocating space for vector of class weights
      class <- vector(mode = "numeric", length = nlevels(y_test))
      
      # Summing weights for each class
      for (i in 1:nlevels(y_test)) {
        class[i] <- sum(weights[which(y_train == levels(y_test)[i])])
      }
      
      # Assigning factor with highest class weight to yhat
      yhat[row] <- levels(y_test)[which(class == max(class))]
    }
    
    # Handling regression
    else {
      
      # Calculating yhat
      yhat[row] <- sum(weights*y_train/tot_weight)
    }
  }
  
  # Handling classification
  if (is.factor(y_train)) {
    
    yhat <- factor(yhat)
    
    # Setting number of correct predictions to 0
    num_correct <- 0
    
    # Getting number of correct predictions
    for (j in 1:length(yhat)) {
      
      # Checking if classes match
      if (y_test[j] == yhat[j]) {
        
        num_correct <- num_correct + 1
      }
    }
    
    # Calculating accuracy and error rate
    accuracy <- num_correct / length(yhat)
    error_rate <- 1 - accuracy
    
    # Getting confusion matrix
    confmat <- confusionMatrix(yhat, y_test)
    
    # Returning named list
    return(invisible(list(yhat = yhat, accuracy = accuracy,
                          error_rate = error_rate, confmat = confmat, k = k)))
  }
  
  # Handling regression
  else {
    
    # Calculating residuals, SSE, MSE, and residual standard error
    residuals <- y_test - yhat
    SSE <- sum(residuals^{2})
    MSE <- SSE/(nrow(test)-2)
    RSE <- sqrt(MSE)
    
    # Returning named list
    return(invisible(list(yhat = yhat, residuals = residuals, SSE = SSE, k = k)))
  }
}