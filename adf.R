######################
ADF <- function(y_t,lags=floor(12*((length(y_t)/100)^(0.25))),trend=T){
  if(lags>=length(y_t)){stop("Number of lags cannot exceed N!")} #need sufficient time span
  #https://papyrus.bib.umontreal.ca/xmlui/bitstream/handle/1866/2114/9423.pdf
  y <- timeSeries::as.timeSeries(as.vector(y_t)) #lags and such only work on this class
  A <- list() -> R #placeholders needed for later
  DF <- do.call(what = "rbind.data.frame", #these are the Dickey-Fuller bootstrapped t-statistics 
                list(c(25,-3.00,-3.60),    #for testing unit roots at the 5% significance level
                     c(50,-2.93,-3.50),
                     c(100,-2.89,-3.45),
                     c(250,-2.88,-3.43),
                     c(500,-2.87,-3.42),
                     c(1000,-2.86,-3.41)))
  colnames(DF) <- c("n","trend.5","notrend.5") #they differ on whether or not trend is assumed
  #################
  for(j in 1:lags){ #default number of lags from rule of thumb in Schwert (1989)
    tam <- list() #another placeholder for immediately below
    ##############
    for(i in 1:j){
      tam[[i]] <- paste0("diff(lag(y,",i,"))") #setting up the ADF functional form
    }#
    R[[j]] <- as.formula(ifelse(trend==T, #take the user's inputs for trend inclusion and create functional forms
                                paste0("diff(y)~c(1:length(y_t))+lag(y,1)+",
                                       paste(unlist(tam),collapse = "+")),
                                paste0("diff(y)~lag(y,1)+",paste(unlist(tam),collapse = "+"))))
    A[[j]] <- AIC((lm(R[[j]]))) #store the AIC values in a list
  }
  temp <- R[[which(unlist(A)==min(unlist(A)))]] #select the number of lags with minimum AIC
  booboo <- ifelse(trend==F, #take the "best" number of lags model and test the null of non-stationarity
                   DF[(as.numeric(length(y_t)>DF$n)[1]),"notrend.5"],
                   DF[(as.numeric(length(y_t)>DF$n)[1]),"trend.5"])
  roo <- ifelse(coef(summary(lm(temp)))["lag(y, 1)",3]<booboo, #print output
                (list(paste(c("Reject Presence of Unit Root: Time Series is Stationary",
                              paste("rho=",round(coef(summary(lm(temp)))["lag(y, 1)",1],4)+1),
                              paste("t-value=",round(coef(summary(lm(temp)))["lag(y, 1)",3],4)))))),
                (list(c("Fail to Reject Unit Root: Time Series is Non-Stationary",
                        paste("rho=",round(coef(summary(lm(temp)))["lag(y, 1)",1],4)+1),
                        paste("t-value=",round(coef(summary(lm(temp)))["lag(y, 1)",3],4))))))
  return(roo)
}

############################Simple example: two time series
set.seed(1234567)
eps  <- rnorm(1000)
eps2 <- rnorm(1000)
y <- eps #initialize the stationary time series
x <- eps2 #initialize the non-stationary time series
for (t in 2:1000) {
  y[t] <- y[t-1] + eps[t] #coefficient is exactly one
  x[t] <- 1.01*x[t-1] + eps2[t] #coefficient is slightly larger than 1 (explosive)
}
plot(y) #plot the stationary series
ADF(y) #test for nonstationarity (should reject)
plot(x) #plot the non-stationary series
ADF(x) #test for nonstationarity (should fail to reject)
