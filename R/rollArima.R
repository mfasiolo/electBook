rollArima<-function (arima.model, ynew, horizon = 1) 
{
  prevARIMA <- array(0, dim = length(ynew))
  prevARIMA[1] <- forecast(arima.model, h = horizon)$mean[horizon]
  for (i in 1:(length(ynew) - 1)) {
    ts2 <- c(arima.model$x, ynew[1:i])
    refit <- Arima(ts2, model = arima.model)
    prevARIMA[i + 1] <- forecast(refit, h = horizon)$mean[horizon]
  }
  return(prevARIMA)
}
