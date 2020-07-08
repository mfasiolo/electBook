#'
#' Rolling Arima model fitting
#'
#' @param arima.model a time series or time series model for which forecasts are required.
#' @param ynew vector of response variables.
#' @param horizon number of periods for forecasting.
#' @return A vector of forecasts from a rolling Arima model.
#' @importFrom forecast forecast Arima
#' @name rollArima
#' @rdname rollArima
#' @export rollArima
#' @export
#'
rollArima <- function(arima.model, ynew, horizon = 1) 
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
