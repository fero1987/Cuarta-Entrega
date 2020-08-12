############# CUARTA ENTREGA - ANÁLISIS DE SERIES TEMPORALES #############
# Alumnos: Martínez-De Elia-Díaz Ugalde

suppressPackageStartupMessages({
  library(tseries)
  library(forecast)
  library(ggplot2)
  library(gridExtra)
  library(car)
  library(nortest)
  library(AdequacyModel)
  library(lmtest)
  library(quantmod)
  library(dygraphs)
  library(lessR)
  library(PASWR2)
  library(dplyr)
  library(psych)
  library(pastecs)
  library(astsa)
  library(tseries)
  library(zoo)
  library(xts)
  library(fma)
  library(expsmooth)
  library(Quandl)
  library(fpp)
  library(urca)
  library(AER)
  library(fUnitRoots)
  library(CADFtest)
  library(fpp2)
  library(datasets)
  library(tidyverse)
  library(magrittr)
  library(ggfortify)
})

# Limpio la memoria
rm( list=ls() )
gc()

# Cargo la base AirPassengers: Total de pasajeros mensuales de 1949 a 1960
ts <- AirPassengers
head(ts)
tail(ts)

# Grafico la serie
autoplot(ts)+ ggtitle("Pasajeros") + ylab("")

# Conjunto de entrenamiento (1949-1959)
train <- window(ts,
                start = 1949, 
                end = c(1959,12))
autoplot(train)+ ggtitle("Train") + ylab("")

# Conjunto de prueba (1960)
test <- window(ts, 
               start = 1960)

# Grafico FAS, FAC y FACP del conjunto de entrenamiento
acf(train,type = "covariance",plot = T)
ggAcf(train) + ggtitle("FAC Train") # La serie decrece linealmente, lo que indica que no es estacionaria
ggAcf(train,type = "partial") + ggtitle("FACP Train")

# Planteo un modelo ingenuo y un modelo ingenuo con estacionalidad
mi <- rwf(train,h = 12)
mie <- snaive(train,h = 12)

########## Analizamos el modelo ingenuo ########## 
summary(mi)
autoplot(mi)
checkresiduals(mi)

# Verificamos los residuos
residuos1 <- resid(mi)

# Cargo la siguiente función que realiza test de normalidad
Normality_Test <- function(ts,type = c("JB", "AD", "SW")){
  require(tseries)
  require(nortest)
  if(type == "JB"){
    p_val = jarque.bera.test(ts)$p.value
    stat  = jarque.bera.test(ts)$statistic
  } else if(type == "AD"){
    p_val = ad.test(ts)$p.value
    stat  = ad.test(ts)$statistic
  } else {
    p_val = shapiro.test(ts)$p.value
    stat  = shapiro.test(ts)$statistic
  }
  
  table = data.frame(P_Value = p_val,
                     Statistic = stat)
  return(table)
}

# Verifico la normalidad de los residuos
Normality_Test(na.omit(residuos1),type = "JB") # p-value: 0.1067532
Normality_Test(na.omit(residuos1),type = "AD") # p-value: 0.3653374
Normality_Test(na.omit(residuos1),type = "SW") # p-value: 0.2315904

# Con ninguno de los tres test rechazo el supuesto de normalidad. Pruebo la aleatoriedad de los residuos

# Cargo la siguiente función de incorrelación que realiza un test de Ljung-Box o Box-Pierce para distintos lags
Incorrelation <- function(ts, type = c("Ljung-Box","Box-Pierce"), fitdf = 0){
  p_ljung_box = NULL
  s_ljung_box = NULL
  for(i in 0:(length(ts)/4)){
    p_ljung_box[i] = Box.test(ts,lag = i,type = type,fitdf = fitdf)$p.value
    s_ljung_box[i] = Box.test(ts,lag = i,type = type,fitdf = fitdf)$statistic
  }
  table = data.frame(j = 1:(length(ts)/4),
                     P_Value = p_ljung_box,
                     Statistic = s_ljung_box)
  return(table)
}

# Planteo el test de Ljung-Box. Si rechazo H0 significa que hay coeficientes de autocorrelación distintos a cero
Incorrelation(residuos1,"Ljung-Box")
inco_wn = Incorrelation(residuos1,"Ljung-Box")

# Grafico los p-value para disntitos lags
autoplot(ts(inco_wn$P_Value)) +
  ggtitle("Test de Ljung-Box", subtitle = "P-Value") +
  ylab("")

# todos los p-value son menores a 0.05. Rechazo el supuesto de incorrelación (independencia)

########## Analizamos el modelo ingenuo con estacionalidad ########## 
summary(mie)
autoplot(mie)
checkresiduals(mie)

# Verificamos los residuos
residuos2 <- resid(mie)

# Verifico la normalidad de los residuos
Normality_Test(na.omit(residuos2),type = "JB") # p-value: 0.3497492
Normality_Test(na.omit(residuos2),type = "AD") # p-value: 0.3658841
Normality_Test(na.omit(residuos2),type = "SW") # p-value: 0.3162318

# Con ninguno de los tres test rechazo el supuesto de normalidad. Pruebo la aleatoriedad de los residuos

# Planteo el test de Ljung-Box. Si rechazo H0 significa que hay coeficientes de autocorrelación distintos a cero
Incorrelation(residuos2,"Ljung-Box")
inco_wn = Incorrelation(residuos2,"Ljung-Box")

# Grafico los p-value para disntitos lags
autoplot(ts(inco_wn$P_Value)) +
  ggtitle("Test de Ljung-Box", subtitle = "P-Value") +
  ylab("")

# todos los p-value son menores a 0.05. Rechazo el supuesto de incorrelación (independencia)

# Verifico la precisión
accuracy(mi, test)
accuracy(mie, test)

# Realizo un CV con el modelo snaive que fue el que mejor resultados arrojó
CV <- tsCV(ts, snaive, drift=TRUE, h=8)

# Calculo el MSE
mse <- colMeans(CV^2, na.rm = T)

# Plotear los valores MSE contra el horizonte de predicciÃ³n
data.frame(h = 1:8, MSE = mse) %>%
  ggplot(aes(x = h, y = MSE)) + geom_point()
