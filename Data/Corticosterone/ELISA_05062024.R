
# Packages ----------------------------------------------------------------

library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)



# Data --------------------------------------------------------------------

setwd("/Users/luisangeltrujillovillarreal/Documents/ELISA_corticosterona/")

database <- read_csv("database_cort_jalil.csv")
database <- as.data.frame(database)

databaseB5 <- database %>% 
  filter(Batch == "B5")
databaseB8 <- database %>% 
  filter(Batch == "B8")

calibration <- read_csv("Calibration_B5_B8.csv")
calibration <- as.data.frame(calibration)
calibrationB5 <- calibration %>% 
                  filter(Batch == "B5")
calibrationB8 <- calibration %>% 
  filter(Batch == "B8")


# 4pl functions ------------------------------------------------------------

model4pl <- function(Concentration, Background, Mid, Slope, Bmax) {
  
  Bmax + ((Background - Bmax) / (1 + ((Concentration/Mid)^Slope)))
  
}

#Calculate concentration of samples

CalcConc <- function(Background, Mid, Slope, Bmax, y) {
  as.numeric(Mid * ((Background - Bmax)/(y - Bmax) - 1)^(1/Slope))
}

# plot standards -----------------------------------------------------------

ggplot(data = calibrationB5) +
  geom_point( aes(Concentration,
                  Abs )) +
  scale_x_log10()

ggplot(data = calibrationB8) +
  geom_point( aes(Concentration,
                  Abs )) +
  scale_x_log10()


# Batch 5 -----------------------------------------------------------------


#MODEL FITTING
#Adjust values by plot above

fit <- nls(Abs ~ model4pl(Concentration, Background, Mid, Slope, Bmax),
           data = calibrationB5,
           start = c(Background=0, Mid=100, Slope=1, Bmax=4),
           control = nls.control(maxiter=1000, warnOnly=TRUE) )
cor(calibrationB5$Abs, predict(fit))




databaseB5$Concentration <- CalcConc(
  coef(fit)["Background"],
  coef(fit)["Mid"],
  coef(fit)["Slope"],
  coef(fit)["Bmax"],
  y = databaseB5$Abs
)


#PLOT

ggplot(data = calibrationB5) +
  geom_point(aes(Concentration, Abs), color = "blue")+
  scale_x_log10() +
  stat_function(data = calibrationB5, fun  = model4pl,
                args = list(Mid = coef(fit)["Mid"],
                            Background = coef(fit)["Background"],
                            Slope = coef(fit)["Slope"],
                            Bmax = coef(fit)["Bmax"])) +
  geom_point(data = databaseB5, aes(x = Concentration, y = Abs), 
             shape = 17, size = 3, color = "green") +
  labs(x = "Concentration (ng/mL)", y = "Abs") +  # Etiquetas de ejes
  theme_minimal()  # Estilo de gráfica minimalista




# Batch 8 -----------------------------------------------------------------


#MODEL FITTING
#Adjust values by plot above

fit <- nls(Abs ~ model4pl(Concentration, Background, Mid, Slope, Bmax),
           data = calibrationB8,
           start = c(Background=0, Mid=100, Slope=1, Bmax=4),
           control = nls.control(maxiter=1000, warnOnly=TRUE) )
cor(calibrationB8$Abs, predict(fit))


databaseB8$Concentration <- CalcConc(
  coef(fit)["Background"],
  coef(fit)["Mid"],
  coef(fit)["Slope"],
  coef(fit)["Bmax"],
  y = databaseB8$Abs
)


#PLOT

ggplot(data = calibrationB8) +
  geom_point(aes(Concentration, Abs), color = "blue")+
  scale_x_log10() +
  stat_function(data = calibrationB8, fun  = model4pl,
                args = list(Mid = coef(fit)["Mid"],
                            Background = coef(fit)["Background"],
                            Slope = coef(fit)["Slope"],
                            Bmax = coef(fit)["Bmax"])) +
  geom_point(data = databaseB8, aes(x = Concentration, y = Abs), 
             shape = 17, size = 3, color = "green") +
  labs(x = "Concentration (ng/mL)", y = "Abs") +  # Etiquetas de ejes
  theme_minimal()  # Estilo de gráfica minimalista




write_csv(databaseB5, "database_B5.csv")
write_csv(databaseB5, "database_B8.csv")
write_csv(calibrationB5, "calibration_B5.csv")
write_csv(calibrationB5, "calibration_B8.csv")


