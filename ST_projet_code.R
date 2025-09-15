# install.packages("zoo")
# install.packages("tseries")
# install.packages("fUnitRoots")
# install.packages("ellipsis")
# install.packages("carData")
# install.packages("car")

# Chargement des bibliothèques nécessaires
library(tseries)     
library(stats)       
library(forecast)    
library(fUnitRoots)  
library(ggplot2)     
library(readr)       
library(zoo)         # Gestion facile des séries temporelles
library(ellipse)     # Fonction pour dessiner des ellipses

# https://www.insee.fr/fr/statistiques/serie/010537812#Tableau
### Importation des données ###

path_data <- "/Users/maximecoppa/Desktop/TP noté 2"
setwd(path_data)
datafile <- "monthly_values.csv" # Définir le fichier de données

# Importer un fichier .csv dans un objet de classe data.frame
data <- read.csv(datafile, sep = ";") 
data.source <- zoo(data[[2]])
head(data)

### PARTIE 1 : les données ###

# Prétraitement des données #
# QUESTION 1 - 2 - 3 #

# Transformation des données afin de pouvoir les étudier
data <- data[4:nrow(data),] # Suppression des trois premières lignes (en-têtes ou informations non pertinentes)
data <- data[,-3] # Suppression de la dernière colonne (non utilisée)
colnames(data) <- c("date", "index_production") # Renommer les colonnes pour plus de clarté

head(data)

# Suppression des 4 dernières lignes (pour comparer plus tard nos prévisions aux valeurs réelles)
data <- data[1:(nrow(data)-4), ]
data

# Conversion des dates
dates <- as.yearmon(seq(from = 2013, to = 2024 + 04/12, by = 1/12)) 
index <- zoo(as.numeric(data$index_production), order.by = dates)

# Différenciation de la série pour la rendre stationnaire
dindex <- diff(index, 1)



# Tracé de la série index
plot.ts(index, xlab = "Dates", ylab = "Index de production") 

# Tracé de la série différenciée dindex 
plot.ts(dindex, xlab = "Dates", ylab = "Index différencié de production de tapis et moquettes")

### Tests de stationnarité
### Nous allons maintenant effectuer trois tests (test de Dickey-Fuller augmenté, test de Phillips-Perron et test KPSS) afin de vérifier la stationnarité des séries index et dindex.

## Tests sur la série index (série non différenciée)

# Test de Dickey-Fuller augmenté

# Avant de réaliser les tests de racine unitaire, nous devons vérifier s'il y a un intercept et/ou une tendance linéaire non nulle.
# Régressons index sur ses dates pour vérifier
summary(lm(index ~ dates))


# Test ADF sans constante et sans tendance
adf <- adfTest(index, lag=0, type="nc") 
print(adf)

# Fonction pour les tests d'autocorrélation des résidus
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l <= fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l, "pval"=pval))
  })
  return(t(pvals))
}

# Tests d'autocorrélation des résidus de l'ADF
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))

# Nous rejetons l'absence d'autocorrélation des résidus au moins une fois, invalidant ainsi le test ADF sans décalages.
# Ajoutons des décalages de ∆Xt jusqu'à ce que les résidus ne soient plus autocorrélés.

adfTest_valid <- function(series, kmax=24, adftype="ct") {
  k <- 0
  noautocorr <- 0
  
  while (noautocorr == 0 && k <= kmax) {
    cat(paste0("ADF avec ", k, " décalages : résidus OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    
    residuals <- adf@test$lm$residuals
    
    pvals <- Qtests(residuals, 24, fitdf = length(adf@test$lm$coefficients))[, "pval"]
    
    if (all(pvals < 0.05, na.rm=TRUE)) {
      noautocorr <- 1
      cat("OK \n")
    } else {
      cat("non \n")
      k <- k + 1
    }
  }
  
  return(adf)
}

adf <- adfTest_valid(index, 24, adftype="nc")

# Nous avons dû considérer N décalages sur le test ADF pour éliminer l'autocorrélation des résidus.
# La racine unitaire n'est pas rejetée au niveau de 95% pour la série en niveaux, la série est donc au moins I(1).

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
adf
# La racine unitaire n'est pas rejetée au niveau de 95% pour la série en niveaux, la série est donc au moins I(1).

# Test de Phillips-Perron : l'hypothèse nulle n'est pas rejetée au niveau de 95%, donc la série index n'est pas stationnaire.
pp.test(index)

# Test KPSS : l'hypothèse nulle est rejetée au niveau de 95%, ce qui confirme que la série index n'est pas stationnaire.
kpss.test(index)

# Test adf : nous testons pour vérifier nos résultats 
adf.test(index)

## Tests sur la série 'dindex' (la différenciée d'ordre 1)

# La représentation graphique précédente semble montrer l'absence de constante et de tendance.
# Vérifions avec une régression :

# Régression pour vérifier la série différenciée
summary(lm(dindex ~ dates[-1]))

adf <- adfTest_valid(dindex, 24, "nc")
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
adf
# Le test rejette l'hypothèse de racine unitaire (p-value < 0.05), nous dirons donc que la série différenciée est "stationnaire".
# Dindex est donc I(1).

pp.test(dindex)
# Le test rejette l'hypothèse de racine unitaire (p-value < 0.05), confirmant ainsi que la série dindex est stationnaire.

kpss.test(dindex)
# Le test accepte l'hypothèse nulle (p-value > 0.10), confirmant ainsi que la série dindex est stationnaire.

### PARTIE 2 : Modèles ARMA ###

#### IDENTIFICATION des ARMA(p,q) ####
# Questions 4 et 5 #

# Regardons les autocorrelogrammes complets et partiels de dindex pour trouver p et q.

par(mfrow = c(1, 2))
pacf(dindex, 24)
acf(dindex, 24) # Nous regardons les décalages de deux ans

pmax <- 4
qmax <- 2

# Fonction de test de la signification statistique des coefficients
signif <- function(estim) {
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef / se
  pval <- (1 - pnorm(abs(t))) * 2
  return(rbind(coef, se, pval))
}

# Le modèle est valide si les résidus ne sont pas autocorrélés. Nous pouvons tester cela en utilisant le test de Ljung-Box de l'hypothèse nulle de nullité conjointe des autocorrélations jusqu'à un ordre donné k.

# Fonction pour estimer un arima et vérifier son ajustement et sa validité (avec tests de signification et de Ljung-Box)
modelchoice <- function(p, q, data = dindex, k = 24) {
  estim <- try(arima(data, c(p, 0, q), optim.control = list(maxit = 20000)))
  if (class(estim) == "try-error") return(c("p" = p, "q" = q, "arsignif" = NA, "masignif" = NA, "resnocorr" = NA, "ok" = NA))
  arsignif <- if (p == 0) NA else signif(estim)[3, p] <= 0.05
  masignif <- if (q == 0) NA else signif(estim)[3, p + q] <= 0.05
  resnocorr <- sum(Qtests(estim$residuals, 24, length(estim$coef) - 1)[,2] <= 0.05, na.rm = TRUE) == 0
  checks <- c(arsignif, masignif, resnocorr)
  ok <- as.numeric(sum(checks, na.rm = TRUE) == (3 - sum(is.na(checks))))
  return(c("p" = p, "q" = q, "arsignif" = arsignif, "masignif" = masignif, "resnocorr" = resnocorr, "ok" = ok))
}

# Fonction pour estimer et vérifier tous les arima(p,q) avec p <= pmax et q <= qmax
armamodelchoice <- function(pmax, qmax) {
  pqs <- expand.grid(0:pmax, 0:qmax)
  t(apply(matrix(1:dim(pqs)[1]), 1, function(row) {
    p <- pqs[row, 1]
    q <- pqs[row, 2]
    cat(paste0("Calcul du ARMA(", p, ",", q, ") \n"))
    modelchoice(p, q)
  }))
}

armamodels <- armamodelchoice(pmax, qmax) # Estimation de tous les arima

selec <- armamodels[armamodels[,"ok"] == 1 & !is.na(armamodels[,"ok"]),] # Modèles bien ajustés et validés
selec
# Nous avons 3 modèles bien ajustés et valides : ARMA (2,0), ARMA(1,1) et ARMA(0, 2)

pqs <- apply(selec, 1, function(row) list("p" = as.numeric(row[1]), "q" = as.numeric(row[2]))) # Création d'une liste des ordres p et q des modèles candidats
names(pqs) <- paste0("arma(", selec[, 1], ",", selec[, 2], ")") # Renommer les éléments de la liste
models <- lapply(pqs, function(pq) arima(dindex, c(pq[["p"]], 0, pq[["q"]]))) # Création d'une liste des modèles candidats estimés
vapply(models, FUN.VALUE = numeric(2), function(m) c("AIC" = AIC(m), "BIC" = BIC(m))) # Calcul des AIC et BIC des modèles candidats
# L'ARIMA(1,1,1) minimise les critères d'information AIC et BIC.

arima111 <- arima(index, c(1, 1, 1), include.mean = FALSE)
model_arima111 <- arima(index, order = c(1, 1, 1))

# Afficher le modèle complet pour plus de détails
summary(model_arima111)

### PARTIE 3 : Prévision ###

# Prévision de 4 valeurs futures

# Questions 6, 7 et 8 #

# Nous traçons la prévision grâce à la méthode forecast
prev <- forecast(arima111, h = 2, level = 95)
plot(prev, xlim = c(2018, 2024 + 7/12))
prev

# Extraction des coefficients du modèle et de la variance des résidus
arima111

# Extraction des prévisions et de leurs variances
XT1 <- prev$mean[1]
XT2 <- prev$mean[2]
sigma2 <- arima111$sigma2

# La commande coef donne d'abord les coefficients AR puis les coefficients MA
phi_1 <- as.numeric(arima111$coef[1])
psi_1 <- as.numeric(arima111$coef[2])

# Calcul de la matrice de covariance Sigma
Sigma <- matrix(c(sigma2, (phi_1 - psi_1) * sigma2, (phi_1 - psi_1) * sigma2, (phi_1 - psi_1)^2 * sigma2 + sigma2), ncol = 2)

# Ensuite, nous représentons la région de confiance bivariée à 95%

require(ellipse)
require(ellipsis)
require(car)
library(ellipse) 


plot(XT1, XT2, xlim = c(70,160), ylim = c(50,200), 
     xlab = "Prévision de X(T+1)", 
     ylab = "Prévision de X(T+2)", 
     main = "Région de confiance bivariée à 95%")

# Ajouter l'ellipse de confiance
lines(ellipse(Sigma, centre=c(XT1, XT2)), type="l", col="red")
abline(h=XT2, v=XT1)

