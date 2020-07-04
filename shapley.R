library(ggplot2)
library(hrbrthemes)
extrafont::loadfonts()
library(viridis)

##### Example ######
#This code replicates calculations in the Excel example
# https://www.real-statistics.com/multiple-regression/shapley-owen-decomposition/
# also see
# https://wernerantweiler.ca/blog.php?item=2014-10-10

rm(list = ls())

vars = letters[1:3]

n = length(vars) #total number of variables
s = 0:n #subset sizes of the variables
weights = factorial(s)*factorial(n-s-1) / factorial(n)

#list of combinations
lcomb = lapply(s, function(i) combn(vars, i, simplify = FALSE))

#list of R2
lr2 = list(0, #intercept-only model
           list(0.844964, 0.703041, 0.733757),
           list(0.882822, 0.849617, 0.763974),
           0.986125 #full model
)

shapley = rep(0, n)
for(v in 1:n){ # v = 1; i = 1; j = 1
    var = vars[v]
    for(i in 1:(length(lr2) - 1)){
        for(j in 1:length(lcomb[[i]])){
            tmp = which(sapply(lcomb[[i + 1]], function(X) setequal(X, c(lcomb[[i]][[j]], var))))
            if (length(tmp) > 0) {
                shapley[v] = shapley[v] + weights[i] * (lr2[[i + 1]][[tmp]] - lr2[[i]][[j]])
            }
        }
    }
}
shapley
sum(shapley)


##### CPUE north ######
rm(list = ls())
load("./dataderived/tmp_cpueN_2020-04-13.RData")
library(mgcv)

K = 5
f = as.formula(paste0(RESPONSE, " ~ s(", paste0(gs_all, collapse = ", k = K) + s("), ", k = K)"))
gam_1 = gamfit = mgcv::gamm(f
                            , select = TRUE
                            , bs = "cr"
                            , method = "REML"
                            , correlation = corARMA(form = ~ 1 | AREA, p = 1, q = 0)
                            , data = DATAnoNA)
fittedGAM <- fitted(gamfit$lme)
SSres <- sum((DATAnoNA[RESPONSE] - fittedGAM)^2)
# SStot <- sum( (DATAnoNA[RESPONSE] - mean(DATAnoNA[RESPONSE]))^2 ) #do not know why gives an error
SStot <- sum( (DATAnoNA[RESPONSE] - mean(DATAnoNA$sqrtCPUE))^2 )
R2tot = 1 - SSres/SStot #non-adjusted R2

#Continue by adopting the example code from above
vars = gs_all
n = length(vars) #total number of variables
s = 0:n #subset sizes of the variables
weights = factorial(s)*factorial(n-s-1) / factorial(n) #In gamma(x + 1) : NaNs produced #is OK for the last element

#list of combinations
lcomb = lapply(s, function(i) combn(vars, i, simplify = FALSE))

#list of R2
# lr2 = as.list(rep(NA, length(lcomb)))
lr2 = lcomb #copy the list structure
lr2[[1]] = 0 #intercept-only model
lr2[[length(lcomb)]] = R2tot #full model
for(i in 2:(length(lcomb) - 1)) {
    for(j in 1:length(lcomb[[i]])) {
        f = as.formula(paste0(RESPONSE, " ~ s(", paste0(lcomb[[i]][[j]], collapse = ", k = K) + s("), ", k = K)"))
        gamfit = mgcv::gamm(f
                            , select = TRUE
                            , bs = "cr"
                            , method = "REML"
                            , correlation = corARMA(form = ~ 1 | AREA, p = 1, q = 0)
                            , data = DATAnoNA)
        fittedGAM <- fitted(gamfit$lme)
        SSres <- sum((DATAnoNA[RESPONSE] - fittedGAM)^2)
        #SStot same as before
        lr2[[i]][[j]] = 1 - SSres/SStot #non-adjusted R2
    }
}

shapley = rep(0, n); names(shapley) = vars
for(v in 1:n){ # v = 1; i = 1; j = 1
    var = vars[v]
    for(i in 1:(length(lr2) - 1)){
        for(j in 1:length(lcomb[[i]])){
            tmp = which(sapply(lcomb[[i + 1]], function(X) setequal(X, c(lcomb[[i]][[j]], var))))
            if (length(tmp) > 0) {
                shapley[v] = shapley[v] + weights[i] * (lr2[[i + 1]][[tmp]] - lr2[[i]][[j]])
            }
        }
    }
}
shapley
# AMO_DJFMA_l7y AMO_annual_l6y    Track_NE191       GSI_l12q
# 0.15955278     0.17063452     0.08348735     0.10115614
sum(shapley)
# [1] 0.5148308
R2tot
# [1] 0.5148308

###### plot ######
rm(list = ls())

shapley = c(0.15955278, 0.17063452, 0.08348735, 0.10115614)
names(shapley) = c("AMO_DJFMA_l7y", "AMO_annual_l6y", "Track_NE191", "GSI_l12q")

shapley = sort(shapley, decreasing = TRUE)
shapley = as.data.frame(shapley)
shapley$Variable = paste0(format(shapley$shapley/sum(shapley$shapley) * 100, digits = 2)
                          ,"% "
                          ,rownames(shapley))

pdf("./figures/shapley_cpueN.pdf", width = 4, height = 4)
ggplot(shapley,
       aes(x = 1, y = shapley, fill = Variable)) +
    geom_bar(position = "fill", stat = "identity", color='NA', width=0.1) +
    scale_y_continuous(labels = scales::percent) +
    labs(x="", y = "") +
    scale_fill_viridis(discrete = T) +
    theme_ipsum_pub(ticks = FALSE, grid = "Y")
dev.off()

shapley$shapley/sum(shapley$shapley) * 100
# [1] 33.14381 30.99131 19.64842 16.21646

##### CPUE south ######
rm(list = ls())
load("./dataderived/tmp_cpueS_2020-04-13.RData")
library(mgcv)

K = 5
f = as.formula(paste0(RESPONSE, " ~ s(", paste0(gs_all, collapse = ", k = K) + s("), ", k = K)"
                      , "+ Time_Block"
                      ))
gam_1 = gamfit = mgcv::gamm(f
                            , select = TRUE
                            , bs = "cr"
                            , method = "REML"
                            , correlation = corARMA(form = ~ 1 | AREA, p = 1, q = 0)
                            , data = DATAnoNA)
fittedGAM <- fitted(gamfit$lme)
SSres <- sum((DATAnoNA[RESPONSE] - fittedGAM)^2)
# SStot <- sum( (DATAnoNA[RESPONSE] - mean(DATAnoNA[RESPONSE]))^2 ) #do not know why gives an error
SStot <- sum( (DATAnoNA[RESPONSE] - mean(DATAnoNA$sqrtCPUE))^2 )
R2tot = 1 - SSres/SStot #non-adjusted R2

#Continue by adopting the example code from above
vars = c(gs_all, "Time_Block")
n = length(vars) #total number of variables
s = 0:n #subset sizes of the variables
weights = factorial(s)*factorial(n-s-1) / factorial(n) #In gamma(x + 1) : NaNs produced #is OK for the last element

#list of combinations
lcomb = lapply(s, function(i) combn(vars, i, simplify = FALSE))
length(unlist(lcomb)) #[1] 11264

#list of R2
# lr2 = as.list(rep(NA, length(lcomb)))
lr2 = lcomb #copy the list structure
lr2[[1]] = 0 #intercept-only model
lr2[[length(lcomb)]] = R2tot #full model
onlytimeblock = FALSE
for(i in 2:(length(lcomb) - 1)) {
    print(Sys.time())
    print(i)
    for(j in 1:length(lcomb[[i]])) {
        if (is.element("Time_Block", lcomb[[i]][[j]])) {
            if (length(lcomb[[i]][[j]]) > 1) { #if there are other variables except time block
                f = as.formula(paste0(RESPONSE, " ~ s(", paste0(setdiff(lcomb[[i]][[j]], "Time_Block"),
                                                                collapse = ", k = K) + s("), ", k = K)"
                                      , "+ Time_Block"
                ))
            } else { #if there is only time block
                onlytimeblock = TRUE
                f = as.formula(paste0(RESPONSE, " ~  Time_Block"))
            }
        } else { #no timeblock among the variables
            f = as.formula(paste0(RESPONSE, " ~ s(", paste0(lcomb[[i]][[j]], collapse = ", k = K) + s("), ", k = K)"))
        }
        if (onlytimeblock) {
            gamfit = gls(f
                         ,correlation = corARMA(form = ~ 1 | AREA, p = 1, q = 0)
                         ,data = DATAnoNA)
            fittedGAM = fitted(gamfit)
        } else {
            gamfit = mgcv::gamm(f
                                , select = TRUE
                                , bs = "cr"
                                , method = "REML"
                                , correlation = corARMA(form = ~ 1 | AREA, p = 1, q = 0)
                                , data = DATAnoNA)
            fittedGAM <- fitted(gamfit$lme)
        }
        SSres <- sum((DATAnoNA[RESPONSE] - fittedGAM)^2)
        #SStot same as before
        lr2[[i]][[j]] = 1 - SSres/SStot #non-adjusted R2
        onlytimeblock = FALSE
    }
}

shapley = rep(0, n); names(shapley) = vars
for(v in 1:n){ # v = 1; i = 1; j = 1
    var = vars[v]
    for(i in 1:(length(lr2) - 1)){
        for(j in 1:length(lcomb[[i]])){
            tmp = which(sapply(lcomb[[i + 1]], function(X) setequal(X, c(lcomb[[i]][[j]], var))))
            if (length(tmp) > 0) {
                shapley[v] = shapley[v] + weights[i] * (lr2[[i + 1]][[tmp]] - lr2[[i]][[j]])
            }
        }
    }
}
shapley
# AREA_CENT_LAT     AMO_DJFMA_l7y    AMO_annual_l2y  FC_Transport_l3m    AMO_annual_l4y  FC_Transport_l2m            avgSST  FC_Transport_l4m
# 0.050069097       0.062007774       0.048408587       0.027572503       0.049859597       0.012864142       0.007139520       0.007346823
# FC_Transport_l7m FC_Transport_l11m        Time_Block
# 0.001518552       0.001750884       0.082851837
sum(shapley)
# [1] 0.3513893
R2tot
# [1] 0.3513893
# save(shapley, file = "./dataderived/shapley_cpueS.RData")

###### plot ######
rm(list = ls())
load("./dataderived/shapley_cpueS.RData")

shapley = sort(shapley, decreasing = TRUE)
shapley = as.data.frame(shapley)
shapley$Variable = paste0(format(shapley$shapley/sum(shapley$shapley) * 100, digits = 2)
                          ,"% "
                          ,rownames(shapley))

pdf("./figures/shapley_cpueS.pdf", width = 4, height = 4)
ggplot(shapley,
       aes(x = 1, y = shapley, fill = Variable)) +
    geom_bar(position = "fill", stat = "identity", color='NA', width=0.1) +
    scale_y_continuous(labels = scales::percent) +
    labs(x="", y = "") +
    scale_fill_viridis(discrete = T) +
    theme_ipsum_pub(ticks = FALSE, grid = "Y")
dev.off()

shapley$shapley/sum(shapley$shapley) * 100
# [1] 23.5783597 17.6464596 14.2488957 14.1892751 13.7763400  7.8467107  3.6609373  2.0907930  2.0317978  0.4982747  0.4321564


##### Landings north ######
rm(list = ls())
load("./dataderived/tmp_2020-04-11.RData")
library(mgcv)

K = 5
gs_all = c("AMO_DJFMA_l7", "NAO_DJF_st_l3", "NAO_DJF_st_l4") #backward-selected
f = as.formula(paste0(RESPONSE, " ~ s(", paste0(gs_all, collapse = ", k = K) + s("), ", k = K)"
                      # , "+ Time_Block"
))
gam_1 = gamfit = mgcv::gamm(f
                            , select = TRUE
                            , bs = "cr"
                            , method = "REML"
                            , correlation = corARMA(form = ~ 1, p = 1, q = 0)
                            , data = DATAnoNA)
fittedGAM <- fitted(gamfit$lme)
SSres <- sum((DATAnoNA[RESPONSE] - fittedGAM)^2)
# SStot <- sum( (DATAnoNA[RESPONSE] - mean(DATAnoNA[RESPONSE]))^2 ) #do not know why gives an error
SStot <- sum( (DATAnoNA[RESPONSE] - mean(DATAnoNA$sqrtLandings))^2 )
R2tot = 1 - SSres/SStot #non-adjusted R2

#Continue by adopting the example code from above
vars = c(gs_all) #, "Time_Block"
n = length(vars) #total number of variables
s = 0:n #subset sizes of the variables
weights = factorial(s)*factorial(n-s-1) / factorial(n) #In gamma(x + 1) : NaNs produced #is OK for the last element

#list of combinations
lcomb = lapply(s, function(i) combn(vars, i, simplify = FALSE))
length(unlist(lcomb)) #[1] 12

#list of R2
# lr2 = as.list(rep(NA, length(lcomb)))
lr2 = lcomb #copy the list structure
lr2[[1]] = 0 #intercept-only model
lr2[[length(lcomb)]] = R2tot #full model
onlytimeblock = FALSE
for(i in 2:(length(lcomb) - 1)) {
    print(Sys.time())
    print(i)
    for(j in 1:length(lcomb[[i]])) {
        if (is.element("Time_Block", lcomb[[i]][[j]])) {
            if (length(lcomb[[i]][[j]]) > 1) { #if there are other variables except time block
                f = as.formula(paste0(RESPONSE, " ~ s(", paste0(setdiff(lcomb[[i]][[j]], "Time_Block"),
                                                                collapse = ", k = K) + s("), ", k = K)"
                                      , "+ Time_Block"
                ))
            } else { #if there is only time block
                onlytimeblock = TRUE
                f = as.formula(paste0(RESPONSE, " ~  Time_Block"))
            }
        } else { #no timeblock among the variables
            f = as.formula(paste0(RESPONSE, " ~ s(", paste0(lcomb[[i]][[j]], collapse = ", k = K) + s("), ", k = K)"))
        }
        if (onlytimeblock) {
            gamfit = gls(f
                         ,correlation = corARMA(form = ~ 1, p = 1, q = 0)
                         ,data = DATAnoNA)
            fittedGAM = fitted(gamfit)
        } else {
            gamfit = mgcv::gamm(f
                                , select = TRUE
                                , bs = "cr"
                                , method = "REML"
                                , correlation = corARMA(form = ~ 1, p = 1, q = 0)
                                , data = DATAnoNA)
            fittedGAM <- fitted(gamfit$lme)
        }
        SSres <- sum((DATAnoNA[RESPONSE] - fittedGAM)^2)
        #SStot same as before
        lr2[[i]][[j]] = 1 - SSres/SStot #non-adjusted R2
        onlytimeblock = FALSE
    }
}

shapley = rep(0, n); names(shapley) = vars
for(v in 1:n){ # v = 1; i = 1; j = 1
    var = vars[v]
    for(i in 1:(length(lr2) - 1)){
        for(j in 1:length(lcomb[[i]])){
            tmp = which(sapply(lcomb[[i + 1]], function(X) setequal(X, c(lcomb[[i]][[j]], var))))
            if (length(tmp) > 0) {
                shapley[v] = shapley[v] + weights[i] * (lr2[[i + 1]][[tmp]] - lr2[[i]][[j]])
            }
        }
    }
}
shapley
# AMO_DJFMA_l7 NAO_DJF_st_l3 NAO_DJF_st_l4
# 0.16616301    0.06777632    0.08260306
sum(shapley)
# [1] 0.3165424
R2tot
# [1] 0.3165424
# save(shapley, file = "./dataderived/shapley_landN.RData")

###### plot ######
rm(list = ls())
load("./dataderived/shapley_landN.RData")

shapley = sort(shapley, decreasing = TRUE)
shapley = as.data.frame(shapley)
shapley$Variable = paste0(format(shapley$shapley/sum(shapley$shapley) * 100, digits = 2)
                          ,"% "
                          ,rownames(shapley))

pdf("./figures/shapley_landN.pdf", width = 4, height = 4)
ggplot(shapley,
       aes(x = 1, y = shapley, fill = Variable)) +
    geom_bar(position = "fill", stat = "identity", color='NA', width=0.1) +
    scale_y_continuous(labels = scales::percent) +
    labs(x="", y = "") +
    scale_fill_viridis(discrete = T) +
    theme_ipsum_pub(ticks = FALSE, grid = "Y")
dev.off()

shapley$shapley/sum(shapley$shapley) * 100
# [1] 52.49313 26.09542 21.41145

