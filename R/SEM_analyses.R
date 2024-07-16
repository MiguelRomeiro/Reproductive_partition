# Description of this script

# This script aims to:
# Calculate effects between multiple variables and ant queen reproductive output using bivariate and piecewise statistical modelling.

# Developed by Miguel Pereira-Romeiro (Details omitted for double anonymised reviewing)
# Authors: Miguel Pereira-Romeiro, Marianne Azevedo-Silva, Henrique Silva Florindo, Gustavo Maruyama Mori, Paulo Silva Oliveira, Anette


# install.packages("raster")
# install.packages("rgeos")
# install.packages("piecewiseSEM")
# install.packages("terra")
# installed.packages("ggplot2")
# install.packages("dplyr")
# install.packages("car")
# install.packages("DHARMa")
# install.packages("performance")
# install.packages("gtools")
# install.packages("corrplot")
# install.packages("effectsize")
# install.packages("renv")

##
#Image generation
# plot(landuse)
# points(coord_rufi, col = as.factor(data_reng$nq), pch = 16)
# points(coord_reng, col = as.factor(data_reng$nq), pch = 17)
# table(values(landuse))

#### Bivariate Analysis ####  POSSÍVELMENTE REMOVER

# summary(glm(nfp ~ nq, family=binomial , weight = nw , data = data_rufi))
# summary(glm(nfp ~ nm, family=binomial , weight = nw , data = data_rufi))
# summary(glm(nfp ~ PHt, family=binomial , weight = nw , data = data_rufi))
# summary(glm(nfp ~ cert, family=binomial , weight = nw , data = data_rufi))
# 
# summary(glm(nfp ~ nq, family=binomial , weight = nw , data = data_reng))
# summary(glm(nfp ~ nm, family=binomial , weight = nw , data = data_reng))
# summary(glm(nfp ~ PHt, family=binomial , weight = nw , data = data_reng))
# summary(glm(nfp ~ cert, family=binomial , weight = nw , data = data_reng))
# 
# ## Proportion of workers x Number of queens per nest ##
# 
# #rufipes
# 
# plot_worker_polygyny_rufi <- ggplot() +
#   geom_smooth(data = data_rufi, aes(x = nq, y = nfp), method = "glm", method.args = list(family="binomial"), se = TRUE, alpha = 0.20, linetype = 1, linewidth = 0.3, fill = "#ce4b4b", col="#660000") +
#   geom_jitter(data = data_rufi, aes(x = nq, y = nfp), shape = 21, size = 2,  color = "#660000", fill = "#bf0d0d", width = 0.07, height = 0) +
#   labs(y= "% of workers per queen", x = "Number of queens per nest") +
#   theme_article()+
#   ylim(0,1)
# 
# plot_worker_polygyny_rufi
# 
# dev.off()
# 
# #renggeri
# 
# svg("nfp_nq_reng.svg")
# 
# plot_worker_polygyny_reng <- ggplot() +
#   geom_smooth(data = data_reng, aes(x = nq, y = nfp), method = "glm", method.args = list(family="binomial"), se = TRUE, alpha = 0.20, linetype = 1, linewidth = 0.3, fill = "#a65c04", col="#a65c04") +
#   geom_jitter(data = data_reng, aes(x = nq, y = nfp), shape = 21, size = 2,  color = "#502d02", fill = "#a65c04", width = 0.07, height = 0) +
#   theme_article()+
#   labs(y= "% of workers per queen", x = "Number of queens per nest") +
#   ylim(0,1)
# 
# plot_worker_polygyny_reng
# 
# dev.off()
# 
# ## Proportion of workers x Number of male mates ##
# 
# #rufipes
# 
# svg("nfp_nm_rufi.svg")
# 
# plot_worker_polyandry_rufi <- ggplot() +
#   geom_smooth(data = data_rufi, aes(x = nm, y = nfp), method = "glm", method.args = list(family="binomial"), se = TRUE, alpha = 0.20, linetype = 1, linewidth = 0.3, fill = "#ce4b4b", col="#660000") +
#   geom_jitter(data = data_rufi, aes(x = nm, y = nfp), shape = 21, size = 2,  color = "#660000", fill = "#bf0d0d", width = 0.07, height = 0) +
#   labs(y= "% of workers per queen", x = "Number of male mates per queen") +
#   theme_article() +
#   ylim(0,1)
# 
# plot_worker_polyandry_rufi
# 
# dev.off()
# 
# #renggeri
# 
# svg("nfp_nm_reng.svg")
# 
# plot_worker_polyandry_reng <- ggplot() +
#   geom_smooth(data = data_reng, aes(x = nm, y = nfp), method = "glm", method.args = list(family="binomial"), se = TRUE, alpha = 0.20, linetype = 1, linewidth = 0.3, fill = "#a65c04", col="#a65c04") +
#   geom_jitter(data = data_reng, aes(x = nm, y = nfp), shape = 21, size = 2,  color = "#502d02", fill = "#a65c04", width = 0.07, height = 0) +
#   labs(y= "% of workers per queen", x = "Number of male mates per queen") +
#   theme_article()+
#   ylim(0,1)
# 
# plot_worker_polyandry_reng
# 
# dev.off()
# 
# ## Proportion of workers x Heterozigosity ##
# 
# #rufipes
# 
# svg("nfp_PHt_rufi.svg")
# 
# plot_worker_het_rufi <- ggplot() +
#   geom_smooth(data = data_rufi, aes(x = PHt, y = nfp), method = "glm", method.args = list(family="binomial"), se = TRUE, alpha = 0.20, linetype = 1, linewidth = 0.3, fill = "#ce4b4b", col="#660000") +
#   geom_jitter(data = data_rufi, aes(x = PHt, y = nfp), shape = 21, size = 2,  color = "#660000", fill = "#bf0d0d", width = 0.007, height = 0) +
#   labs(y= "% of workers per queen", x = "Proportion of heterozygous loci per queen") +
#   theme_article() +
#   ylim(0,1)
# 
# plot_worker_het_rufi
# 
# dev.off()
# 
# #renggeri
# 
# svg("nfp_PHt_reng.svg")
# 
# plot_worker_het_reng <- ggplot() +
#   geom_smooth(data = data_reng, aes(x = PHt, y = nfp), method = "glm", method.args = list(family="binomial"), se = TRUE, alpha = 0.20, linetype = 1, linewidth = 0.3, fill = "#a65c04", col="#a65c04") +
#   geom_jitter(data = data_reng, aes(x = PHt, y = nfp), shape = 21, size = 2,  color = "#502d02", fill = "#a65c04", width = 0.007, height = 0) +
#   labs(y= "% of workers per queen", x = "Proportion of heterozygous loci per queen") +
#   theme_article() +
#   ylim(0,1)
# 
# plot_worker_het_reng
# 
# dev.off()
# 
# ## Proportion of workers x Proportion of Cerrado ##
# 
# #rufipes
# 
# svg("nfp_cert_rufi.svg")
# 
# plot_worker_cert_rufi <- ggplot() +
#   geom_jitter(data = data_rufi, aes(x = cert, y = nfp), shape = 21, size = 2,  color = "#660000", fill = "#bf0d0d", width = 0.007, height = 0) +
#   labs(y= "% of workers per queen", x = "Proportion of Cerrado around nest") +
#   theme_article() +
#   ylim(0,1)
# 
# plot_worker_cert_rufi
# 
# dev.off()
# 
# #renggeri
# 
# svg("nfp_cert_reng.svg")
# 
# plot_worker_cert_reng <- ggplot() +
#   geom_smooth(data = data_reng, aes(x = cert, y = nfp), method = "glm", method.args = list(family="binomial"), se = TRUE, alpha = 0.20, linetype = 1, linewidth = 0.3, fill = "#a65c04", col="#a65c04") +
#   geom_jitter(data = data_reng, aes(x = cert, y = nfp), shape = 21, size = 2,  color = "#502d02", fill = "#a65c04", width = 0.007, height = 0) +
#   labs(y= "% of workers per queen", x = "Proportion of Cerrado around nest") +
#   theme_article() +
#   ylim(0,1)
# 
# plot_worker_cert_reng
# 
# dev.off()


#### PieceWiseSEM ####

#C. rufipes

nq_ru <- glm(nq ~ cert + nm, family = poisson, data=data_rufi)
nfp_ru <- glm(nfp ~ nq + nm + PHt, family=binomial , weight = nw , data= data_rufi)

sem_rufi <- psem(nq_ru, nfp_ru)
summary(sem_rufi, .progressBar = F)
summary(sem_rufi, .progressBar = F)$AIC

# Remove non-significant terms stepwise

# Polyandry (nm) removed from polygyny (nq) model
nq_ru2 <- glm(nq ~ cert , data=data_rufi, family = poisson)
nfp_ru2 <- glm(nfp ~ nq + nm + PHt, family=binomial , weight = nw, data=data_rufi)

sem_rufi2 <- psem(nq_ru2, nfp_ru2)
summary(sem_rufi2, .progressBar = F)
summary(sem_rufi2, .progressBar = F)$AIC

# Polygyny (nq) model removed altogether
nfp_ru3 <- glm(nfp ~ nq + nm + PHt, family=binomial , weight = nw, data=data_rufi)

sem_rufi3 <- psem(nfp_ru3)
summary(sem_rufi3, .progressBar = F)
summary(sem_rufi3, .progressBar = F)$AIC
# Only the linear regression
summary(nfp_ru3, .progressBar = F)

# C. renggeri

nq_re <- glm(nq ~ cert + nm, family = poisson, data = data_reng)
nfp_re <- glm(nfp ~ nq + nm + PHt, family=binomial , weight = nw , data= data_reng)

sem_reng <- psem(nq_re, nfp_re)
summary(sem_reng, .progressBar = F)
summary(sem_reng, .progressBar = F)$AIC

# Remove non-significant terms stepwise

# Polyandry (nm) removed from polygny (nq) model
nq_re2 <- glm(nq ~ cert,  family = poisson, data = data_reng)
nfp_re2 <- glm(nfp ~ nq + nm + PHt, family=binomial , weight = nw , data= data_reng)

sem_reng2 <- psem(nq_re2, nfp_re2)
summary(sem_reng2, .progressBar = F)
summary(sem_reng2, .progressBar = F)$AIC

# Polyandry (nm) removed from reproductive partitioning (nfp) model
nq_re3 <- glm(nq ~ cert,  family = poisson, data = data_reng)
nfp_re3 <- glm(nfp ~ nq + PHt, family=binomial , weight = nw , data= data_reng)

sem_reng3 <- psem(nq_re3, nfp_re3)
summary(sem_reng3, .progressBar = F)
summary(sem_reng3, .progressBar = F)$AIC


#### Diagnostics ####

## Collinearity ##
vif(nq_ru)
vif(nfp_ru)

vif(nq_re)
vif(nfp_re)

check_collinearity(nq_ru)		
check_collinearity(nfp_ru)		
simulateResiduals(nq_ru,plot=T)
simulateResiduals(nfp_ru,plot=T)

check_collinearity(nq_re)		
check_collinearity(nfp_re)		
simulateResiduals(nq_re,plot=T)
simulateResiduals(nfp_re,plot=T)