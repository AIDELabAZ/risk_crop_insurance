############################
#Load Data & Variable Defs
############################

#Load Librariies
library(foreign)
library(arm)
library(stargazer)
library(R2WinBUGS)


#Load Data
data <- read.dta("rci_data.dta")
attach(data)

#Center variables in dataset
y <- lny - mean (lny)
rice_labor <- croXlnl_r - mean(croXlnl_r)
rice_fert <- croXlnf_r - mean(croXlnf_r)
rice_mech <- croXlnm_r - mean(croXlnm_r)
rice_pest <- croXlnp_r - mean(croXlnp_r)
sorghum_labor <- croXlnl_s - mean(croXlnl_s)
sorghum_fert <- croXlnf_s - mean(croXlnf_s)
sorghum_mech <- croXlnm_s - mean(croXlnm_s)
sorghum_pest <- croXlnp_s - mean(croXlnp_s)
wheat_labor <- croXlnl_w - mean(croXlnl_w)
wheat_fert <- croXlnf_w - mean(croXlnf_w)
wheat_mech <- croXlnm_w - mean(croXlnm_w)
wheat_pest <- croXlnp_w - mean(croXlnp_w)
maize_labor <- croXlnl_m - mean(croXlnl_m)
maize_fert <- croXlnf_m - mean(croXlnf_m)
maize_mech <- croXlnm_m - mean(croXlnm_m)
maize_pest <- croXlnp_m - mean(croXlnp_m)
cotton_labor <- croXlnl_c - mean(croXlnl_c)
cotton_fert <- croXlnf_c - mean(croXlnf_c)
cotton_mech <- croXlnm_c - mean(croXlnm_c)
cotton_pest <- croXlnp_c - mean(croXlnp_c)

#Define the dataset
n <- length(y)
x <-cbind(rice_labor, rice_fert, rice_mech, rice_pest, sorghum_labor, sorghum_fert, sorghum_mech, sorghum_pest, wheat_labor, wheat_fert, wheat_mech,
          wheat_pest, maize_labor, maize_fert, maize_mech, maize_pest, cotton_labor, cotton_fert, cotton_mech, cotton_pest)
L <-ncol(x)



# get parcel index variable
parcel.name <- as.vector(parcel_id)
uniq.p <- unique(parcel.name)
J <- length(uniq.p)
parcelid <- rep (NA, J)
for (i in 1:J){
  parcel_id[parcel.name==uniq.p[i]] <- i
}

# get household index variable
house.name <- as.vector(VDSA_HH_id)
uniq.h <- unique(house.name)
K <- length(uniq.h)
houseid <- rep (NA, K)
for (i in 1:K){
  VDSA_HH_id[house.name==uniq.h[i]] <- i
}

# get time index variable
time.name <- as.vector(tindex)
uniq.t <- unique(time.name)
H <- length(uniq.t)
tindex <- rep (NA, H)
for (i in 1:H){
  tindex[time.name==uniq.t[i]] <- i
}

# get village index variable
vil.name <- as.vector(vil_id)
uniq.v <- unique(vil.name)
V <- length(uniq.v)
vil_id <- rep (NA, V)
for (i in 1:V){
  vil_id[vil.name==uniq.v[i]] <- i
}

# get village-time variable
viltime.name <- as.vector(vil_t)
uniq.c <- unique(viltime.name)
C <- length(uniq.c)
viltime <- rep (NA, C)
for (i in 1:C){
  vil_t[viltime.name==uniq.c[i]] <- i
}


######################################################
#Run multilevel regressions with Rainfall as covariate
#####################################################

## Classical complete pooling regression
fit.3 <- lm(lny ~  x + factor(crop) + factor(vil_t) + factor(vil_id) + factor(tindex) -1)
display(fit.3, digits=3, detail=TRUE)

## Null model with village, season, and village-season
M0 <- lmer(lny ~ factor(crop) -1 + (1 | vil_id) + (1 | tindex) + (1 | vil_t) )
display(M0, digits=3, detail=TRUE)

## Null model with village, season, village-season, and household
M1 <- lmer(lny ~ factor(crop) -1 + (1 | VDSA_HH_id) + (1 | vil_id) + (1 | tindex) + (1 | vil_t) )
display(M1, digits=3, detail=TRUE)

## Null model with village, season, village-season, household, and parcel
M2 <- lmer(lny ~ factor(crop) -1 + (1 | parcel_id) + (1 | VDSA_HH_id) + (1 | vil_id) + (1 | tindex) + (1 | vil_t) )
display(M2, digits=3, detail=TRUE)

## Model with covariates, village, season, village-season, household, and parcel
M3 <- lmer(lny ~ x + factor(crop) -1 + (1 | parcel_id) + (1 | VDSA_HH_id)  + (1 | vil_id) + (1 | tindex) + (1 | vil_t) )
display(M3, digits=3, detail=TRUE)

stargazer(fit.3, M3, M3, title="Results of Classical Estimation of Production Function", covariate.labels=c("log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "countryIndia"), align=TRUE, no.space=TRUE)

stargazer(M0, M1, M2, M3, title="Results of Classical Estimation of Production Function", covariate.labels=c("log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "countryIndia"), align=TRUE, no.space=TRUE)


############################
#Graph Results
############################

#Define residuals
M3.resid <- residuals(M3)
M3.sd.resid <- sd(M3.resid)


#Zeta Plots
library(lattice)
library(robustlmm)

pr3 <- profile(M3, which =  1:26)

z=c(expression(tau[gamma]),expression(tau[delta]),expression(tau[theta]),expression(tau[xi]),
    expression(tau[omega]),expression(sigma[y]),"r.labr","r.fert","r.mech","r.pest",
    "s.labr","s.fert","s.mech","s.pest","w.labr","w.fert","w.mech","w.pest","m.labr",
    "m.fert","m.mech","m.pest","c.lab","c.fert","c.mech","c.pest")

xyplot(pr3, absVal = TRUE, aspect = 0.7, strip= FALSE, strip.left = strip.custom(factor.levels=z), 
       layout = c(7,4))

confint(pr3)

getME(M3, "beta")

#Graph Conditional Modes of Groups
ranef(M3)
str(ranef(M3))
print(qqmath(ranef(M3, condVar=TRUE), strip = FALSE))



############################
#BUGS multilevel regressions
############################


##
## Bugs model with village, season, and village-season
##
prod0.data <-  list ("y", "vil_t", "tindex", "vil_id", "n", "C", "V", "H")
prod0.inits <- function(){
  list (mu.a=rnorm(1), a=rnorm(C), mu.d=rnorm(1), d=rnorm(V), mu.e=rnorm(1), e=rnorm(H),
        sigma.y=runif(1), sigma.a=runif(1), sigma.d=runif(1), sigma.e=runif(1))
}

prod0.parameters <- c ("a", "d", "e", "sigma.a", "sigma.d", "sigma.e", "sigma.y")

B.0 <- bugs (prod0.data, prod0.inits, prod0.parameters, "india.0.bug", n.chains=4, 
             n.iter=10000, working.directory=NULL, clearWD=TRUE, 
             debug=TRUE )

plot(B.0)
options(max.print=1000000)
print(B.0, digits.summary = 3)


##
## Bugs model with village, season, village-season, and household
##
prod1.data <-  list ("y", "vil_t", "tindex", "vil_id", "VDSA_HH_id", "n", "C", "V", "H", "K")
prod1.inits <- function(){
  list (mu.a=rnorm(1), a=rnorm(C), mu.d=rnorm(1), d=rnorm(V), mu.e=rnorm(1), e=rnorm(H),
        mu.f=rnorm(1), f=rnorm(K), sigma.y=runif(1), sigma.a=runif(1), 
        sigma.d=runif(1), sigma.e=runif(1), sigma.f=runif(1))
}

prod1.parameters <- c ("a", "d", "e", "f", "sigma.f", "sigma.a", "sigma.d", "sigma.e", "sigma.y")

B.1 <- bugs (prod1.data, prod1.inits, prod1.parameters, "india.1.bug", n.chains=4, 
             n.iter=10000, working.directory=NULL, clearWD=TRUE, 
             debug=TRUE )

plot(B.1)
options(max.print=1000000)
print(B.1, digits.summary = 3)


##
## Bugs model with village, season, village-season, household, and parcel
##
prod2.data <-  list ("y", "vil_t", "tindex", "vil_id", "VDSA_HH_id", "parcel_id", "n", "C", "V", "H", "K", "J")
prod2.inits <- function(){
  list (mu.a=rnorm(1), a=rnorm(C), mu.d=rnorm(1), d=rnorm(V), mu.e=rnorm(1), e=rnorm(H),
        mu.f=rnorm(1), f=rnorm(K), mu.g=rnorm(1), g=rnorm(J), sigma.y=runif(1), sigma.a=runif(1), 
        sigma.d=runif(1), sigma.e=runif(1), sigma.f=runif(1), sigma.g=runif(1))
}

prod2.parameters <- c ("a", "d", "e", "f", "g", "sigma.g", "sigma.f", "sigma.a", "sigma.d", "sigma.e", "sigma.y")

B.2 <- bugs (prod2.data, prod2.inits, prod2.parameters, "india.2.bug", n.chains=4, 
             n.iter=10000, working.directory=NULL, clearWD=TRUE, 
             debug=TRUE )

plot(B.2, display.parallel = TRUE)
#options(max.print=1000000)
print(B.2, digits.summary = 3)


##
## Bugs model with covariates, village, season, village-season, household, and parcel
##
prod3.data <-  list ("y", "x", "vil_t", "tindex", "vil_id", "VDSA_HH_id", "parcel_id", "n", "L","C", "V", "H", "K", "J")
prod3.inits <- function(){
  list (b=rnorm(1), mu.a=rnorm(1), a=rnorm(C), mu.d=rnorm(1), d=rnorm(V), mu.e=rnorm(1), e=rnorm(H),
        mu.f=rnorm(1), f=rnorm(K), mu.g=rnorm(1), g=rnorm(J), sigma.y=runif(1), sigma.a=runif(1), 
        sigma.d=runif(1), sigma.e=runif(1), sigma.f=runif(1), sigma.g=runif(1))
}

prod3.parameters <- c ("a", "d", "e", "f", "g", "b", "sigma.g", "sigma.f", "sigma.a", "sigma.d", "sigma.e", "sigma.y")

B.3 <- bugs (prod3.data, prod3.inits, prod3.parameters, "india.3.bug", n.chains=4, 
               n.iter=10000, working.directory=NULL, clearWD=TRUE, 
               debug=TRUE )

plot(B.3)
options(max.print=1000000)
print(B.3, digits.summary = 3)

##
## histograms
##
attach.all(B.3$sims.list)

sigma.y <- sigma.y^2
sigma.a <- sigma.a^2
sigma.e <- sigma.e^2
sigma.d <- sigma.d^2
sigma.f <- sigma.f^2
sigma.g <- sigma.g^2

hist(sigma.y, xlab=expression(sigma[y]^2), breaks=20, ylab="Frequency", 
    main=paste("Histogram of Idiosyncratic Error Term"), cex.main=2, col="burlywood1")
hist(sigma.a, xlab=expression(tau[theta]^2), breaks=20, ylab="Frequency", 
    main=paste("Histogram of Weather Error Term"), cex.main=2, col="burlywood1")
hist(sigma.e, xlab=expression(tau[omega]^2), breaks=20, ylab="Frequency", 
    main=paste("Histogram of Time Error Term"), cex.main=2, col="burlywood1")
hist(sigma.d, xlab=expression(tau[zeta]^2), breaks=20, ylab="Frequency", 
    main=paste("Histogram of Village Error Term"), cex.main=2, col="burlywood1")
hist(sigma.f, xlab=expression(tau[delta]^2), breaks=20, ylab="Frequency", 
    main=paste("Histogram of Household Error Term"), cex.main=2, col="burlywood1")
hist(sigma.g, xlab=expression(tau[gamma]^2), breaks=20, ylab="Frequency", 
    main=paste("Histogram of Parcel Error Term"), cex.main=2, col="burlywood1")

##
## mixture plots
##
A <-as.mcmc.list(B.3)

## idiosyncratic variance 
x <- as.vector(time(A))

s.y <- which(varnames(A)=="sigma.y")
s.yy <- A[ , s.y, drop=TRUE]
s.yy <- do.call("cbind", s.yy)
s.yy <- s.yy^2

s.a <- which(varnames(A)=="sigma.a")
s.ay <- A[ , s.a, drop=TRUE]
s.ay <- do.call("cbind", s.ay)
s.ay <- s.ay^2

s.e <- which(varnames(A)=="sigma.e")
s.ey <- A[ , s.e, drop=TRUE]
s.ey <- do.call("cbind", s.ey)
s.ey <- s.ey^2

s.d <- which(varnames(A)=="sigma.d")
s.dy <- A[ , s.d, drop=TRUE]
s.dy <- do.call("cbind", s.dy)
s.dy <- s.dy^2

s.f <- which(varnames(A)=="sigma.f")
s.fy <- A[ , s.f, drop=TRUE]
s.fy <- do.call("cbind", s.fy)
s.fy <- s.fy^2

s.g <- which(varnames(A)=="sigma.g")
s.gy <- A[ , s.g, drop=TRUE]
s.gy <- do.call("cbind", s.gy)
s.gy <- s.gy^2
  
matplot(x,s.yy, type="l", col= c("darksalmon", "darkviolet", "cyan", "black"), xlab = "Iterations", 
        ylab = expression(sigma[y]^2), main = expression("Trace of "*sigma[y]^2))
matplot(x,s.ay, type="l", col= c("darksalmon", "darkviolet", "cyan", "black"), xlab = "Iterations", 
        ylab = expression(tau[theta]^2), main = expression("Trace of "*tau[theta]^2))
matplot(x,s.ey, type="l", col= c("darksalmon", "darkviolet", "cyan", "black"), xlab = "Iterations", 
        ylab = expression(tau[omega]^2), main = expression("Trace of "*tau[omega]^2))
matplot(x,s.dy, type="l", col= c("darksalmon", "darkviolet", "cyan", "black"), xlab = "Iterations", 
        ylab = expression(tau[xi]^2), main = expression("Trace of "*tau[xi]^2))
matplot(x,s.fy, type="l", col= c("darksalmon", "darkviolet", "cyan", "black"), xlab = "Iterations", 
        ylab = expression(tau[delta]^2), main = expression("Trace of "*tau[delta]^2))
matplot(x,s.gy, type="l", col= c("darksalmon", "darkviolet", "cyan", "black"), xlab = "Iterations", 
        ylab = expression(tau[gamma]^2), main = expression("Trace of "*tau[gamma]^2))



###################
#Robustness Checks
###################

## Season only
M4 <- lmer(lny ~ x + factor(crop) -1 + (1 | parcel_id) + (1 | VDSA_HH_id) + (1 | tindex) )
display(M4, digits=3, detail=TRUE)

## Village and season only
M5 <- lmer(lny ~ x + factor(crop) -1 + (1 | parcel_id) + (1 | VDSA_HH_id) + (1 | vil_id) )
display(M5, digits=3, detail=TRUE)

## Village-season only
M6 <- lmer(lny ~ x + factor(crop) -1 + (1 | parcel_id) + (1 | VDSA_HH_id) + (1 | vil_t) )
display(M6, digits=3, detail=TRUE)

## Village-season and season only
M7 <- lmer(lny ~ x + factor(crop) -1 + (1 | parcel_id) + (1 | VDSA_HH_id) + (1 | tindex) + (1 | vil_t) )
display(M7, digits=3, detail=TRUE)

## Village-season and village only
M8 <- lmer(lny ~ x + factor(crop) -1 + (1 | parcel_id) + (1 | VDSA_HH_id) + (1 | vil_id) + (1 | vil_t) )
display(M8, digits=3, detail=TRUE)

stargazer(M4, M5, M6, M7, M8, M3, title="Results of Robustness Checks", covariate.labels=c("log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "log labor", "log fertilizer", "log mechanization", "log pesticides", "countryIndia"), align=TRUE, no.space=TRUE)

