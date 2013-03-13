require(lavaan)

##############################################################################
##
## Table 3.2
## (nothing lavaan here)
##
##############################################################################

# sample covariance matrix
S <- matrix(c(4.0, 2.4, 2.8,
		2.4, 4.0, 2.0,
		2.8, 2.0, 4.0), nrow=3)

# predicted covariance matrix
SIGM <- matrix(c(4.0, 2.4, 1.2,
		2.4, 4.0, 2.0,
		1.2, 2.0, 4.0), nrow=3)

RES <- S - SIGM             # RESIDUAL COVARIANCE MATRIX
SDET  <- det(S)             # DETERMINANT OF THE SAMPLE COVARIANCE MATRIX
SIGMDET  <- det(SIGM)       # DETERMINANT OF THE PREDICTED COVARIANCE MATRIX
LOGS  <- log(SDET)          # NATURAL LOG OF SAMPLE MATRIX DETERMINANT
LOGSIGM  <- log(SIGMDET)    # NATURAL LOG OF PREDICTED MATRIX DETERMINANT
SIGMINV = solve(SIGM)       # INVERSE OF PREDICTED MATRIX
SDIV = S %*% SIGMINV        # MULTIPLICATION OF SAMPLE MATRIX AND PREDICTED INVERSE
STRACE = sum(diag(SDIV))    # TRACE OF THE RESULTING SDIV MATRIX
SORDER = dim(S)[1]          # ORDER OF SAMPLE MATRIX = NUMBER OF INDICATORS

#  calculation of Fml
FML = abs((LOGS - LOGSIGM) + STRACE - SORDER);


# print results
S
SIGM
RES
SDET
SIGMDET
LOGS 
LOGSIGM
SDIV
STRACE
SORDER
FML


##############################################################################
##
## Table 4.1 - 4.5
## Two-factor cfa model of neuroticism and extraversion
##
##############################################################################
sds <- '5.7  5.6  6.4  5.7  6.0  6.2  5.7  5.6'

cors <- '
 1.000
 0.767  1.000 
 0.731  0.709  1.000 
 0.778  0.738  0.762  1.000 
-0.351  -0.302  -0.356  -0.318  1.000 
-0.316  -0.280  -0.300  -0.267  0.675  1.000 
-0.296  -0.289  -0.297  -0.296  0.634  0.651  1.000 
-0.282  -0.254  -0.292  -0.245  0.534  0.593  0.566  1.000'

covs <- getCov(cors, sds = sds, names = c("n1", "n2", "n3", "n4", "e1", "e2", "e3", "e4"))

model <- '

  neuroticism  =~ n1 + n2 + n3 + n4
  extraversion =~ e1 + e2 + e3 + e4

'

fit <- cfa(model, sample.cov = covs, sample.nobs = 250, mimic = "mplus")


##############################################################################
##
## Table 4.2 
## (continued from Table 4.1)
##
##############################################################################
# Sample Variance-Covariance Matrix (S)
covs

# Model-Implied Variance-Covariance Matrix (Sigma)
fitted(fit)$cov

# Fitted Residual Matrix
resid(fit)$cov

# Standardized Residual Matrix 
# Brown reports LISREL 8.72 output which is different from the lavaan result
resid(fit, type = "standardized")$cov


##############################################################################
##
## Table 4.3
## (continued from Table 4.1)
## Modification indices and EPC's for the two-factor cfa model of neuroticism and extraversion
##
##############################################################################
modindices(fit)

# mi: modification indices
# sepc.all: completely standardized expected change

##############################################################################
##
## Table 4.4
## (continued from Table 4.1)
## standardized parameter estimates and r-square
##
##############################################################################
summary(fit, standardized = TRUE, rsquare = TRUE)

##############################################################################
##
## Table 4.5
## Measurement model of health status involving latent variables and single indicators
##
##############################################################################
Data <- read.table("http://people.bu.edu/tabrown/Ch4/fig4.3.dat")
names(Data) <- c("subjid", "activ", "soma", "pain", "menth", "socf", "vital", "genhlth", "age")

model <- '
  
  # factors
  physicalf =~ activ + soma + pain
  mentalf   =~ menth + socf + vital
  
  # pseudofactors for single indicators
  gwb       =~ genhlth
  agef      =~ age
  
  # fix residual variance of pseudofactor indicators
  age       ~~ 0 * age
  genhlth   ~~ 7.88 * genhlth
  
  # residual covariances
  activ ~~ soma

'

fit <- cfa(model, data = Data)
summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

##############################################################################
##
## Table 5.1
## Standardized residuals and modification indices for a one-factor cfa model
## of indicators of neuroticism and extraversion
##
##############################################################################

sds <- '5.7  5.6  6.4  5.7  6.0  6.2  5.7  5.6'

cors <- '
 1.000
 0.767  1.000 
 0.731  0.709  1.000 
 0.778  0.738  0.762  1.000 
-0.351  -0.302  -0.356  -0.318  1.000 
-0.316  -0.280  -0.300  -0.267  0.675  1.000 
-0.296  -0.289  -0.297  -0.296  0.634  0.651  1.000 
-0.282  -0.254  -0.292  -0.245  0.534  0.593  0.566  1.000'

covs <- getCov(cors, sds = sds, names = c("n1", "n2", "n3", "n4", "e1", "e2", "e3", "e4"))

model <- 'onefactor  =~ n1 + n2 + n3 + n4 + e1 + e2 + e3 + e4'

fit <- cfa(model, sample.cov = covs, sample.nobs = 250, mimic = "eqs")
summary(fit, fit.measures = TRUE)
## RMSEA, CFI, TLI differ

# Standardized Residual Matrix
# Brown reports LISREL 8.72 output which is different from the lavaan result
resid(fit, type = "standardized")$cov


# modification indices
modindices(fit)


##############################################################################
##
## Table 5.2
## Base model
##
##############################################################################
sds <- '2.06  1.52  1.92  1.41  1.73  1.77  2.49  2.27  2.68  1.75  2.57  2.66'

cors <- '
  1.000 
  0.300  1.000 
  0.229  0.261  1.000 
  0.411  0.406  0.429  1.000 
  0.172  0.252  0.218  0.481  1.000 
  0.214  0.268  0.267  0.579  0.484  1.000 
  0.200  0.214  0.241  0.543  0.426  0.492  1.000 
  0.185  0.230  0.185  0.545  0.463  0.548  0.522  1.000 
  0.134  0.146  0.108  0.186  0.122  0.131  0.108  0.151  1.000 
  0.134  0.099  0.061  0.223  0.133  0.188  0.105  0.170  0.448  1.000 
  0.160  0.131  0.158  0.161  0.044  0.124  0.066  0.061  0.370  0.350  1.000 
  0.087  0.088  0.101  0.198  0.077  0.177  0.128  0.112  0.356  0.359  0.507  1.000'

covs <- getCov(cors, sds = sds, names = paste("x", 1:12, sep = ""))

model <- '

  copingm  =~ x1 + x2 + x3 + x4
  socialm  =~ x4 + x5 + x6 + x7 + x8
  enhancem =~ x9 + x10 + x11 + x12

  x11 ~~ x12 

'

fit <- cfa(model, sample.cov = covs, sample.nobs = 500, mimic = "mplus")
summary(fit, fit.measures = TRUE)

##############################################################################
##
## Table 5.3
## no cross-loading x4
##
##############################################################################
model <- '

  copingm  =~ x1 + x2 + x3 + x4
  socialm  =~ x5 + x6 + x7 + x8
  enhancem =~ x9 + x10 + x11 + x12

  x11 ~~ x12 '

fit <- cfa(model, sample.cov = covs, sample.nobs = 500, mimic = "mplus")

# Standardized Residuals
# Brown reports LISREL 8.72 output which is different from the lavaan result
resid(fit, type = "standardized")$cov


# Modification Indices
modindices(fit)

# Completely Standardized Solution
standardizedSolution(fit)


##############################################################################
##
## Table 5.4
## x12 on wrong factor
##
##############################################################################
model <- '

  copingm  =~ x1 + x2 + x3 + x4
  socialm  =~ x4 +x5 + x6 + x7 + x8 + x12
  enhancem =~ x9 + x10 + x11

  x11 ~~ x12 '

fit <- cfa(model, sample.cov = covs, sample.nobs = 500, mimic = "mplus")
## WARNING: some fit measures are a bit different

# Standardized Residuals
# Brown reports LISREL 8.72 output which is different from the lavaan result
resid(fit, type = "standardized")$cov

# Modification Indices
modindices(fit)

# Completely Standardized Solution
standardizedSolution(fit)


##############################################################################
##
## Table 5.5
## without correlated error of item 11 and item 12
##
##############################################################################
model <- '

  coping  =~ x1 + x2 + x3 + x4
  social  =~ x4 +x5 + x6 + x7 + x8
  enhance =~ x9 + x10 + x11 + x12
'

fit <- cfa(model, sample.cov = covs, sample.nobs = 500, mimic = "mplus")

# Standardized Residuals
# Brown reports LISREL 8.72 output which is different from the lavaan result
resid(fit, type = "standardized")$cov

# Modification Indices
modindices(fit)

# Completely Standardized Solution
standardizedSolution(fit)


##############################################################################
##
## Table 5.6
## efa
## not in lavaan but e.g. with the factanal package
##
##############################################################################
Data <- read.table("http://people.bu.edu/tabrown/Ch5/efa.dat")
Data <- Data[,-13]
names(Data) <- paste("x", 1:12, sep = "")

fit <- factanal(Data, 3, rotation = "promax")
print(fit)

# eigenvalues
eigen(cor(Data))$values


##############################################################################
##
## Table 5.7 + 5.8
## Drinking motives e/cfa
##
##############################################################################
Data <- read.table("http://people.bu.edu/tabrown/Ch5/efa.dat")
names(Data) <- paste("x", 1:13, sep = "") # x13 is not used in the model


model <- '
  coping  =~ NA*x1 + x2 + x3 + x4 + x5 + x6 + x7 + 0*x8 + x9 + x10 + x11 + 0*x12
  social  =~ NA*x8 + 0*x1 + x2 + x3 + x4 + x5 + x6 + x7 + x9 + x10 + x11 + 0*x12
  enhance =~ NA*x12 + 0*x1 + x2 + x3 + x4 + x5 + x6 + x7 + 0*x8 + x9 + x10 + x11
  
  coping  ~~ 1*coping
  social  ~~ 1*social
  enhance ~~ 1*enhance
'

fit <- cfa(model, data = Data)
summary(fit, fit.measures = TRUE)

mis <- modindices(fit)
mis[mis$mi >= 10 & !is.na(mis$mi),]

##############################################################################
##
## Table 6.2
## Correlated methods CFA
##
##############################################################################

# no results provided in book

##############################################################################
##
## Table 6.3 + 6.4
## Correlated uniqueness CFA of the MTMM matrix of cluster A personality disorders
##
##############################################################################
sds <- '3.61  3.66  3.59  2.94  3.03  2.85  2.22  2.42  2.04'

cors <- '
  1.000 
  0.290  1.000 
  0.372  0.478  1.000 
  0.587  0.238  0.209  1.000 
  0.201  0.586  0.126  0.213  1.000 
  0.218  0.281  0.681  0.195  0.096  1.000 
  0.557  0.228  0.195  0.664  0.242  0.232  1.000 
  0.196  0.644  0.146  0.261  0.641  0.248  0.383  1.000 
  0.219  0.241  0.676  0.290  0.168  0.749  0.361  0.342  1.000'

covs <- getCov(cors, sds = sds, names = c("pari", "szti", "szdi", "parc", "sztc", "szdc", "paro", "szto", "szdo"))

model <- '
  paranoid    	=~ pari + parc + paro
  schizotypal 	=~ szti + sztc + szto
  schizoid    	=~ szdi + szdc + szdo

  pari ~~ szti + szdi
  szti ~~ szdi
  parc ~~ sztc + szdc
  sztc ~~ szdc
  paro ~~ szto + szdo
  szto ~~ szdo
'  

fit <- cfa(model, sample.cov = covs, sample.nobs = 500, std.lv = TRUE)
summary(fit, fit.measures = TRUE, standardized = TRUE)


##############################################################################
##
## Table 7.1
## Two-factor model of memory - congeneric
##
##############################################################################
sds <- '2.610  2.660  2.590  1.940  2.030  2.050'

cors <-'
  1.000
  0.661  1.000
  0.630  0.643  1.000
  0.270  0.300  0.268  1.000
  0.297  0.265  0.225  0.805  1.000
  0.290  0.287  0.248  0.796  0.779  1.000'

covs <- getCov(cors, sds = sds, names = paste("x", 1:6, sep = ""))

model.congeneric <- '

  auditorymemory =~ x1 + x2 + x3
  visualmemory   =~ x4 + x5 + x6

'

fit.congeneric <- cfa(model.congeneric, sample.cov = covs, sample.nobs = 200, std.lv = TRUE)
summary(fit.congeneric, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)


##############################################################################
##
## Table 7.3
## Congeneric, Tau-equivalent and parallel solutions of a
## two-factor model  of memory (N=200)
##
##############################################################################

# tau equivalent: auditory memory
model.tau.a <- '
  auditorymemory =~ x1 + v1*x1 + v1*x2 + v1*x3
  visualmemory   =~ x4 + x5 + x6
'

fit.tau.a <- cfa(model.tau.a, sample.cov = covs, sample.nobs = 200, std.lv = TRUE)

# tau equivalent: auditory & visual memory
model.tau.av <- '
  auditorymemory =~ x1 + v1*x1 + v1*x2 + v1*x3
  visualmemory   =~ x4 + v2*x4 + v2*x5 + v2*x6
'

fit.tau.av <- cfa(model.tau.av, sample.cov = covs, sample.nobs = 200, std.lv = TRUE)

# parallel: auditory memory
model.parallel.a <- '
  auditorymemory =~ x1 + v1*x1 + v1*x2 + v1*x3
  visualmemory   =~ x4 + v2*x4 + v2*x5 + v2*x6

  x1 ~~ v3 * x1
  x2 ~~ v3 * x2
  x3 ~~ v3 * x3
'

fit.parallel.a <- cfa(model.parallel.a, sample.cov = covs, sample.nobs = 200, std.lv = TRUE)

# parallel: auditory & visual memory
model.parallel.av <- '
  auditorymemory =~ x1 + v1*x1 + v1*x2 + v1*x3
  visualmemory   =~ x4 + v2*x4 + v2*x5 + v2*x6

  x1 ~~ v3 * x1
  x2 ~~ v3 * x2
  x3 ~~ v3 * x3
 
  x4 ~~ v4 * x4
  x5 ~~ v4 * x5
  x6 ~~ v4 * x6
'
fit.parallel.av <- cfa(model.parallel.av, sample.cov = covs, sample.nobs = 200, std.lv = TRUE)

# Chi square difference tests
anova(fit.congeneric, fit.tau.a, fit.tau.av, fit.parallel.a, fit.parallel.av, test = "chisq")


##############################################################################
##
## Table 7.4
## Two-factor model of memory - parallel
## (continued from Table 7.3)
##
##############################################################################
summary(fit.parallel.av, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

##############################################################################
##
## Table 7.6 + 7.7
## Longitudinal model of job satisfaction
##
##############################################################################
sds <- '1.940  2.030  2.050  1.990  2.610  2.660  2.590  2.550'

cors <- '
  1.000
  0.736  1.000
  0.731  0.648  1.000
  0.771  0.694  0.700  1.000
  0.685  0.512  0.496  0.508  1.000
  0.481  0.638  0.431  0.449  0.726  1.000
  0.485  0.442  0.635  0.456  0.743  0.672  1.000
  0.508  0.469  0.453  0.627  0.759  0.689  0.695  1.000'

covs <- getCov(cors, sds = sds, names = c("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2"))

ms <- c(1.500,  1.320,  1.450,  1.410,  6.600,  6.420,  6.560,  6.310) # this should be a numeric vector, not a character string

model.equalform <- '

    satis1 =~ A1 + B1 + C1 + D1
    satis2 =~ A2 + B2 + C2 + D2

    A1 ~~ A2
    B1 ~~ B2
    C1 ~~ C2
    D1 ~~ D2
    
    # fix indicator intercepts to 0
    A1 ~ 0*1
    A2 ~ 0*1
    
    # free factor intercepts
    satis1 ~ 1
    satis2 ~ 1
'

fit.equalforms <- cfa(model.equalform, sample.cov = covs, sample.nobs = 250, sample.mean = ms, meanstructure = TRUE)
summary(fit.equalforms, standardized = TRUE, fit.measures = TRUE)

##############################################################################
##
## Table 7.7
## Longitudinal model of job satisfaction
## Equality of factor loadings over the two assessment points
## ..continues from previous table
##
##############################################################################

model.equalfl <- '

    # equality of factor loadings
    satis1 =~ v1*A1 + v2*B1 + v3*C1 + v4*D1
    satis2 =~ v1*A2 + v2*B2 + v3*C2 + v4*D2

    A1 ~~ A2
    B1 ~~ B2
    C1 ~~ C2
    D1 ~~ D2
    
    # fix indicator intercepts to 0
    A1 ~ 0*1
    A2 ~ 0*1
    
    # free factor intercepts
    satis1 ~ 1
    satis2 ~ 1
'

fit.equalfl <- cfa(model.equalfl, sample.cov = covs, sample.nobs = 250, sample.mean = ms, meanstructure = TRUE)
summary(fit.equalfl, standardized = TRUE, fit.measures = TRUE)
anova(fit.equalforms, fit.equalfl, test = "chisq")

##############################################################################
##
## Table 7.7
## Longitudinal model of job satisfaction
## Equality of indicator intercepts over the two assessment points
## ..continued from previous table
##
##############################################################################
model.equali <- '

    # equality of factor loadings
    satis1 =~ v1*A1 + v2*B1 + v3*C1 + v4*D1
    satis2 =~ v1*A2 + v2*B2 + v3*C2 + v4*D2

    A1 ~~ A2
    B1 ~~ B2
    C1 ~~ C2
    D1 ~~ D2
    
    # fix indicator intercepts to 0
    A1 ~ 0*1
    A2 ~ 0*1

    # free factor intercepts
    satis1 ~ 1
    satis2 ~ 1

    # equal indicator intercepts
    B1 ~ equal("B2~1")*1
    C1 ~ equal("C2~1")*1
    D1 ~ equal("D2~1")*1
'

fit.equali <- cfa(model.equali, sample.cov = covs, sample.nobs = 250, sample.mean = ms, meanstructure = TRUE)
summary(fit.equali, standardized = TRUE, fit.measures = TRUE)
anova(fit.equalfl, fit.equali,  test = "chisq")

##############################################################################
##
## Table 7.7
## Longitudinal model of job satisfaction
## Equality of indicator error variances over the two assessment points
## ..continued from previous table
##
##############################################################################
model.equalrv <- '

    # equality of factor loadings
    satis1 =~ v1*A1 + v2*B1 + v3*C1 + v4*D1
    satis2 =~ v1*A2 + v2*B2 + v3*C2 + v4*D2

    A1 ~~ A2
    B1 ~~ B2
    C1 ~~ C2
    D1 ~~ D2
    
    # fix indicator intercepts to 0
    A1 ~ 0*1
    A2 ~ 0*1

    # free factor intercepts
    satis1 ~ 1
    satis2 ~ 1

    # equal indicator intercepts
    B1 ~ equal("B2~1")*1
    C1 ~ equal("C2~1")*1
    D1 ~ equal("D2~1")*1

    # equal residual variances
    A1 ~~ v5*A1
    B1 ~~ v6*B1
    C1 ~~ v7*C1
    D1 ~~ v8*D1
    A2 ~~ v5*A2
    B2 ~~ v6*B2
    C2 ~~ v7*C2
    D2 ~~ v8*D2
'

fit.equalrv <- cfa(model.equalrv, sample.cov = covs, sample.nobs = 250, sample.mean = ms, meanstructure = TRUE)
summary(fit.equalrv, standardized = TRUE, fit.measures = TRUE)
anova(fit.equali, fit.equalrv,  test = "chisq")


##############################################################################
##
## Table 7.9
## 
## Tests of measurement invariance and population heterogeneity of DSM-IV
## major depressive disorder in men and women
##
##############################################################################
Data <- read.table("http://people.bu.edu/tabrown/Ch7/MDDALL.dat")
names(Data) <- c("sex", paste("mdd", 1:9, sep = ""))
Data$sex <- factor(Data$sex, levels = c(0, 1), labels = c("female", "male"))

model.mdd <- '
  MDD =~ mdd1 + mdd2 + mdd3 + mdd4 + mdd5 + mdd6 + mdd7 + mdd8 + mdd9
  mdd1 ~~ mdd2
'
# Single group solution (men)
fit.men <- cfa(model.mdd, data = Data[Data$sex == "male",])
summary(fit.men, fit.measures = TRUE)

# Single group solution (women)
fit.women <- cfa(model.mdd, data = Data[Data$sex == "female",])
summary(fit.women, fit.measures = TRUE)

# Measurement Invariance using semTools
require(semTools)
measurementInvariance(model.mdd, data = Data, group = "sex", strict = TRUE)

# using lavaan
# measurementInvariance doesn't do equal factor variance. But, this can be accomplished as follows
fit.ef <- cfa(model.mdd, data = Data, group = "sex", meanstructure = TRUE) # equal form
fit.efl <- update(fit.ef, group.equal = c("loadings")) # equal factor laodings
fit.eii <- update(fit.efl, group.equal = c("loadings", "intercepts")) # equal indicator intercepts
fit.eir <- update(fit.eii, group.equal = c("loadings", "intercepts", "residuals")) # equal indicator error variances
fit.fv <- update(fit.eir, group.equal = c("loadings", "intercepts", "residuals", "lv.variances")) # equal factor variances
fit.fm <- update(fit.fv, group.equal = c("loadings", "intercepts", "residuals", "lv.variances", "means")) # equal latent means

# chi-squared diff tests
anova(fit.ef, fit.efl, fit.eii, fit.eir, fit.fv, fit.fm, test = "chisq")

##############################################################################
##
## Table 7.11
## 
## Parameter estimates from the equal form measurement model of major 
## depression in men and woman
## .. continued from previous chunk 
##############################################################################
summary(fit.ef, standardized = TRUE, rsquare = TRUE)


##############################################################################
##
## Table 7.14-7.15
## 
## MIMIC model of Social Phobia and Agoraphobia 
##
##############################################################################
sds <- '2.26 2.73 2.11 2.32 2.61 2.44 0.50'

cors <- '
  1.000
  0.705   1.000
  0.724   0.646   1.000
  0.213   0.195   0.190   1.000
  0.149   0.142   0.128   0.521   1.000
  0.155   0.162   0.135   0.557   0.479   1.000
  -0.019  -0.024  -0.029  -0.110  -0.074  -0.291   1.000'

covs <- getCov(cors, sds = sds, names = c("S1", "S2", "S3", "A1", "A2", "A3", "sex"))

model <- '
  socialphobia =~ S1 + S2 + S3 
  agoraphobia  =~ A1 + A2 + A3

  socialphobia + agoraphobia + A3 ~ sex

  socialphobia ~~ agoraphobia
'

fit <- cfa(model, sample.cov = covs, sample.nobs = 730)
summary(fit, standardized = TRUE, fit.measures = TRUE)

##############################################################################
##
## Table 8.1
## Four-factor (first-order) measurement model of coping
##
## Table 8.3
## Higher-order model of coping
##
##############################################################################
sds <- '1.40  2.10  1.30  2.30  2.40  1.90  2.00  2.90  2.30  3.10  2.20  1.20'

cors <- '
  1.00
  0.78  1.00
  0.80  0.77  1.00
  0.56  0.51  0.48  1.00
  0.52  0.51  0.46  0.78  1.00
  0.59  0.51  0.51  0.80  0.79  1.00
  0.16  0.15  0.17  0.14  0.18  0.16  1.00
  0.19  0.13  0.18  0.14  0.16  0.16  0.81  1.00
  0.12  0.17  0.17  0.17  0.20  0.16  0.75  0.80  1.00
  0.16  0.13  0.17  0.15  0.16  0.18  0.56  0.52  0.50  1.00
  0.16  0.14  0.18  0.15  0.16  0.18  0.51  0.58  0.51  0.81  1.00
  0.16  0.15  0.14  0.16  0.16  0.14  0.52  0.57  0.52  0.80  0.79  1.00'

covs <- getCov(cors, sds = sds, names = c("P1", "P2", "P3", "C1", "C2", "C3", "E1", "E2", "E3", "S1", "S2", "S3"))

model8.1 <- '
  probslv         =~ P1 + P2 + P3
  cogres          =~ C1 + C2 + C3
  expremot        =~ E1 + E2 + E3
  socspt          =~ S1 + S2 + S3
'

fit8.1 <- cfa(model8.1, sample.cov = covs, sample.nobs = 275)
inspect(fit8.1, "std.coef")

model8.3 <- '
  probslv         =~ P1 + P2 + P3
  cogres          =~ C1 + C2 + C3
  expremot        =~ E1 + E2 + E3
  socspt          =~ S1 + S2 + S3
  
  probfoc           =~ probslv + cogres
  emotfoc           =~ expremot + socspt
'

fit8.3 <- cfa(model8.3, sample.cov = covs, sample.nobs = 275, std.lv = TRUE)
inspect(fit8.3, "std.coef")

##############################################################################
##
## Table 8.5
##
## Selected unstandardized parameter estimates of measurement model of the 
## reactions to traumatic stress questionnaire
## 
##############################################################################
sds <- '1.150  1.200  1.570  2.820  1.310  1.240  1.330  1.290'

cors <- '
  1.000
  0.594  1.000
  0.607  0.613  1.000
  0.736  0.765  0.717  1.000
  0.378  0.321  0.360  0.414  1.000
  0.314  0.301  0.345  0.363  0.732  1.000
  0.310  0.262  0.323  0.337  0.665  0.583  1.000
  0.317  0.235  0.276  0.302  0.632  0.557  0.796  1.000'

covs <- getCov(cors, sds = sds, names = paste("Y", 1:8, sep = ""))

model <- '
  intrus =~ Y1 + Y2 + Y3 + Y4
  avoid  =~ Y5 + Y6 + Y7 + Y8

  Y7 ~~ Y8
'
fit <- cfa(model, sample.cov = covs, sample.nobs = 500, std.lv = TRUE, mimic = "EQS")
inspect(fit, "coef")

##############################################################################
##
## Table 8.6, 8.7 & 8.8
## 
## Syntax for computing point estimates for scale reliability and 95% CI
## 
##############################################################################
sds <- '1.15 1.20 1.57 2.82 1.31 1.24 1.33 1.29'

cors <- '
  1.000
  0.594  1.000
  0.607  0.613  1.000
  0.736  0.765  0.717  1.000
  0.378  0.321  0.360  0.414  1.000
  0.314  0.301  0.345  0.363  0.732  1.000
  0.310  0.262  0.323  0.337  0.665  0.583  1.000
  0.317  0.235  0.276  0.302  0.632  0.557  0.796  1.000'

covs <- getCov(cors, sds = sds, names = paste("Y", 1:8, sep = ""))


model <- '
  # main model.
  intrus =~ Y1 + l1*Y1 + l2*Y2 + l3*Y3 + l4*Y4
  avoid  =~ Y5 + l5*Y5 + l6*Y6 + l7*Y7 + l8*Y8

  # label the residual variances
  Y1 ~~ e1*Y1
  Y2 ~~ e2*Y2
  Y3 ~~ e3*Y3
  Y4 ~~ e4*Y4
  Y5 ~~ e5*Y5
  Y6 ~~ e6*Y6
  Y7 ~~ e7*Y7
  Y8 ~~ e8*Y8

  # covariance between Y7 and Y8
  Y7 ~~ e78*Y8

  # defined parameters
  intr.tru := (l1 + l2 + l3 + l4)^2
  intr.tot := (l1 + l2 + l3 + l4)^2 + e1 + e2 + e3 + e4
  intr.rel := intr.tru/intr.tot  

  avoid.tru := (l5 + l6 + l7 + l8)^2
  avoid.tot := (l5 + l6 + l7 + l8)^2 + e5 + e6 + e7 + e8 + 2*e78
  avoid.rel := avoid.tru/avoid.tot

'

fit <- cfa(model, sample.cov = covs, sample.nobs = 500, mimic = "EQS", std.lv = TRUE)
summary(fit, standardized = TRUE, fit.measures = TRUE)

# 95% CI
parameterEstimates(fit)


##############################################################################
##
## Table 8.9
## Model with a formative construct
## 
##############################################################################
sds <- '2.5  2.1  3.0  4.1  3.9  4.4  1.2  1.0  1.2'

cors <- '
  1.000
  0.700  1.000
  0.713  0.636  1.000
  0.079  0.066  0.076  1.000
  0.088  0.058  0.070  0.681  1.000
  0.084  0.056  0.074  0.712  0.633  1.000
  0.279  0.248  0.240  0.177  0.155  0.170  1.000
  0.250  0.214  0.222  0.157  0.143  0.152  0.373  1.000
  0.280  0.236  0.251  0.173  0.178  0.171  0.448  0.344  1.000 '

covs <- getCov(cors, sds = sds, names = c(paste("Y", 1:6, sep = ""), paste("X", 1:3, sep = "")))

model <- '
  satis =~ Y1 + Y2 + Y3
  optim =~ Y4 + Y5 + Y6
  stress =~ NA*satis + optim

  stress ~ 1*X1 + X2 + X3

  X1 ~~ X2 + X3
  X2 ~~ X3
'

fit <- cfa(model, sample.cov = covs, sample.nobs = 500)
summary(fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

##############################################################################
##
## Table 9.1
## Estimation of a CFA model with missing data using direct ML
## 
##############################################################################
Data <- read.table("http://people.bu.edu/tabrown/Ch9/cfamiss.dat", na.strings = "9")
names(Data) <- c("subject", "s1", "s2", "s3", "s4")

model <- '
  esteem =~ s1 + s2 + s3 + s4
  s2 ~~ s4
'

fit <- cfa(model, data = Data, missing = "ml")
summary(fit, fit.measures = TRUE)

## Question: could the following be combined into an "inspect function"?
# number of patterns
fit@Data@Mp[[1]]$npatterns

# patterns and their frequency
pats <- fit@Data@Mp[[1]]$pat * 1L
colnames(pats) <- fit@Data@ov.names[[1]]
print(pats)

# covariance coverage
coverage <- fit@Data@Mp[[1]]$coverage
colnames(coverage) <- rownames(coverage) <- fit@Data@ov.names[[1]]
print(coverage)


##############################################################################
##
## Table 9.2
## 
## Estimation of CFA model with missing data using multiple imputation
##
##############################################################################
data9.2 <- read.table("http://people.bu.edu/tabrown/Ch9/cfamiss.dat", na.strings = "9")
names(data9.2) <- c("subject", "s1", "s2", "s3", "s4")

model9.2 <- '
  esteem =~ s1 + s2 + s3 + s4
  s2 ~~ s4
'
require(semTools)
fit9.2 <- cfa.mi(model9.2, data = data9.2, m = 5, miPackage = "mice", seed = 44176)
inspect(fit9.2, "fit")
inspect(fit9.2, "impute")



##############################################################################
##
## Table 9.5 (syntax) + Table 9.7 (results)
## lavaan syntax for conducting CFA with non-normal, continuous data using
## robust maximum likelihood 
##
##############################################################################
Data <- read.table("http://people.bu.edu/tabrown/Ch9/NONML.DAT", nrows = 870)
names(Data) <- c("x1", "x2", "x3", "x4", "x5")

model <- '
  f1 =~ x1 + x2 + x3 + x4 + x5
  x1 ~~ x3 # added in the second run
'

fit <- cfa(model, data = Data, mimic = "EQS", estimator = "MLM")
summary(fit, fit.measures = TRUE)

##############################################################################
##
## Table 9.9 (syntax) + Table 9.10 (results)
## lavaan syntax for conducting CFA with categorical indicators
##
##############################################################################
Data <- read.fwf("http://people.bu.edu/tabrown/Ch9/BINARY.DAT", width = c(1,1,1,1,1,1), n = 750)
names(Data) <- c(paste("y", 1:6, sep = ""))

model1 <- '
  etoh =~ y1 + y2 + y3 + y4 + y5 + y6
'

fit1 <- cfa(model1, data = Data, ordered = names(Data), estimator = "WLSMVS", mimic = "mplus")
summary(fit1, fit.measures = TRUE)


##############################################################################
##
## Table 9.11
## lavaan syntax for conducting nested model comparison with WLSMV: One factor
## CFA of alcohol dependence with binary indicators (factor loadings constrained
## to equality)
##
##############################################################################

# continued from previous example

model2 <- '
  etoh =~ y1 + l1*y2 + l1*y3 + l1*y4 + l1*y5 + l1*y6
'

fit2 <- cfa(model2, data = Data, ordered = names(Data), estimator = "WLSMVS", mimic = "mplus")
summary(fit2, fit.measures = TRUE)

# diff test
anova(fit1, fit2)

# The output does not correspond exactly to what is reported by mplus.
# see: https://groups.google.com/forum/?fromgroups=#!topic/lavaan/LxqIagOPmRU

##############################################################################
##
## Table 9.14
## lavaan syntax for bootstrapping a one-factor CFA with non-normal, continuous 
## indicators
## 
##############################################################################
Data <- read.table("http://people.bu.edu/tabrown/Ch9/NONML.DAT", nrows = 870)
names(Data) <- c("x1", "x2", "x3", "x4", "x5")

model <- '
  f1 =~ x1 + x2 + x3 + x4 + x5
  x1 ~~ x3 
'

fit <- cfa(model, data = Data, se = "bootstrap", bootstrap = 500)
summary(fit, fit.measures = TRUE) 
parameterEstimates(fit, ci = TRUE, level = .90, boot.ci.type = "bca.simple")


##############################################################################
##
## Table 10.1
## Satorra-Saris method of determining power to detect that factor covariances
## is significantly different from zero
## 
##############################################################################

## Step 1: generate population covariance matrix from H1 model

covs <- matrix(c(
  1, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0,
  0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 1), nrow = 6)

colnames(covs) <- rownames(covs) <- paste("x", 1:6, sep = "")

model <- '
  esteem  =~ .65*x1 + .70*x2 + .72*x3
  depress =~ .60*x4 + .70*x5 + .65*x6
  
  esteem ~~ 1*esteem
  depress ~~ 1*depress
  
  x1 ~~ .5775*x1
  x2 ~~ .51*x2
  x3 ~~ .4816*x3
  x4 ~~ .64*x4
  x5 ~~ .51*x5
  x6 ~~ .5775*x6
  
  esteem ~~ .35*depress
'
fit <- cfa(model, sample.cov = covs, sample.nobs = 500)

covs.pop <- fitted(fit)$cov

## Step 2: analyze residual covariance matrix to ensure that population values are recovered

model <- '
  esteem  =~ x1 + x2 + x3
  depress =~ x4 + x5 + x6
'

fit <- cfa(model, sample.cov = covs.pop, sample.nobs = 500)
residuals(fit)$cov

## Step 3: fit H0 model that contains the misspecified parameter and target sample size

model <- '
  esteem  =~ x1 + x2 + x3
  depress =~ x4 + x5 + x6
  
  esteem ~~ 0*depress
'

fit <- cfa(model, sample.cov = covs.pop, sample.nobs = 100)


## Step 4: use X2 from step 3 as noncentrality parameter to estimate power at targeted sample sizes
dfs <- 1
alfa <- 0.05
crit <- qchisq(1-alfa,dfs)
lamba <- inspect(fit, "fit")[1]
power <- 1-pchisq(crit, 1, lamba)
power





