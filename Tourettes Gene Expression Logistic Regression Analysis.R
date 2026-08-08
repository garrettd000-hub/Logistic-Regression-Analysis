

#####################
#### Week 8 code ####
#####################

##########################################
####  Week 8 Notes Part 1 2-way ANOVA ####
##########################################
######################
#### Code Chunk 1 ####
######################
library(ALL); data(ALL)
fus.prot <- pData(ALL)$"fusion protein" # extract data…
boxplot(exprs(ALL)["1970_s_at",]~fus.prot,main="1970_s_at")
lm.1970.s.at.fp <- lm(exprs(ALL)["1970_s_at",]~fus.prot) # fit model
summary(lm.1970.s.at.fp)$coef
anova(lm.1970.s.at.fp) # get anova p-value, and we are done!
df.tmp <- data.frame(model.matrix(lm.1970.s.at.fp),
                     fus.prot[!is.na(fus.prot)])
colnames(df.tmp) <- c("Intercept","p190.p210","p210","FP")
df.tmp[1:10,]

######################
#### Code Chunk 2 ####
######################
# how many measurements (patients) do we have in each group defined
# by combinations of the two variables? : 
table(pData(ALL)[,c("sex","relapse")])
df.tmp <- data.frame(
  expr.158.at=exprs(ALL)["158_at",],
  sex=pData(ALL)$sex,
  relapse=pData(ALL)$relapse)
# run linear model with interaction term and apply anova() –
# this will give us two-way ANOVA, just that simple:
anova(lm(expr.158.at~sex*relapse,df.tmp))
boxplot(expr.158.at~sex+relapse,df.tmp) # stratify by sex/relapse
boxplot(expr.158.at~relapse,df.tmp) # stratify by relapse only

######################
#### Code Chunk 3 ####
######################
model.matrix(lm(expr.158.at~sex*relapse,df.tmp))[1:5,]
model.matrix(lm(expr.158.at~sex*relapse,df.tmp))[40:45,]
summary(lm(expr.158.at~sex*relapse,df.tmp))$coef
mean(df.tmp$expr.158.at[df.tmp$sex=="F" & df.tmp$relapse==F],na.rm=T)
mean(df.tmp$expr.158.at[df.tmp$sex=="M" & df.tmp$relapse==F],na.rm=T)

######################
#### Code Chunk 4 ####
######################
library(limma)
b.mask <- !is.na(pData(ALL)$sex)&!is.na(pData(ALL)$relapse)
sex.wo.na <- pData(ALL)$sex[b.mask]
relapse.wo.na <- pData(ALL)$relapse[b.mask]
exprs.wo.na <- exprs(ALL)[,b.mask]
design <- model.matrix(~sex.wo.na*relapse.wo.na)
design[1:5,]
limma.fit.e <- eBayes(lmFit(exprs.wo.na,design))
#note that we ask for most significant cross-terms below:
topTable(limma.fit.e,"sex.wo.naM:relapse.wo.naTRUE",5)
# we have looked at 158_at already, let’s check other hits from topTable:
boxplot(exprs(ALL)["35808_at",]~sex+relapse,df.tmp)
boxplot(exprs(ALL)["33700_at",]~sex+relapse,df.tmp)

######################
#### Code Chunk 5 ####
######################
df.tmp <- data.frame(
  expr.1803.at=exprs(ALL)["1803_at",],
  sex=pData(ALL)$sex,
  relapse=pData(ALL)$relapse)
anova(lm(expr.1803.at~sex*relapse,df.tmp))
boxplot(expr.1803.at~sex+relapse,df.tmp)

######################
#### Code Chunk 6 ####
######################
#### load data (shown in previous notes) ####
ALL.pdat <- pData(ALL)
date.cr.chr <- as.character(ALL.pdat$date.cr)
diag.chr <- as.character(ALL.pdat$diagnosis)
date.cr.t <- strptime(date.cr.chr,"%m/%d/%Y")
diag.t <- strptime(diag.chr,"%m/%d/%Y")
days2remiss <- as.numeric(date.cr.t - diag.t)
x.d2r <- as.numeric(days2remiss)
exprs.34852 <- exprs(ALL)["34852_g_at",]
d2r.34852 <- data.frame(G=exprs.34852,D2R=x.d2r)
################################################
df.41690.fp <- data.frame(D2R=x.d2r,
                          G=exprs(ALL)["41690_at",],
                          FP=as.character(fus.prot))
summary(lm(D2R~G+FP,df.41690.fp))$coef # fit without interaction
anova(lm(D2R~G+FP,df.41690.fp))
summary(lm(D2R~G*FP,df.41690.fp))$coef # fit with interaction
anova(lm(D2R~G*FP,df.41690.fp))

#####################################
####  Week 8 Notes Part 2 Extras ####
#####################################
######################
#### Code Chunk 1 ####
######################
set.seed(1234)
n.batch <- 7 # number of batches
n.reps <- 3 # n replicated measurements in each batch
sd.batch <- 1.0 # sd of random effect term (sigma_r) 
sd.rep <- 0.2 # sd of measurement error, epsilon
mu <- 1.0 # grand mean
# now we will generate 21 measurements: 3 per batch, 7 batches:
# sample 7 distinct values for alpha_i; replicate 3 times: 
batch.eff <- rep(rnorm(n.batch,sd=sd.batch),n.reps) # alpha-term 
rep.eff <- rnorm(n.reps*n.batch,sd=sd.rep) # sample measurement err
re.df <- data.frame(batch=rep(paste("b",1:n.batch),n.reps),
                    y=mu+batch.eff+rep.eff) # make sure everything in matching order
boxplot(y~batch,re.df)
library(lme4)
lmer(y~(1|batch),re.df)

######################
#### Code Chunk 2 ####
######################
library(ALL); data(ALL)
b.mask <- !is.na(pData(ALL)$relapse)
relapse.wo.na <- pData(ALL)$relapse[b.mask]
exprs.wo.na <- exprs(ALL)[,b.mask]
# save expr. levels into variable with a short name (pure convenience):
ge <- exprs.wo.na["1584_at",] # we pick a gene here, just an example
br <- seq(min(ge),max(ge),by=(max(ge)-min(ge))/10)
cnt.1 <- hist(ge[relapse.wo.na],breaks=br)$counts
cnt.0 <- hist(ge[!relapse.wo.na],breaks=br)$counts
glm.1584.at <- glm(relapse.wo.na~ge, family="binomial")
summary(glm.1584.at)
old.par<-par(ps=16)
plot(ge,relapse.wo.na, xlab="1584_at",ylab="Relapse")
points((br[1:10]+br[2:11])/2, 
       cnt.1/(cnt.0+cnt.1),pch=19,col='red',cex=1.2)
segments(br[1:10], cnt.1/(cnt.0+cnt.1),
         br[2:11], cnt.1/(cnt.0+cnt.1),col='red')
y.pred <- predict(glm.1584.at)
points(ge, exp(y.pred)/(1+exp(y.pred)),col="blue",pch='+')
par(old.par)

######################
#### Code Chunk 3 ####
######################
library(survival) # load survival analysis library
# copying old code here: extract days-to-remission:
ALL.pdat <- pData(ALL)
date.cr.chr <- as.character(ALL.pdat$date.cr)
diag.chr <- as.character(ALL.pdat$diagnosis)
date.cr.t <- strptime(date.cr.chr,"%m/%d/%Y")
diag.t <- strptime(diag.chr,"%m/%d/%Y")
d2r <- as.numeric(date.cr.t - diag.t)
# start the survival analysis:
d2r.ind <-as.numeric(!is.na(d2r)) # T at positions where d2r has data
d2r.surv <- Surv(d2r,d2r.ind)
plot(survfit(d2r.surv~sex,data=pData(ALL)),
     col=c(2,4),xlab="Days",
     ylab="1-p(remission)")
legend(100,0.8,legend=c('M','F'),lwd=1,col=c('blue','red'))

######################
#### Code Chunk 4 ####
######################
coxph(Surv(d2r,as.numeric(!is.na(d2r)))~sex, data=pData(ALL))





















