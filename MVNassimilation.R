# Generate and run Gibbs for bivariate normal linear regression model
# Note 1: Could generalize to p-variate normal
# Note 2: HERE, THE PRIOR FOR BETA IS CONDITIONAL ON SIGMA -- 4 x 1 beta=(beta1',beta2')' given 2 x 2 sigma ~ N_4(beta0, T0^{-1} %x% sigma) 
#   IN OTHER WORDS, BETA AND SIGMA ARE NOT A PRIORI INDEPENDENT
#   IN THE CASE OF NON-INFORMATIVE PRIORS, HOWEVER, THIS PRIOR SPECIFICATION IS SIMILAR TO ASSUMING PRIOR INDEPENDENCE
#   FOR MORE DETAILS, SEE FRUHWIRTH-SCHNATTER AND PYNE (2010) WEB SUPPLEMENT EQ (30)
# May, 2014
###########################################

###
# Plot leaf temp vs air temp 
###


rm(list = ls()) 
# Libraries
library(mvtnorm)
library(MCMCpack)   # iwish distribution used in updating Sigma
library(spam) #sparse matrices
library(RcppArmadillo)
library(Rcpp)
library(mvnfast)


setwd('C:/Users/mck14/Dropbox')
setwd('~/Dropbox')
sourceCpp("dMvn.cpp")
sourceCpp("rMvn.cpp")
source('mckFUNCTIONS.r')
source('clarkFunctions.R')
source( 'WarmGrow/clarkRcppFunctions.r' )

#style
par(cex=.7,pch=20,bty = 'n',bg = 'white',col = nColor()[1],
    col.axis = "#696969", col.lab = "#696969", col.main = "#696969")
palette(nColor(10,distinct = T))

#all input data
dDir = 'warming/Li6400/licorPROCESSED'
eDir = 'warming/ENVDATA/Read_Env/processedData'


points = read.csv(file = 'warming/Li6400/licorPROCESSED/all.rates.points.csv',stringsAsFactors = F)

#add instantaneous temperature

whr = which(points[,'Photo'] > 0 & 
            points[,'species'] == 'quru' &
            points[,'Light'] == 'Gap')
#plot(points[whr,'JD2009'],points[whr,'Photo'])
plot(points[whr,'JD2009'],points[whr,'Am'])

dayL = read.csv(file = 'DayLength.csv',stringsAsFactors = F)
dayL[1:which.max(dayL[,'Hours']),'Hours'] = max(dayL[,'Hours'])

each(mean,sd,range)(points$sm)
each(mean,sd,range)(points$sm)
each(mean,sd,range)(points$par)
each(mean,sd,range)(points$meanPar)
each(mean,sd,range)(points$Am,na.rm=T)
quantile(points$Am,c(.025,.5,.975),na.rm=T)
length(unique(points$Uindiv))

range(points$t)
points[which.min(points$sm),]
gapSM = mean(points$sm[grep('gap',all_treatment(points$Plot))],na.rm=T)
shdSM = mean(points$sm[grep('shade',all_treatment(points$Plot))],na.rm=T)
hotSM = mean(points$sm[grep('hot',all_treatment(points$Plot))],na.rm=T)
notSM = mean(points$sm[grep('zero',all_treatment(points$Plot))],na.rm=T)
mean(points$VpdA[grep('hot',all_treatment(points$Plot))],na.rm=T)
mean(points$VpdA[grep('zero',all_treatment(points$Plot))],na.rm=T)
mean(points$VpdL[grep('hot',all_treatment(points$Plot))],na.rm=T)
mean(points$VpdL[grep('zero',all_treatment(points$Plot))],na.rm=T)
mean(points$par[grep('gap',all_treatment(points$Plot))],na.rm=T)
mean(points$par[grep('shade',all_treatment(points$Plot))],na.rm=T)


(gapSM-shdSM) / shdSM
(hotSM-notSM) / notSM

species = sort(unique(points[,'species']))

gStom.fit = list()

for(nSp in 1:length(species)){
 # nSp=2

points = read.csv(file = 'warming/Li6400/licorPROCESSED/all.rates.points.csv',stringsAsFactors = F)
points = points[!is.na(points[,'R']),]
points = points[points[,'R'] < 0,]
points = points[points[,'species'] == species[nSp],]
#speciesLight = paste(points[,'species'],points[,'Light'],sep='')
#points = points[speciesLight == 'qualShade',]


#add instantaneous temperature
dfTemp = read.csv(paste(eDir,'gapfill-data-Duke-AT.csv',sep='/'))
times = points[,'JD2009'] +timeToFrac(points[,'HHMMSS'])
jdR = julian.date(c("01/01/2012","12/31/2013"),origin = c(month = 12, day = 31, year = 2008)) 
doys = times - julian.date(c("01/01/2012"),origin = c(month = 12, day = 31, year = 2008)) 
doys[doys > 366] = doys[doys > 366] - 366

dayL = read.csv(file = 'DayLength.csv',stringsAsFactors = F)
dayLength = dayL[match(floor(doys),dayL[,1]),2]
theDOY = floor(doys)

whr2 = apply(matrix(times,ncol=1),1, 
             function(x) nearWHICH(dfTemp[,'times'],x))

whr = match(points[,'Plot'],tolower(str_sub(colnames(dfTemp), start = 1, end = 3)))

points[,ncol(points)+1] = dfTemp[cbind(whr2,whr)]
colnames(points)[ncol(points)] = 'Tempi'

dfTemp = read.csv(paste(eDir,'gapfill-data-Duke-Q.csv',sep='/'))
times = points[,'JD2009'] +timeToFrac(points[,'HHMMSS'])

whr2 = apply(matrix(times,ncol=1),1, 
             function(x) nearWHICH(dfTemp[,'times'],x))

whr = match(points[,'Plot'],tolower(str_sub(colnames(dfTemp), start = 1, end = 3)))

points[,ncol(points)+1] = dfTemp[cbind(whr2,whr)]
colnames(points)[ncol(points)] = 'Pari'


### leaf temperature vs sensor air temperature
plot(points[,'Tleaf'],points[,'Tempi'],xlim = c(20,50))
abline(0,1,col=2)

fit = lm(Tleaf ~ t + dayFRAC + par, data = points)
plot(points[,'Tleaf'],predict(fit),xlim = c(20,50))
abline(0,1,col=2)

### gap = 1
fit = lm(points$Tleaf ~ points$Tempi + points$dayFRAC + points$Light, data = points)
plot(points[,'Tleaf'],predict(fit),xlim = c(20,50),col=as.factor(points[,'Light']))
abline(0,1,col=2)
fit = coefficients(fit)
fit = c(8.7580453,0.6793779,9.6999956,-1.7867873)
names(fit) = c('int','temp','dayFrac','light')

###make predictive design matrix xPred
c('ZeroXT','HotXT','dTemp','Sm')


dfTmp = read.csv(paste(eDir,'gapfill-data-Duke-AT.csv',sep='/'))

#grab dates
whr = which(dfTmp[,'JD2009'] >= jdR[1] & dfTmp[,'JD2009'] <= jdR[2] )
days = dfTmp[whr,'JD2009']
length(whr)
#indicators for temp treatment repeted for gap or not
xPred = rbind(
          cbind(rep(0,length(whr)),rep(0,length(whr)),rep(0,length(whr))),#000 n,z,s
          cbind(rep(0,length(whr)),rep(0,length(whr)),rep(1,length(whr))),#001 n,z,g
          cbind(rep(1,length(whr)),rep(0,length(whr)),rep(0,length(whr))),#100 c,z,s
          cbind(rep(1,length(whr)),rep(0,length(whr)),rep(1,length(whr))),#101 c,z,g
          cbind(rep(1,length(whr)),rep(1,length(whr)),rep(0,length(whr))),#110 c,h,s
          cbind(rep(1,length(whr)),rep(1,length(whr)),rep(1,length(whr)))#111 c,h,g
        )
xPred = cbind(rep(1,nrow(xPred)),xPred)
colnames(xPred) = c('Int','Cham','Hot','Light')

tmpX = rbind(
  cbind(rep(1,length(whr)),rep(0,length(whr)),rep(1,length(whr)),rep(1,length(whr))),#1011
  cbind(rep(1,length(whr)),rep(0,length(whr)),rep(1,length(whr)),rep(1,length(whr))),#1011
  cbind(rep(1,length(whr)),rep(0,length(whr)),rep(1,length(whr)),rep(1,length(whr))),#1011
  cbind(rep(1,length(whr)),rep(0,length(whr)),rep(1,length(whr)),rep(1,length(whr))),#1011
  cbind(rep(1,length(whr)),rep(1,length(whr)),rep(1,length(whr)),rep(1,length(whr))),#0111
  cbind(rep(1,length(whr)),rep(1,length(whr)),rep(1,length(whr)),rep(1,length(whr))) #0111
)
colnames(tmpX) = c('At','HotXT','dTemp','Sm')
xPred = cbind(xPred,tmpX)
dim(xPred)
#treat X temp interactions
dfTmp = read.csv(paste(eDir,'daily-Duke-AT.csv',sep='/'))
whr2 = match(days,dfTmp[,'JD2009'])
length(whr2)
dfTmpM = Tempvpd = dfTmp

##make DOY and day Length
predDOY = dfTmp[whr2,'JD2009'] - 365 - 365 - 365
predDOY[predDOY > 366] = predDOY[predDOY > 366] - 366

predDayL = dayL[match(predDOY,dayL[,1]),2]
XpredDayL = rep(predDayL,6)
XpredDOY = rep(predDOY,6)

dfTmpM[,10:33] = apply(dfTmp[,10:33],2,rollprev,roll = 14, include = T)

whrC = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Gapcontrol')
whr0 = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Gapzero')
whr3 = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Gapthree')
whr5 = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Gapfive')
whrCs = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Shadecontrol')
whr0s = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Shadezero')
whr3s = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Shadethree')
whr5s = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Shadefive')


tmpMat = rep(c(apply(dfTmpM[whr2,c(whr3s,whr5s)],1,mean),apply(dfTmpM[whr2,c(whr3,whr5)],1,mean)),3)
xPred[,'HotXT'] = tmpMat * xPred[,'HotXT']
dim(xPred)

tmpMat = c(apply(dfTmpM[whr2,c(whrCs)],1,mean),
           apply(dfTmpM[whr2,c(whrC)],1,mean),
           apply(dfTmpM[whr2,c(whr0s)],1,mean),
           apply(dfTmpM[whr2,c(whr0)],1,mean),
           apply(dfTmpM[whr2,c(whr3s,whr5s)],1,mean),
           apply(dfTmpM[whr2,c(whr3,whr5)],1,mean)
)
xPred[,'At'] = tmpMat * xPred[,'At']



dfTmp = read.csv(paste(eDir,'daily-Duke-SM.csv',sep='/'))
whr2 = match(days,dfTmp[,'JD2009'])
length(whr2)
dfTmpM = dfTmp
dfTmpM[,10:33] = apply(dfTmp[,10:33],2,rollprev,roll = 2, include = T)
whrC = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Gapcontrol')
whr0 = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Gapzero')
whr3 = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Gapthree')
whr5 = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Gapfive')
whrCs = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Shadecontrol')
whr0s = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Shadezero')
whr3s = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Shadethree')
whr5s = which(treatment(str_sub(colnames(dfTmp), start = 1, end = 3),combine =F,light=T) == 'Shadefive')

tmpMat = c(apply(dfTmpM[whr2,c(whrCs)],1,mean),
           apply(dfTmpM[whr2,c(whrC)],1,mean),
           apply(dfTmpM[whr2,c(whr0s)],1,mean),
           apply(dfTmpM[whr2,c(whr0)],1,mean),
           apply(dfTmpM[whr2,c(whr3s,whr5s)],1,mean),
           apply(dfTmpM[whr2,c(whr3,whr5)],1,mean)
           )
xPred[,'Sm'] = tmpMat * xPred[,'Sm']
dim(xPred)

dfTmp = read.csv(paste(eDir,'daily-Duke-AT.csv',sep='/'))
whr2 = match(days,dfTmp[,'JD2009'])
length(whr2)
dfTemp = read.csv(paste(eDir,'gapfill-data-Duke-AT.csv',sep='/'))
whr = which(dfTemp[,'JD2009'] >= jdR[1] & dfTemp[,'JD2009'] <= jdR[2] )
length(whr)
whrC = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Gapcontrol')
whr0 = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Gapzero')
whr3 = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Gapthree')
whr5 = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Gapfive')
whrCs = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Shadecontrol')
whr0s = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Shadezero')
whr3s = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Shadethree')
whr5s = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Shadefive')
tmpMat =  -xPred[,'At'] + 
  c(apply(dfTemp[whr,whrCs],1,mean,na.rm=T),
    apply(dfTemp[whr,whrC],1,mean,na.rm=T),
    apply(dfTemp[whr,whr0s],1,mean,na.rm=T),
    apply(dfTemp[whr,whr0],1,mean,na.rm=T),
    apply(dfTemp[whr,c(whr3s,whr5s)],1,mean,na.rm=T),
    apply(dfTemp[whr,c(whr3,whr5)],1,mean,na.rm=T)) 
  
  fitTemp = -xPred[,'At'] + 
  (fit[1] +
  fit[2] * 
    c(apply(dfTemp[whr,whrCs],1,mean,na.rm=T),
      apply(dfTemp[whr,whrC],1,mean,na.rm=T),
      apply(dfTemp[whr,whr0s],1,mean,na.rm=T),
      apply(dfTemp[whr,whr0],1,mean,na.rm=T),
      apply(dfTemp[whr,c(whr3s,whr5s)],1,mean,na.rm=T),
      apply(dfTemp[whr,c(whr3,whr5)],1,mean,na.rm=T)) +
  fit[3] * rep(dfTemp[whr,'dayFraction'],6) +
  fit[4] * xPred[,'Light']) 

xPred[,'dTemp'] = fitTemp * xPred[,'dTemp']

dfTemp = read.csv(paste(eDir,'gapfill-data-Duke-Q.csv',sep='/'))

whr = which(dfTemp[,'JD2009'] >= jdR[1] & dfTemp[,'JD2009'] <= jdR[2] )
length(whr)
whrC = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Gapcontrol')
whr0 = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Gapzero')
whr3 = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Gapthree')
whr5 = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Gapfive')
whrCs = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Shadecontrol')
whr0s = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Shadezero')
whr3s = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Shadethree')
whr5s = which(treatment(str_sub(colnames(dfTemp), start = 1, end = 3),combine =F,light=T) == 'Shadefive')

iPred =
  c(apply(dfTemp[whr,whrCs],1,mean,na.rm=T),
    apply(dfTemp[whr,whrC],1,mean,na.rm=T),
    apply(dfTemp[whr,whr0s],1,mean,na.rm=T),
    apply(dfTemp[whr,whr0],1,mean,na.rm=T),
    apply(dfTemp[whr,c(whr3s,whr5s)],1,mean,na.rm=T),
    apply(dfTemp[whr,c(whr3,whr5)],1,mean,na.rm=T)) 
iPred[iPred<0] = 0

photoSyn = n.rect.h(iPred,Amax = 5,Rd = -1,Q = .05,Theta = .1)

plot(rep(days,6),iPred,type='l')

dfTemp = read.csv('warming/ENVDATA/processedData/Treatment-RH-DF.csv')
whr = which(dfTemp[,'JD2009'] >= jdR[1] & dfTemp[,'JD2009'] <= jdR[2] )

XpredRH = c(dfTemp[whr,'Shadecontrol'],
            dfTemp[whr,'Gapcontrol'],
            dfTemp[whr,'Shadezero'],
            dfTemp[whr,'Gapzero'],
            dfTemp[whr,'ShadeH'],
            dfTemp[whr,'GapH'])
tmpMat = c(apply(Tempvpd[whr2,c(whrCs)],1,mean),
           apply(Tempvpd[whr2,c(whrC)],1,mean),
           apply(Tempvpd[whr2,c(whr0s)],1,mean),
           apply(Tempvpd[whr2,c(whr0)],1,mean),
           apply(Tempvpd[whr2,c(whr3s,whr5s)],1,mean),
           apply(Tempvpd[whr2,c(whr3,whr5)],1,mean)
)


XpredVPD = rhTOvpd(XpredRH,tmpMat)
            
          

xPred = cbind(xPred,XpredDayL) #Daylength
xPred = cbind(xPred,XpredDOY) #Day of year
xPred = cbind(xPred,XpredVPD) #Vpd
xPred = cbind(xPred,photoSyn) #Photosynthesis squared
xPred = cbind(xPred,tmpMat ) #photo X instantaneous temperature
xPred = cbind(xPred,photoSyn) #Photosynthesis squared
xPred = cbind(xPred,XpredVPD) #Photo X vpd
xPred = cbind(xPred,tmpMat ) #photo X instantaneous temperature
xPred = cbind(xPred,photoSyn) #photo X Vpd X instantaneous temp


#####
#Make xmat
write.csv(points, file = 'warming/Li6400/licorPROCESSED/allPoints.csv',row.names=F)

# Treatments
# treats = paste(points[,'Light'],points[,'Temp'],sep="") #with out PAR included
treats = points[,'Temp'] #with PAR included
chamber = points[,'Temp']
  chamber[chamber == 'Control'] = 0
  chamber[chamber != 0] = 1
  chamber = as.numeric(chamber)
hot = points[,'Temp']
  hot[hot == 'Control'] = 0
  hot[hot == 'Ambient'] = 0
  hot[hot != 0] = 1
  hot = as.numeric(hot)
at = points[,'meanT']
treatZ = (hot-1) * -1 * points[,'meanT']
treatH = hot * points[,'meanT']
#treatSm = treatD * points[,'sm']
#treatPar = treatD * points[,'meanPar']
treatL = rep(0,nrow(points))
treatL[points[,'Light'] == 'Gap'] = 1
int = rep(1,nrow(points))
sm =  points[,'sm']
par = points[,'meanPar']
I = points[,'PARi'] #incident light
P = points[,'Photo']
iAt = points[,'t']

doyTable =read.csv('warming/DayOfYearTable.csv')
library(lubridate)
dat = as.Date(points[,'Date'], format = "%m/%e/%Y")
datTable = as.Date(apply(doyTable[,c(3,4,5)],1,function(x)paste(x,collapse='-')),format = "%b-%e-%Y")

pointDOY = doyTable[match(dat,datTable),2]
pointDayL = dayL[match(pointDOY,dayL[,1]),2]

#differnce between two week mean and leaf temp
#thetaM = vaggregate(theta[,1],points[,'ID'],mean)
leaf = tapply(points[,'Tleaf'],points[,'ID'],mean)
whr = match(points[,'ID'],names(leaf))

deltaT = (leaf[whr] - points[,'meanT'])


tmp = read.csv('warming/ENVDATA/processedData/Treatment-RH-DF.csv')
pp = points[,'JD2009'] + points[,'dayFRAC']
pt = tmp[,'JD2009'] + tmp[,'dayFraction']
whr = fuzWHICH(pp,pt)
whr2 = match(treatment(points[,'Plot'],combine =T,light=T),colnames(tmp))
Rh = tmp[cbind(whr,whr2)]
points = cbind(points,Rh)
Vpd = rhTOvpd(points[,'Rh'],points[,'t'])
points = cbind(points,Vpd)

#basic x martrix
xmat = cbind(int,chamber,hot,treatL,at,hot*at,deltaT,sm,pointDayL,pointDOY,Vpd,P,iAt,P,Vpd,P,iAt)


colnames(xmat) = colnames(xPred) = c('Int','Chamb','Hot','Light','At','HXTemp','DTemp','Sm','DayL','DOY','Vpd','Photo','iAt',
                                     'P2','PXVpd','PXiAt','PXVpdXiAt')


each(mean,sd)(xmat[,'At'])
each(mean,sd)(xmat[,'HXTemp'])
each(mean,sd)(xmat[,'Sm'])
each(mean,sd)(xmat[,'DTemp'])
each(mean,sd)(xmat[,'DayL'])
each(mean,sd)(xmat[,'DOY'])
each(mean,sd)(xmat[,'Vpd'])
each(mean,sd)(xmat[,'Photo'])
each(mean,sd)(xmat[,'iAt'])

stdMS = rbind(c(24.257083,2.737504),c(0.12579359,0.04777856) ,c(8.625303,4.540125),c(13.9828759 ,0.6635037 ),c(195.1239,35.81265),c(0.9831747,0.6612403 ))
colnames(stdMS) = c('mean','std')
rownames(stdMS) = c('temp','sm','dtemp','DayL','DOY','Vpd')

xmat[,'At'] = (xmat[,'At'] - 24.257083) / 2.737504
xmat[,'Sm'] = (xmat[,'Sm'] - 0.17579359) / 0.04777856
xmat[,'DTemp'] = (xmat[,'DTemp'] - 8.625303) / 4.540125
xmat[,'DayL'] = (xmat[,'DayL'] - 13.9828759) / 0.6635037
xmat[,'DOY'] = (xmat[,'DOY'] - 195.1239) / 35.81265
xmat[,'Vpd'] = (xmat[,'Vpd'] - 0.941747 ) /0.6605477 
xmat[,'Photo'] = (xmat[,'Photo'] - 2.510538 ) /3.931949 
xmat[,'iAt'] = (xmat[,'iAt'] - 24.507388 ) / 3.502204 
xmat[,'HXTemp'] = xmat[,'At'] * xmat[,'Hot']
xmat[,'PXiAt'] = xmat[,'Photo'] * xmat[,'iAt'] 
xmat[,'P2'] = xmat[,'Photo'] * xmat[,'Photo']
xmat[,'PXVpd'] = xmat[,'Photo'] * xmat[,'Vpd'] 
xmat[,'PXiAt'] = xmat[,'Photo'] * xmat[,'iAt'] 
xmat[,'PXVpdXiAt'] = xmat[,'Photo'] * xmat[,'Vpd'] * xmat[,'iAt'] 

xPred[,'At'] = (xPred[,'At'] - 24.257083) / 2.737504
xPred[,'Sm'] = (xPred[,'Sm'] - 0.17579359) / 0.04777856
xPred[,'DTemp'] = (xPred[,'DTemp'] - 8.625303) / 4.540125
xPred[,'DayL'] = (xPred[,'DayL'] - 13.9828759) / 0.6635037
xPred[,'DOY'] = (xPred[,'DOY'] - 195.1239) / 35.81265
xPred[,'Vpd'] = (xPred[,'Vpd'] - 0.9831747 ) /0.6612403 
xPred[,'Photo'] = (xPred[,'Photo'] - 2.510538 ) /3.931949 
xPred[,'iAt'] = (xPred[,'iAt'] - 24.507388 ) / 3.502204 
xPred[,'HXTemp'] = xPred[,'At'] * xPred[,'Hot']
xPred[,'PXiAt'] = xPred[,'Photo'] * xPred[,'iAt'] 
xPred[,'P2'] =  xPred[,'Photo'] * xPred[,'Photo']
xPred[,'PXVpd'] = xPred[,'Photo'] * xPred[,'Vpd'] 
xPred[,'PXiAt'] = xPred[,'Photo'] * xPred[,'iAt'] 
xPred[,'PXVpdXiAt'] = xPred[,'Photo'] * xPred[,'Vpd'] * xPred[,'iAt'] 

each(mean,sd)(xmat[,'At'])
each(mean,sd)(xmat[,'HXTemp'])
each(mean,sd)(xmat[,'Sm'])
each(mean,sd)(xmat[,'DTemp'])

#remove interaction
xmatHold = xmat
xpredHold = xPred
xmatHold -> xmat
xpredHold -> xPred

whr.na = match(c( "Int","Chamb","Hot","Photo","P2","PXVpd","PXiAt","PXVpdXiAt"),colnames(xmat))
xPred2 = xPred[,whr.na] 
xmat2 = xmat[,whr.na] 

whr.na = match(c("Int","Chamb","Hot","At","HXTemp","DTemp","Sm","DayL"),colnames(xmat))
xPred = xPred[,whr.na] 
xmat = xmat[,whr.na] 

save(fit,stdMS,file = 'warming/forAssPredNoX.Rdata')

#make response variable
Y1 = points[,'Am']
Y2 = points[,'R']
Y3 = points[,'Q']
Y4 = points[,'Th']
Y5 = points[,'Cond']
P = points[,'Photo']

miss.photo = which(is.na(P)) #no observation 
obs.photo = which(!is.na(P)) #observation

Y = cbind(Y1,Y2,Y3,Y4) #wide
miss = which(is.na(Y),arr.ind=T) #no observation wide
obs = which(!is.na(Y),arr.ind=T) #observation wide
y = c(Y1,Y2,Y3,Y4) #long
miss.long = which(is.na(y)) #no observation 
obs.long = which(!is.na(y)) #observation

Y1[is.na(Y1)] = mean(Y1,na.rm=T) #fill gaps
Y2[is.na(Y2)] = mean(Y2,na.rm=T) #fill gaps
Y3[is.na(Y3)] = mean(Y3,na.rm=T) #fill gaps
Y4[is.na(Y4)] = mean(Y4,na.rm=T) #fill gaps
Y5[is.na(Y5)] = mean(Y5,na.rm=T) #fill gaps
P[is.na(P)] = mean(P,na.rm=T) #fill gaps

Y = cbind(Y1,Y2,Y3,Y4) #wide
y = c(Y1,Y2,Y3,Y4) #long

library(MASS)


whr = which(points[,'PARi'] > 2) # & points[,'Cond'] < .01*points[,'Photo'] + .001*points[,'Photo']*points[,'Photo'] +.08)
xmat2 = xmat2[whr,]  #| grepl('LC',points[,'instance'])
Y5 = Y5[whr ]
fit = lm(Y5~xmat2-1)
plot(P[whr],Y5)
pp = predict(fit,se.fit = F)
points(P[whr],pp,pch=1.3,col=4)


plot(P[whr],Y5,pch=1.3,col=4)
Yp  = Y5; Yp[P[whr] < 0 ] = 0
fit = smooth.spline(x = P[whr],y = Yp,nknots = 4)
pp = predict(fit,P[whr])
sap = pp[[2]]
sap[sap < 0] = quantile(sap[-0.5 <= P[whr] & P[whr] < .5],.05)
points(P[whr],sap)
Yp[Yp > (2*sap)] = sap[Yp > (2*sap)]
fit = smooth.spline(x = P[whr],y = Yp,nknots = 4)
pp = predict(fit,P[whr])
Yp[P[whr] < 0] = quantile(sap[-0.5 <= P[whr] & P[whr] < .5],.025)
fit = smooth.spline(x = P[whr],y = Yp,nknots = 4)
pp = predict(fit,P[whr])
plot(P[whr],Yp,pch=1.3,col=4)
points(pp)

predict(fit,0)
Ym = Y5[P[whr] < 0 ]
plot(P[whr][P[whr] < 0],Ym)

plot(P[whr],Yp,pch=1.3,col=4)
points(pp)
abline(predict(fit,0)$y,coefficients(fitP)[2])

#gStom.fit[[nSp]] = fitP 
#gStom.fit[[nSp+4]] = fit 
#} 
#saveRDS(gStom.fit,file = 'gStomFit.rds')


whr =grep('g',rownames(xmat))
#light = rep(0,nrow(xmat))
#light[whr] = 1
#xmat = cbind(xmat,light)

#xmat[xmat[,'Par'] < 0,'Par'] = 0

Xstar = bdiag.spam(xmat,xmat,xmat,xmat) #wide
X = xmat #long

# Data
n <- length(Y1) #number of 
n.o = length(obs.long)
p<-ncol(xmat)  # Number of covariates
k<-ncol(Y)           # Number of outcomes

pG = ncol(xmat2) #covariates in gs
nG = nrow(xmat2)

# Priors
nu0<-4		    # prior df for Sigma
c0<-diag(k)	# Prior precision (INVERSE scale) for Sigma
# So if you want non-diag scale, D, then c0=solve(D)

b0<-matrix(0,k,p)	  # Prior mean of (b11 b12
#                b21 b22)
T0<-diag(.01,p) 	    # p x p Prior Precision beta1 and beta2, where p denotes number of Parameters
# Note that prior for 4 x 1 beta is actually 4 x 4 T0 %x% sigma^{-1} -- this is like the usual "conjugate" MVN-W prior
# Analogous to N-IG prior in the univariate case
# Inits
beta<-rep(0,k*p)
sigma<-diag(k)      # Covariance of Y_1,..Y_k

nu0<-3
sig0 = .5
beta0 <- c(0,3,2,1)

#process on y priors
pr.mu = 1
s1 = 10
s2 = pr.mu*(s1 -1)

pr.muG = 1
s1G = 10
s2G = pr.mu*(s1 -1)

bG0<-rep(0,pG)    # Prior mean for Gs
TG0<-diag(.01,pG)    # Prior mean for Gs

#observation error priors
ob.mu = .5
so1 = 10
so2 = ob.mu*(so1 -1)

#theta truncations
lo.t = c(0.001,-10,0.0001,0.001)
hi.t = c(30,-.001,.99,.99)

#betaG truncations
loBG = rep(-Inf,pG)
hiBG = rep(Inf,pG)

#######
#Inits#
#######
theta<-c(5,-1,0.1,.1)
Tau<-diag(4)  # Sigma^{-1}
sig = 1
tg = .5
beta = c(1,4,1,1)
counter = 0

# Store
nsim<-500
pb <- txtProgressBar(min = 1, max = nsim, style = 3)
burnin = round(nsim/4)

Sigmas<-matrix(0,nsim,k*k)	# Store Samples
Beta<-matrix(0,nsim,p*k)

Sig = rep(0,nsim); Sig[1] <- sig
tgibbs = rep(0,nsim); tgibbs[1] <- tg
Ynew = matrix(0,n,4)

sigG = 1
SigG = rep(0,nsim); SigG[1] <- sigG
BetaG<-matrix(0,nsim,pG)

Y.o = thetaold = Y
y.o =  y
Pold = rnorm(len(P),P,2)

Asum = A2sum = Rsum = R2sum = Qsum = Q2sum = Thsum = Th2sum = Psum = P2sum = rep(0,n)
Acount = Rcount = Qcount = Thcount = rep(0,n)
allTotsum = allTot2sum = allcount = rep(0,length(iPred))
allPsum =  allP2sum = totcount = rep(0,length(iPred))
allGsum = allG2sum = allGcount = rep(0,length(iPred))

Sigmamulti = rep(1,n)
accept = acceptI= ii =  0
########
# MCMC #
########

ptm<-proc.time()
for (i in 2:nsim) {
  
  # Update beta (update V and m using Matrix Normal form, then vectorize to form mstar)
  V<-solve(T0+crossprod(X))     # p x p post. covariance = 2 x 2 here
  m<-V%*%(T0%*%t(b0) + t(X)%*%Y)   # p x k = 2 x 2, where p is number of covariates (including int.) and k is number of outcomes
  mstar<-c(m)                   # p*k x 1 = (b11, b12, b21, b22)
  Beta[i,]<-beta<-c(rmvnorm(1,mstar,sigma%x%V))
  betastar<-matrix(beta,p,k)    # Return to matrix normal form for updating Sigma
  
  # Update Sigma
  resid<-matrix(y-Xstar%*%beta,n,k,byrow=F)
  sigma<-riwish(n+nu0+p,c0+crossprod(resid,resid)+t(betastar-t(b0))%*%T0%*%(betastar-t(b0)))    # Note: need to account for mat-norm prior for beta|sigma
  #       in both df and in posterior scale
  Sigmas[i,]<-c(sigma)
  
  #predict theta(Ys) and share across response  #note should balance with observations
  #theta = X%*%betastar
    #tnorm.mvtRcpp
  #theta = t(apply(theta,1,function(x) rmvnorm(1,x,sigma)))
  #theta = t(apply(theta,1,function(x) tnorm.mvtRcpp(x,x,sigma,lo = lo.t,hi=hi.t,times=5)))
  #thetaM = vaggregate(theta[,1],points[,'ID'],mean)
  #thetaM = tapply(theta[,1],points[,'ID'],mean)
  #whr = match(points[,'ID'],names(thetaM))
  #theta[,1] = thetaM[whr]
  #theta[,2] = tapply(theta[,2],points[,'ID'],mean)[whr]
  #theta[,3] = tapply(theta[,3],points[,'ID'],mean)[whr]
  #theta[,4] = tapply(theta[,4],points[,'ID'],mean)[whr]
  
  #Y[miss] = theta[miss] #predict missing thetas
  #(theta*bg)/sg + sum(y[which(timeI==time[1])],na.rm=T)/tg 3 balance with observations  
  
  #y = c(Y[,1],Y[,2],Y[,3],Y[,4])
  
  #observation error on y
  #w1 <- so1 + n.o/2
  #w2 <- so2 + .5*sum( (Pold[obs.photo] - n.rect.h(I[obs.photo],thetaold[obs.photo,1],thetaold[obs.photo,2],thetaold[obs.photo,3],thetaold[obs.photo,4]))^2 )
  #tgibbs[i] <- tg <-  1/rgamma(1,w1,w2)
  
  #process error
  u1 <- s1 + n/2
  u2 <- s2 + .5*sum( (Pold - n.rect.h(I,thetaold[,1],thetaold[,2],thetaold[,3],thetaold[,4]))^2 )
  Sig[i] = sig = 1/rgamma(1,u1,u2)
  
  #updateY
  if(T){
    
    #propose new thetas
    thetanew = X%*%betastar
    tempTnew = cbind(thetanew,Sigmamulti)
    thetanew = t(apply(tempTnew,1,function(x) tnorm.mvtRcpp(x[1:4],x[1:4],sigma*x[5],lo = lo.t,hi=hi.t,times=5)))
    thetaM = tapply(thetanew[,1],points[,'ID'],last)
    whr = match(points[,'ID'],names(thetaM)) #share parameters across indiv samples
    thetanew[,1] = thetaM[whr]
    thetanew[,2] = tapply(thetanew[,2],points[,'ID'],last)[whr]
    thetanew[,3] = tapply(thetanew[,3],points[,'ID'],last)[whr]
    thetanew[,4] = tapply(thetanew[,4],points[,'ID'],last)[whr]
    #near observations
    #thetanew = t(apply(theta,1,function(x) tnorm.mvtRcpp(x,x,sigma/2,lo = lo.t,hi=hi.t,times=5)))
    #thetanewhold = rnorm(nrow(obs),Y.o[obs[,2]],sigma/2)
    #thetanew = t(apply(thetanew,1,function(x) tnorm.mvtRcpp(x,x,sigma/2,lo = lo.t,hi=hi.t,times=5)))
    
    
    #P from proposed theta 
    Pnew = n.rect.h(I,thetanew[,1],thetanew[,2],thetanew[,3],thetanew[,4])  
    #Pnew[obs.photo] = rnorm(length(obs.photo),P[obs.photo],sqrt(tg)) #near observation
    
    whr = which(is.na(Pnew)) #remove any bad theta combinations
    if(length(whr) > 0){
      thetanew2 = t(apply(tempTnew[whr,],1,function(x) tnorm.mvtRcpp(x,x,sigma/100,lo = lo.t,hi=hi.t,times=5)))
      thetaM = tapply(thetanew2[,1],points[,'ID'],last)
      whr2 = match(points[,'ID'],names(thetaM)) #share parameters across indiv samples
      thetanew2[whr2,1] = thetaM[whr2]
      thetanew2[whr2,2] = tapply(thetanew2[,2],points[,'ID'],last)[whr2]
      thetanew2[whr2,3] = tapply(thetanew2[,3],points[,'ID'],last)[whr2]
      thetanew2[whr2,4] = tapply(thetanew2[,4],points[,'ID'],last)[whr2]
      thetanew[whr] =  thetanew2
      Pnew[whr] = n.rect.h(I,thetanew2[whr,1],thetanew2[whr,2],thetanew[whr,3],thetanew[whr,4])  
    }
    whr = which(is.na(Pnew)) #remove any bad theta combinations
    if(length(whr) > 0){
      Pnew[whr] = Pold[whr]
      if(i < burnin) Pold[whr] = rnorm(len(whr),P[whr],2)
      thetanew[whr,] = thetaold[whr,]
    }
    
    pnow <- rep(0,n)
    pnew <- pnow
    
    #observation error
    #pnow[obs.photo] <- dnorm(P[obs.photo],Pold[obs.photo],sqrt(tg),log=T)
    #pnew[obs.photo] <- dnorm(P[obs.photo],Pnew[obs.photo],sqrt(tg),log=T)
    
    #process error on p
    pnow <- pnow + dnorm(P, Pold , sqrt(sig),log=T)
    pnew <- pnew + dnorm(P, Pnew , sqrt(sig),log=T)
    
    #parameter error
    #pnow = pnow + dmvnorm(Y.o[miss],thetaold[miss], sigma[1:2,1:2], log = T)
    #pnew = pnew + dmvnorm(matrix(Y.o[miss[,1],],ncol=4)[,c(1,2)],matrix(thetanew[miss[,1],],ncol=4)[,c(1,2)], sigma[1:2,1:2], log = T)
    pnow[miss[,1]] = pnow[miss[,1]] + 
          dmvn_arma(matrix(Y.o[miss[,1],],ncol=4)[,c(1,2)],
                      matrix(thetaold[miss[,1],],ncol=4)[,c(1,2)], 
                      sigma[1:2,1:2], logd = T)    
    pnew[miss[,1]] = pnew[miss[,1]] + 
          dmvn_arma(matrix(Y.o[miss[,1],],ncol=4)[,c(1,2)],
                      matrix(thetanew[miss[,1],],ncol=4)[,c(1,2)], 
                      sigma[1:2,1:2], logd = T)
    
    #pnow = pnow + dmvnorm(Y.o[obs],thetaold[obs], sigma, log = T)
    #pnew = pnew + dmvnorm(Y.o[obs],thetanew[obs[,2]], sigma, log = T)
    
    pnow = pnow + dmvn_arma(Y.o,thetaold,sigma, logd = T)    
    pnew = pnew + dmvn_arma(Y.o,thetanew,sigma, logd = T)

    anow <- pnow
    anew <- pnew
    
    a <- exp(anew - anow)
    z <- runif(1,0,1)
    
    w = which(z < a)
    accept = accept + (z < a)
    acceptI = acceptI + (z < a) #accept over last interval
    ii = ii + 1
    if (i%%100==0){
      #plot(accept/(i-1),main = mean(accept/(i-1)))
      Sigmamulti[acceptI/ii < .15] = Sigmamulti[acceptI/ii < .15] *.75
      Sigmamulti[acceptI/ii > .4] = Sigmamulti[acceptI/ii < .4] * 1.25
      Sigmamulti[Sigmamulti > 1.75] = 1.75
      acceptI = 0
      ii = 0
    }
    
    thetaold[w,]  <- thetanew[w,]

  }
  
  ##Add assimilation estimates for each treatment
  thetaP = xPred%*%betastar
  #thetaP = apply(thetaP,2,rollprev,168)
  thetaP[thetaP[,1] < lo.t[1],1] = lo.t[1];thetaP[thetaP[,1] > hi.t[1],1] = hi.t[1]
  thetaP[thetaP[,2] < lo.t[2],2] = lo.t[2];thetaP[thetaP[,2] > hi.t[2],2] = hi.t[2]
  thetaP[thetaP[,3] < lo.t[3],3] = lo.t[3];thetaP[thetaP[,3] > hi.t[3],3] = hi.t[3]
  thetaP[thetaP[,4] < lo.t[4],4] = lo.t[4];thetaP[thetaP[,4] > hi.t[4],4] = hi.t[4]
  thetaP = t(apply(thetaP,1,function(x) tnorm.mvtRcpp(x,x,sigma,lo = lo.t,hi=hi.t,times=2)))
  allP = rnorm(length(iPred),n.rect.h(iPred,thetaP[,1],thetaP[,2],thetaP[,3],thetaP[,4]),sqrt(sig))  
  allTot = cumsum(allP)
  
  
  xPred2[,'Photo'] = (allP - 2.510538) / 3.931949 
  xPred2[,'P2'] = xPred2[,'Photo'] * xPred2[,'Photo']
  xPred2[,'PXVpd'] = xPred2[,'Photo'] * (xpredHold[,'Vpd'] - 0.9831747 ) /0.6612403 
  xPred2[,'PXiAt'] = xPred2[,'Photo'] * (xpredHold[,'iAt'] - 24.507388 ) / 3.502204 
  xPred2[,'PXVpdXiAt'] = xPred2[,'Photo'] * (xpredHold[,'Vpd'] - 0.9831747 ) /0.6612403 * (xpredHold[,'iAt'] - 24.507388 ) / 3.502204 
  
  ##update beta for conductance 
    
  IV <- crossprod(xmat2)/(sigG) + TG0
  V  <- chol2inv(chol(IV))
  v  <- crossprod(xmat2,Y5)/(sigG) + TG0%*%bG0
  mu <- V%*%v
  betaG = BetaG[i,] <- matrix(tnorm.mvt(mu,mu,V,loBG,hiBG),ncol=1)
  
  ## update sig for conductane
  
  u1 <- s1G + nG/2
  u2 <- s2G + .5*sum( (Y5 - xmat2%*%betaG)^2 )
  sigG = SigG[i] = 1/rgamma(1,u1,u2)
  
  tmpG =  xPred2%*%betaG
  tmpG[tmpG <= 0] = 0.00001
  tmpG[which(is.na(xPred2%*%betaG))] = mean( tmpG,na.rm=T)
  allGs =  tnorm(nrow(xPred2),0,1,tmpG, sigG)
  
  ##Predict conductance
  
  if(i > burnin){
    counter <- counter + 1

    Asum  <- rowSums(cbind(Asum,thetaold[,1]),na.rm=T)
    A2sum <- rowSums(cbind(A2sum,thetaold[,1]^2),na.rm=T)
    Acount <- rowSums(cbind(Acount,!is.na(thetaold[,1])),na.rm=T)
    
    Rsum  <- rowSums(cbind(Rsum,thetaold[,2]),na.rm=T)
    R2sum <- rowSums(cbind(R2sum,thetaold[,2]^2),na.rm=T)
    Rcount <- rowSums(cbind(Rcount,!is.na(thetaold[,2])),na.rm=T)
    
    Qsum  <- rowSums(cbind(Qsum,thetaold[,3]),na.rm=T)
    Q2sum <- rowSums(cbind(Q2sum,thetaold[,3]^2),na.rm=T)
    Qcount <- rowSums(cbind(Qcount,!is.na(thetaold[,3])),na.rm=T)
    
    Thsum  <- rowSums(cbind(Thsum,thetaold[,4]),na.rm=T)
    Th2sum <- rowSums(cbind(Th2sum,thetaold[,4]^2),na.rm=T)
    Thcount <- rowSums(cbind(Thcount,!is.na(thetaold[,4])),na.rm=T)
    
    allPsum  <- rowSums(cbind(allPsum,allP),na.rm=T)
    allP2sum <- rowSums(cbind(allP2sum,allP^2),na.rm=T)
    allcount <- rowSums(cbind(allcount,!is.na(allP)),na.rm=T)
    
    allTotsum  <- rowSums(cbind(allTotsum,allTot),na.rm=T)
    allTot2sum <- rowSums(cbind(allTot2sum,allTot^2),na.rm=T)
    totcount <- rowSums(cbind(totcount,!is.na(allP)),na.rm=T)
    
    allGsum  <- rowSums(cbind(allGsum,allGs),na.rm=T)
    allG2sum <- rowSums(cbind(allG2sum,allGs^2),na.rm=T)
    allGcount <- rowSums(cbind(allGcount,!is.na(allGs)),na.rm=T)
    
  }
  
  progress(i,nsim,100)
  
}
close(pb)
proc.time()-ptm

###########
# Summary #
###########
m.beta<-matrix(colMeans(Beta[burnin:nsim,]),nrow=4,byrow=T)    # True vals: 5, -2; -5, 3
matrix(apply(Beta[burnin:nsim,],2,quantile,.975),nrow=4,byrow=T)
matrix(apply(Beta[burnin:nsim,],2,quantile,.0275),nrow=4,byrow=T)
colnames(m.beta) = colnames(xmat)
m.beta

m.betaG = apply(BetaG[burnin:nsim,],2,mean)

m.sigma<-colMeans(Sigmas[burnin:nsim,]) # True vals: 1, -1 ,-1 3
plot(burnin:nsim, Beta[burnin:nsim,3],type="l",col="lightgreen")
abline(h=m.beta[3],col="blue")
matrix(m.sigma,4,4)
m.sig = mean(Sig[burnin:nsim])

thetaP = xPred%*%t(m.beta)
thetaP = n.rect.h(iPred,thetaP[1],thetaP[2],thetaP[3],thetaP[4])
thetaG = xPred2%*%m.betaG

plot(thetaP,thetaG)

Amean <- Asum/Acount
Asd   <- sqrt(A2sum/Acount - Amean^2)  
Asd[is.na(Asd)] = 0
Rmean <- Rsum/Rcount
Rsd   <- sqrt(R2sum/Rcount - Rmean^2) 
Rsd[is.na(Rsd)] = 0
Qmean <- Qsum/Qcount
Qsd   <- sqrt(Q2sum/Qcount - Qmean^2) 
Qsd[is.na(Qsd)] = 0
Thmean <- Thsum/Qcount
Thsd   <- sqrt(Th2sum/Thcount - Thmean^2) 
Thsd[is.na(Thsd)] = 0

allPmean <- allPsum/allcount
allPsd   <- sqrt(allP2sum/allcount - allPmean^2) 
allTotmean <- allTotsum/totcount
allTotsd   <- sqrt(allTot2sum/totcount - allTotmean^2) 
allGmean <- allGsum/allGcount
allGsd   <- sqrt(allG2sum/allGcount - allGmean^2) 

par(mfrow=c(2,2))
plot(allPmean,allGmean)
plot(allPmean[10000:12000],allGmean[10000:12000])
plot(allGmean[10000:12000])
plot(allPmean[10000:12000])
par(mfrow=c(3,2))
wh = which(1246 < days & days < 1396)
plot(allPmean[wh],main= 'Shade Control');abline(h=0,col=2)
plot(allPmean[wh+17489],main= 'Gap Control');abline(h=0,col=2)
plot(allPmean[wh+17489*2],main= 'Shade Zero');abline(h=0,col=2)
plot(allPmean[wh+17489*3],main= 'Gap Zero');abline(h=0,col=2)
plot(allPmean[wh+17489*4],main= 'Shade Hot');abline(h=0,col=2)
plot(allPmean[wh+17489*5],main= 'Gap Hot');abline(h=0,col=2)
par(mfrow=c(1,1))

save.image(file=paste('warming/',species[nSp],'MVNassimilation11.Rdata',sep=""))

}

load(file='warming/acruMVNassimilation.Rdata')

sam = round(runif(1,1,17489))
plot(allPmean[sam:(sam+500)],type='l',ylim = c(-20,20))
lines(allPmean[sam:(sam+500)]+allPsd[sam:(sam+500)],col=2)
lines(allPmean[sam:(sam+500)]-allPsd[sam:(sam+500)],col=2)
plot(iPred[1:(17489*1)],allPmean[(1):(17489*1)])
plot(iPred[(1+17489*1):(17489*2)],allPmean[(1+17489*1):(17489*2)])
plot(iPred[(1+17489*2):(17489*3)],allPmean[(1+17489*2):(17489*3)])
plot(iPred[(1+17489*3):(17489*4)],allPmean[(1+17489*3):(17489*4)])
plot(iPred[(1+17489*4):(17489*5)],allPmean[(1+17489*4):(17489*5)])
plot(iPred[(1+17489*5):(17489*6)],allPmean[(1+17489*5):(17489*6)])

i=1
whr = points[,'ID'] == unique(points[,'ID'])[i]
plot(points[whr,'PARi'],points[whr,'Photo'], ylim = c(min(points[whr,'Photo'])-2,max(points[whr,'Photo'])+5))
unique(points[,'ID'])[i]
points(allPmean[11230:11330])
abline(h=points[6,'Photo'])
  
head(beta)
allIndiv = unique(points[,'ID'])

theta = X%*%t(m.beta)
#theta = t(apply(tt,1,function(x) tnorm.mvtRcpp(x,x,sigma,lo = lo.t,hi=hi.t)))
thetaM = tapply(theta[,1],points[,'ID'],mean)
whr = match(points[,'ID'],names(thetaM)) #share parameters across indiv samples
theta[,1] = thetaM[whr]
theta[,2] = tapply(theta[,2],points[,'ID'],mean)[whr]
theta[,3] = tapply(theta[,3],points[,'ID'],mean)[whr]
theta[,4] = tapply(theta[,4],points[,'ID'],mean)[whr]


whr = which(points[,'ID'] == allIndiv[i])

plot(points[whr,'PARi'],points[whr,'Photo'], ylim = c(min(points[whr,'Photo'])-2,max(points[whr,'Photo'])+5))
lines(0:2000,n.rect.h(0:2000,theta[whr,1],theta[whr,2],theta[whr,3],theta[whr,4]),col=2)
i = i+1


pheno = read.csv('/home/mckwit/Google Drive/Budburst/indivPHENO.csv',stringsAsFactors=F)
pheno = pheno[pheno[,'site'] == 'DF',]

whr = match(pheno[,'Species'],species)
whr = which(!is.na(whr))
pheno = pheno[whr,]

cham = gsub("_[a-zA-Z0-9]",'',pheno[,'Chamber'])
head(pheno)


lig = rep('S',nrow(pheno))
lig[grep('G',cham, ignore.case=T)] = 'G'
tem = treatment(cham,combine=F)
tem[grep('control',tem, ignore.case=T)] = 'con'
tem[grep('zero',tem, ignore.case=T)] = 'zero'
tem[grep('three',tem, ignore.case=T)] = 'warm'
tem[grep('five',tem, ignore.case=T)] = 'warm'
temLig = paste(pheno[,'Species'],tem,lig,sep="")

devel3 = rbind(matrix(tapply(pheno[,'Y3D3'],temLig,median,na.rm=T),4,6,byrow=T),
               matrix(tapply(pheno[,'Y3D5'],temLig,median,na.rm=T),4,6,byrow=T))
devel5 = rbind(matrix(tapply(pheno[,'Y4D5'],temLig,median,na.rm=T),4,6,byrow=T),
               matrix(tapply(pheno[,'Y4D5'],temLig,median,na.rm=T),4,6,byrow=T))
colnames(devel3) = colnames(devel5) = c('warmG','zeroG','conG','warmS','zeroS','conS')
rownames(devel3) = rownames(devel5) = paste(c(sort(species),sort(species)),c(3,3,3,3,5,5,5,5),sep='/')

load('/home/mckwit/Google Drive/Budburst/bbModel.RData')
