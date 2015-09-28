#####
#Build plots for mvn assimilation model
#
#MCKwit March 2015
####
path = 'C:/Users/mck14/Dropbox/'
path = '/home/mckwit/Dropbox/'

#style
par(cex=.7,pch=20,bty = 'n',bg = 'white',col = nColor()[1],
    col.axis = "#696969", col.lab = "#696969", col.main = "#696969")
palette(nColor(10,distinct = T))
palette(c("#00A0B0","#6A4A3C","#CC333F","#EB6841","#EDC951"))
thecol=c("#00A0B0","#6A4A3C","#CC333F","#EB6841","#EDC951")
thecol=c("#00A0B036","#6A4A3C36","#CC333F36","#EB684136","#EDC95136")

load(file=paste(path,'warming/qualMVNassimilation.Rdata',sep=''))

allA = allR = allQ = allTh = m.beta
allA.lo = allR.lo = allQ.lo = allTh.lo = m.beta
allA.hi = allR.hi = allQ.hi = allTh.hi = m.beta
species = sort(species)
rownames(allA) = rownames(allR) = rownames(allQ) = rownames(allTh) = sort(species)
allTab = as.character()
for( ii in 1:length(species)){
  iii = ii
  load(file=paste(path,'warming/',sort(species)[ii],'MVNassimilation.Rdata',sep = ""))
  ii = iii
  m.beta.lo = matrix(apply(Beta[burnin:nsim,],2,quantile,.05),ncol=8,byrow=T)
  m.beta.hi = matrix(apply(Beta[burnin:nsim,],2,quantile,.95),ncol=8,byrow=T)
  allA[ii,] = m.beta[1,]  ; allA.lo[ii,]  = m.beta.lo[1,]; allA.hi[ii,]   = m.beta.hi[1,];
  allR[ii,] = m.beta[2,]  ; allR.lo[ii,]  = m.beta.lo[2,]; allR.hi[ii,]   = m.beta.hi[2,];
  allQ[ii,] = m.beta[3,]  ; allQ.lo[ii,]  = m.beta.lo[3,]; allQ.hi[ii,]   = m.beta.hi[3,];
  allTh[ii,] = m.beta[4,] ; allTh.lo[ii,] = m.beta.lo[4,]; allTh.hi[ii,]  = m.beta.hi[4,];
  allTab = rbind(allTab,
                 c(species[ii], 'Am',paste(round(m.beta[1,],2),"(", round(m.beta.lo[1,],2),"," ,round(m.beta.hi[1,],2),")",sep="")),
                 c(species[ii], 'R',paste(round(m.beta[2,],2),"(", round(m.beta.lo[2,],2),"," ,round(m.beta.hi[2,],2),")",sep="")),
                 c(species[ii], 'Q',paste(round(m.beta[3,],3),"(", round(m.beta.lo[3,],3),"," ,round(m.beta.hi[3,],3),")",sep="")),
                 c(species[ii], 'Th',paste(round(m.beta[4,],3),"(", round(m.beta.lo[4,],3),"," ,round(m.beta.hi[4,],3),")",sep=""))
                 )
}
colnames(allTab) = c("Spec","Param",colnames(m.beta))
write.csv(allTab,file = paste(path,'warming/allTab.csv',sep=''),row.names=F)
save(allA,allR,allQ,allTh,
     allA.lo,allR.lo,allQ.lo,allTh.lo,
     allA.hi,allR.hi,allQ.hi,allTh.hi,
     file = paste(path,'warming/allBetas.RData',sep=''))


par(mfrow=c(1,3),mar=c(3, 1, 2, 1), oma=c(0,0,0,0),xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
kk=1
for(k in c('allA','allR','allQ','allTh')){

co = get(k); co.lo = get(paste(k,'.lo',sep='')); co.hi = get(paste(k,'.hi',sep=''))
co = co[,c(-1,-2)];co.hi = co.hi[,c(-1,-2)];co.lo = co.lo[,c(-1,-2)]
num = 4*ncol(co) + ncol(co)
step = 1/num
ys = step
minX = min(c(co.lo,co.hi))
maxX = max(c(co.lo,co.hi))
minY = 0; maxY = 1
xTit = ''; yTit=''

Title = c('Max Photosynthesis','Respiration','Light Use Efficiency','Curvature')[kk]
kk=1+kk
plot(0, xlim=c(minX, maxX), ylim=c(minY, maxY), type="n", bty="n", las=1, 
     main='', xlab=bquote(bold(.(xTit))), ylab=bquote(bold(.(yTit))), 
     family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8, xaxt="n", yaxt="n")
for(i in seq(0,ncol(co),by=2)){
  rect(minX-10,i*5*step,
       maxX+10,i*5*step +  5*step,
       col = rgb(.95,.95,.95), border=NA
  )
  
}
xs = rep(seq(step,4*step,by = step),ncol(co)) + rep(seq(0,1-step*5,step*5),each=nrow(co))
title(bquote(bold(.(Title))),line=.9)
if(k == 'allQ') axis(4,at = xs,labels=rep(c('ACRU','LITU','QUAL','QURU'),ncol(co)),tick=F,cex.axis=.55,hadj = .6)
if(k == 'allQ'){
  axis(1,at = c(round(minX,2),0,round(maxX,2)),tick=F,cex = .9);#axis(3,tick=F,cex = .9)
}else{
  axis(1,at = c(round(minX,1),0,round(maxX,1)),tick=F,cex = .9);#axis(3,tick=F,cex = .9)
}
axis(2,at = seq(5*step,33.5*step,by = 5*step),
     labels=c('Heated','Light',expression(bar('AT')),expression(paste('HeatedX',bar('AT'))),expression(paste(Delta,' AT')),expression(bar('SM'))),
     tick=F,cex.axis=.9,hadj = 0,padj=1.1,adj = 0,line=-.5)

abline(h = seq(0,1,5*step),lty=3)

abline(v=0,lty=2,col='grey')
pty1 = c(0,1,2,5)
pty2 = c(22,21,24,23)
for(i in 1:ncol(co)){
  if(i %% 2 == 0){b.col = 'white'}else{b.col =  rgb(.95,.95,.95)}
  for(j in 1:4){
    if((co.lo[j,i] < 0 & co.hi[j,i] > 0) | (co.lo[j,i] > 0 & co.hi[j,i] < 0 )){
      lines(c(co.lo[j,i],co.hi[j,i]),c(ys,ys),col= b.col,lwd = 8)
      lines(c(co.lo[j,i],co.hi[j,i]),c(ys,ys),col=j,lwd = .7)
      points(co[j,i],ys,col=j,cex=.75,pch=pty1[j])
    }else{
      lines(c(co.lo[j,i],co.hi[j,i]),c(ys,ys),col= b.col,lwd = 8)
      lines(c(co.lo[j,i],co.hi[j,i]),c(ys,ys),col=j,lwd = 1.5)
      points(co[j,i],ys,col=j,cex=.75,pch=pty2[j],bg=j)
    }
    ys = ys + step  
  }
  ys = ys + step  
}
}

all = all.lo = all.hi = numeric()
for(i in 1:4){
  all = rbind(all,allA[i,-c(1,2)],allR[i,-c(1,2)],allQ[i,-c(1,2)])
  all.lo = rbind(all.lo,allA.lo[i,-c(1,2)],allR.lo[i,-c(1,2)],allQ.lo[i,-c(1,2)])
  all.hi = rbind(all.hi,allA.hi[i,-c(1,2)],allR.hi[i,-c(1,2)],allQ.hi[i,-c(1,2)])
}

radarPlot(all,all.lo,all.hi,scale.within = T,legend = F,color='')


library(DT)
allTab = read.csv(file =  paste(path,'warming/allTab.csv',sep=''),stringsAsFactors=F)
datatable(allTab,rownames = FALSE, 
          options = list(pageLength = 8, dom = 'tip'),
          caption = htmltools::tags$caption(
            style = 'caption-side: top; text-align: left;',
            'Table 1: ', htmltools::em('Parameters of the light response model.'))
)

#print by treatments
load(file=paste(path,'warming/qualMVNassimilation.Rdata',sep=''))
species = sort(species)
npred = length(days)
days = days-(3*365)
yr = days
yr[days <= 366] = 1
yr[days > 366] = 2
#days[days > 366] = days[days > 366] - 366
on = days
on[days < 91] = 0  
on[days > 287] = 0
on[on !=0] = 1

trts = treatment(pheno$Plot,combine=T,light = TRUE,combineAandC = FALSE)
s1 = tapply(pheno[,'SP12'],paste0(pheno[,'Spec'],trts),mean) - 3*365 - 1
e1 = tapply(pheno[,'FA12'],paste0(pheno[,'Spec'],trts),mean) - 3*365 - 1
s2 = tapply(pheno[,'SP13'],paste0(pheno[,'Spec'],trts),mean) - 4*365 - 1
e2 = tapply(pheno[,'FA13'],paste0(pheno[,'Spec'],trts),mean) - 4*365 - 1

jdR = julian.date(c("01/01/2012","12/31/2013"),origin = c(month = 12, day = 31, year = 2008)) 
whr = which(dfTmp[,'JD2009'] >= jdR[1] & dfTmp[,'JD2009'] <= jdR[2] )
days = dfTmp[whr,'JD2009']

allAssim = totAssim = allAssSd = totAssSd = numeric()
for( ii in 1:length(species)){
  iii = ii
  load(file=paste(path,'warming/',sort(species)[ii],'MVNassimilation.Rdata',sep = ""))
  ii = iii
  species = sort(species)
  
  #on = rep(0,len(days))
  #on[(s1[] <= on & on <= e1[]) | (s2[] <= on & on <= e2[])] = 1
  #on[days < 91] = 0  
  #on[days > 287] = 0
  #on[on !=0] = 1
  
  allAssim = rbind(allAssim, rbind(allPmean[1:npred]*on ,
                                   allPmean[(npred+1):(2*npred)]*on ,
                                   allPmean[(2*npred+1):(3*npred)]*on,
                                   allPmean[(3*npred+1):(4*npred)]*on,
                                   allPmean[(4*npred+1):(5*npred)]*on,
                                   allPmean[(5*npred+1):(6*npred)]*on ))
  allAssSd = rbind(allAssSd, rbind(allPsd[1:npred]*on ,
                                   allPsd[(npred+1):(2*npred)]*on ,
                                   allPsd[(2*npred+1):(3*npred)]*on ,
                                   allPsd[(3*npred+1):(4*npred)]*on,
                                   allPsd[(4*npred+1):(5*npred)]*on,
                                   allPsd[(5*npred+1):(6*npred)]*on ))
  
  cs = allTotmean[1:npred] * on
  cs[cs > 0] = cs[cs > 0] - min(cs[cs > 0])
  cs[yr == 2 & cs > 0] = cs[yr == 2 & cs > 0] - min(cs[yr == 2 & cs > 0])
  cg = allTotmean[(npred+1):(2*npred)] * on
  cg[cg > 0] = cg[cg > 0] - min(cg[cg > 0])
  cg[yr == 2 & cg > 0] = cg[yr == 2 & cg > 0] - min(cg[yr == 2 & cg > 0])
  zs =  allTotmean[(2*npred+1):(3*npred)] * on 
  zs[zs > 0] = zs[zs > 0] - min(zs[zs > 0])
  zs[yr == 2 & zs > 0] = zs[yr == 2 & zs > 0] - min(zs[yr == 2 & zs > 0])
  zg =  allTotmean[(3*npred+1):(4*npred)] * on 
  zg[zg > 0] = zg[zg > 0] - min(zg[zg > 0])
  zg[yr == 2 & zg > 0] = zg[yr == 2 & zg > 0] - min(zg[yr == 2 & zg > 0])
  hs = allTotmean[(4*npred+1):(5*npred)] * on
  hs[hs > 0] = hs[hs > 0] - min(hs[hs > 0])
  hs[yr == 2 & hs > 0] = hs[yr == 2 & hs > 0] - min(hs[yr == 2 & hs > 0])
  hg  = allTotmean[(5*npred+1):(6*npred)] * on 
  hg[hg > 0] = hg[hg > 0] - min(hg[hg > 0])
  hg[yr == 2 & hg > 0] = hg[yr == 2 & hg > 0] - min(hg[yr == 2 & hg > 0])
     
  totAssim = rbind(totAssim, rbind(cs*on ,cg*on ,zs*on, zg*on,hs*on,hg*on ))
  
  
  cs = allTotsd[1:npred] * on
  cs[cs > 0] = cs[cs > 0] - min(cs[cs > 0])
  cs[yr == 2 & cs > 0] = cs[yr == 2 & cs > 0] - min(cs[yr == 2 & cs > 0])
  cg = allTotsd[(npred+1):(2*npred)] * on
  cg[cg > 0] = cg[cg > 0] - min(cg[cg > 0])
  cg[yr == 2 & cg > 0] = cg[yr == 2 & cg > 0] - min(cg[yr == 2 & cg > 0])
  zs =  allTotsd[(2*npred+1):(3*npred)] * on 
  zs[zs > 0] = zs[zs > 0] - min(zs[zs > 0])
  zs[yr == 2 & zs > 0] = zs[yr == 2 & zs > 0] - min(zs[yr == 2 & zs > 0])
  zg =  allTotsd[(3*npred+1):(4*npred)] * on 
  zg[zg > 0] = zg[zg > 0] - min(zg[zg > 0])
  zg[yr == 2 & zg > 0] = zg[yr == 2 & zg > 0] - min(zg[yr == 2 & zg > 0])
  hs = allTotsd[(4*npred+1):(5*npred)] * on
  hs[hs > 0] = hs[hs > 0] - min(hs[hs > 0])
  hs[yr == 2 & hs > 0] = hs[yr == 2 & hs > 0] - min(hs[yr == 2 & hs > 0])
  hg  = allTotsd[(5*npred+1):(6*npred)] * on 
  hg[hg > 0] = hg[hg > 0] - min(hg[hg > 0])
  hg[yr == 2 & hg > 0] = hg[yr == 2 & hg > 0] - min(hg[yr == 2 & hg > 0])
  
  totAssSd = rbind(totAssSd, rbind(cs*on ,cg*on ,zs*on, zg*on,hs*on,hg*on ))

}
whr2 = apply(matrix(times,ncol=1),1,function(x) nearWHICH(dfTemp[,'times'],x))
whr = which(dfTemp[,'JD2009'] >= jdR[1] & dfTemp[,'JD2009'] <= jdR[2] )
colnames(allAssim) = colnames(totAssim) = colnames(allAssSd) = colnames(totAssSd) = dfTmp[whr,'times']
rownames(allAssim) = rownames(totAssim) = rownames(allAssSd) = rownames(totAssSd) = c('AcruShadeControl','AcruGapControl','AcruShadeZero','AcruGapZero','AcruShadeHot','AcruGapHot',
                                            'LituShadeControl','LituGapControl','LituShadeZero','LituGapZero','LituShadeHot','LituGapHot',
                                            'QualShadeControl','QualGapControl','QualShadeZero','QualGapZero','QualShadeHot','QualGapHot',
                                            'QuruShadeControl','QuruGapControl','QuruShadeZero','QuruGapZero','QuruShadeHot','QuruGapHot')
save(allAssim,totAssim,allAssSd,totAssSd,file = paste(path,'warming/Assim.RData',sep=''))

#######
load(file=paste(path,'warming/qualMVNassimilation.Rdata',sep=''))
pheno = read.csv(file=paste(path,'DissertationEssentials/Biomass/BiomassPhen.csv',sep=''),stringsAsFactors = F)
species = sort(species)


allAssim = totAssim = allAssSd = totAssSd = matrix(NA,nrow=0,ncol=npred)
for( ii in 1:nrow(pheno)){
  iii = ii
  load(file=paste(path,'warming/',pheno[ii,'Spec'],'MVNassimilation.Rdata',sep = ""))
  ii = iii
  npred = length(days)
  days = days-(2*365)
  yr = days
  yr[days <= (366+365)] = 1
  yr[days > 366+365] = 2
  #on = rep(1,len(days))
  on = days
  on[days < pheno[ii,'SP10']] = 0  
  on[days > pheno[ii,'FA10'] & days < pheno[ii,'SP11']] = 0  
  on[days > pheno[ii,'FA11']] = 0
  on[on !=0] = 1
  trt = treatment(pheno[ii,'Plot'],combine =FALSE,light = TRUE)
  if(trt == 'Shadecontrol'){
    allAssim = rbind(allAssim,allPmean[1:npred]*on)
    allAssSd = rbind(allAssSd, allPsd[1:npred]*on)
    cs = allTotmean[1:npred] * on
    cs[cs > 0] = cs[cs > 0] - min(cs[cs > 0])
    cs[yr == 2 & cs > 0] = cs[yr == 2 & cs > 0] - min(cs[yr == 2 & cs > 0])
    totAssim = rbind(totAssim,cs*on)
    cs = allTotsd[1:npred] * on
    cs[cs > 0] = cs[cs > 0] - min(cs[cs > 0])
    cs[yr == 2 & cs > 0] = cs[yr == 2 & cs > 0] - min(cs[yr == 2 & cs > 0])
    totAssSd = rbind(totAssSd,cs*on) 
  }
  if(trt == 'Gapcontrol'){
    allAssim = rbind(allAssim,allPmean[(npred+1):(2*npred)]*on)
    allAssSd = rbind(allAssSd, allPsd[(npred+1):(2*npred)]*on)
    cg = allTotmean[(npred+1):(2*npred)] * on
    cg[cg > 0] = cg[cg > 0] - min(cg[cg > 0])
    cg[yr == 2 & cg > 0] = cg[yr == 2 & cg > 0] - min(cg[yr == 2 & cg > 0])
    totAssim = rbind(totAssim,cg*on)
    cg = allTotsd[(npred+1):(2*npred)] * on
    cg[cg > 0] = cg[cg > 0] - min(cg[cg > 0])
    cg[yr == 2 & cg > 0] = cg[yr == 2 & cg > 0] - min(cg[yr == 2 & cg > 0])
    totAssSd = rbind(totAssSd,cg*on) 
  }
  if(trt == 'Shadezero'){
    allAssim = rbind(allAssim,allPmean[(2*npred+1):(3*npred)]*on)
    allAssSd = rbind(allAssSd, allPsd[(2*npred+1):(3*npred)]*on)
    zs =  allTotmean[(2*npred+1):(3*npred)] * on 
    zs[zs > 0] = zs[zs > 0] - min(zs[zs > 0])
    zs[yr == 2 & zs > 0] = zs[yr == 2 & zs > 0] - min(zs[yr == 2 & zs > 0])
    totAssim = rbind(totAssim,zs*on)
    zs =  allTotsd[(2*npred+1):(3*npred)] * on 
    zs[zs > 0] = zs[zs > 0] - min(zs[zs > 0])
    zs[yr == 2 & zs > 0] = zs[yr == 2 & zs > 0] - min(zs[yr == 2 & zs > 0])
    totAssSd = rbind(totAssSd,zs*on) 
  }
  if(trt == 'Gapzero'){
    allAssim = rbind(allAssim,allPmean[(3*npred+1):(4*npred)]*on)
    allAssSd = rbind(allAssSd, allPsd[(3*npred+1):(4*npred)]*on)
    zg =  allTotmean[(3*npred+1):(4*npred)] * on 
    zg[zg > 0] = zg[zg > 0] - min(zg[zg > 0])
    zg[yr == 2 & zg > 0] = zg[yr == 2 & zg > 0] - min(zg[yr == 2 & zg > 0])
    totAssim = rbind(totAssim,zg*on)
    zg =  allTotsd[(3*npred+1):(4*npred)] * on 
    zg[zg > 0] = zg[zg > 0] - min(zg[zg > 0])
    zg[yr == 2 & zg > 0] = zg[yr == 2 & zg > 0] - min(zg[yr == 2 & zg > 0])
    totAssSd = rbind(totAssSd,zg*on) 
  }
  if(trt == 'Shadethree' | trt == 'Shadefive'){
    allAssim = rbind(allAssim,allPmean[(4*npred+1):(5*npred)]*on)
    allAssSd = rbind(allAssSd, allPsd[(4*npred+1):(5*npred)]*on)
    hs = allTotmean[(4*npred+1):(5*npred)] * on
    hs[hs > 0] = hs[hs > 0] - min(hs[hs > 0])
    hs[yr == 2 & hs > 0] = hs[yr == 2 & hs > 0] - min(hs[yr == 2 & hs > 0])
    totAssim = rbind(totAssim,hs*on)
    hs = allTotsd[(4*npred+1):(5*npred)] * on
    hs[hs > 0] = hs[hs > 0] - min(hs[hs > 0])
    hs[yr == 2 & hs > 0] = hs[yr == 2 & hs > 0] - min(hs[yr == 2 & hs > 0])
    totAssSd = rbind(totAssSd,hs*on) 
  }
  if(trt == 'Gapthree' | trt == 'Gapfive'){
    allAssim = rbind(allAssim,allPmean[(5*npred+1):(6*npred)]*on )
    allAssSd = rbind(allAssSd, allPsd[(5*npred+1):(6*npred)]*on)
    hg  = allTotmean[(5*npred+1):(6*npred)] * on 
    hg[hg > 0] = hg[hg > 0] - min(hg[hg > 0])
    hg[yr == 2 & hg > 0] = hg[yr == 2 & hg > 0] - min(hg[yr == 2 & hg > 0])
    totAssim = rbind(totAssim,hg*on)
    hg  = allTotsd[(5*npred+1):(6*npred)] * on 
    hg[hg > 0] = hg[hg > 0] - min(hg[hg > 0])
    hg[yr == 2 & hg > 0] = hg[yr == 2 & hg > 0] - min(hg[yr == 2 & hg > 0])
    totAssSd = rbind(totAssSd,hg*on) 
  }
  print(ii)
}
colnames(allAssim) = colnames(totAssim) = colnames(allAssSd) = colnames(totAssSd) = dfTmp[whr,'times']
rownames(allAssim) = rownames(totAssim) = rownames(allAssSd) = rownames(totAssSd) = pheno[,'Indiv']
save(allAssim,totAssim,allAssSd,totAssSd,file = paste(path,'warming/AssimIndiv.RData',sep=''))

load(file = paste(path,'warming/AssimIndiv.RData',sep=''))

##
#Scale All Ass to total leaf area
stom = function(spec,A){
  whr = match(spec,c("acru","litu","qual","quru"))
  g = A
  g[A < 0]  =  coefficients(gStom.fit[[whr]])[2] * A[A < 0] + predict(gStom.fit[[whr+4]],0)$y
  g[A >= 0] =  predict(gStom.fit[[whr+4]],A[A >= 0])$y
  g
}


gStom.fit = readRDS(file = paste(path,'gStomFit.rds',sep=''))


allCond = allAssim
for( i in 1:dim(allAssim)[1]){
  allCond[i,] = stom(pheno[i,'Spec'],allAssim[i,])
  print(i)
}

saveRDS(allCond,file='warming/ENVDATA/processedData/Conduct.rds')
#allCond mol h20 per m2 per s
#allCond mol h20 per m2 per s * leaf Area cm2 / 10000 m2 per cm2
#allCond mol h20 per m2 per s * leaf Area cm2 / 10000 m2 per cm2 * 60 sec per min * 60 min per hour
#allCond mol h20 per m2 per s * leaf Area cm2 / 10000 m2 per cm2 * 60 sec per min * 60 min per hour * 18.01528 g / mol

par(mfrow=c(2,1),mar=c(3,2,0,0))


whr = which(days == 536 | days == 537)
pp = allCond*(pheno[,'T.Area.cm.2.']/10000)*18.01528*60*60 #grams per plant per hour
plot(pheno[,'T.Area.cm.2.']*.01,apply(pp[,whr],1,sum))
points(c(30,45,60,75),c(190,250,325,400),col='blue')

allAssLA = allAssim*(pheno[,'T.Area.cm.2.']/10000)  #per plant 
allConLA = allCond*(pheno[,'T.Area.cm.2.']/10000)*18.015280 #grams per plant 
 #mol h20 per m2 per s

plot(allAssim[120,]) #per second
plot(allCond[13,]*18.01528)

plot(allAssLA[1,]) #per hour
plot(allConLA[4,])

plot(allAssLA[1,]*3600) #per hour
plot(allConLA[1,]*3600)

#days = JD
DayAssLA = t(tapply.matrix(allAssLA*3600,1,days,sum) )
DayConLA = t(tapply.matrix(allConLA*3600,1,days,sum) )

barplot(DayAssLA[2,])
barplot(DayConLA[2,])

#0.0555084350618  mole per gram
#18.01528 gram per mole


##
load(file = paste(path,'warming/Assim.RData',sep=''))
par(bg='white')
totAssim[totAssim == 0] = NA
top = totAssim + (.4 + totAssSd)*10
bot = totAssim - (.4 + totAssSd)*10
top[which(is.na(totAssim))] = 0
bot[which(is.na(totAssim))] = 0
whrC = grep('GapZero',rownames(totAssim))
whrH = grep('GapHot',rownames(totAssim))
plot(totAssim[whrC[1],]/10000,
     ylim=c(-0,1.5),type='n')
for(i in 1:len(whrC)){
  bCol = c("#00A0B080","#6A4A3C80","#CC333F80","#EB684180","#EDC95180")
  polygon(c(1:ncol(bot),(ncol(bot):1)),c(top[whrC[i],]/10000,rev(bot[whrC[i],]/10000)),col=bCol[i],border=NA)
  polygon(c(1:ncol(bot),(ncol(bot):1)),c(top[whrH[i],]/10000,rev(bot[whrH[i],]/10000)),col=bCol[i],border=NA)
}

for(i in 1:len(whrC)){
  lines(totAssim[whrC[i],]/10000,col=i,lty=1,lwd = 3)
  lines(totAssim[whrH[i],]/10000,col=i,lty=2, lwd= 3)
}
legend(0,1.5,species,1:4)

par(bg='white')
totAssim[totAssim == 0] = NA
Msd = apply(totAssSd,1,mean)
top = totAssim + (.4 + totAssSd)*10
bot = totAssim - (.4 + totAssSd)*10
top[which(is.na(totAssim))] = 0
bot[which(is.na(totAssim))] = 0
whrC = grep('ShadeZero',rownames(totAssim))
whrH = grep('ShadeHot',rownames(totAssim))
plot(totAssim[whrC[1],]/10000,
     ylim=c(-0,.4),type='n')
for(i in 1:len(whrC)){
  bCol = c("#00A0B080","#6A4A3C80","#CC333F80","#EB684180","#EDC95180")
  polygon(c(1:ncol(bot),(ncol(bot):1)),c(top[whrC[i],]/10000,rev(bot[whrC[i],]/10000)),col=bCol[i],border=NA)
  polygon(c(1:ncol(bot),(ncol(bot):1)),c(top[whrH[i],]/10000,rev(bot[whrH[i],]/10000)),col=bCol[i],border=NA)
}

for(i in 1:len(whrC)){
  lines(totAssim[whrC[i],]/10000,col=i,lty=1,lwd = 3)
  lines(totAssim[whrH[i],]/10000,col=i,lty=2, lwd= 3)
}
legend(0,.4,species,1:4)

lines(top[whrC[i],]/10000,)




par(bg='white')
totAssim[totAssim == 0] = NA
plot(totAssim[10,]*pheno[10,'T.Area.cm.2.']/10000,
     ylim=c(-0,10000),type='n')
for(i in 1:nrow(pheno)){
lines(totAssim[i,]*pheno[i,'T.Area.cm.2.']/10000,col=match(pheno[i,'Spec'],species))
}
legend(0,10000,species,1:4)


totAssim[totAssim == 0] = NA
plot(totAssim[10,]*pheno[10,'T.Area.cm.2.']/10000/pheno[10,'MassG'],
    type='n')
for(i in 1:nrow(pheno)){
  lines(totAssim[i,]*pheno[i,'T.Area.cm.2.']/10000/pheno[i,'MassG'],col=match(pheno[i,'Spec'],species))
}
legend(0,50000,species,1:4)

par(mfrow=c(2,2))
totAssim[totAssim == 0] = NA
plot(totAssim[10,]*pheno[10,'T.Area.cm.2.']/10000,
     type='n',ylab = 'Gap per Plant Assimilation', ylim=c(-0,10000))
whr = grep('g',pheno[,'Plot'])
for(i in whr){
  lines(totAssim[i,]*pheno[i,'T.Area.cm.2.']/10000,col=match(pheno[i,'Spec'],species))
}
legend(0,10000,species,1:4,text.col='black')


totAssim[totAssim == 0] = NA
plot(totAssim[10,]*pheno[10,'T.Area.cm.2.']/10000,
     type='n',ylab = 'Shade per Plant Assimilation', ylim=c(-0,400))
whr = grep('s',pheno[,'Plot'])
for(i in whr){
  lines(totAssim[i,]*pheno[i,'T.Area.cm.2.']/10000,col=match(pheno[i,'Spec'],species))
}
legend(0,400,species,1:4)



totAssim[totAssim == 0] = NA
plot(totAssim[10,]*pheno[10,'T.Area.cm.2.']/10000/pheno[10,'MassG'],
    type='n',ylab = 'Gap per Plant/Mass Assimilation', ylim=c(-0,3500))
whr = grep('g',pheno[,'Plot'])
for(i in whr){
  lines(totAssim[i,]*pheno[i,'T.Area.cm.2.']/10000/pheno[i,'MassG'],col=match(pheno[i,'Spec'],species))
}
legend(0,3500,species,1:4)


totAssim[totAssim == 0] = NA
plot(totAssim[10,]*pheno[10,'T.Area.cm.2.']/10000/pheno[10,'MassG'],
     type='n',ylab = 'Shade per Plant/Mass Assimilation', ylim=c(-0,7000))
whr = grep('s',pheno[,'Plot'])
for(i in whr){
  lines(totAssim[i,]*pheno[i,'T.Area.cm.2.']/10000/pheno[i,'MassG'],col=match(pheno[i,'Spec'],species))
}
legend(0,7000,species,1:4)

tot1 = apply(allAssim[,1:(ncol(allAssim)/2)]/10000 ,1,sum,na.rm=T)
tot2 = apply(allAssim[,(ncol(allAssim)/2):ncol(allAssim)]/10000 ,1,sum,na.rm=T)

totper1 = apply(allAssim[,1:(ncol(allAssim)/2)]*pheno[,'T.Area.cm.2.']/10000 ,1,sum,na.rm=T)
totper2 = apply(allAssim[,(ncol(allAssim)/2):ncol(allAssim)]*pheno[,'T.Area.cm.2.']/10000 ,1,sum,na.rm=T)
par(mfrow=c(2,2),mar=c(3,3,2,0))
whr = grep('g',pheno[,'Plot'])
boxplot(tot1[whr]~ pheno[whr,'Spec'] + pheno[whr,'Temp'],col =  1:4)
boxplot(tot2[whr]~ pheno[whr,'Spec'] + pheno[whr,'Temp'],outline=F ,col =  1:4)

whr = grep('s',pheno[,'Plot'])
boxplot(tot1[whr]~ pheno[whr,'Spec'] + pheno[whr,'Temp'] )
boxplot(tot2[whr]~ pheno[whr,'Spec'] + pheno[whr,'Temp'],outline=F  ,col =  1:4)

whr = grep('g',pheno[,'Plot'])
boxplot(totper1[whr]~ pheno[whr,'Spec'] + pheno[whr,'Temp'] )
boxplot(totper2[whr]~ pheno[whr,'Spec'] + pheno[whr,'Temp'],outline=F ,col =  1:4 )

whr = grep('s',pheno[,'Plot'])
boxplot(totper1[whr]~ pheno[whr,'Spec'] + pheno[whr,'Temp'] )
boxplot(totper2[whr]~ pheno[whr,'Spec'] + pheno[whr,'Temp'] ,outline=F ,col =  1:4)

whr = grep('s',pheno[,'Plot'])
boxplot(totper1[whr]~ pheno[whr,'Spec'] + pheno[whr,'Temp'] )
boxplot(totper2[whr]~ pheno[whr,'Spec'] + pheno[whr,'Temp'] ,outline=F ,col =  1:4)

par(mfrow=c(1,1))
boxplot(c(totper1[whr],totper2[whr])~ c(pheno[whr,'Temp'],pheno[whr,'Temp']) + c(pheno[whr,'Spec'],pheno[whr,'Spec'])  ,outline=F ,col =  rep(1:4,each=2),
        ylab="Total assimialtion per plant",names = c("acru","heated", "litu","heated", "qual","heated", "quru","heated"))

whr = grep('g',pheno[,'Plot'])
boxplot( pheno[whr,'FA11'] - pheno[whr,'SP11']~ pheno[whr,'Spec'] + pheno[whr,'Temp'],outline=F  ,col =  1:3,horizontal =T)
whr = grep('s',pheno[,'Plot'])
boxplot( pheno[whr,'FA11'] - pheno[whr,'SP11']~ pheno[whr,'Spec'] + pheno[whr,'Temp'],outline=F  ,col =  1:4,horizontal =T)
whr = grep('g',pheno[,'Plot'])
boxplot(pheno[whr,'SP11']~ pheno[whr,'Spec'] + pheno[whr,'Temp'],outline=F  ,col =  1:3,horizontal =T,xlim=c(790,1090))
boxplot(pheno[whr,'FA11']~ pheno[whr,'Spec'] + pheno[whr,'Temp'],outline=F  ,col =  1:3,horizontal =T)
whr = grep('s',pheno[,'Plot'])
boxplot( pheno[whr,'FA11'] - pheno[whr,'SP11']~ pheno[whr,'Spec'] + pheno[whr,'Temp'],outline=F  ,col =  1:4,horizontal =T)

stripchart(pheno[whr,'SP11']~ pheno[whr,'Spec'] + pheno[whr,'Temp'],outline=F  ,col =  1:3,horizontal =T,xlim=c(790,1090))
stripchart(pheno[whr,'FA11']~ pheno[whr,'Spec'] + pheno[whr,'Temp'],outline=F  ,col =  1:3,horizontal =T,add=T)


par(mfrow=c(1,1),mar=c(4,2,1,1),bg = 'white')
split.screen(c(1,2))
colors = c("#00A0B096","#00A0B0","#6A4A3C96","#6A4A3C","#CC333F96","#CC333F","#EB684196","#EB6841","#EDC951")
colors=c("black","grey","black","grey","black","grey","black","grey")
shape = c(15,16,17,18)
plot(0,0,ylim=c(0,1),xlim=c(790,1090),yaxt='n',xaxt='n',xlab='Day of Year',ylab='')
axis(1,seq(800,1100,50),seq(800,1100,50)-2*365)
whr = order(pheno[,'SP11'])
lwd = 1/nrow(pheno)
widths = table(paste(pheno[,'Spec'],pheno[,'Light'],sep=""))*lwd
st = c(0,cumsum(widths)[-7])
sp = sort(unique(paste(pheno[,'Spec'],pheno[,'Light'],sep="")))
lig = sort(unique(pheno[,'Light']))
numac = st[1];numacs=st[2];
numli= st[3];numlis=st[4];
numqa= st[5];numqas=st[6];
numqrs=st[7];
for(i in 1:nrow(pheno)){
  id = match(paste(pheno[i,'Spec'],pheno[i,'Light'],sep=""),sp)
  if(id == 1){num = numac}
  if(id == 2){num = numacs}
  if(id == 3){num = numli}
  if(id == 4){num = numlis}
  if(id == 5){num = numqa}
  if(id == 6){num = numqas}
  if(id == 7){num = numqrs}
  print(num)
  rect(pheno[whr[i],'FA11'],num,pheno[whr[i],'SP11'],num+lwd,col=colors[id],border="#F8F8FF")
  if(id == 1){numac = lwd + numac}
  if(id == 2){numacs = lwd + numacs}
  if(id == 3){numli = lwd + numli}
  if(id == 4){numlis = lwd + numlis}
  if(id == 5){numqa = lwd + numqa}
  if(id == 6){numqas = lwd + numqas}
  if(id == 7){numqrs = lwd + numqrs}
}
widths = table(paste(pheno[,'Spec'],pheno[,'Light'],sep=""))*lwd
st = c(0,cumsum(widths))
points(tapply(pheno[,'SP11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),median),((st[-1]+st[-8])/2), col= 'black',cex=2)
points(tapply(pheno[,'SP11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),quantile,.05),((st[-1]+st[-8])/2),pch = '|', col= 'black',cex=2)
points(tapply(pheno[,'SP11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),quantile,.95),((st[-1]+st[-8])/2),pch = '|', col= 'black',cex=2)
points(tapply(pheno[,'FA11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),median),((st[-1]+st[-8])/2), col= 'black',cex=2)
points(tapply(pheno[,'FA11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),quantile,.05),((st[-1]+st[-8])/2),pch = '|', col= 'black',cex=2)
points(tapply(pheno[,'FA11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),quantile,.95),((st[-1]+st[-8])/2),pch = '|', col= 'black',cex=2)
widths = table(pheno[,'Spec'])*lwd
st = c(0,cumsum(widths))
sp = sort(unique(pheno[,'Spec']))
text(x = rep(195+365*2,4),y = ((st[-1]+st[-5])/2),sp,col='black',cex=2)

#assim per meter squared correlation to season length
par(mar=c(0,1,4,2))
split.screen(c(2,1),2)
colors = c("#00A0B0","#6A4A3C","#CC333F","#EB6841","#EDC951")
screen(3)
whr = grep('g',pheno[,'Plot'])
plot(pheno[whr,'FA13'] - pheno[whr,'SP13'],tot2[whr],col = colors[match( pheno[whr,'Spec'],sp)],xlab = 'Growing Season Length',ylab = 'Assimilation' )
points(pheno[whr,'FA12'] - pheno[whr,'SP12'],tot1[whr],col = colors[match( pheno[whr,'Spec'],sp)])
cors = numeric()
for( i in 1:4){
  whr2 = which(pheno[whr,'Spec'] == sp[i])
  if(length(whr2) == 0)next
  fit = lm( c(tot2[whr][whr2],tot1[whr][whr2])~c(pheno[whr[whr2],'FA13'] - pheno[whr[whr2],'SP13'],pheno[whr[whr2],'FA12'] - pheno[whr[whr2],'SP12']))
  abline(coefficients(fit),col = colors[i])
  cors = c(cors,round(cor(c(tot2[whr][whr2],tot1[whr][whr2]),c(pheno[whr[whr2],'FA13'] - pheno[whr[whr2],'SP13'],pheno[whr[whr2],'FA12'] - pheno[whr[whr2],'SP12'])),2))  
}
legend(170,.22,paste(sp," (",cors,")",sep=""),fill = colors,title = 'Shade',text.col = 'black')
legend(190,2,paste(sp[1:3]," (",cors,")",sep=""),fill = colors,title = 'Gap',text.col = 'black')

#assim per plant correlation to season length
screen(4)
colors = c("#00A0B0","#6A4A3C","#CC333F","#EB6841","#EDC951")
colors = c("black","black","black","black","black")
colors = c("grey","grey","grey","grey","grey")
whr = grep('g',pheno[,'Plot'])
plot(pheno[whr,'FA13'] - pheno[whr,'SP13'],totper2[whr],col = colors[match( pheno[whr,'Spec'],sp)],xlab = '',ylab = 'Assimilation',
     ylim = c(0,300),pch = shape[match( pheno[whr,'Spec'],sp)],type='n',xaxt='n',yaxt='n') #5000, 300
grid(lty=3)
points(pheno[whr,'FA13'] - pheno[whr,'SP13'],totper2[whr],col = colors[match( pheno[whr,'Spec'],sp)],pch = shape[match( pheno[whr,'Spec'],sp)])
points(pheno[whr,'FA12'] - pheno[whr,'SP12'],totper1[whr],col = colors[match( pheno[whr,'Spec'],sp)],pch = shape[match( pheno[whr,'Spec'],sp)])
axis(1,tick = F)
axis(4,tick = F)

cors = numeric()
for( i in 1:4){
  whr2 = which(pheno[whr,'Spec'] == sp[i])
  if(length(whr2) == 0)next
  fit = lm( c(totper2[whr][whr2],totper1[whr][whr2])~c(pheno[whr[whr2],'FA13'] - pheno[whr[whr2],'SP13'],pheno[whr[whr2],'FA12'] - pheno[whr[whr2],'SP12']))
  abline(coefficients(fit),col = colors[i],lty=i)
  cors = c(cors,round(cor(c(totper2[whr][whr2],totper1[whr][whr2]),c(pheno[whr[whr2],'FA13'] - pheno[whr[whr2],'SP13'],pheno[whr[whr2],'FA12'] - pheno[whr[whr2],'SP12'])),2))  
}
legend(170,280,paste(sp," (",cors,")",sep=""),col=colors,title = 'Shade',text.col = 'black',pch=15:18)
legend(170,4000,paste(sp[1:3]," (",cors,")",sep=""),col=colors,title = 'Gap',text.col = 'black',pch=15:18)

close.screen(,all.screens=T)

library(dygraphs)
load(file = paste(path,'warming/Assim.RData',sep=''))
samp = floor(seq(1,length(totAssim[2,]),length.out = 1000))
acruGapC =    ts(totAssim[7,samp],frequency=1,start =  as.Date("2000-05-06"))
acruGapC.hi = ts(totAssim[7,samp] + totAssSd[1,samp],frequency=1,start =  as.Date("2000-05-06"))
acruGapC.lo = ts(totAssim[7,samp] - totAssSd[1,samp],frequency=1,start =  as.Date("2000-05-06"))
acruGapZ =    ts(totAssim[9,samp],frequency=1,start =  as.Date("2000-05-06"))
acruGapZ.hi = ts(totAssim[9,samp] + totAssSd[9,samp],frequency=1,start =  as.Date("2000-05-06"))
acruGapZ.lo = ts(totAssim[9,samp] - totAssSd[9,samp],frequency=1,start =  as.Date("2000-05-06"))
acruGapH =    ts(totAssim[11,samp],frequency=1,start =  as.Date("2000-05-06"))
acruGapH.hi = ts(totAssim[11],samp] + totAssSd[11,samp],frequency=1,start =  as.Date("2000-05-06"))
acruGapH.lo = ts(totAssim[11,samp] - totAssSd[11,samp],frequency=1,start =  as.Date("2000-05-06"))
pr = ts.union(acruGapC.lo,acruGapC,acruGapC.hi,
              acruGapZ.lo,acruGapZ,acruGapZ.hi,
              acruGapH.lo,acruGapH,acruGapH.hi)
dygraph(pr, main = "ACRU Gap Carbon Accumulation") %>%
  dyRangeSelector() %>%
  dySeries(c("acruGapC.lo","acruGapC","acruGapC.hi"),label = "Control") %>%
  dySeries(c("acruGapZ.lo","acruGapZ","acruGapZ.hi"),label = "Zero") %>%
  dySeries(c("acruGapH.lo","acruGapH","acruGapH.hi"),label = "Hot")


mCol = nColor(number = 4,base.col = '#AF2F03',trans = 1)
bCol = nColor(number = 4,base.col = '#AF2F03',trans = .25)
iii = 0
samp = floor(seq(1,length(totAssim[2,]),length.out = 1000))
plot(totAssim[7,samp],type='n',ylim=c(-40,4500),ylab= c('Shade Carbon Accumulation'),xlab='Time')
for(ii in c(3,9,15,21)){
  iii = iii +1
  lines(totAssim[ii,samp],col=mCol[iii],lty=2)
  lines(totAssim[ii+2,samp],col=mCol[iii])
  polygon(c(1:length(samp),length(samp):1),c(totAssim[ii,samp] + totAssSd[ii,samp],rev(totAssim[ii,samp] - totAssSd[ii,samp])),col=bCol[iii],border=NA)
  polygon(c(1:length(samp),length(samp):1),c(totAssim[ii+2,samp] + totAssSd[ii+2,samp],rev(totAssim[ii+2,samp] - totAssSd[ii+2,samp])),col=bCol[iii],border=NA)
}
legend(1,4600,fill = mCol, legend = c("ACRU","LITU","QUAL","QURU"),bty='n')

par(mfrow = c(1,1))
plot(allAssim[12,],ylim = c(-5,10),type='n')
for(mk in c(24,18,12,6)){
  lines(allAssim[mk,],col=mk);lines(allAssim[mk,]+allAssSd[mk,],col=mk,lty=3);lines(allAssim[mk,]-allAssSd[mk,],col=mk,lty=3)
}

par(mfrow = c(1,1))
plot(totAssim[1,],type='n',ylim=c(0,12000))
for(mk in c(13,14,15,16)){
  lines(totAssim[mk,],col=mk);lines(totAssim[mk,]+totAssSd[mk,],col=mk,lty=3);lines(totAssim[mk,]-totAssSd[mk,],col=mk,lty=3)
}

par(mfrow = c(1,1))
plot(totAssim[2,],type='l',ylim=c(0,10000))
lines(totAssim[8,],col=1)
lines(totAssim[14,],col=3)
lines(totAssim[20,],col=4)
lines(totAssim[5,],col=5)
lines(totAssim[6,],col=6)



radarPlot(allA[,-1],allA.lo[,-1],allA.hi[,-1],scale.within = F,legend = F)
radarPlot(allA[,-1],allA.lo[,-1],allA.hi[,-1],scale.within = T,legend = F)
radarPlot(allR[,-1],allR.lo[,-1],allR.hi[,-1],scale.within = F,legend = F)
radarPlot(allR[,-1],allR.lo[,-1],allR.hi[,-1],scale.within = T,legend = F)
radarPlot(allQ[,-1],allQ.lo[,-1],allQ.hi[,-1],scale.within = F,legend = F)
radarPlot(allQ[,-1],allQ.lo[,-1],allQ.hi[,-1],scale.within = T,legend = F)
radarPlot(allTh,allTh.lo,allTh.hi,scale.within = F,legend = F)
radarPlot(allTh,allTh.lo,allTh.hi,scale.within = T,legend = F)


allPoints = read.csv(file = paste(path,'warming/Li6400/licorPROCESSED/all.rates.points.csv',sep=''),stringsAsFactors = F)
colors = c("#00A0B0","#6A4A3C","#CC333F","#EB6841","#EDC951")
colorsT = c("#00A0B080","#6A4A3C80","#CC333F80","#EB684180","#EDC95180")
plot(allPoints[,'PARi'],allPoints[,'Photo'])

par(mfrow=c(2,2),mar=c(0,4,4,0),asp=1)
whr = which(allPoints[,'species'] == 'acru' & allPoints[,'Light'] == 'Gap')
plot(allPoints[whr,'PARi'],allPoints[whr,'Photo'],pch=20,xlim=c(0,2000),ylim=c(-5,20),col=colors[1],cex=.75 ,xaxt='n',ylab='Photosynthesis')
whr = which(allPoints[,'species'] == 'acru' & allPoints[,'Light'] == 'Shade' )
points(allPoints[whr,'PARi'],allPoints[whr,'Photo'],col=rgb(.2,.2,.2,.5),cex=.6,pch=3)
title(main='acru',line=-1)
abline(h = 0,lty=3)

par(mar=c(0,0,4,4))
whr = which(allPoints[,'species'] == 'litu' & allPoints[,'Light'] == 'Gap' )
plot(allPoints[whr,'PARi'],allPoints[whr,'Photo'],pch=20,xlim=c(0,2000),ylim=c(-5,20),col=colors[2],cex=.75 ,xaxt='n',yaxt='n')
whr = which(allPoints[,'species'] == 'litu' & allPoints[,'Light'] == 'Shade' )
points(allPoints[whr,'PARi'],allPoints[whr,'Photo'],col=rgb(.2,.2,.2,.5),cex=.6,pch=3)
abline(h = 0,lty=3)
title(main='litu',line=-1)
axis(4)

par(mar=c(4,4,0,0))
whr = which(allPoints[,'species'] == 'qual' & allPoints[,'Light'] == 'Gap' )
plot(allPoints[whr,'PARi'],allPoints[whr,'Photo'],pch=20,xlim=c(0,2000),ylim=c(-5,20),col=colors[3],cex=.75 ,ylab='',xlab='PPFD')
whr = which(allPoints[,'species'] == 'qual' & allPoints[,'Light'] == 'Shade' )
points(allPoints[whr,'PARi'],allPoints[whr,'Photo'],col=rgb(.2,.2,.2,.5),cex=.6,pch=3)
title(main='qual',line=-1)
abline(h = 0,lty=3)

par(mar=c(4,0,0,4))
whr = which(allPoints[,'species'] == 'quru' & allPoints[,'Light'] == 'Gap' )
plot(allPoints[whr,'PARi'],allPoints[whr,'Photo'],pch=20,xlim=c(0,2000),ylim=c(-5,20),col=colors[4],cex=.75, yaxt='n',xlab='')
whr = which(allPoints[,'species'] == 'quru' & allPoints[,'Light'] == 'Shade' )
points(allPoints[whr,'PARi'],allPoints[whr,'Photo'],col=rgb(.2,.2,.2,.5),cex=.6,pch=3)
abline(h = 0,lty=3)
title(main='quru',line=-1)
axis(4)





###
#
###
#3pplots
###
#
###

#firstplot species data
#second plot individual response
#third plot growth respons
load(file = paste(path,'warming/Assim.RData',sep=''))

coldAss = matrix(NA,4,10000)
hotAss = matrix(NA,4,10000)
for(k in 1:4){
  for(j in 1:10000){
    tmp = rnorm(len((ncol(allAssim)/2):ncol(allAssim)),allAssim[whrC[k],(ncol(allAssim)/2):ncol(allAssim)]*3600/1000000,allAssSd[whrC[k],(ncol(allAssim)/2):ncol(allAssim)]*3600/1000000)
    tmp[allAssim[whrC[k],(ncol(allAssim)/2):ncol(allAssim)]*3600/1000000 == 0] = 0
    coldAss[k,j] <- sum(tmp)
    tmp = rnorm(len((ncol(allAssim)/2):ncol(allAssim)),allAssim[whrH[k],(ncol(allAssim)/2):ncol(allAssim)]*3600/1000000,allAssSd[whrH[k],(ncol(allAssim)/2):ncol(allAssim)]*3600/1000000)
    tmp[allAssim[whrH[k],(ncol(allAssim)/2):ncol(allAssim)]*3600/1000000 == 0] = 0
    hotAss[k,j] <- sum(tmp)
  }
  print(k)
}
cold = apply(coldAss,1,quantile,c(.025,.5,.975))
hot = apply(hotAss,1,quantile,c(.025,.5,.975))

load(file = paste(path,'warming/AssimIndiv.RData',sep=''))
par(mfrow = c(2,1),mar=c(2,4,1,3))
par(mar=c(3, 4, 3, 2), oma=c(0,0,0,0), bg="#F0F0F0", xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
library(lubridate)
dts = as.Date(days, origin = mdy("12/31/2008"))
mns = months(dts,abbreviate=T)
par(bg='white')
totAssim[totAssim == 0] = NA
Msd = apply(totAssSd,1,mean)
top = totAssim + (.4 + totAssSd)*10
bot = totAssim - (.4 + totAssSd)*10
top[which(is.na(totAssim))] = 0
bot[which(is.na(totAssim))] = 0
whrC = grep('ShadeZero',rownames(totAssim))
whrH = grep('ShadeHot',rownames(totAssim))
whrYr2 = (ncol(allAssim)/2):ncol(allAssim)
labs = seq(1,len(whrYr2),length=12)
labs2 = seq(1,len(whrYr2),length=6)[c(-1,-6)]
rects = diff(labs2)[1]/2
plot(totAssim[whrC[1],whrYr2]*3600/1000000,ylim=c(-1,17),xlim=c(1,(2*len(whrYr2))),type='n',xaxt='n',xlab='',
     ylab=expression("Total Assimiation mol CO"[2] * " m"^-2))
#rect(-1000,0,(2*len(whrYr2))+1000,5,col=rgb(.2,.2,.2,.2),border=NA)
#rect(-1000,10,(2*len(whrYr2))+1000,15,col=rgb(.2,.2,.2,.2),border=NA)
axis(1,at = labs,labels = mns[labs])
axis(1,at = labs2+len(whrYr2),labels = c('ACRU','LITU','QUAL','QURU'))
axis(4)
abline(h=0,col='black')
abline(h=c(5,10,15),lty=3,col=rgb(.5,.5,.5))

#for(i in 1:len(whrC)){
#  bCol = c("#00A0B080","#6A4A3C80","#CC333F80","#EB684180","#EDC95180")
#  polygon(c(1:ncol(bot),(ncol(bot):1)),c(top[whrC[i],]/10000,rev(bot[whrC[i],]/10000)),col=bCol[i],border=NA)
#  polygon(c(1:ncol(bot),(ncol(bot):1)),c(top[whrH[i],]/10000,rev(bot[whrH[i],]/10000)),col=bCol[i],border=NA)
#}

for(i in 1:len(whrC)){
  lines(totAssim[whrC[i],(ncol(allAssim)/2):ncol(allAssim)]*3600/1000000,col=i,lty=1,lwd = 3)
    lines(c(6874,(labs2+len(whrYr2))[i]-50),c(totAssim[whrC[i],15617]*3600/1000000,totAssim[whrC[i],15617]*3600/1000000),col=thecol[i],lwd=2)
  lines(totAssim[whrH[i],(ncol(allAssim)/2):ncol(allAssim)]*3600/1000000,col=i,lty=2, lwd= 3)
    lines(c(6874,(labs2+len(whrYr2))[i] +rects - 50),c(totAssim[whrH[i],15617]*3600/1000000,totAssim[whrH[i],15617]*3600/1000000),col=thecol[i],lty=2,lwd=2)
}  
for(i in 1:len(whrC)){
  rect((labs2+len(whrYr2))[i]- rects +50 ,cold[1,i],
       -50+(labs2+len(whrYr2))[i]        ,cold[3,i],col = thecol[i],lty=1,lwd=3,border=i)
  rect((labs2+len(whrYr2))[i]+50         ,hot[1,i],
       (labs2+len(whrYr2))[i] +rects -50 ,hot[3,i], col = thecol[i],lty=2,lwd=3,border=i)
}
text(6874,totAssim[whrH[4],15617]*3600/1000000,adj=c(0,-.1),labels='Heated',col=4)
text(6874,totAssim[whrC[4],15617]*3600/1000000,adj=c(0,-.1),labels='Ambient',col=4)
text(6874,totAssim[whrH[1],15617]*3600/1000000,adj=c(0,-.1),labels='Heated',col=1)
text(6874,totAssim[whrC[1],15617]*3600/1000000,adj=c(0,-.1),labels='Ambient',col=1)

legend(0,15,species,1:4,text.col='black')
#legend((labs2+len(whrYr2)-1.4*rects)[4],4,c('Ambient','Heated'),col='black',lty=c(1,2),text.col='black')


for(i in 1:4){
lines(c(7155+i*200,7155+i*200),c(bot[whrC[i],15617]*3600/1000000,top[whrC[i],15617]*3600/1000000),col=i,lwd=4)
lines(c(7155+100+i*200,7155+100+i*200),c(bot[whrH[i],15617]*3600/1000000,top[whrH[i],15617]*3600/1000000),col=i,lwd=4,lty=1)
}

bs = apply(Beta[burnin:nsim,],2,mean)
bCov = cov(Beta[burnin:nsim,])

library(mvtnorm)


tmpBet = rmvnorm(1,bs,bCov)
tnorm.mvtRcpp(,lo =lo.t,hi=hi.t)
thetanew = t(apply(thetanew,1,function(x) tnorm.mvtRcpp(x,x,sigma/4,lo = lo.t,hi=hi.t,times=5)))
norm(n.rect.h(230,tmpBet))

thetaP = xPred%*%betastar
allP = rnorm(length(iPred),n.rect.h(iPred,thetaP[,1],thetaP[,2],thetaP[,3],thetaP[,4]),sqrt(sig))  
allTot = cumsum(allP)

Title=''
xTit = 'Day of Year'


par(mar=c(3, 4, 3, 2), oma=c(0,0,0,0), bg="white", xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
rng = 4880:5500
ind = sample(1:100,1)
plot(allAssim[ind,rng],type='l', bty="n", las=1, main=Title, xlab=bquote(bold(.(xTit))), ylab=expression("Assimiation mol CO"[2] * " m"^-2 * " s"^-1), family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8, xaxt="n", yaxt="n")
 tmp = allAssim[ind,rng]
 tmp[allAssim[ind,rng] < 0] = 0
 polygon(x=c(1:len(rng),rev(1:len(rng))),y=c(tmp,rep(0,len(rng))),col=rgb(.2,.2,.7,.5),border=NA)
  tmp = allAssim[ind,rng]
  tmp[allAssim[ind,rng] >= 0] = 0
  polygon(x=c(1:len(rng),rev(1:len(rng))),y=c(tmp,rep(0,len(rng))),col=rgb(.7,.2,.5,.5),border=NA)
 abline(h=0,col='black',type=3)
 grid()
 at=seq(1,len(rng),length = 10)
 axis(1,at= at,days[rng][at]-3*365,tick=F)
 axis(2,tick=F)

load(file='warming/ENVDATA/processedData/sapENV.RData') 
plots = TtoC(rownames(allAssim),"\\|")[,1]
whrEnv = match(plots,STDplots(gsub('Q',"",colnames(SM)[grep('Q',colnames(SM))]))) + 8
Q[,whrEnv[ind]]
days
n = table(days)
dayFrac= numeric()
for(i in 1:len(n)){
  dayFrac = c(dayFrac,seq(0,1-(1/24),length= n[i]))
}
DayF = days+dayFrac

par(mfrow=c(3,1),mar=c(0, 4, 3, 1), oma=c(0,0,0,0), bg="white", xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
rng = 4900:5500
ind = sample(1:100,1)
plot(DayF[rng],allAssim[ind,rng],type='l', bty="n", las=1, main=Title, xlab=bquote(bold(.(xTit))), ylab=expression("Assimiation mol CO"[2] * " m"^-2 * " s"^-1), family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8, xaxt="n", yaxt="n")
tmp = allAssim[ind,rng]
tmp[allAssim[ind,rng] < 0] = 0
polygon(x=c(DayF[rng],rev(DayF[rng])),y=c(tmp,rep(0,len(rng))),col=rgb(.2,.2,.7,.25),border=NA)
tmp = allAssim[ind,rng]
tmp[allAssim[ind,rng] >= 0] = 0
polygon(x=c(DayF[rng],rev(DayF[rng])),y=c(tmp,rep(0,len(rng))),col=rgb(.7,.2,.5,.25),border=NA)
abline(h=0,col='black',type=3)
grid()

axis(2,tick=F)

par(mar=c(1.5, 4, 1.5, 1))
whr = which(DayF[rng[1]] <= Q[,'times'] & Q[,'times'] <= DayF[last(rng)])
tmp = Q[whr,whrEnv[ind]] 
tmp[tmp<0] = 0
tmp = rollmean(tmp,14)
plot(Q[whr,'times'],tmp,type='l', bty="n", las=1, main=Title, xlab=bquote(bold(.(xTit))), ylab=expression("PPFD mol CO"[2] * " m"^-2 * " s"^-1), family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8, xaxt="n", yaxt="n")
polygon(x=c(Q[whr,'times'],rev(Q[whr,'times'])),y=c(tmp,rep(0,len(whr))),col=rgb(.2,.2,.7,.25),border=NA)
abline(h=0,col='black',lty=3)
axis(2,tick=F)
grid()
par(mar=c(3, 4, 0, 1))
plot(SM[whr,'times'],SM[whr,whrEnv[ind]],type='l', bty="n", las=1, main=Title, xlab=bquote(bold(.(xTit))), ylab=expression("Assimiation mol CO"[2] * " m"^-2 * " s"^-1), family="Helvetica", cex.main=1.5, cex.axis=0.8, cex.lab=0.8, xaxt="n", yaxt="n")
axis(2,tick=F)
abline(h=0,col='black',lty=3)
at=seq(1,len(whr),length = 10)
axis(1,at= Q[whr[at],'times'],floor(Q[whr[at],'times'])-3*365,tick=F)
grid() 

###
# GS plot
##


par(mfrow=c(1,1),mar=c(3,2,1,1),bg = 'white')
split.screen(c(1,2))
colors = c("#00A0B096","#00A0B0","#6A4A3C96","#6A4A3C","#CC333F96","#CC333F","#EB684196","#EB6841","#EDC951")
colors=c("black","grey","black","grey","black","grey","black","grey")
shape = c(15,16,17,18)
plot(0,0,ylim=c(0,1),xlim=c(790,1090),yaxt='n',xaxt='n',xlab='Day of Year',ylab='')
axis(1,seq(800,1100,50),seq(800,1100,50)-2*365)
whr = order(pheno[,'SP11'])
lwd = 1/nrow(pheno)
widths = table(paste(pheno[,'Spec'],pheno[,'Light'],sep=""))*lwd
st = c(0,cumsum(widths)[-7])
sp = sort(unique(paste(pheno[,'Spec'],pheno[,'Light'],sep="")))
lig = sort(unique(pheno[,'Light']))
numac = st[1];numacs=st[2];
numli= st[3];numlis=st[4];
numqa= st[5];numqas=st[6];
numqrs=st[7];
for(i in 1:nrow(pheno)){
  id = match(paste(pheno[i,'Spec'],pheno[i,'Light'],sep=""),sp)
  if(id == 1){num = numac}
  if(id == 2){num = numacs}
  if(id == 3){num = numli}
  if(id == 4){num = numlis}
  if(id == 5){num = numqa}
  if(id == 6){num = numqas}
  if(id == 7){num = numqrs}
  print(num)
  rect(pheno[whr[i],'FA11'],num,pheno[whr[i],'SP11'],num+lwd,col=colors[id],border="#F8F8FF")
  if(id == 1){numac = lwd + numac}
  if(id == 2){numacs = lwd + numacs}
  if(id == 3){numli = lwd + numli}
  if(id == 4){numlis = lwd + numlis}
  if(id == 5){numqa = lwd + numqa}
  if(id == 6){numqas = lwd + numqas}
  if(id == 7){numqrs = lwd + numqrs}
}
widths = table(paste(pheno[,'Spec'],pheno[,'Light'],sep=""))*lwd
st = c(0,cumsum(widths))
colss = c('grey','black','grey','black','grey','black','grey','black')
ltys = c(1,1,2,2,3,3,4,4)
for(j in 1:7){
lines(c(tapply(pheno[,'SP11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),median)[j],tapply(pheno[,'FA11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),median)[j]),
      c(((st[-1]+st[-8])/2)[j],((st[-1]+st[-8])/2)[j]),lwd=3,lty=ltys[j],col=colss[j])
}
points(tapply(pheno[,'SP11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),median),((st[-1]+st[-8])/2), col= c('grey','black'),cex=1.5,pch=c(15,15,16,16,17,17,18))
points(tapply(pheno[,'SP11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),quantile,.05),((st[-1]+st[-8])/2),pch = '|', col=  c('grey','black'),cex=2)
points(tapply(pheno[,'SP11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),quantile,.95),((st[-1]+st[-8])/2),pch = '|', col=  c('grey','black'),cex=2)
points(tapply(pheno[,'FA11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),median),((st[-1]+st[-8])/2), col=  c('grey','black'),cex=1.5,pch=c(15,15,16,16,17,17,18))
points(tapply(pheno[,'FA11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),quantile,.05),((st[-1]+st[-8])/2),pch = '|', col=  c('grey','black'),cex=2)
points(tapply(pheno[,'FA11'],paste(pheno[,'Spec'],pheno[,'Light'],sep=""),quantile,.95),((st[-1]+st[-8])/2),pch = '|', col=  c('grey','black'),cex=2)
widths = table(pheno[,'Spec'])*lwd
st = c(0,cumsum(widths))
sp = sort(unique(pheno[,'Spec']))
#text(x = rep(195+365*2,4),y = ((st[-1]+st[-5])/2),sp,col= c('grey'),cex=2)

#assim per meter squared correlation to season length
par(mar=c(0,0,0,0))
split.screen(c(2,1),2)
par(mar=c(1,1,3,3))
screen(3)
colors = c("black","black","black","black","black")
whr = grep('s',pheno[,'Plot'])
plot(pheno[whr,'FA13'] - pheno[whr,'SP13'],totper2[whr],col = colors[match( pheno[whr,'Spec'],sp)],xlab = '',ylab = '',
     xlim = c(170,270),ylim = c(0,300),pch = shape[match( pheno[whr,'Spec'],sp)],type='n',xaxt='n',yaxt='n') #5000, 300
grid(lty=3)
points(pheno[whr,'FA13'] - pheno[whr,'SP13'],totper2[whr],col = colors[match( pheno[whr,'Spec'],sp)],pch = shape[match( pheno[whr,'Spec'],sp)])
points(pheno[whr,'FA12'] - pheno[whr,'SP12'],totper1[whr],col = colors[match( pheno[whr,'Spec'],sp)],pch = shape[match( pheno[whr,'Spec'],sp)])
axis(4,tick = F)

cors = numeric()
for( i in 1:4){
  whr2 = which(pheno[whr,'Spec'] == sp[i])
  if(length(whr2) == 0)next
  fit = lm( c(totper2[whr][whr2],totper1[whr][whr2])~c(pheno[whr[whr2],'FA13'] - pheno[whr[whr2],'SP13'],pheno[whr[whr2],'FA12'] - pheno[whr[whr2],'SP12']))
  abline(coefficients(fit),col = colors[i],lty=i,lwd=3)
  cors = c(cors,round(cor(c(totper2[whr][whr2],totper1[whr][whr2]),c(pheno[whr[whr2],'FA13'] - pheno[whr[whr2],'SP13'],pheno[whr[whr2],'FA12'] - pheno[whr[whr2],'SP12'])),2))  
}
title("Shade")
legend(170,300,sp,col=colors,
       lty = 1:4,text.col = 'black',pch=15:18,
       bty='n')

#assim per plant correlation to season length
screen(4)
par(mar=c(3,1,1,3))
colors = c("grey","grey","grey","grey","grey")
whr = grep('g',pheno[,'Plot'])
plot(pheno[whr,'FA13'] - pheno[whr,'SP13'],totper2[whr],col = colors[match( pheno[whr,'Spec'],sp)],xlab = 'Season Length (Days)',ylab = '',
     xlim = c(170,270),ylim = c(0,5000),pch = shape[match( pheno[whr,'Spec'],sp)],type='n',xaxt='n',yaxt='n') #5000, 300
grid(lty=3)
points(pheno[whr,'FA13'] - pheno[whr,'SP13'],totper2[whr],col = colors[match( pheno[whr,'Spec'],sp)],pch = shape[match( pheno[whr,'Spec'],sp)])
points(pheno[whr,'FA12'] - pheno[whr,'SP12'],totper1[whr],col = colors[match( pheno[whr,'Spec'],sp)],pch = shape[match( pheno[whr,'Spec'],sp)])
axis(1,tick = F)
axis(4,tick = F)

cors = numeric()
for( i in 1:4){
  whr2 = which(pheno[whr,'Spec'] == sp[i])
  if(length(whr2) == 0)next
  fit = lm( c(totper2[whr][whr2],totper1[whr][whr2])~c(pheno[whr[whr2],'FA13'] - pheno[whr[whr2],'SP13'],pheno[whr[whr2],'FA12'] - pheno[whr[whr2],'SP12']))
  abline(coefficients(fit),col = colors[i],lty=i,lwd=3)
  cors = c(cors,round(cor(c(totper2[whr][whr2],totper1[whr][whr2]),c(pheno[whr[whr2],'FA13'] - pheno[whr[whr2],'SP13'],pheno[whr[whr2],'FA12'] - pheno[whr[whr2],'SP12'])),2))  
}
title("Gap",col='grey')
#legend(170,4000,paste(sp[1:3]," (",cors,")",sep=""),col=colors,
#       lty = 1:4,title = 'Gap',text.col = 'black',pch=15:18)

close.screen(,all.screens=T)

