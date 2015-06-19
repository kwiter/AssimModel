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

load(file=paste(path,'warming/qualMVNassimilation.Rdata',sep=''))

allA = allR = allQ = allTh = m.beta
allA.lo = allR.lo = allQ.lo = allTh.lo = m.beta
allA.hi = allR.hi = allQ.hi = allTh.hi = m.beta

rownames(allA) = rownames(allR) = rownames(allQ) = rownames(allTh) = sort(species)
allTab = as.character()
for( ii in 1:length(species)){
  
  load(file=paste(path,'warming/',sort(species)[ii],'MVNassimilation.Rdata',sep = ""))
  species = sort(species)
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


par(mfrow=c(2,2),mar=c(3, 4, 2, 2), oma=c(0,0,0,0),xpd=FALSE, xaxs="r", yaxs="i", mgp=c(2.1,.3,0), las=1, col.axis="#434343", col.main="#343434", tck=0, lend=1)
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

Title = c('Maximum Photosynthesis','Respiration','Light Use Efficiency','Curvature')[kk]
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
axis(4,at = xs,labels=rep(c('ACRU','LITU','QUAL','QURU'),ncol(co)),tick=F,cex.axis=.55,hadj = .6)
axis(1,tick=F,cex = .9);#axis(3,tick=F,cex = .9)
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
days[days > 366] = days[days > 366] - 366
on = days
on[days < 91] = 0  
on[days > 287] = 0
on[on !=0] = 1

allAssim = totAssim = allAssSd = totAssSd = numeric()
for( ii in 1:length(species)){
  
  load(file=paste(path,'warming/',sort(species)[ii],'MVNassimilation.Rdata',sep = ""))
  species = sort(species)
  
  allAssim = rbind(allAssim, rbind(allPmean[1:npred]*on ,
                                   allPmean[(npred+1):(2*npred)]*on ,
                                   allPmean[(2*npred+1):(3*npred)]*on ,
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
colnames(allAssim) = colnames(totAssim) = colnames(allAssSd) = colnames(totAssSd) = dfTmp[whr2,'times']
rownames(allAssim) = rownames(totAssim) = rownames(allAssSd) = rownames(totAssSd) = c('AcruShadeControl','AcruGapControl','AcruShadeZero','AcruGapZero','AcruShadeHot','AcruGapHot',
                                            'LituShadeControl','LituGapControl','LituShadeZero','LituGapZero','LituShadeHot','LituGapHot',
                                            'QualShadeControl','QualGapControl','QualShadeZero','QualGapZero','QualShadeHot','QualGapHot',
                                            'QuruShadeControl','QuruGapControl','QuruShadeZero','QuruGapZero','QuruShadeHot','QuruGapHot')
save(allAssim,totAssim,allAssSd,totAssSd,file = paste(path,'warming/Assim.RData',sep=''))

#######
load(file=paste(path,'warming/qualMVNassimilation.Rdata',sep=''))
pheno = read.csv(file=paste(path,'DissertationEssentials/Biomass/BiomassPhen.csv',sep=''),stringsAsFactors = F)
species = sort(species)


allAssim = totAssim = allAssSd = totAssSd = matrix(NA,nrow=0,ncol=len(whr2))
for( ii in 1:nrow(pheno)){
  load(file=paste(path,'warming/',pheno[ii,'Spec'],'MVNassimilation.Rdata',sep = ""))
  npred = length(days)
  days = days-(2*365)
  yr = days
  yr[days <= (366+365)] = 1
  yr[days > 366+365] = 2
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
colnames(allAssim) = colnames(totAssim) = colnames(allAssSd) = colnames(totAssSd) = dfTmp[whr2,'times']
rownames(allAssim) = rownames(totAssim) = rownames(allAssSd) = rownames(totAssSd) = pheno[,'Indiv']
save(allAssim,totAssim,allAssSd,totAssSd,file = paste(path,'warming/AssimIndiv.RData',sep=''))

par(bg='white')
totAssim[totAssim == 0] = NA
plot(totAssim[10,]*pheno[10,'T.Area.cm.2.']/10000,
     ylim=c(-100,50000),type='n')
for(i in 1:nrow(pheno)){
lines(totAssim[i,]*pheno[i,'T.Area.cm.2.']/10000,col=match(pheno[i,'Spec'],species))
}
legend(0,50000,species,1:4)


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
     type='n',ylab = 'Gap per Plant Assimilation')
whr = grep('g',pheno[,'Plot'])
for(i in whr){
  lines(totAssim[i,]*pheno[i,'T.Area.cm.2.']/10000,col=match(pheno[i,'Spec'],species))
}
legend(0,50000,species,1:4,text.col='black')


totAssim[totAssim == 0] = NA
plot(totAssim[10,]*pheno[10,'T.Area.cm.2.']/10000,
     type='n',ylab = 'Shade per Plant Assimilation')
whr = grep('s',pheno[,'Plot'])
for(i in whr){
  lines(totAssim[i,]*pheno[i,'T.Area.cm.2.']/10000,col=match(pheno[i,'Spec'],species))
}
legend(0,500,species,1:4)



totAssim[totAssim == 0] = NA
plot(totAssim[10,]*pheno[10,'T.Area.cm.2.']/10000/pheno[10,'MassG'],
    type='n',ylab = 'Gap per Plant/Mass Assimilation')
whr = grep('g',pheno[,'Plot'])
for(i in whr){
  lines(totAssim[i,]*pheno[i,'T.Area.cm.2.']/10000/pheno[i,'MassG'],col=match(pheno[i,'Spec'],species))
}
legend(0,6000,species,1:4)


totAssim[totAssim == 0] = NA
plot(totAssim[10,]*pheno[10,'T.Area.cm.2.']/10000/pheno[10,'MassG'],
     type='n',ylab = 'Shade per Plant/Mass Assimilation')
whr = grep('s',pheno[,'Plot'])
for(i in whr){
  lines(totAssim[i,]*pheno[i,'T.Area.cm.2.']/10000/pheno[i,'MassG'],col=match(pheno[i,'Spec'],species))
}
legend(0,20000,species,1:4)

tot1 = apply(allAssim[,1:(ncol(allAssim)/2)]/10000 ,1,sum,na.rm=T)
tot2 = apply(allAssim[,(ncol(allAssim)/2):ncol(allAssim)]/10000 ,1,sum,na.rm=T)

totper1 = apply(allAssim[,1:(ncol(allAssim)/2)]*pheno[,'MassG']/10000 ,1,sum,na.rm=T)
totper2 = apply(allAssim[,(ncol(allAssim)/2):ncol(allAssim)]*pheno[,'MassG']/10000 ,1,sum,na.rm=T)
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


par(mfrow=c(1,1),mar=c(4,1,0,1),bg = '#F8F8FF')
split.screen(c(1,2))
colors = c("#00A0B096","#00A0B0","#6A4A3C96","#6A4A3C","#CC333F96","#CC333F","#EB684196","#EB6841","#EDC951")
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
par(mar=c(4,0,0,0))
split.screen(c(2,1),2)
colors = c("#00A0B0","#6A4A3C","#CC333F","#EB6841","#EDC951")
screen(3)
whr = grep('s',pheno[,'Plot'])
plot(pheno[whr,'FA11'] - pheno[whr,'SP11'],tot2[whr],col = colors[match( pheno[whr,'Spec'],sp)],xlab = 'Growing Season Length',ylab = 'Assimilation' )
points(pheno[whr,'FA10'] - pheno[whr,'SP10'],tot1[whr],col = colors[match( pheno[whr,'Spec'],sp)])
cors = numeric()
for( i in 1:4){
  whr2 = which(pheno[whr,'Spec'] == sp[i])
  if(length(whr2) == 0)next
  fit = lm( c(tot2[whr][whr2],tot1[whr][whr2])~c(pheno[whr[whr2],'FA11'] - pheno[whr[whr2],'SP11'],pheno[whr[whr2],'FA10'] - pheno[whr[whr2],'SP10']))
  abline(coefficients(fit),col = colors[i])
  cors = c(cors,round(cor(c(tot2[whr][whr2],tot1[whr][whr2]),c(pheno[whr[whr2],'FA11'] - pheno[whr[whr2],'SP11'],pheno[whr[whr2],'FA10'] - pheno[whr[whr2],'SP10'])),2))  
}
legend(130,.7,paste(sp," (",cors,")",sep=""),fill = colors,title = 'Shade',text.col = 'black')
legend(190,2,paste(sp[1:3]," (",cors,")",sep=""),fill = colors,title = 'Gap',text.col = 'black')

#assim per plant correlation to season length
screen(4)
colors = c("#00A0B0","#6A4A3C","#CC333F","#EB6841","#EDC951")
whr = grep('s',pheno[,'Plot'])
plot(pheno[whr,'FA11'] - pheno[whr,'SP11'],totper2[whr],col = colors[match( pheno[whr,'Spec'],sp)],xlab = 'Growing Season Length',ylab = 'Assimilation' )
points(pheno[whr,'FA10'] - pheno[whr,'SP10'],totper1[whr],col = colors[match( pheno[whr,'Spec'],sp)])
cors = numeric()
for( i in 1:4){
  whr2 = which(pheno[whr,'Spec'] == sp[i])
  if(length(whr2) == 0)next
  fit = lm( c(totper2[whr][whr2],totper1[whr][whr2])~c(pheno[whr[whr2],'FA11'] - pheno[whr[whr2],'SP11'],pheno[whr[whr2],'FA10'] - pheno[whr[whr2],'SP10']))
  abline(coefficients(fit),col = colors[i])
  cors = c(cors,round(cor(c(totper2[whr][whr2],totper1[whr][whr2]),c(pheno[whr[whr2],'FA11'] - pheno[whr[whr2],'SP11'],pheno[whr[whr2],'FA10'] - pheno[whr[whr2],'SP10'])),2))  
}
legend(130,.2,paste(sp," (",cors,")",sep=""),fill = colors,title = 'Shade',text.col = 'black')
legend(190,60,paste(sp[1:3]," (",cors,")",sep=""),fill = colors,title = 'Gap',text.col = 'black')



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



radarPlot(allA,allA.lo,allA.hi,scale.within = F,legend = F)
radarPlot(allA,allA.lo,allA.hi,scale.within = T,legend = F)
radarPlot(allR,allR.lo,allR.hi,scale.within = F,legend = F)
radarPlot(allR,allR.lo,allR.hi,scale.within = T,legend = F)
radarPlot(allQ,allQ.lo,allQ.hi,scale.within = F,legend = F)
radarPlot(allQ,allQ.lo,allQ.hi,scale.within = T,legend = F)
radarPlot(allTh,allTh.lo,allTh.hi,scale.within = F,legend = F)
radarPlot(allTh,allTh.lo,allTh.hi,scale.within = T,legend = F)


allPoints = read.csv(file = paste(path,'warming/Li6400/licorPROCESSED/all.rates.points.csv',sep=''),stringsAsFactors = F)
colors = c("#00A0B0","#6A4A3C","#CC333F","#EB6841","#EDC951")
colorsT = c("#00A0B080","#6A4A3C80","#CC333F80","#EB684180","#EDC95180")
plot(allPoints[,'PARi'],allPoints[,'Photo'])

par(mfrow=c(2,2),mar=c(0,4,4,0))
whr = which(allPoints[,'species'] == 'acru' & allPoints[,'Light'] == 'Gap')
plot(allPoints[whr,'PARi'],allPoints[whr,'Photo'],pch=20,xlim=c(0,2000),ylim=c(-5,20),col=colors[1],cex=.75 ,xaxt='n',ylab='Photosynthesis')
whr = which(allPoints[,'species'] == 'acru' & allPoints[,'Light'] == 'Shade' )
points(allPoints[whr,'PARi'],allPoints[whr,'Photo'],col=rgb(.2,.2,.2,.7),cex=.75,pch=3)
title(main='acru',line=-1)
abline(h = 0,lty=3)

par(mar=c(0,0,4,4))
whr = which(allPoints[,'species'] == 'litu' & allPoints[,'Light'] == 'Gap' )
plot(allPoints[whr,'PARi'],allPoints[whr,'Photo'],pch=20,xlim=c(0,2000),ylim=c(-5,20),col=colors[2],cex=.75 ,xaxt='n',yaxt='n')
whr = which(allPoints[,'species'] == 'litu' & allPoints[,'Light'] == 'Shade' )
points(allPoints[whr,'PARi'],allPoints[whr,'Photo'],col=rgb(.2,.2,.2,.7),cex=.75,pch=3)
abline(h = 0,lty=3)
title(main='litu',line=-1)
axis(4)

par(mar=c(4,4,0,0))
whr = which(allPoints[,'species'] == 'qual' & allPoints[,'Light'] == 'Gap' )
plot(allPoints[whr,'PARi'],allPoints[whr,'Photo'],pch=20,xlim=c(0,2000),ylim=c(-5,20),col=colors[3],cex=.75 ,ylab='',xlab='PPFD')
whr = which(allPoints[,'species'] == 'qual' & allPoints[,'Light'] == 'Shade' )
points(allPoints[whr,'PARi'],allPoints[whr,'Photo'],col=rgb(.2,.2,.2,.7),cex=.75,pch=3)
title(main='qual',line=-1)
abline(h = 0,lty=3)

par(mar=c(4,0,0,4))
whr = which(allPoints[,'species'] == 'quru' & allPoints[,'Light'] == 'Gap' )
plot(allPoints[whr,'PARi'],allPoints[whr,'Photo'],pch=20,xlim=c(0,2000),ylim=c(-5,20),col=colors[4],cex=.75, yaxt='n',xlab='')
whr = which(allPoints[,'species'] == 'quru' & allPoints[,'Light'] == 'Shade' )
points(allPoints[whr,'PARi'],allPoints[whr,'Photo'],col=rgb(.2,.2,.2,.7),cex=.75,pch=3)
abline(h = 0,lty=3)
title(main='quru',line=-1)
axis(4)
