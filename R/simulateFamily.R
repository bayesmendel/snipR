makeSimDataset <- function(nFams,pDup, singleP, stdP, parents, gparents, sibs,ants,cuz, errRt,numUnique,birthrt) {
    ip.list <- simIPs(numUnique=numUnique)
    fams <- list()
      for (i in 1:nFams) {
        id <- paste(sample(seq(0,9,1),10,replace=TRUE),collapse="")
        tfam <- simulateFamily(id, pDup, id, ip.list,birthrt)
        tfam$matchID<-id
        tfam$nErrors<-0
        trimfam <- trimFamily(tfam,singleP=singleP,stdP=stdP,parents=parents,gparents=gparents,sibs=sibs,ants=ants,cuz=cuz)
        #if(sum(is.na(trimfam))) {browser()}
        fams[[id]] <- trimfam
        if (fams[[id]][1,"Duplicate"]>0) {
            fams[[id]]$Duptype <- "original"
            for (j in 1:fams[[id]][1,"nDuplicates"]) {
                  tfam2 <- simulateErrors(trimfam, errRt=errRt)
                  errs <- tfam2[[2]]
                 tfam2 <- tfam2[[1]]
                  tid <- paste(sample(seq(0,9,1),10,replace=TRUE),collapse="")
                  tfam2$requestId <- tid
                  tfam2$Duptype <- "isDup"
                  tfam2$nErrors <- errs
                fams[[tid]] <- tfam2
            }
        } else {fams[[id]]$Duptype <- "noDup"}
      }
    #fams <- do.call(rbind.data.frame,fams)
    fams <- plyr::ldply(fams, data.frame)
    return(fams)
}

simulateFamily <- function(famID, pDup, uniqueID, ip.list, birthrt=c(2,2)) {
  data(BRCApenet.metaDSL.2008,death.othercauses, compriskSurv, CBRCApenet.2012, BrOvJointDsn.2014,package="BayesMendel5")
  data(panc.fam,pancpenet.2008,package="BayesMendel5")
  data(MMR.fam,MMRpenet.2008,package="BayesMendel5")
  data(mela.fam,melapenet.HBI.2009,package="BayesMendel5")

  ###### Simulation parameters ######
  # pbrca1,pbrca2's: various q, allele frequency (from Danielle)
  pbrca1 <- .01
  pbrca2 <- .011
  pmlh1 <- .008
  pmsh2 <- .014
  pmsh6 <- .006
  ppanc <- .5
  pp16 <- .008
  q<-c(pbrca1,pbrca2,pmlh1,pmsh2,pmsh6,ppanc,pp16)

  ## generate size of each generation
  f1 <- rpois(1,birthrt[1]) #daughters
  if (f1==0){
    f1<-1
    f1.flag<-1
  } else {f1.flag<-0}
  m1 <- rpois(1,birthrt[2]) #sons
  if (m1==0){
    m1<-1
    m1.flag<-1
  } else {m1.flag<-0}
  g1 <- f1 + m1
  f2 <- rpois(g1,birthrt[1]) #daughters of daughters/sons
  if (any(f2==0)){
    i<-which(f2==0)
    f2[i]<-1
    f2.flag<-i
  } else {f2.flag<-0}
  m2 <- rpois(g1,birthrt[2]) #sons of daughters/sons
  if (any(m2==0)){
    i <- which(m2==0)
    m2[i]<-1
    m2.flag<-i
  } else {m2.flag<-0}
  g2 <- sum(f2)+sum(m2)
  fam.size <- 2*g1+g2+2

  ###### Simulate allele sharing over generations ######
  ## generate parent (generation 0) genotypes
  geno.f=geno.m<-matrix(0,1,length(q)*2)
  geno.fF=geno.mF<-matrix(0, 1, length(q))
  for (j in 1:length(q)){
    geno.f[,(2*j-1):(2*j)] = cbind(rbinom(1,size=1, prob=q[j]), rbinom(1,size=1,prob=q[j]))
    geno.m[,(2*j-1):(2*j)] = cbind(rbinom(1,size=1, prob=q[j]), rbinom(1,size=1,prob=q[j]))
    geno.fF[,j]<-geno.f[,(2*j-1)]+geno.f[,(2*j)]
    geno.fF[,j]<-ifelse(geno.fF[,j]>1,1,geno.fF[,j])
    geno.mF[,j]<-geno.m[,(2*j-1)]+geno.m[,(2*j)]
    geno.mF[,j]<-ifelse(geno.mF[,j]>1,1,geno.mF[,j])
  }
  # which allele to pass to daughters
  geno.d=geno.dd=geno.ds<-matrix(0,1,length(q)*2)
  geno.h<-matrix(0,f1,length(q)*2)
  geno.dF=geno.ddF=geno.dsF<-matrix(0,1,length(q))
  geno.hF<-matrix(0,f1,length(q))
  GenD=GenH<-matrix(0,f1,length(q))
  GenDD=GenDS<-list(f1)
  for (i in 1:f1){
    GenDD[[i]] <- matrix(0,f2[i],length(q))
    GenDS[[i]] <- matrix(0,m2[i],length(q))
    for (j in 1:length(q)){
      indx.f0 = rbinom(1,size=1,prob=0.5)
      indx.m0 = rbinom(1,size=1,prob=0.5)
      geno.d[,(2*j-1)] = geno.f[,(2*j-1)] *indx.f0+ geno.f[,(2*j)]* (1-indx.f0)
      geno.d[,(2*j)] = geno.m[,(2*j-1)] *indx.m0+ geno.m[,(2*j)]*(1-indx.m0)
      geno.dF[,j]<-geno.d[,(2*j-1)]+geno.d[,(2*j)]
      geno.dF[,j]<-ifelse(geno.dF[,j]>1,1,geno.dF[,j])
      GenD[i,j]<-geno.dF[,j]
      # generate husbands' genotypes
      geno.h[i,(2*j-1):(2*j)]=cbind(rbinom(1,size=1, prob=q[j]), rbinom(1,size=1,prob=q[j]))
      geno.hF[i,j]<-geno.h[i,(2*j-1)]+geno.h[i,(2*j)]
      geno.hF[i,j]<-ifelse(geno.hF[i,j]>1,1,geno.hF[i,j])
      GenH[i,j]<-geno.hF[i,j]
      # which alleles to pass to their sons
      if (m2[i]>0) {
        for (k1 in 1:m2[i]){
          indx.f0 = rbinom(1,size=1,prob=0.5)
          indx.m0 = rbinom(1,size=1,prob=0.5)
          geno.ds[,(2*j-1)] = geno.d[,(2*j-1)]*indx.f0+geno.d[,(2*j)]*(1-indx.f0)
          geno.ds[,(2*j)] = geno.h[i,(2*j-1)]*indx.m0+geno.h[i,(2*j)]*(1-indx.m0)
          geno.dsF[,j]<-geno.ds[,(2*j-1)]+geno.ds[,(2*j)]
          geno.dsF[,j]<-ifelse(geno.dsF[,j]>1,1,geno.dsF[,j])
          GenDS[[i]][k1,j]<-geno.dsF[,j]
        }
      }
      # which alleles to pass to their daughters
      if (f2[i]>0) {
        for (k2 in 1:f2[i]){
          indx.f0 = rbinom(1,size=1,prob=0.5)
          indx.m0 = rbinom(1,size=1,prob=0.5)
          geno.dd[,(2*j-1)] = geno.d[,(2*j-1)]*indx.f0+geno.d[,(2*j)]*(1-indx.f0)
          geno.dd[,(2*j)] = geno.h[i,(2*j-1)]*indx.m0+geno.h[i,(2*j)]*(1-indx.m0)
          geno.ddF[,j]<-geno.dd[,(2*j-1)]+geno.dd[,(2*j)]
          geno.ddF[,j]<-ifelse(geno.ddF[,j]>1,1,geno.ddF[,j])
          GenDD[[i]][k2,j]<-geno.ddF[,j]
        }
      }
    }
  }
  # which allele to pass to sons
  geno.s=geno.sd=geno.ss<-matrix(0,1,length(q)*2)
  geno.w<-matrix(0,m1,length(q)*2)
  geno.sF=geno.sdF=geno.ssF<-matrix(0,1,length(q))
  geno.wF<-matrix(0,m1,length(q))
  GenS=GenW<-matrix(0,m1,length(q))
  GenSD=GenSS<-list(m1)
  for (i in 1:m1){
    GenSD[[i]] <- matrix(0,f2[i+f1],length(q))
    GenSS[[i]] <- matrix(0,m2[i+f1],length(q))
    for (j in 1:length(q)){
      indx.f0 = rbinom(1,size=1,prob=0.5)
      indx.m0 = rbinom(1,size=1,prob=0.5)
      geno.s[,(2*j-1)] = geno.f[,(2*j-1)]*indx.f0 + geno.f[,(2*j)]*(1-indx.f0)
      geno.s[,(2*j)] = geno.m[,(2*j-1)]*indx.m0 + geno.m[,(2*j)]*(1-indx.m0)
      geno.sF[,j]<-geno.s[,(2*j-1)]+geno.s[,(2*j)]
      geno.sF[,j]<-ifelse(geno.sF[,j]>1,1,geno.sF[,j])
      GenS[i,j]<-geno.sF[,j]
      # generate wives' genotypes
      geno.w[i,(2*j-1):(2*j)]=cbind(rbinom(1,size=1, prob=q[j]), rbinom(1,size=1,prob=q[j]))
      geno.wF[i,j]<-geno.w[i,(2*j-1)]+geno.w[i,(2*j)]
      geno.wF[i,j]<-ifelse(geno.wF[i,j]>1,1,geno.wF[i,j])
      GenW[i,j]<-geno.wF[i,j]
      #which alleles to pass to their sons
       if (m2[i+f1]>0) {
        for (k1 in 1:m2[i+f1]){
          indx.f0 = rbinom(1,size=1,prob=0.5)
          indx.m0 = rbinom(1,size=1,prob=0.5)
          geno.ss[,(2*j-1)] = geno.s[,(2*j-1)]*indx.f0+geno.s[,(2*j)]*(1-indx.f0)
          geno.ss[,(2*j)] = geno.w[i,(2*j-1)]*indx.m0+geno.w[i,(2*j)]*(1-indx.m0)
          geno.ssF[,j]<-geno.ss[,(2*j-1)]+geno.ss[,(2*j)]
          geno.ssF[,j]<-ifelse(geno.ssF[,j]>1,1,geno.ssF[,j])
          GenSS[[i]][k1,j]<-geno.ssF[,j]
        }
      }
      #which alleles to pass to their daughters
      if (f2[i+f1]>0) {
        for (k2 in 1:f2[i+f1]){
          indx.f0 = rbinom(1,size=1,prob=0.5)
          indx.m0 = rbinom(1,size=1,prob=0.5)
          geno.sd[,(2*j-1)] = geno.s[,(2*j-1)]*indx.f0+geno.s[,(2*j)]*(1-indx.f0)
          geno.sd[,(2*j)] = geno.w[i,(2*j-1)]*indx.m0+geno.w[i,(2*j)]*(1-indx.m0)
          geno.sdF[,j]<-geno.sd[,(2*j-1)]+geno.sd[,(2*j)]
          geno.sdF[,j]<-ifelse(geno.sdF[,j]>1,1,geno.sdF[,j])
          GenSD[[i]][k2,j]<-geno.dsF[,j]
        }
      }
    }
  }

  ###### Fit models ######
  #begin with "brcapro" model, let j=1 correspond to brca1, j=2 correspond to brca2.
  #generate affectedbreast for all male relatives.
  tot=fam.size
  #GenOverall<-array(rep(0,length(q)*tot),dim=c(n,length(q), tot))
  GenOverall<-matrix(0,tot,length(q))
  GenOverall[1,]<-geno.mF
  GenOverall[2,]<-geno.fF
  cD<-2+f1
  GenOverall[3:cD,]<-GenD
  cH<-cD+f1
  GenOverall[(cD+1):cH,]<-GenH
  cS<-cH+m1
  GenOverall[(cH+1):cS,]<-GenS
  cW<-cS+m1
  GenOverall[(cS+1):cW,]<-GenW
  #cDD<-cW+(f1*f2)
  #kids of daughter
  tem<-cW
  #daughters of daughter
  for (ll in 1:f1){
    if (length(GenDD[[ll]])>0) {
      GenOverall[(tem+1):(tem+f2[ll]),]<-GenDD[[ll]]
      tem<-tem+f2[ll]
    }
  }
  cDD<-tem
  #sons of daughter
  for (ll in 1:f1){
    if (length(GenDS[[ll]])>0) {
      GenOverall[(tem+1):(tem+m2[ll]),]<-GenDS[[ll]]
      tem<-tem+m2[ll]
    }
  }
  #kids of son
  cDS<-tem
  #daughters of son
  for (ll in 1:m1){
    if (length(GenSD[[ll]])>0) {
      GenOverall[(tem+1):(tem+f2[f1+ll]),]<-GenSD[[ll]]
      tem<-tem+f2[f1+ll]
    }
  }
  cSD<-tem
  #sons of son
  for (ll in 1:m1){
    #sons of son
    if (length(GenSS[[ll]])>0) {
      GenOverall[(tem+1):(tem+m2[f1+ll]),]<-GenSS[[ll]]
      tem<-tem+m2[f1+ll]
    }
  }
  cSS<-tem

  ## Calculate failure times
  failBC=failOC=failCC=failEC=failPC=failMC<- matrix(0,1,tot)
  timeBC=timeOC=timeCC=timeEC=timePC=timeMC<-c()

  tempM<-c(1,(cD+1):cH,(cH+1):cS,(cDD+1):cDS,(cSD+1):cSS )
  tempF<-c(2,3:cD,(cS+1):cW, (cW+1):cDD, (cDS+1):cSD)

  #all males
  for (j in tempM){
    if (GenOverall[j,1]==1 & GenOverall[j,2]==1) {
      w=BRCApenet.metaDSL.2008$fMX[,5]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,1]==1 & GenOverall[j,2]==0) {
      w=BRCApenet.metaDSL.2008$fMX[,2]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,1]==0 & GenOverall[j,2]==1) {
      w=BRCApenet.metaDSL.2008$fMX[,4]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,1]==0 & GenOverall[j,2]==0) {
      w=BRCApenet.metaDSL.2008$fMX[,1]
      ws=append(w, 1-sum(w))
    }
    new<-rmultinom(1,1,prob=ws)
    timeBC[j]<-which(new[,1]==1)
    timeOC[j]<-111

    if (GenOverall[j,3]==1 & GenOverall[j,4]==1 & GenOverall[j,5]==1) {
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M111")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==1 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M110")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==0 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M100")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==0 & GenOverall[i,5]==0){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M000")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==0 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M001")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==1 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M010")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==1 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M011")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==0 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M101")]
      ws=append(w, 1-sum(w))
    }

    new1<-rmultinom(1,1,prob=ws)
    timeCC[j]<-which(new1[,1]==1)
    timeEC[j]<-111

    if (GenOverall[j,6]==1){
      w=pancpenet.2008$fMX[,2]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,6]==0){
      w=pancpenet.2008$fMX[,1]
      ws=append(w, 1-sum(w))
    }

    new<-rmultinom(1,1,prob=ws)
    timePC[j]<-which(new[,1]==1)

    if (GenOverall[j,7]==1){
      w=melapenet.HBI.2009$fMX[,2]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,7]==0){
      w=melapenet.HBI.2009$fMX[,1]
      ws=append(w, 1-sum(w))
    }

    new<-rmultinom(1,1,prob=ws)
    timeMC[j]<-which(new[,1]==1)
  }

  #all females
  for (j in tempF){
    if (GenOverall[j,1]==1 & GenOverall[j,2]==1){
      w=BRCApenet.metaDSL.2008$fFX[,5]
      ws=append(w, 1-sum(w))
      wo=BRCApenet.metaDSL.2008$fFY[,5]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,1]==1 & GenOverall[j,2]==0){
      w=BRCApenet.metaDSL.2008$fFX[,2]
      ws=append(w, 1-sum(w))
      wo=BRCApenet.metaDSL.2008$fFY[,2]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,1]==0 & GenOverall[j,2]==1){
      w=BRCApenet.metaDSL.2008$fFX[,4]
      ws=append(w, 1-sum(w))
      wo=BRCApenet.metaDSL.2008$fFY[,4]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,1]==0 & GenOverall[j,2]==0){
      w=BRCApenet.metaDSL.2008$fFX[,1]
      ws=append(w, 1-sum(w))
      wo=BRCApenet.metaDSL.2008$fFY[,1]
      wos<-append(wo, 1-sum(wo))
    }
    
    new1<-rmultinom(1,1,prob=ws)
    new2<-rmultinom(1,1,prob=wos)
    timeBC[j]<-which(new1[,1]==1)
    timeOC[j]<-which(new2[,1]==1)

    if (GenOverall[j,3]==1 & GenOverall[j,4]==1 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M111")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M111")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==1 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M110")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M110")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==0 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M100")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M100")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==0 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M000")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M000")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==0 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M001")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M001")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==1 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M010")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M010")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==1 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M011")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M011")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==0 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M101")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M101")]
      wos<-append(wo, 1-sum(wo))
    }

    new1<-rmultinom(1,1,prob=ws)
    new2<-rmultinom(1,1,prob=wos)
    timeCC[j]<-which(new1[,1]==1)
    timeEC[j]<-which(new2[,1]==1)

    if (GenOverall[j,6]==1){
      w=pancpenet.2008$fFX[,2]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,6]==0){
      w=pancpenet.2008$fFX[,1]
      ws=append(w, 1-sum(w))
    }

    new<-rmultinom(1,1,prob=ws)
    timePC[j]<-which(new[,1]==1)

    if (GenOverall[j,7]==1){
      w=melapenet.HBI.2009$fFX[,2]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,7]==0){
      w=melapenet.HBI.2009$fFX[,1]
      ws=append(w, 1-sum(w))
    }
    new<-rmultinom(1,1,prob=ws)
    timeMC[j]<-which(new[,1]==1)
  }
  failBC = timeBC
  failOC = timeOC
  failCC = timeCC
  failEC = timeEC
  failPC = timePC
  failMC = timeMC
  cen<-round(rnorm(tot,55,10))

  failBC<-as.vector(t(failBC))
  indBC = ifelse(((failBC <= cen) & (failBC!=111)),1,0)
  obsBC = failBC
  obsBC[indBC==0] = cen[indBC==0]

  failOC<-as.vector(t(failOC))
  indOC = ifelse(((failOC <= cen) & (failOC!=111)),1,0)
  obsOC = failOC
  obsOC[indOC==0] = cen[indOC==0]

  failCC<-as.vector(t(failCC))
  indCC = ifelse(((failCC <= cen) & (failCC!=111)),1,0)
  obsCC = failCC
  obsCC[indCC==0] = cen[indCC==0]

  failEC<-as.vector(t(failEC))
  indEC = ifelse(((failEC <= cen) & (failEC!=111)),1,0)
  obsEC = failEC
  obsEC[indEC==0] = cen[indEC==0]

  failPC<-as.vector(t(failPC))
  indPC = ifelse(((failPC <= cen) & (failPC!=111)),1,0)
  obsPC = failPC
  obsPC[indPC==0] = cen[indPC==0]

  failMC<-as.vector(t(failMC))
  indMC = ifelse(((failMC <= cen) & (failMC!=111)),1,0)
  obsMC = failMC
  obsMC[indMC==0] = cen[indMC==0]

  ## output the data as the brcapro
  ID=rep(1:tot) 
  gen<-c(1:tot)
  gen[tempM]<-1
  gen[tempF]<-0
  Gender <- gen
  # husbands as fathers
  tf<-c()
  tf2<-c()
  for (jj in 1:f1){
    for (kk in 1:f2[jj]) {
      tf<-c(tf,2+f1+jj)
    }
    for (kk in 1:m2[jj]) {
      tf2<-c(tf2,2+f1+jj)
    }
  }
  # sons as fathers
  tf1<-c()
  tf3<-c()
  for (jj in 1:m1){
    for (kk in 1:f2[f1+jj]) {
      tf1<-c(tf1,2+2*f1+jj)
    }
    for (kk in 1:m2[f1+jj]) {
      tf3<-c(tf3,2+2*f1+jj)
    }
  }
  # assignment of fathers
  ftem<-c(0,0, rep(c(1,0),each=f1), rep(c(1,0),each=m1), tf, tf2, tf1, tf3)
  FatherID=ftem
  
  # dauthers as mothers
  tf<-c()
  tf2<-c()
  for (jj in 1:f1){
    for (kk in 1:f2[jj]) {
      tf<-c(tf,2+jj)
    }
    for (kk in 1:m2[jj]) {
      tf2<-c(tf2,2+jj)
    }
  }
  # wives as mothers
  tf1<-c()
  tf3<-c()
  for (jj in 1:m1){
    for (kk in 1:f2[f1+jj]) {
      tf1<-c(tf1,2+2*f1+m1+jj)
    }
    for (kk in 1:m2[f1+jj]) {
      tf3<-c(tf3,2+2*f1+m1+jj)
    }
  }
  # assignment of mothers
  mtem<-c(0,0, rep(c(2,0),each=f1), rep(c(2,0),each=m1), tf, tf2, tf1, tf3)
  MotherID=mtem

  AffectedBreast = indBC
  AffectedOvary = indOC
  AffectedColon = indCC
  AffectedEndometrium = indEC
  AffectedPancreas = indPC
  AffectedSkin = indMC

  AgeBreast =obsBC
  AgeOvary=obsOC
  AgeColon=obsCC
  AgeEndometrium=obsEC
  AgePancreas=obsPC
  AgeSkin=obsMC

  BRCA1 = as.vector(t(GenOverall[,1]))
  BRCA2 = as.vector(t(GenOverall[,2]))
  MLH1 = as.vector(t(GenOverall[,3]))
  MSH2 = as.vector(t(GenOverall[,4]))
  MSH6 = as.vector(t(GenOverall[,5]))
  PANC = as.vector(t(GenOverall[,6]))
  P16 = as.vector(t(GenOverall[,7]))

  FamilyID=famID

  Twins= rep(0,tot)
  Ethnic=rep(sample(c("nonAJ","AJ"),1,prob=c(.92,.08)), tot)
  dat = data.frame(cbind(ID,Gender, FatherID, MotherID, AffectedBreast,
                       AffectedOvary, AffectedColon,AffectedEndometrium, AffectedPancreas, AffectedSkin, AgeBreast, AgeOvary, AgeColon, AgeEndometrium, AgePancreas, AgeSkin,Twins, Ethnic,
                       BRCA1, BRCA2, MLH1, MSH2, MSH6, PANC, P16))

  ##  Remove generation 1 and generation 2 placeholders (f1.flag, m1.flag, etc.)
  # Remove generation 2 daughters
  if (any(f2.flag>0)) {
    for (i in 1:length(f2.flag)) {
      if (f2.flag[i]<=f1) {
        # Daughters of daughters
        mom.id <- as.numeric(toString(dat[2+f2.flag[i],"ID"]))
        idx <- which(dat$MotherID==mom.id & dat$Gender==0)
        dat <- dat[-idx,]
      } else {
        # Daughters of sons
        dad.id <- as.numeric(toString(dat[2+f1+f2.flag[i],"ID"]))
        idx <- which(dat$FatherID==dad.id & dat$Gender==0)
        dat <- dat[-idx,]
      }
    }
  }
  # Remove generation 2 sons
  if (any(m2.flag>0)) {
    for (i in 1:length(m2.flag)) {
      if (m2.flag[i]<=f1) {
        # Sons of daughters
        mom.id <- as.numeric(toString(dat[2+m2.flag[i],"ID"]))
        idx <- which(dat$MotherID==mom.id & dat$Gender==1)
        dat <- dat[-idx,]
      } else {
        # Sons of sons
        dad.id <- as.numeric(toString(dat[2+f1+m2.flag[i],"ID"]))
        idx <- which(dat$FatherID==dad.id & dat$Gender==1)
        dat <- dat[-idx,]
      }
    }
  }
  # Remove generation 1 sons and wives
  if (m1.flag>0) {
    dad.id <- 2+2*f1+1
    dat <- dat[-c(dad.id,dad.id+1),]
    # and subsequent generation 2 offspring
    dat <- dat[dat$FatherID!=dad.id,]
  }
  # Remove generation 1 daughters and husbands
  if (f1.flag>0) {
    dat <- dat[-c(3,4),]
    # and subsequent generation 2 offspring
    dat <- dat[dat$MotherID!=3,]
  }
  # randomly select generation 0 individual if no gen 1 individuals
  if (f1.flag>0 & m1.flag>0) {
    dat <- dat[dat$Gender==1,]
    dups <- runif(1)
    Duplicate <- as.numeric(dups<pDup)
    nDuplicates <- Duplicate
    if (Duplicate[1]==1) {
      dups <- runif(1)
      while (dups<pDup) {
        nDuplicates <- nDuplicates+1
        dups <- runif(1)
      }
    }
    dat$Duplicate<-Duplicate
    dat$nDuplicates<-nDuplicates
    dat$FamilyID<-FamilyID
    dat$requestId<-uniqueID
    dat$senderIP <- ip.list[[sample(length(ip.list),1,prob=c(rep(.18,3),rep(.03,5), runif(length(ip.list)-8,0,.01)))]]
    dat$ID <- 1
    return(dat)
  # randomly select generation 1 individual if no gen 2 individuals
  } else if (max(as.numeric(dat$FatherID)-1)==1) {
    g1 <- which(dat$FatherID==1)
    if (length(g1)==1) {
      ind <- g1
    } else {
      ttt <- 0
      while (ttt == 0) {
        ind <- sample(g1,1)
        if (dat[ind,"Gender"]==0 | sum(dat[,"Gender"]==0)==0 | runif(1)<.01) {
          ttt <- 1
        }
      }
    }
    dat <- rbind(dat[ind,],dat[-ind,])
    fam.size <- nrow(dat)
    zID <- dat$ID
    tID <- 1:fam.size
    oldID <- as.numeric(as.vector(rownames(dat)))
    tfID <- as.numeric(as.vector(dat$FatherID))
    tfID[tfID!=0] <- tfID[tfID!=0]+1
    tmID <- as.numeric(as.vector(dat$MotherID))
    tmID[tmID!=0] <- tmID[tmID!=0]+1
    dat$ID <- tID
    dat$FatherID<-tfID
    dat$MotherID<-tmID
    preorder <- c(zID[ind],zID[-ind])
    rownames(dat) <- as.character(preorder)
  
    dups <- runif(1)
    Duplicate <- as.numeric(dups<pDup)
    nDuplicates <- Duplicate
    if (Duplicate[1]==1) {
      dups <- runif(1)
      while (dups<pDup) {
        nDuplicates <- nDuplicates+1
        dups <- runif(1)
      }
    }
    dat$Duplicate<-Duplicate
    dat$nDuplicates<-nDuplicates
    dat$FamilyID<-FamilyID
    dat$requestId<-uniqueID
    dat$senderIP <- ip.list[[sample(length(ip.list),1,prob=c(rep(.18,3),rep(.03,5), runif(length(ip.list)-8,0,.01)))]]
    return(dat)
  ## Randomly define self in generation 2
  } else {
    #dat2 <- dat
    fam.size <- nrow(dat)
    if (f1.flag) {
      old.f1 <- f1
      f1<-0
    }
    if (m1.flag) {
      old.m1 <- m1
      m1<-0
    }
    zID <- dat$ID
    g2 <- fam.size-2-2*f1-2*m1
    if (g2==1) {
      ind <- fam.size
    } else {
      ttt <- 0
      while (ttt == 0) {
        ind <- sample((2+2*f1+2*m1+1):fam.size,length(g2))  
        #if (is.na((dat[ind,"Gender"]==0 | runif(1)<.005 | sum(dat[,"Gender"]==0)==0))) {browser()}
        if (dat[ind,"Gender"]==0 | runif(1)<.005 | sum(dat[,"Gender"]==0)==0) {
          ttt <- 1
        }
      }
    }
    dat <- rbind(dat[ind,],dat[-ind,])
    tID <- 1:fam.size
    oldID <- as.numeric(as.vector(rownames(dat)))
    tfID <- as.numeric(as.vector(dat$FatherID))
    if (f1.flag) {
      tfID[tfID!=0 & tfID>2] <- tfID[tfID!=0 & tfID>2]-1
      tfID[tfID!=0 & tfID==1] <- tfID[tfID!=0 & tfID==1]+1
    }
    if (!f1.flag) {tfID[tfID!=0] <- tfID[tfID!=0]+1}
    tmID <- as.numeric(as.vector(dat$MotherID))
    if (f1.flag) {
      tmID[tmID!=0 & tmID>2] <- tmID[tmID!=0 & tmID>2]-1
      tmID[tmID!=0 & tmID==2] <- tmID[tmID!=0 & tmID==2]+1
    }
    if (!f1.flag) {tmID[tmID!=0] <- tmID[tmID!=0]+1}
    dat$ID <- tID
    dat$FatherID<-tfID
    dat$MotherID<-tmID
    if (f1.flag) {f1<-old.f1}
    if (m1.flag) {m1<-old.m1}
    preorder <- c(zID[ind],zID[-ind])
    rownames(dat) <- as.character(preorder)

    ## Identify and build out other side of family
    fID <- dat[1,"FatherID"]
    mID <- dat[1,"MotherID"]
    #if (length(dat[fID,"FatherID"])==0) {browser()}
    if (dat[fID,"FatherID"]==0) {
      allls <- geno.h[as.numeric(dat[mID,"ID"])-3,]
      idy<-fID
    } else {
      if (f1.flag) {
        allls <- geno.w[as.numeric(dat[fID,"ID"])-3,]  
      } else {
        allls <- geno.w[as.numeric(dat[fID,"ID"])-3-2*f1,]
      }
      idy<-mID
    } 
    #if (is.na(allls[1])) {browser()}
    fam2 <- buildInLaws(allls, max(preorder))
    dat[idy,"FatherID"]<-max(preorder)+1
    dat[idy,"MotherID"]<-max(preorder)+2
    dat <- rbind(dat,fam2)
    dups <- runif(1)
    Duplicate <- as.numeric(dups<pDup)
    nDuplicates <- Duplicate
    if (Duplicate[1]==1) {
      dups <- runif(1)
      while (dups<pDup) {
        nDuplicates <- nDuplicates+1
        dups <- runif(1)
      }
    }
    dat$Duplicate<-Duplicate
    dat$nDuplicates<-nDuplicates
    dat$FamilyID<-FamilyID
    dat$requestId<-uniqueID
    dat$senderIP <- ip.list[[sample(length(ip.list),1,prob=c(rep(.18,3),rep(.03,5), runif(length(ip.list)-8,0,.01)))]]
    return(dat)
  }
}

buildInLaws <- function(alleles, add.num, birthrt=c(2,2)) {
###### Simulation parameters ######
    # pbrca1,pbrca2's: various q, allele frequency (from Danielle)
    pbrca1 <- .01
    pbrca2 <- .011
    pmlh1 <- .008
    pmsh2 <- .014
    pmsh6 <- .006
    ppanc <- .5
    pp16 <- .008
    q<-c(pbrca1,pbrca2,pmlh1,pmsh2,pmsh6,ppanc,pp16)

    ## generate size of each generation
  f1 <- rpois(1,birthrt[1]) #daughters
  if (f1==0){
    f1<-1
    f1.flag<-1
  } else {f1.flag<-0}
  m1 <- rpois(1,birthrt[2]) #sons
  if (m1==0){
    m1<-1
    m1.flag<-1
  } else {m1.flag<-0}
  g1 <- f1 + m1
  f2 <- rpois(g1,birthrt[1]) #daughters of daughters/sons
  if (any(f2==0)){
    i<-which(f2==0)
    f2[i]<-1
    f2.flag<-i
  } else {f2.flag<-0}
  m2 <- rpois(g1,birthrt[2]) #sons of daughters/sons
  if (any(m2==0)){
    i <- which(m2==0)
    m2[i]<-1
    m2.flag<-i
  } else {m2.flag<-0}
  g2 <- sum(f2)+sum(m2)
  fam.size <- 2*g1+g2+2

    ###### Simulate allele sharing over generations ######
    ## generate parent (generation 0) genotypes
    geno.f=geno.m<-matrix(0,1,length(q)*2)
    geno.fF=geno.mF<-matrix(0, 1, length(q))
    for (j in 1:length(q)){
        ff <- sample(c(0,1),1)
        geno.f[,(2*j-1):(2*j)] = cbind(rbinom(1,size=1, prob=q[j]), ff*alleles[j-1]+(1-ff)*alleles[j])
        geno.m[,(2*j-1):(2*j)] = cbind(rbinom(1,size=1, prob=q[j]), (1-ff)*alleles[j-1]+ff*alleles[j])
        geno.fF[,j]<-geno.f[,(2*j-1)]+geno.f[,(2*j)]
        geno.fF[,j]<-ifelse(geno.fF[,j]>1,1,geno.fF[,j])
        geno.mF[,j]<-geno.m[,(2*j-1)]+geno.m[,(2*j)]
        geno.mF[,j]<-ifelse(geno.mF[,j]>1,1,geno.mF[,j])
    }
    # which allele to pass to daughters
  geno.d=geno.dd=geno.ds<-matrix(0,1,length(q)*2)
  geno.h<-matrix(0,f1,length(q)*2)
  geno.dF=geno.ddF=geno.dsF<-matrix(0,1,length(q))
  geno.hF<-matrix(0,f1,length(q))
  GenD=GenH<-matrix(0,f1,length(q))
  GenDD=GenDS<-list(f1)
  for (i in 1:f1){
    GenDD[[i]] <- matrix(0,f2[i],length(q))
    GenDS[[i]] <- matrix(0,m2[i],length(q))
    for (j in 1:length(q)){
      indx.f0 = rbinom(1,size=1,prob=0.5)
      indx.m0 = rbinom(1,size=1,prob=0.5)
      geno.d[,(2*j-1)] = geno.f[,(2*j-1)] *indx.f0+ geno.f[,(2*j)]* (1-indx.f0)
      geno.d[,(2*j)] = geno.m[,(2*j-1)] *indx.m0+ geno.m[,(2*j)]*(1-indx.m0)
      geno.dF[,j]<-geno.d[,(2*j-1)]+geno.d[,(2*j)]
      geno.dF[,j]<-ifelse(geno.dF[,j]>1,1,geno.dF[,j])
      GenD[i,j]<-geno.dF[,j]
      # generate husbands' genotypes
      geno.h[i,(2*j-1):(2*j)]=cbind(rbinom(1,size=1, prob=q[j]), rbinom(1,size=1,prob=q[j]))
      geno.hF[i,j]<-geno.h[i,(2*j-1)]+geno.h[i,(2*j)]
      geno.hF[i,j]<-ifelse(geno.hF[i,j]>1,1,geno.hF[i,j])
      GenH[i,j]<-geno.hF[i,j]
      # which alleles to pass to their sons
      if (m2[i]>0) {
        for (k1 in 1:m2[i]){
          indx.f0 = rbinom(1,size=1,prob=0.5)
          indx.m0 = rbinom(1,size=1,prob=0.5)
          geno.ds[,(2*j-1)] = geno.d[,(2*j-1)]*indx.f0+geno.d[,(2*j)]*(1-indx.f0)
          geno.ds[,(2*j)] = geno.h[i,(2*j-1)]*indx.m0+geno.h[i,(2*j)]*(1-indx.m0)
          geno.dsF[,j]<-geno.ds[,(2*j-1)]+geno.ds[,(2*j)]
          geno.dsF[,j]<-ifelse(geno.dsF[,j]>1,1,geno.dsF[,j])
          GenDS[[i]][k1,j]<-geno.dsF[,j]
        }
      }
      # which alleles to pass to their daughters
      if (f2[i]>0) {
        for (k2 in 1:f2[i]){
          indx.f0 = rbinom(1,size=1,prob=0.5)
          indx.m0 = rbinom(1,size=1,prob=0.5)
          geno.dd[,(2*j-1)] = geno.d[,(2*j-1)]*indx.f0+geno.d[,(2*j)]*(1-indx.f0)
          geno.dd[,(2*j)] = geno.h[i,(2*j-1)]*indx.m0+geno.h[i,(2*j)]*(1-indx.m0)
          geno.ddF[,j]<-geno.dd[,(2*j-1)]+geno.dd[,(2*j)]
          geno.ddF[,j]<-ifelse(geno.ddF[,j]>1,1,geno.ddF[,j])
          GenDD[[i]][k2,j]<-geno.ddF[,j]
        }
      }
    }
  }
  # which allele to pass to sons
  geno.s=geno.sd=geno.ss<-matrix(0,1,length(q)*2)
  geno.w<-matrix(0,m1,length(q)*2)
  geno.sF=geno.sdF=geno.ssF<-matrix(0,1,length(q))
  geno.wF<-matrix(0,m1,length(q))
  GenS=GenW<-matrix(0,m1,length(q))
  GenSD=GenSS<-list(m1)
  for (i in 1:m1){
    GenSD[[i]] <- matrix(0,f2[i+f1],length(q))
    GenSS[[i]] <- matrix(0,m2[i+f1],length(q))
    for (j in 1:length(q)){
      indx.f0 = rbinom(1,size=1,prob=0.5)
      indx.m0 = rbinom(1,size=1,prob=0.5)
      geno.s[,(2*j-1)] = geno.f[,(2*j-1)]*indx.f0 + geno.f[,(2*j)]*(1-indx.f0)
      geno.s[,(2*j)] = geno.m[,(2*j-1)]*indx.m0 + geno.m[,(2*j)]*(1-indx.m0)
      geno.sF[,j]<-geno.s[,(2*j-1)]+geno.s[,(2*j)]
      geno.sF[,j]<-ifelse(geno.sF[,j]>1,1,geno.sF[,j])
      GenS[i,j]<-geno.sF[,j]
      # generate wives' genotypes
      geno.w[i,(2*j-1):(2*j)]=cbind(rbinom(1,size=1, prob=q[j]), rbinom(1,size=1,prob=q[j]))
      geno.wF[i,j]<-geno.w[i,(2*j-1)]+geno.w[i,(2*j)]
      geno.wF[i,j]<-ifelse(geno.wF[i,j]>1,1,geno.wF[i,j])
      GenW[i,j]<-geno.wF[i,j]
      #which alleles to pass to their sons
       if (m2[i+f1]>0) {
        for (k1 in 1:m2[i+f1]){
          indx.f0 = rbinom(1,size=1,prob=0.5)
          indx.m0 = rbinom(1,size=1,prob=0.5)
          geno.ss[,(2*j-1)] = geno.s[,(2*j-1)]*indx.f0+geno.s[,(2*j)]*(1-indx.f0)
          geno.ss[,(2*j)] = geno.w[i,(2*j-1)]*indx.m0+geno.w[i,(2*j)]*(1-indx.m0)
          geno.ssF[,j]<-geno.ss[,(2*j-1)]+geno.ss[,(2*j)]
          geno.ssF[,j]<-ifelse(geno.ssF[,j]>1,1,geno.ssF[,j])
          GenSS[[i]][k1,j]<-geno.ssF[,j]
        }
      }
      #which alleles to pass to their daughters
      if (f2[i+f1]>0) {
        for (k2 in 1:f2[i+f1]){
          indx.f0 = rbinom(1,size=1,prob=0.5)
          indx.m0 = rbinom(1,size=1,prob=0.5)
          geno.sd[,(2*j-1)] = geno.s[,(2*j-1)]*indx.f0+geno.s[,(2*j)]*(1-indx.f0)
          geno.sd[,(2*j)] = geno.w[i,(2*j-1)]*indx.m0+geno.w[i,(2*j)]*(1-indx.m0)
          geno.sdF[,j]<-geno.sd[,(2*j-1)]+geno.sd[,(2*j)]
          geno.sdF[,j]<-ifelse(geno.sdF[,j]>1,1,geno.sdF[,j])
          GenSD[[i]][k2,j]<-geno.dsF[,j]
        }
      }
    }
  }

  ###### Fit models ######
  #begin with "brcapro" model, let j=1 correspond to brca1, j=2 correspond to brca2.
  #generate affectedbreast for all male relatives.
  tot=fam.size
  #GenOverall<-array(rep(0,length(q)*tot),dim=c(n,length(q), tot))
  GenOverall<-matrix(0,tot,length(q))
  GenOverall[1,]<-geno.mF
  GenOverall[2,]<-geno.fF
  cD<-2+f1
  GenOverall[3:cD,]<-GenD
  cH<-cD+f1
  GenOverall[(cD+1):cH,]<-GenH
  cS<-cH+m1
  GenOverall[(cH+1):cS,]<-GenS
  cW<-cS+m1
  GenOverall[(cS+1):cW,]<-GenW
  #cDD<-cW+(f1*f2)
  #kids of daughter
  tem<-cW
  #daughters of daughter
  for (ll in 1:f1){
    if (length(GenDD[[ll]])>0) {
      GenOverall[(tem+1):(tem+f2[ll]),]<-GenDD[[ll]]
      tem<-tem+f2[ll]
    }
  }
  cDD<-tem
  #sons of daughter
  for (ll in 1:f1){
    if (length(GenDS[[ll]])>0) {
      GenOverall[(tem+1):(tem+m2[ll]),]<-GenDS[[ll]]
      tem<-tem+m2[ll]
    }
  }
  #kids of son
  cDS<-tem
  #daughters of son
  for (ll in 1:m1){
    if (length(GenSD[[ll]])>0) {
      GenOverall[(tem+1):(tem+f2[f1+ll]),]<-GenSD[[ll]]
      tem<-tem+f2[f1+ll]
    }
  }
  cSD<-tem
  #sons of son
  for (ll in 1:m1){
    #sons of son
    if (length(GenSS[[ll]])>0) {
      GenOverall[(tem+1):(tem+m2[f1+ll]),]<-GenSS[[ll]]
      tem<-tem+m2[f1+ll]
    }
  }
  cSS<-tem

  ## Calculate failure times
  failBC=failOC=failCC=failEC=failPC=failMC<- matrix(0,1,tot)
  timeBC=timeOC=timeCC=timeEC=timePC=timeMC<-c()

  tempM<-c(1,(cD+1):cH,(cH+1):cS,(cDD+1):cDS,(cSD+1):cSS )
  tempF<-c(2,3:cD,(cS+1):cW, (cW+1):cDD, (cDS+1):cSD)

  #all males
  for (j in tempM){
    if (GenOverall[j,1]==1 & GenOverall[j,2]==1) {
      w=BRCApenet.metaDSL.2008$fMX[,5]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,1]==1 & GenOverall[j,2]==0) {
      w=BRCApenet.metaDSL.2008$fMX[,2]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,1]==0 & GenOverall[j,2]==1) {
      w=BRCApenet.metaDSL.2008$fMX[,4]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,1]==0 & GenOverall[j,2]==0) {
      w=BRCApenet.metaDSL.2008$fMX[,1]
      ws=append(w, 1-sum(w))
    }
    new<-rmultinom(1,1,prob=ws)
    timeBC[j]<-which(new[,1]==1)
    timeOC[j]<-111

    if (GenOverall[j,3]==1 & GenOverall[j,4]==1 & GenOverall[j,5]==1) {
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M111")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==1 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M110")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==0 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M100")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==0 & GenOverall[i,5]==0){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M000")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==0 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M001")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==1 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M010")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==1 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M011")]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==0 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fMX[,which(names(MMRpenet.2008$fMX[1,])=="M101")]
      ws=append(w, 1-sum(w))
    }

    new1<-rmultinom(1,1,prob=ws)
    timeCC[j]<-which(new1[,1]==1)
    timeEC[j]<-111

    if (GenOverall[j,6]==1){
      w=pancpenet.2008$fMX[,2]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,6]==0){
      w=pancpenet.2008$fMX[,1]
      ws=append(w, 1-sum(w))
    }

    new<-rmultinom(1,1,prob=ws)
    timePC[j]<-which(new[,1]==1)

    if (GenOverall[j,7]==1){
      w=melapenet.HBI.2009$fMX[,2]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,7]==0){
      w=melapenet.HBI.2009$fMX[,1]
      ws=append(w, 1-sum(w))
    }

    new<-rmultinom(1,1,prob=ws)
    timeMC[j]<-which(new[,1]==1)
  }

  #all females
  for (j in tempF){
    if (GenOverall[j,1]==1 & GenOverall[j,2]==1){
      w=BRCApenet.metaDSL.2008$fFX[,5]
      ws=append(w, 1-sum(w))
      wo=BRCApenet.metaDSL.2008$fFY[,5]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,1]==1 & GenOverall[j,2]==0){
      w=BRCApenet.metaDSL.2008$fFX[,2]
      ws=append(w, 1-sum(w))
      wo=BRCApenet.metaDSL.2008$fFY[,2]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,1]==0 & GenOverall[j,2]==1){
      w=BRCApenet.metaDSL.2008$fFX[,4]
      ws=append(w, 1-sum(w))
      wo=BRCApenet.metaDSL.2008$fFY[,4]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,1]==0 & GenOverall[j,2]==0){
      w=BRCApenet.metaDSL.2008$fFX[,1]
      ws=append(w, 1-sum(w))
      wo=BRCApenet.metaDSL.2008$fFY[,1]
      wos<-append(wo, 1-sum(wo))
    }
    
    new1<-rmultinom(1,1,prob=ws)
    new2<-rmultinom(1,1,prob=wos)
    timeBC[j]<-which(new1[,1]==1)
    timeOC[j]<-which(new2[,1]==1)

    if (GenOverall[j,3]==1 & GenOverall[j,4]==1 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M111")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M111")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==1 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M110")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M110")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==0 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M100")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M100")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==0 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M000")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M000")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==0 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M001")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M001")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==1 & GenOverall[j,5]==0){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M010")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M010")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==0 & GenOverall[j,4]==1 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M011")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M011")]
      wos<-append(wo, 1-sum(wo))
    } else if (GenOverall[j,3]==1 & GenOverall[j,4]==0 & GenOverall[j,5]==1){
      w=MMRpenet.2008$fFX[,which(names(MMRpenet.2008$fFX[1,])=="M101")]
      ws=append(w, 1-sum(w))
      wo=MMRpenet.2008$fFY[,which(names(MMRpenet.2008$fFY[1,])=="M101")]
      wos<-append(wo, 1-sum(wo))
    }

    new1<-rmultinom(1,1,prob=ws)
    new2<-rmultinom(1,1,prob=wos)
    timeCC[j]<-which(new1[,1]==1)
    timeEC[j]<-which(new2[,1]==1)

    if (GenOverall[j,6]==1){
      w=pancpenet.2008$fFX[,2]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,6]==0){
      w=pancpenet.2008$fFX[,1]
      ws=append(w, 1-sum(w))
    }

    new<-rmultinom(1,1,prob=ws)
    timePC[j]<-which(new[,1]==1)

    if (GenOverall[j,7]==1){
      w=melapenet.HBI.2009$fFX[,2]
      ws=append(w, 1-sum(w))
    } else if (GenOverall[j,7]==0){
      w=melapenet.HBI.2009$fFX[,1]
      ws=append(w, 1-sum(w))
    }
    new<-rmultinom(1,1,prob=ws)
    timeMC[j]<-which(new[,1]==1)
  }
  failBC = timeBC
  failOC = timeOC
  failCC = timeCC
  failEC = timeEC
  failPC = timePC
  failMC = timeMC
  cen<-round(rnorm(tot,55,10))

  failBC<-as.vector(t(failBC))
  indBC = ifelse(((failBC <= cen) & (failBC!=111)),1,0)
  obsBC = failBC
  obsBC[indBC==0] = cen[indBC==0]

  failOC<-as.vector(t(failOC))
  indOC = ifelse(((failOC <= cen) & (failOC!=111)),1,0)
  obsOC = failOC
  obsOC[indOC==0] = cen[indOC==0]

  failCC<-as.vector(t(failCC))
  indCC = ifelse(((failCC <= cen) & (failCC!=111)),1,0)
  obsCC = failCC
  obsCC[indCC==0] = cen[indCC==0]

  failEC<-as.vector(t(failEC))
  indEC = ifelse(((failEC <= cen) & (failEC!=111)),1,0)
  obsEC = failEC
  obsEC[indEC==0] = cen[indEC==0]

  failPC<-as.vector(t(failPC))
  indPC = ifelse(((failPC <= cen) & (failPC!=111)),1,0)
  obsPC = failPC
  obsPC[indPC==0] = cen[indPC==0]

  failMC<-as.vector(t(failMC))
  indMC = ifelse(((failMC <= cen) & (failMC!=111)),1,0)
  obsMC = failMC
  obsMC[indMC==0] = cen[indMC==0]

  ## output the data as the brcapro
  ID=rep(1:tot) 
  gen<-c(1:tot)
  gen[tempM]<-1
  gen[tempF]<-0
  Gender <- gen
  # husbands as fathers
  tf<-c()
  tf2<-c()
  for (jj in 1:f1){
    for (kk in 1:f2[jj]) {
      tf<-c(tf,2+f1+jj)
    }
    for (kk in 1:m2[jj]) {
      tf2<-c(tf2,2+f1+jj)
    }
  }
  # sons as fathers
  tf1<-c()
  tf3<-c()
  for (jj in 1:m1){
    for (kk in 1:f2[f1+jj]) {
      tf1<-c(tf1,2+2*f1+jj)
    }
    for (kk in 1:m2[f1+jj]) {
      tf3<-c(tf3,2+2*f1+jj)
    }
  }
  # assignment of fathers
  ftem<-c(0,0, rep(c(1,0),each=f1), rep(c(1,0),each=m1), tf, tf2, tf1, tf3)
  FatherID=ftem
  
  # dauthers as mothers
  tf<-c()
  tf2<-c()
  for (jj in 1:f1){
    for (kk in 1:f2[jj]) {
      tf<-c(tf,2+jj)
    }
    for (kk in 1:m2[jj]) {
      tf2<-c(tf2,2+jj)
    }
  }
  # wives as mothers
  tf1<-c()
  tf3<-c()
  for (jj in 1:m1){
    for (kk in 1:f2[f1+jj]) {
      tf1<-c(tf1,2+2*f1+m1+jj)
    }
    for (kk in 1:m2[f1+jj]) {
      tf3<-c(tf3,2+2*f1+m1+jj)
    }
  }
  # assignment of mothers
  mtem<-c(0,0, rep(c(2,0),each=f1), rep(c(2,0),each=m1), tf, tf2, tf1, tf3)
  MotherID=mtem

  AffectedBreast = indBC
  AffectedOvary = indOC
  AffectedColon = indCC
  AffectedEndometrium = indEC
  AffectedPancreas = indPC
  AffectedSkin = indMC

  AgeBreast =obsBC
  AgeOvary=obsOC
  AgeColon=obsCC
  AgeEndometrium=obsEC
  AgePancreas=obsPC
  AgeSkin=obsMC

  BRCA1 = as.vector(t(GenOverall[,1]))
  BRCA2 = as.vector(t(GenOverall[,2]))
  MLH1 = as.vector(t(GenOverall[,3]))
  MSH2 = as.vector(t(GenOverall[,4]))
  MSH6 = as.vector(t(GenOverall[,5]))
  PANC = as.vector(t(GenOverall[,6]))
  P16 = as.vector(t(GenOverall[,7]))

  Twins= rep(0,tot)
  Ethnic=rep(sample(c("nonAJ","AJ"),1,prob=c(.92,.08)), tot)
  dat = data.frame(cbind(ID,Gender, FatherID, MotherID, AffectedBreast,
                       AffectedOvary, AffectedColon,AffectedEndometrium, AffectedPancreas, AffectedSkin, AgeBreast, AgeOvary, AgeColon, AgeEndometrium, AgePancreas, AgeSkin,Twins, Ethnic,
                       BRCA1, BRCA2, MLH1, MSH2, MSH6, PANC, P16),stringsAsFactors = FALSE)

  ##  Remove generation 1 and generation 2 placeholders (f1.flag, m1.flag, etc.)
  # Remove generation 2 daughters
  if (any(f2.flag>0)) {
    for (i in 1:length(f2.flag)) {
      if (f2.flag[i]<=f1) {
        # Daughters of daughters
        mom.id <- as.numeric(toString(dat[2+f2.flag[i],"ID"]))
        idx <- which(dat$MotherID==mom.id & dat$Gender==0)
        dat <- dat[-idx,]
      } else {
        # Daughters of sons
        dad.id <- as.numeric(toString(dat[2+f1+f2.flag[i],"ID"]))
        idx <- which(dat$FatherID==dad.id & dat$Gender==0)
        dat <- dat[-idx,]
      }
    }
  }
  # Remove generation 2 sons
  if (any(m2.flag>0)) {
    for (i in 1:length(m2.flag)) {
      if (m2.flag[i]<=f1) {
        # Sons of daughters
        mom.id <- as.numeric(toString(dat[2+m2.flag[i],"ID"]))
        idx <- which(dat$MotherID==mom.id & dat$Gender==1)
        dat <- dat[-idx,]
      } else {
        # Sons of sons
        dad.id <- as.numeric(toString(dat[2+f1+m2.flag[i],"ID"]))
        idx <- which(dat$FatherID==dad.id & dat$Gender==1)
        dat <- dat[-idx,]
      }
    }
  }
  # Remove generation 1 sons and wives
  if (m1.flag>0) {
    dad.id <- 2+2*f1+1
    dat <- dat[-c(dad.id,dad.id+1),]
    # and subsequent generation 2 offspring
    dat <- dat[dat$FatherID!=dad.id,]
  }
  # Remove generation 1 daughters and husbands
  if (f1.flag>0) {
    dat <- dat[-c(3,4),]
    # and subsequent generation 2 offspring
    dat <- dat[dat$MotherID!=3,]
  }
    dat$ID <- as.numeric(as.vector(dat$ID))+add.num
    tfID <- as.numeric(as.vector(dat$FatherID))
    tfID[tfID!=0] <- tfID[tfID!=0]+add.num
    tmID <- as.numeric(as.vector(dat$MotherID))
    tmID[tmID!=0] <- tmID[tmID!=0]+add.num
    dat$FatherID<-tfID
    dat$MotherID<-tmID
    rownames(dat) <- as.character(dat$ID)

    return(dat)
}

simulateErrors <- function(family, errRt = 1.5, keyVar.bin, keyVar.cont, keyWt) {
    errs <- rpois(1, errRt) # overall errors number coding errors
    # binary.v <- c("Gender","AffectedBreast","AffectedOvary","AffectedColon","AffectedEndometrium","AffectedPancreas",   
    #             "AffectedSkin","Twins",            
    #             "BRCA1","BRCA2","MLH1","MSH2","MSH6","PANC","P16")
    # cont.v <- c("AgeBreast","AgeOvary","AgeColon","AgeEndometrium","AgePancreas","AgeSkin")
    vars <- c(keyVar.bin, keyVar.cont)
    # var.wts <- rep(1,length(vars))
    var.wts <- cumsum(keyWt/sum(keyWt))
    fam.size <- nrow(family)
    fam.wts <- rep(1,fam.size)
    fam.wts <- cumsum(fam.wts/sum(fam.wts))
    if (errs>0) {
        for (i in 1:errs) {
            what <- sample(vars,1,prob=var.wts)
            if (what %in% keyVar.bin) {
                ll <- 0
                while (ll<1) {
                    who <- sample(fam.size,1,prob=fam.wts)
                    if (!is.na(family[who,what])) {
                      if (length(levels(family[,what]))==1) {
                        levels(family[,what]) <- c("0","1")
                      }
                      family[who,what] <- as.character(1-(as.numeric(family[who,what])-1))
                      ll<-1
                    }
                }
            } else {
                ll <- 0
                while (ll<1) {
                    who <- sample(fam.size,1,prob=fam.wts)
                    if(!is.na(family[who,what])) {
                        if (is.factor(family[,what])) {
                            tvec <- as.numeric(levels(family[,what])[family[,what]])
                        } else {
                            tvec <- as.numeric(family[,what])
                        }
                        tv <- as.character(tvec[who])
                        ii <- sample(1:nchar(tv),1)
                        abb <- as.character(as.numeric(substr(tv,ii,ii))+round(runif(1,0,9)))
                        if (nchar(abb)>1) {abb <- substr(abb,2,2)}
                        substr(tv,ii,ii) <- abb
                        tvec[who] <- as.numeric(tv)
                        #if(dim(family)[1]>length(as.factor(as.character(tvec)))){browser()}
                        family[,what] <- as.factor(as.character(tvec))
                        ll<-1
                    }
                }
            }
        }
    }
    return(list(family,errs))
}

simIPs <- function(numUnique=10) {
    IP.list <- lapply(seq(numUnique), function(i) {
                    paste0(toString(round(runif(1,0,9))),toString(round(runif(1,0,9))),toString(round(runif(1,0,9))),".",
                           toString(round(runif(1,0,9))),toString(round(runif(1,0,9))),".",
                           toString(round(runif(1,0,9))),toString(round(runif(1,0,9))),toString(round(runif(1,0,9))),".",
                           toString(round(runif(1,0,9))))
        })
    return(IP.list)
}

trimFamily <- function(family,singleP=.5,stdP=.4,parents=.98,gparents=.9,sibs=.5,ants=.5,cuz=.1) {
    if (nrow(family)==0) {
      keeper <- family
    } else if (runif(1)<singleP | nrow(family)<7) {
      keeper <- family[family$ID==1,]
    } else {
        if (runif(1)<(stdP)) {
            target.size <- 7
        } else {
            target.size <- 7 + round(rexp(1,.1))
        }
        fam.size <- sum(family$FatherID>0 & family$MotherID>0)+4
        if (target.size>=fam.size) {
            keeper <- family
        } else {
            keeper <- family[family$ID==1,]
            # parents
            i <- which(family[,"ID"]==keeper[1,"FatherID"])
            mem <- family[i,]
            if (runif(1)<parents) {
                keeper <- rbind(keeper,mem)
            }
            j <- which(family[,"ID"]==keeper[1,"MotherID"])
            mem <- family[j,]
            if (runif(1)<parents) {
                keeper <- rbind(keeper,mem)
            }
            # grandparents
            k <- which(family[,"ID"]==family[i,"FatherID"])
            mem <- family[k,]
            if (runif(1)<gparents) {
                keeper <- rbind(keeper,mem)
            }
            l <- which(family[,"ID"]==family[i,"MotherID"])
            mem <- family[l,]
            if (runif(1)<gparents) {
                keeper <- rbind(keeper,mem)
            }
            m <- which(family[,"ID"]==family[j,"FatherID"])
            mem <- family[m,]
            if (runif(1)<gparents) {
                keeper <- rbind(keeper,mem)
            }
            n <- which(family[,"ID"]==family[j,"MotherID"])
            mem <- family[n,]
            if (runif(1)<gparents) {
                keeper <- rbind(keeper,mem)
            }
            fam <- family[-c(1,i,j,k,l,m,n),]
            while (nrow(keeper)<target.size) {
                i <- sample(nrow(fam),1)
                mem <- fam[i,]
                if (mem$FatherID==keeper[1,"FatherID"] | mem$MotherID==keeper[1,"MotherID"]) {
                    # siblings
                    if  (runif(1)<sibs) {
                        keeper <- rbind(keeper,mem)
                        fam <- fam[-i,]
                    }
                } else if ((mem$FatherID==family[which(family[,"ID"]==keeper[1,"FatherID"]),"FatherID"]) |
                           (mem$MotherID==family[which(family[,"ID"]==keeper[1,"FatherID"]),"MotherID"]) |
                           (mem$FatherID==family[which(family[,"ID"]==keeper[1,"MotherID"]),"FatherID"]) |
                           (mem$MotherID==family[which(family[,"ID"]==keeper[1,"MotherID"]),"MotherID"])) {
                    # aunts/uncles (excluding who married in)
                    if  (runif(1)<ants) {
                        keeper <- rbind(keeper,mem)
                        fam <- fam[-i,]
                    }
                # cousins
                } else if (mem$FatherID!=0 & mem$MotherID!=0) {
                    if  (runif(1)<cuz) {
                        keeper <- rbind(keeper,mem)
                        fam <- fam[-i,]
                    }
                }
            }
        }
    }
    #if (nrow(keeper) < 7 & nrow(keeper)>1) {
    #  browser()
    #}
    return(keeper)
}

buildFamily <- function(dd, id, sz) {
  if (sz>1) {
    tempfam <- dd[dd$requestId==id,]
    tempself <- tempfam[tempfam$ID==1,]
    tempfam <- tempfam[tempfam$ID!=1,]
    # Find father
    tempfather <- findFather(tempself,tempfam,"Father")
    rez <- cbind(tempself,tempfather)
    colnames(tempfather) <- colnames(tempself)
    # Find mother
    tempmom <- findMother(tempself,tempfam, "Mother")
    rez <- cbind(rez,tempmom)
    colnames(tempmom) <- colnames(tempself)
    # Find fathers father
    tempfatherfather <- findFather(tempfather,tempfam,"FFather")
    rez <- cbind(rez,tempfatherfather)
    # Find fathers mother
    tempfathermother <- findMother(tempfather,tempfam,"FMother")
    rez <- cbind(rez,tempfathermother)
    # Find mothers father
    tempmotherfather <- findFather(tempmom,tempfam,"MFather")
    rez <- cbind(rez,tempmotherfather)
    # Find mothers mother
    tempmothermother <- findMother(tempmom,tempfam,"MMother")
    rez <- cbind(rez,tempmothermother)
    # Find siblings
    sibs <- findSiblings(tempself,tempfam, "Sibling")
    rez <- cbind(rez,sibs)
    # Find fathers siblings
    dubsibs <- findSiblings(tempfather,tempfam, "FSibling")
    rez <- cbind(rez,dubsibs)
    # Find mothers siblings
    dubsibs <- findSiblings(tempmom,tempfam, "MSibling")
    rez <- cbind(rez,dubsibs)
    # Find cousins
    cousins <- findCousins(tempself,tempfam,"Cousin")
    rez <- cbind(rez,cousins)
    # Find children
    children <- findChildren(tempself,tempfam,"Child")
    rez <- cbind(rez,children)
    # Find others
    others <- findOthers(rez,tempfam, "Other")
    rez <- cbind(rez,others)
  } else {
    tempfam <- dd[dd$requestId==id,]
    tempself <- tempfam[tempfam$ID==1,]
    rez <- tempself
  }
  return(rez)
}