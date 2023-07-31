rm(list = ls())
options(stringsAsFactors = F)
setwd()
#unique matrix
readin<-read.csv("0122_GNHS_allsample_to_use.csv",header = T)
gene<-read.csv("sp_protein_gene.csv",header = T)
for (i in 19:456) {
  indsel<-which(colnames(readin)[i]==gene$prot)
  if(length(indsel)>0){
    colnames(readin)[i]<-paste(colnames(readin)[i],gene$gene[indsel],sep = "_")
  }
}
readin[,19:456]<-log2(readin[,19:456])
pool<-readin[grepl("pool",readin$prot),]
QC1<-readin[grepl("NL9998",readin$match),]
QC2<-readin[grepl("NL9999",readin$match),]
QC3<-readin[grepl("NL9997",readin$match),]
QC<-rbind(QC1,QC2,QC3)
pool_QC<-rbind(pool,QC)
sample<-readin[!readin$prot%in%pool_QC$prot,]
nasamp<-sample[is.na(sample$match),]
sample<-sample[!is.na(sample$match),]
# write.csv(nasamp,"0122_GNHS_sample_NA_to_use.csv",row.names = F)
# write.csv(sample,"0122_GNHS_sample_to_analysis_matrix.csv",row.names = F)
# write.csv(pool_QC,"0122_GNHS_pool_qc_of_uessamples.csv",row.names = F)

##technical repetition & biological repetition####
# a<-sample$phase_batch[duplicated(sample$phase_batch)]
protmat<-sample[,19:456]
row.names(protmat)<-sample$prot
TechRep<-unique(sample$phase_batch[duplicated(sample$phase_batch)])

corTechRep<-vector()
m<-1
for(i in 1:length(TechRep)){
  indSel<-which(sample$phase_batch==TechRep[i])
  if(length(indSel)>0){
    matCor<-as.matrix(protmat[indSel,])
    comAll<-combn(nrow(matCor),2)
    for(j in 1:ncol(comAll)){
      #corrlation#
      corTechRep[m]<-cor(matCor[comAll[1,j],],matCor[comAll[2,j],], use = "na.or.complete")
      names(corTechRep)[m]<-paste0(row.names(matCor)[comAll[,j]],collapse = ";")
      m<-m+1
    }
  }else{
    print(i)
    
  }
}
# b<-sample$match[duplicated(sample$match)]
# d<-unique(sample$match)
bioRep<-unique(sample$match[duplicated(sample$match)])
corbioRep<-vector()
q<-1
for(i in 1:length(bioRep)){
  indSel<-which(sample$match==bioRep[i])
  if(length(indSel)>0){
    matCor<-as.matrix(protmat[indSel,])
    comAll<-combn(nrow(matCor), 2)
    for(j in 1:ncol(comAll)){
      corbioRep[q]<-cor(matCor[comAll[1,j],],matCor[comAll[2,j],], use = "na.or.complete")
      names(corbioRep)[q]<-paste0(row.names(matCor)[comAll[,j]],collapse = ";")
      q<-q+1
    }
  }else{
    print(i)
  }
}
BioRep<-corbioRep[!names(corbioRep)%in%names(corTechRep)]
bioname<-unique(unlist(strsplit(names(BioRep),";")))
Biosample<-sample[sample$prot%in%bioname,]
unibio<-unique(Biosample$match)
sum(sample$biological_replicate==1)
# write.csv(corbioRep,"0122_GNHS_techrep_corrlation.csv")
# write.csv(BioRep,"0122_GNHS_biorep_corrlation.csv")
medtech<-median(corTechRep)
medbio<-median(BioRep)
library(vioplot)
library(zoo)
pdf("0122_GNHS_sample_correlation_techrep_biorep.pdf",width =8,height = 6)
listRep<-list(corTechRep,BioRep)
names(listRep)<-c("corTechRep","BioRep")
vioplot(listRep,col="gray",ylim=c(0,1))
legend("bottomleft",bty="n",legend = c("Samples:10668 ","Technical replicates:5331","Biological replicates:518"),cex=1.5)
dev.off()

######
unibatch<-unique(Biosample$phase_batch)
unibio<-data.frame(matrix(NA,nrow = 1,ncol = ncol(Biosample)))
colnames(unibio)<-colnames(Biosample)
m<-1
for (i in 1:length(unibatch)) {
  indsel<-Biosample[which(Biosample$phase_batch==unibatch[i]),]
  unibio[m,]<-indsel[which.min(indsel$time),]
  m<-m+1
  
}
# c<-unibio$match[duplicated(unibio$match)]
# sum(Biosample$biological_replicate==1)
bioRep1<-unique(unibio$match[duplicated(unibio$match)])
protmat1<-unibio[,19:456]
row.names(protmat1)<-unibio$prot
corbioRep1<-vector()
q<-1
for(i in 1:length(bioRep1)){
  indSel<-which(unibio$match==bioRep1[i])
  if(length(indSel)>1){
    matCor<-as.matrix(protmat1[indSel,])
    comAll<-combn(nrow(matCor), 2)
    for(j in 1:ncol(comAll)){
      corbioRep1[q]<-cor(matCor[comAll[1,j],],matCor[comAll[2,j],], use = "na.or.complete")
      names(corbioRep1)[q]<-paste0(row.names(matCor)[comAll[,j]],collapse = ";")
      q<-q+1
    }
  }else{
    print(i)
  }
}
write.csv(corbioRep1,"0122_GNHS_biorep_true_corrlation.csv")

trubio<-median(corbioRep1)


library(vioplot)
library(zoo)
pdf("0122_GNHS_sample_correlation_techrep_biorep_true.pdf",width =8,height = 6)
listRep<-list(corTechRep,corbioRep1)
names(listRep)<-c("corTechRep","BioRep")
vioplot(listRep,col="gray",ylim = c(0,1))
legend("bottomleft",bty="n",legend = c("Samples:10668 ","Technical replicates:5331","Biological replicates:518"),cex=1.5)
dev.off()

MeanNA<-data.frame(matrix(NA,nrow=length(TechRepNA),ncol=ncol(NAprot)))
colnames(MeanNA)<-colnames(NAprot)
row.names(MeanNA)<-TechRepNA
for (i in 1:length(TechRepNA)) {
  indsel<-which(nasamp$phase_batch==TechRepNA[i])
  if(length(indsel)>1){
    MeanNA[i,]<-log2(apply(2^NAprot[indsel,],2,mean,na.rm=T)) 
  }
  
}
NAother<-nasamp[!nasamp$phase_batch%in%row.names(MeanNA),]
NAotherprot<-NAother[,19:456]
row.names(NAotherprot)<-NAother$phase_batch
analysis<-rbind(Meanmatrix,MeanNA,NAotherprot)
write.csv(analysis,"0122_GNHS_analysis_matrix_include_NA.csv",row.names = T)
missinganalysis<-sum(is.na(analysis))/(nrow(analysis)*ncol(analysis))
sampinf<-read.csv("1214_GNHS_phase1-3_sampleinf_check_inculdFH.csv",header = T)
analysis$match<-row.names(analysis)
sampinf<-sampinf[!duplicated(sampinf$match),]
analysis<-merge(analysis,sampinf[,12:13],by="match",all.x = T)
write.csv(analysis,"0122_GNHS_analysis_matrix_include_NA_and_inf.csv",row.names = F)


################################
#linear mixing model
GNHS<-read.csv("0106_GNHS_matrix.csv",header = T)
GNHS$pat<-GNHS$match
GNHS$pat<-gsub("baseline","",GNHS$pat)
GNHS$pat<-gsub("F2","",GNHS$pat)
GNHS$pat<-gsub("F3","",GNHS$pat)
patient<-read.csv("0127_GNHS_analysis_patientinf.csv",header = T)
patBS<-patient[,grepl("bs",colnames(patient))]
patBS$pat<-patient$X
patF2<-patient[,grepl("f2",colnames(patient))]
patF2$pat<-patient$X
patF3<-patient[,grepl("f3",colnames(patient))]
patF3$pat<-patient$X
GNHSBS<-GNHS[grepl("baseline",GNHS$match),]
GNHSF2<-GNHS[grepl("F2",GNHS$match),]
GNHSF3<-GNHS[grepl("F3",GNHS$match),]
BS<-merge(GNHSBS,patBS,by="pat",all.x = T)
colnames(BS)<-gsub("_bs","",colnames(BS))
F2<-merge(GNHSF2,patF2,by="pat",all.x = T)
colnames(F2)<-gsub("_f2","",colnames(F2))
F3<-merge(GNHSF3,patF3,by="pat",all.x = T)
colnames(F3)<-gsub("_f3","",colnames(F3))
GNHSmat<-rbind(BS,F2,F3)
write.csv(GNHSmat,"0122_GNHS_mat_patinf.csv",row.names = F)

#multiple linear regression analysis##
library(lmerTest)
patmat<-GNHSmat[,c(1:284,299,307)]
singmat<-protP<-sexP<-ageP<-vector()
m<-1
q<-1
for (i in 3:(ncol(patmat)-3)) {
  fm1=lmer(mets~patmat[,i]+sex+age+(1|pat),
           data = patmat)
  c=summary(fm1)
  
  if(isSingular(fm1)){
    singmat[q]<-colnames(patmat)[i]
    q<-q+1
  }
  else{
    
    # anova(fm1)
    protP[m]<-summary(fm1)$coefficients[,5][2] 
    sexP[m]<-summary(fm1)$coefficients[,5][3]
    ageP[m]<-summary(fm1)$coefficients[,5][4]
    names(protP)[m]<-colnames(patmat)[i]
    names(sexP)[m]<-colnames(patmat)[i]
    names(ageP)[m]<-colnames(patmat)[i]
    m<-m+1
  }}
# write.csv(protP,"0122_GNHS_xianxinhunhemoxing_p_all.csv")
p<-protP[protP<0.05]
write.csv(p,"0122_GNHS_xianxinhunhemoxing_p0.05_159prot.csv")
samP<-data.frame(cbind(protP,sexP,ageP))
# r<-a[a<0.1]
# write.csv(samP,"0127_GNHS_xianxinhunhemoxing_p_age_sex_prot.csv")

########################################
#Comparative analysis of different underlying disease types
#t-test by mets=1&different classifications 
library(ROSE)
mets<-gnhs[gnhs$mets==1,]
mets$type<-paste(mets$fbg,mets$hdl,mets$sbpdbp,mets$tg,mets$waist,sep="_")
ddd<-data.frame(table(mets$type))
colnames(ddd)<-c("fbg_hdl_sbpdbp_tg_waist","sum")
#DM VS hyper
#1:6
a<-mets[mets$type=="0_0_1_1_1",]
b<-mets[mets$type=="1_0_0_1_1",]
pat1<-a[!a$pat_ID%in%b$pat_ID,]
pat2<-b[!b$pat_ID%in%a$pat_ID,]
# m<-setdiff(pat1$pat_ID,pat2$pat_ID)
# aaa<-unique(pat1$pat_ID)
temp<-rbind(pat1,pat2)
temp<-temp[,c(1:441,452,454,469,470)]
temp[,2:439][is.na(temp[,2:439])]<-66666
# table(temp$test_type)
balb<-ovun.sample(type~.,data = temp,
                  p=0.5,seed=1,method = "over")$data

# table(balb$test_type)
balb[balb==66666]<-NA
pat1<-balb[balb$type=="0_0_1_1_1",2:439]
pat2<-balb[balb$type=="1_0_0_1_1",2:439]

FC<-vector();
pvalue<-vector();
for(i in 1:ncol(pat1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = F,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}

pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)
# write.csv(FC_PAd,"0809_0_x_1_allprot_FC_PAd.csv")


protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]

#4:12
a<-mets[mets$type=="0_1_1_1_0",]
b<-mets[mets$type=="1_1_0_1_0",]
pat1<-a[!a$pat_ID%in%b$pat_ID,]
pat2<-b[!b$pat_ID%in%a$pat_ID,]
# m<-setdiff(pat1$pat_ID,pat2$pat_ID)
# aaa<-unique(pat1$pat_ID)
temp<-rbind(pat1,pat2)
temp<-temp[,c(1:441,452,454,469,470)]
temp[,2:439][is.na(temp[,2:439])]<-66666
# table(temp$test_type)
balb<-ovun.sample(type~.,data = temp,
                  p=0.5,seed=1,method = "over")$data

# table(balb$test_type)
balb[balb==66666]<-NA
pat1<-balb[balb$type=="0_1_1_1_0",2:439]
pat2<-balb[balb$type=="1_1_0_1_0",2:439]

FC<-vector();
pvalue<-vector();
for(i in 1:ncol(pat1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = F,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}

pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)
# write.csv(FC_PAd,"0809_0_x_1_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]

#13:5
a<-mets[mets$type=="1_1_0_1_1",]
b<-mets[mets$type=="0_1_1_1_1",]
pat1<-a[!a$pat_ID%in%b$pat_ID,]
pat2<-b[!b$pat_ID%in%a$pat_ID,]
# m<-setdiff(pat1$pat_ID,pat2$pat_ID)
# aaa<-unique(pat1$pat_ID)
temp<-rbind(pat1,pat2)
temp<-temp[,c(1:441,452,454,469,470)]
temp[,2:439][is.na(temp[,2:439])]<-66666
# table(temp$test_type)
balb<-ovun.sample(type~.,data = temp,
                  p=0.5,seed=1,method = "over")$data

# table(balb$test_type)
balb[balb==66666]<-NA
pat1<-balb[balb$type=="1_1_0_1_1",2:439]
pat2<-balb[balb$type=="0_1_1_1_1",2:439]

FC<-vector();
pvalue<-vector();
for(i in 1:ncol(pat1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = F,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}

pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)
# write.csv(FC_PAd,"0809_0_x_1_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]

#DM VS (TG or HDL)
#7:5
a<-mets[mets$type=="1_0_1_0_1",]
b<-mets[mets$type=="0_1_1_1_1",]
pat1<-a[!a$pat_ID%in%b$pat_ID,]
pat2<-b[!b$pat_ID%in%a$pat_ID,]

pat1<-pat1[,2:439]
pat2<-pat2[,2:439]

FC<-vector();
pvalue<-vector();
for(i in 1:ncol(pat1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = F,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}

pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)
# write.csv(FC_PAd,"0809_0_x_1_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#7:1
a<-mets[mets$type=="1_0_1_0_1",]
b<-mets[mets$type=="0_0_1_1_1",]
pat1<-a[!a$pat_ID%in%b$pat_ID,]
pat2<-b[!b$pat_ID%in%a$pat_ID,]
# m<-setdiff(pat1$pat_ID,pat2$pat_ID)
# aaa<-unique(pat1$pat_ID)
temp<-rbind(pat1,pat2)
temp<-temp[,c(1:441,452,454,469,470)]
temp[,2:439][is.na(temp[,2:439])]<-66666
# table(temp$test_type)
balb<-ovun.sample(type~.,data = temp,
                  p=0.5,seed=1,method = "over")$data

# table(balb$test_type)
balb[balb==66666]<-NA
pat1<-balb[balb$type=="1_0_1_0_1",2:439]
pat2<-balb[balb$type=="0_0_1_1_1",2:439]

FC<-vector();
pvalue<-vector();
for(i in 1:ncol(pat1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = F,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}

pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)
# write.csv(FC_PAd,"0809_0_x_1_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]

#7:3
a<-mets[mets$type=="1_0_1_0_1",]
b<-mets[mets$type=="0_1_1_0_1",]
pat1<-a[!a$pat_ID%in%b$pat_ID,]
pat2<-b[!b$pat_ID%in%a$pat_ID,]
# m<-setdiff(pat1$pat_ID,pat2$pat_ID)
# aaa<-unique(pat1$pat_ID)
temp<-rbind(pat1,pat2)
temp<-temp[,c(1:441,452,454,469,470)]
temp[,2:439][is.na(temp[,2:439])]<-66666
# table(temp$test_type)
balb<-ovun.sample(type~.,data = temp,
                  p=0.5,seed=1,method = "over")$data

# table(balb$test_type)
balb[balb==66666]<-NA
pat1<-balb[balb$type=="1_0_1_0_1",2:439]
pat2<-balb[balb$type=="0_1_1_0_1",2:439]

FC<-vector();
pvalue<-vector();
for(i in 1:ncol(pat1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = F,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}

pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)
# write.csv(FC_PAd,"0809_0_x_1_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#hyper vs (TG/HDL)
#7:13
a<-mets[mets$type=="1_0_1_0_1",]
b<-mets[mets$type=="1_1_0_1_1",]
pat1<-a[!a$pat_ID%in%b$pat_ID,]
pat2<-b[!b$pat_ID%in%a$pat_ID,]
# m<-setdiff(pat1$pat_ID,pat2$pat_ID)
# aaa<-unique(pat1$pat_ID)
temp<-rbind(pat1,pat2)
temp<-temp[,c(1:441,452,454,469,470)]
temp[,2:439][is.na(temp[,2:439])]<-66666
# table(temp$test_type)
balb<-ovun.sample(type~.,data = temp,
                  p=0.5,seed=1,method = "over")$data

# table(balb$test_type)
balb[balb==66666]<-NA
pat1<-balb[balb$type=="1_0_1_0_1",2:439]
pat2<-balb[balb$type=="1_1_0_1_1",2:439]

FC<-vector();
pvalue<-vector();
for(i in 1:ncol(pat1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = F,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}

pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)
# write.csv(FC_PAd,"0809_0_x_1_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]

#7:6
a<-mets[mets$type=="1_0_1_0_1",]
b<-mets[mets$type=="1_0_0_1_1",]
pat1<-a[!a$pat_ID%in%b$pat_ID,]
pat2<-b[!b$pat_ID%in%a$pat_ID,]
# m<-setdiff(pat1$pat_ID,pat2$pat_ID)
# aaa<-unique(pat1$pat_ID)
temp<-rbind(pat1,pat2)
temp<-temp[,c(1:441,452,454,469,470)]
temp[,2:439][is.na(temp[,2:439])]<-66666
# table(temp$test_type)
balb<-ovun.sample(type~.,data = temp,
                  p=0.5,seed=1,method = "over")$data

# table(balb$test_type)
balb[balb==66666]<-NA
pat1<-balb[balb$type=="1_0_1_0_1",2:439]
pat2<-balb[balb$type=="1_0_0_1_1",2:439]

FC<-vector();
pvalue<-vector();
for(i in 1:ncol(pat1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = F,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}

pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)
# write.csv(FC_PAd,"0809_0_x_1_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]

#####################
#Dysregulated proteins overlap
##DM VS hyper
#up
up1<-read.csv("20220314_gnhs1-3_dm_hyper_1_6_protUP.csv",header = T)
up2<-read.csv("20220314_gnhs1-3_dm_hyper_4_12_protUP.csv",header = T)
up3<-read.csv("20220314_gnhs1-3_dm_hyper_13_5_protUP.csv",header = T)
prot<-up1[up1$X%in%up2$X,]
prot<-prot[prot$X%in%up3$X,]

#down
down1<-read.csv("20220314_gnhs1-3_dm_hyper_1_6_protdown.csv",header = T)
down2<-read.csv("20220314_gnhs1-3_dm_hyper_4_12_protdown.csv",header = T)
down3<-read.csv("20220314_gnhs1-3_dm_hyper_13_5_protdown.csv",header = T)
prot<-down1[down1$X%in%down2$X,]
prot<-prot[prot$X%in%down3$X,]

##DM VS tg/hdl
up1<-read.csv("20220314_gnhs1-3_dm_hyper_7_5_protUP.csv",header = T)
up2<-read.csv("20220314_gnhs1-3_dm_hyper_7_1_protUP.csv",header = T)
up3<-read.csv("20220314_gnhs1-3_dm_hyper_7_3_protUP.csv",header = T)
prot<-up1[up1$X%in%up2$X,]
prot<-prot[prot$X%in%up3$X,]
# write.csv(prot,"20220314_gnhs1-3_dm_tg_hdl_overlap_protUP.csv",row.names = F)

#down
down1<-read.csv("20220314_gnhs1-3_dm_hyper_7_5_protdown.csv",header = T)
down2<-read.csv("20220314_gnhs1-3_dm_hyper_7_1_protdown.csv",header = T)
down3<-read.csv("20220314_gnhs1-3_dm_hyper_7_3_protdown.csv",header = T)
prot<-down1[down1$X%in%down2$X,]
prot<-prot[prot$X%in%down3$X,]
# write.csv(prot,"20220314_gnhs1-3_dm_tg_hdl_overlap_protdown.csv",row.names = F)


##hyper VS tg/hdl
up1<-read.csv("20220314_gnhs1-3_dm_hyper_7_13_protUP.csv",header = T)
up2<-read.csv("20220314_gnhs1-3_dm_hyper_7_6_protUP.csv",header = T)
prot<-up1[up1$X%in%up2$X,]
# write.csv(prot,"20220314_gnhs1-3_hyper_tg_hdl_overlap_protUP.csv",row.names = F)

#down
down1<-read.csv("20220314_gnhs1-3_dm_hyper_7_13_protdown.csv",header = T)
down2<-read.csv("20220314_gnhs1-3_dm_hyper_7_6_protdown.csv",header = T)
prot<-down1[down1$X%in%down2$X,]
# write.csv(prot,"20220314_gnhs1-3_hyper_tg_hdl_overlap_protdown.csv",row.names = F)

########################################################################################################
#mets0 and 1 were grouped separately, and then difference analysis by group
readin<-read.csv("20220304_GNHS1-3_uniquesample_inculde_inf.csv",header = T)

readin$pat_ID<-readin$match
readin$pat_ID<-gsub("baseline","",readin$pat_ID)
readin$pat_ID<-gsub("F2","",readin$pat_ID)
readin$pat_ID<-gsub("F3","",readin$pat_ID)
gnhs<-readin[readin$mets==1|readin$mets==0,]
mets<-gnhs[gnhs$mets==1,]
mets$type<-paste(mets$fbg,mets$hdl,mets$sbpdbp,mets$tg,mets$waist,sep="_")
ddd<-data.frame(table(mets$type))
colnames(ddd)<-c("fbg_hdl_sbpdbp_tg_waist","sum")

con<-gnhs[gnhs$mets==0,]
con$type<-paste(con$fbg,con$hdl,con$sbpdbp,con$tg,con$waist,sep="_")
ccc<-data.frame(table(con$type))
colnames(ddd)<-c("fbg_hdl_sbpdbp_tg_waist","sum")

#DM_METS
#A6 VS B5
A6<-mets[mets$type=="1_0_0_1_1",]
B5<-con[con$type=="0_0_0_1_1",]
# m<-A6$pat_ID[A6$pat_ID%in%B5$pat_ID]
# a6<-A6[!A6$pat_ID%in%m,]
# b5<-B5[!B5$pat_ID%in%m,]
b5<-B5[!B5$pat_ID%in%A6$pat_ID,]
temp<-rbind(A6,b5)
# install.packages("MatchIt")
library(MatchIt)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.2)
m.data <- match.data(m.out)

hist(unlist(A6[,2:439]))
plot(density(unlist(A6[,2:439]),na.rm = T))
qqnorm(unlist(A6[,2:439]))
qqline(unlist(A6[,2:439]))

hist(unlist(b5[,2:439]))
plot(density(unlist(b5[,2:439]),na.rm = T))
qqnorm(unlist(b5[,2:439]))
qqline(unlist(b5[,2:439]))


# wilcox.test(test1,test2,paired=TRUE)

pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"0809_0_x_1_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
# protdown<-FC_PAd[down,]



#A7 VS B7
A7<-mets[mets$type=="1_0_1_0_1",]
B7<-con[con$type=="0_0_1_0_1",]
b7<-B7[!B7$pat_ID%in%A7$pat_ID,]
temp<-rbind(A7,b7)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.2)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(A7[,2:439]))
plot(density(unlist(A7[,2:439]),na.rm = T))
qqnorm(unlist(A7[,2:439]))
qqline(unlist(A7[,2:439]))

hist(unlist(b7[,2:439]))
plot(density(unlist(b7[,2:439]),na.rm = T))
qqnorm(unlist(b7[,2:439]))
qqline(unlist(b7[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]

#A8 VS B9
A8<-mets[mets$type=="1_0_1_1_0",]
B9<-con[con$type=="0_0_1_1_0",]
b9<-B9[!B9$pat_ID%in%A8$pat_ID,]
temp<-rbind(A8,b9)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.2)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(A8[,2:439]))
plot(density(unlist(A8[,2:439]),na.rm = T))
qqnorm(unlist(A8[,2:439]))
qqline(unlist(A8[,2:439]))

hist(unlist(b9[,2:439]))
plot(density(unlist(b9[,2:439]),na.rm = T))
qqnorm(unlist(b9[,2:439]))
qqline(unlist(b9[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
# protdown<-FC_PAd[down,]


#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#HDL_METS

#A2 VS B5
A2<-mets[mets$type=="0_1_0_1_1",]
B5<-con[con$type=="0_0_0_1_1",]
b5<-B5[!B5$pat_ID%in%A2$pat_ID,]
temp<-rbind(A2,b5)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.2)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(A2[,2:439]))
plot(density(unlist(A2[,2:439]),na.rm = T))
qqnorm(unlist(A2[,2:439]))
qqline(unlist(A2[,2:439]))

hist(unlist(b5[,2:439]))
plot(density(unlist(b5[,2:439]),na.rm = T))
qqnorm(unlist(b5[,2:439]))
qqline(unlist(b5[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#A3 VS B7
A3<-mets[mets$type=="0_1_1_0_1",]
B7<-con[con$type=="0_0_1_0_1",]
b7<-B7[!B7$pat_ID%in%A3$pat_ID,]
temp<-rbind(A3,b7)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.2)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(A3[,2:439]))
plot(density(unlist(A3[,2:439]),na.rm = T))
qqnorm(unlist(A3[,2:439]))
qqline(unlist(A3[,2:439]))

hist(unlist(b7[,2:439]))
plot(density(unlist(b7[,2:439]),na.rm = T))
qqnorm(unlist(b7[,2:439]))
qqline(unlist(b7[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#A4 VS B9
A4<-mets[mets$type=="0_1_1_1_0",]
B9<-con[con$type=="0_0_1_1_0",]
b9<-B9[!B9$pat_ID%in%A4$pat_ID,]
temp<-rbind(A4,b9)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.2)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(A4[,2:439]))
plot(density(unlist(A4[,2:439]),na.rm = T))
qqnorm(unlist(A4[,2:439]))
qqline(unlist(A4[,2:439]))

hist(unlist(b9[,2:439]))
plot(density(unlist(b9[,2:439]),na.rm = T))
qqnorm(unlist(b9[,2:439]))
qqline(unlist(b9[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")
protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]

#TG METS
#A1 VS B7
A1<-mets[mets$type=="0_0_1_1_1",]
B7<-con[con$type=="0_0_1_0_1",]
b7<-B7[!B7$pat_ID%in%A1$pat_ID,]
temp<-rbind(A1,b7)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.2)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(A1[,2:439]))
plot(density(unlist(A1[,2:439]),na.rm = T))
qqnorm(unlist(A1[,2:439]))
qqline(unlist(A1[,2:439]))

hist(unlist(b7[,2:439]))
plot(density(unlist(b7[,2:439]),na.rm = T))
qqnorm(unlist(b7[,2:439]))
qqline(unlist(b7[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#A2 VS B14
A2<-mets[mets$type=="0_1_0_1_1",]
B14<-con[con$type=="0_1_0_0_1",]
a2<-A2[!A2$pat_ID%in%B14$pat_ID,]
temp<-rbind(a2,B14)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.5)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(a2[,2:439]))
plot(density(unlist(a2[,2:439]),na.rm = T))
qqnorm(unlist(a2[,2:439]))
qqline(unlist(a2[,2:439]))

hist(unlist(B14[,2:439]))
plot(density(unlist(B14[,2:439]),na.rm = T))
qqnorm(unlist(B14[,2:439]))
qqline(unlist(B14[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#A4 VS B16
A4<-mets[mets$type=="0_1_1_1_0",]
B16<-con[con$type=="0_1_1_0_0",]
a4<-A4[!A4$pat_ID%in%B14$pat_ID,]
temp<-rbind(a4,B16)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.2)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(a4[,2:439]))
plot(density(unlist(a4[,2:439]),na.rm = T))
qqnorm(unlist(a4[,2:439]))
qqline(unlist(a4[,2:439]))

hist(unlist(B16[,2:439]))
plot(density(unlist(B16[,2:439]),na.rm = T))
qqnorm(unlist(B16[,2:439]))
qqline(unlist(B16[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#A6 VS B18
A6<-mets[mets$type=="1_0_0_1_1",]
B18<-con[con$type=="1_0_0_0_1",]
b18<-B18[!B18$pat_ID%in%A7$pat_ID,]
temp<-rbind(A6,b18)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.2)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(A6[,2:439]))
plot(density(unlist(A6[,2:439]),na.rm = T))
qqnorm(unlist(A6[,2:439]))
qqline(unlist(A6[,2:439]))

hist(unlist(b18[,2:439]))
plot(density(unlist(b18[,2:439]),na.rm = T))
qqnorm(unlist(b18[,2:439]))
qqline(unlist(b18[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#A8 VS B20
A8<-mets[mets$type=="1_0_1_1_0",]
B20<-con[con$type=="1_0_1_0_0",]
b20<-B20[!B20$pat_ID%in%A8$pat_ID,]
temp<-rbind(A8,b20)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.2)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(A8[,2:439]))
plot(density(unlist(A8[,2:439]),na.rm = T))
qqnorm(unlist(A8[,2:439]))
qqline(unlist(A8[,2:439]))

hist(unlist(b20[,2:439]))
plot(density(unlist(b20[,2:439]),na.rm = T))
qqnorm(unlist(b20[,2:439]))
qqline(unlist(b20[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]



#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]



#A3 VS B16
A3<-mets[mets$type=="0_1_1_0_1",]
B16<-con[con$type=="0_1_1_0_0",]
a3<-A3[!A3$pat_ID%in%B16$pat_ID,]
temp<-rbind(a3,B16)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.5)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(a3[,2:439]))
plot(density(unlist(a3[,2:439]),na.rm = T))
qqnorm(unlist(a3[,2:439]))
qqline(unlist(a3[,2:439]))

hist(unlist(B16[,2:439]))
plot(density(unlist(B16[,2:439]),na.rm = T))
qqnorm(unlist(B16[,2:439]))
qqline(unlist(B16[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]



#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]



#A7 VS B20
A7<-mets[mets$type=="1_0_1_0_1",]
B20<-con[con$type=="1_0_1_0_0",]
a7<-A7[!A7$pat_ID%in%B20$pat_ID,]
temp<-rbind(a7,B20)
m.out <- matchit(mets ~age+sex , data = temp,
                 method = "nearest",ratio=1,caliper=0.5)
m.data <- match.data(m.out)
# shapiro.test(unlist(m.data[,2:439]))
hist(unlist(a7[,2:439]))
plot(density(unlist(a7[,2:439]),na.rm = T))
qqnorm(unlist(a7[,2:439]))
qqline(unlist(a7[,2:439]))

hist(unlist(B20[,2:439]))
plot(density(unlist(B20[,2:439]),na.rm = T))
qqnorm(unlist(B20[,2:439]))
qqline(unlist(B20[,2:439]))


pat1<-m.data[m.data$mets==1,c(2:439,473)]
pat2<-m.data[m.data$mets==0,c(2:439,473)]
pat2<-pat2[match(pat1$subclass,pat2$subclass),]
pat1$subclass
pat2$subclass
#t-test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>10&sum(!is.na(pat2[,i]))>10){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-t.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]


#wilcox.test
FC<-vector();
pvalue<-vector();
for(i in 1:(ncol(pat1)-1)){
  if(sum(!is.na(pat1[,i]))>1&sum(!is.na(pat2[,i]))>1){
    FC[i]<-mean(2^pat1[,i],na.rm=T)/mean(2^pat2[,i],na.rm = T)
    pvalue[i]<-wilcox.test(2^pat1[,i],2^pat2[,i], paired = T,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
}


pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd)
row.names(FC_PAd)<-colnames(pat1)[1:438]
# write.csv(FC_PAd,"20220322_gnhs1-3_DM_mets_A6_B7__wilcox_test_paired_allprot_FC_PAd.csv")

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]









