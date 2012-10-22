bc<-read.csv('GSE1456-MAS5data.csv')
clin<-read.csv('GSE1456-sampleinfo.csv',row.names=1)

sum(colnames(bc)==rownames(clin))
match(colnames(bc),rownames(clin))

###

names(clin)

table(clin$Grade)

attach(clin)

table(Grade)

table(Grade,Rec.event)
prop.table(table(Grade,Rec.event),1)
prop.table(table(Grade,Rec.event),2)
prop.table(table(Grade,Rec.event))

fisher.test(table(Grade,Rec.event))

library(survival)
Surv(Rec.time,Rec.event)

survfit(Surv(Rec.time,Rec.event)~1)

plot(survfit(Surv(Rec.time,Rec.event) ~ 1), xlab="Time (years)",ylab="Proportion relapse-free")

plot(survfit(Surv(Rec.time,Rec.event) ~ Grade), xlab="Time (years)",ylab="Proportion relapse-free", col=3:1)
legend(1,0.5,fill=3:1,c("Grade 1","Grade 2","Grade3"))

#legend(1,0.5,fill=3:1,names(table(Grade)))
#legend(1,0.5,fill=3:1,gsub("G","Grade ",names(table(Grade))))

survdiff(Surv(Rec.time,Rec.event) ~ Grade)

#########################

library(limma)

design<-model.matrix(~ER.status)
fit<-lmFit(bc,design)
fit<-eBayes(fit)
tt<-topTable(fit,coef=2,n=nrow(bc))

library(hgu133a.db)
hgu133a()
sym<-as.list(hgu133aSYMBOL)
gname<-as.list(hgu133aGENENAME)
unlist(sym[1:5])

names(tt)

tt$ID[1]

match(tt$ID[1],names(sym))

sym[match(tt$ID[1],names(sym))]
gname[match(tt$ID[1],names(gname))]

###

library(gplots)

bc.sig<-as.matrix(bc[match(tt$ID[tt$adj.P.Val<0.001],rownames(bc)),])

heatmap.2(bc.sig,trace='none',col=greenred(50),scale='row')

heatmap.2(bc.sig,trace='none',col=greenred(50),scale='row',ColSideColors=ifelse(ER.status=="ER+","blue","red"))

zz<-t(apply(as.matrix(bc.sig),1,function(x) (x-mean(x))/sd(x)))
heatmap.2(zz,trace='none',col=greenred(50),scale='none',ColSideColors=ifelse(ER.status=="ER+","blue","red"))
zz[zz>3]<-3
zz[zz< -3]<- -3
heatmap.2(zz,trace='none',col=greenred(50),scale='none',ColSideColors=ifelse(ER.status=="ER+","blue","red"))

# Correlation distance
heatmap.2(zz,trace='none',col=greenred(50),scale='none',ColSideColors=ifelse(ER.status=="ER+","blue","red"),distfun=function(x) as.dist(1-cor(t(x))))

plot(hclust(dist(t(zz))),hang=-1)

ct<-cutree(hclust(dist(t(zz))),2)

table(ct)
table(ct,ER.status)

#st.col<-c("red","purple","lightblue","blue","green")[as.numeric(as.factor(Subtype))]
#heatmap.2(zz,trace='none',col=greenred(50),scale='none',ColSideColors=st.col)

sym[match(tt$ID[1],names(sym))]

ppic.dat<-as.matrix(bc)[match(tt$ID[1],rownames(bc)),]

boxplot(ppic.dat~Rec.event)

############

library(pamr)

bc.dat<-list(x=as.matrix(bc),y=as.vector(Rec.event))

train<-pamr.train(bc.dat)
cv<-pamr.cv(train,bc.dat)
pamr.plotcv(cv)

