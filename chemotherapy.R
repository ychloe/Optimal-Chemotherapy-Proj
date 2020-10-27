###############
##    EDA    ##
###############
drugdata <- read.delim("C:/Users/patrickverdun/Desktop/drugdata.txt")
summary(drugdata)
sd(drugdata$yn)
sd(drugdata$yc)
min(drugdata$yc) #0.01140725
minyc=drugdata[drugdata$yc==min(drugdata$yc),] #min cancer cell 
minyn=drugdata[drugdata$yn==min(drugdata$yn),]
maxyc=drugdata[drugdata$yc==max(drugdata$yc),] #max cancer cell 
maxyn=drugdata[drugdata$yn==max(drugdata$yn),]
nrow(minyc) #48
minyc[order(minyc$yn),] #none of them greater than 0.5
max(drugdata$yc)
newyn=drugdata[drugdata$yn>=0.5,]
newyn[order(newyn$yc),]
#run:296  A:10.0 B:10.0 C:100.0 yn:0.5014498 yc:0.2159635
Am=aggregate(drugdata$yc, by=list(Dosage=drugdata$A), FUN=mean)
Bm=aggregate(drugdata$yc, by=list(Dosage=drugdata$B), FUN=mean)
Cm=aggregate(drugdata$yc, by=list(Dosage=drugdata$C), FUN=mean)

Am=aggregate(drugdata$yn, by=list(Dosage=drugdata$A), FUN=mean)
Bm=aggregate(drugdata$yn, by=list(Dosage=drugdata$B), FUN=mean)
Cm=aggregate(drugdata$yn, by=list(Dosage=drugdata$C), FUN=mean)


###############
##  Model    ##
###############
f=lm(yn~factor(A)*factor(B)*factor(C),drugdata)
t=summary.aov(f) #no df for error

drugdata$A=factor(drugdata$A)
drugdata$B=factor(drugdata$B)
drugdata$C=factor(drugdata$C)

m1=lm(yc ~ A + B + C + A:C + A:B + B:C,drugdata)
s1=summary(m1)
m2=lm(yc ~ C + B + A + C:A + B:A,drugdata)
s2=summary(m2)
m3=lm(yc ~ C + B + A + B:A + C:B,drugdata)
s3=summary(m3)
m4=lm(yc ~ C + B + A + C:A + C:B,drugdata)
s4=summary(m4)
m5=lm(yc ~ C + B + A + C:A,drugdata)
s5=summary(m5)
m6=lm(yc ~ C + B + A + C:B,drugdata)
s6=summary(m6)
m7=lm(yc ~ C + B + A + B:A,drugdata)
s7=summary(m7)
m8=lm(yc ~ C + B + A,drugdata)
s8=summary(m8)
list(s1$r.squared,s2$r.squared,s3$r.squared,s4$r.squared,s5$r.squared,s6$r.squared,s7$r.squared,s8$r.squared)
list(s1$sigma,s2$sigma,s3$sigma,s4$sigma,s5$sigma,s6$sigma,s7$sigma,s8$sigma)
list(AIC(m1,k = 2),AIC(m2,k = 2),AIC(m3,k = 2),AIC(m4,k = 2),AIC(m5,k = 2),AIC(m6,k = 2),AIC(m7,k = 2),AIC(m8,k = 2))
list(BIC(m1),BIC(m2),BIC(m3),BIC(m4),BIC(m5),BIC(m6),BIC(m7),BIC(m8))
plot(m1)

ABm=aggregate(drugdata$yc, by=list(A=drugdata$A,B=drugdata$B), FUN=mean)
ACm=aggregate(drugdata$yc, by=list(A=drugdata$A,C=drugdata$C), FUN=mean)
BCm=aggregate(drugdata$yc, by=list(B=drugdata$B,C=drugdata$C), FUN=mean)
mu=mean(drugdata$yc)
A=Am$x-mu
B=Bm$x-mu
C=Cm$x-mu
AB=ABm$x-mu+Am$x[sapply(ABm$A,function(x)which(Am$Dosage==x))]-Bm$x[sapply(ABm$B,function(x)which(Bm$Dosage==x))]
AC=ACm$x-mu+Am$x[sapply(ACm$A,function(x)which(Am$Dosage==x))]-Cm$x[sapply(ACm$C,function(x)which(Cm$Dosage==x))]
BC=BCm$x-mu+Am$x[sapply(BCm$B,function(x)which(Bm$Dosage==x))]-Bm$x[sapply(BCm$C,function(x)which(Cm$Dosage==x))]
RMSE=0.03830322

m1=lm(yn ~ A + B + C + A:B + A:C + B:C,drugdata)
s1=summary(m1)
m2=lm(yn ~ C + B + A + C:A + B:A,drugdata)
s2=summary(m2)
m3=lm(yn ~ C + B + A + B:A + C:B,drugdata)
s3=summary(m3)
m4=lm(yn ~ C + B + A + C:A + C:B,drugdata)
s4=summary(m4)
m5=lm(yn ~ C + B + A + C:A,drugdata)
s5=summary(m5)
m6=lm(yn ~ C + B + A + C:B,drugdata)
s6=summary(m6)
m7=lm(yn ~ C + B + A + B:A,drugdata)
s7=summary(m7)
m8=lm(yn ~ C + B + A,drugdata)
s8=summary(m8)
list(s1$r.squared,s2$r.squared,s3$r.squared,s4$r.squared,s5$r.squared,s6$r.squared,s7$r.squared,s8$r.squared)
list(s1$sigma,s2$sigma,s3$sigma,s4$sigma,s5$sigma,s6$sigma,s7$sigma,s8$sigma)
list(AIC(m1,k = 2),AIC(m2,k = 2),AIC(m3,k = 2),AIC(m4,k = 2),AIC(m5,k = 2),AIC(m6,k = 2),AIC(m7,k = 2),AIC(m8,k = 2))
list(BIC(m1),BIC(m2),BIC(m3),BIC(m4),BIC(m5),BIC(m6),BIC(m7),BIC(m8))
t=summary.aov(m1)
ABm=aggregate(drugdata$yn, by=list(A=drugdata$A,B=drugdata$B), FUN=mean)
ACm=aggregate(drugdata$yn, by=list(A=drugdata$A,C=drugdata$C), FUN=mean)
BCm=aggregate(drugdata$yn, by=list(B=drugdata$B,C=drugdata$C), FUN=mean)
mu=mean(drugdata$yn)
A=Am$x-mu
B=Bm$x-mu
C=Cm$x-mu
AB=ABm$x-mu+Am$x[sapply(ABm$A,function(x)which(Am$Dosage==x))]-Bm$x[sapply(ABm$B,function(x)which(Bm$Dosage==x))]
AC=ACm$x-mu+Am$x[sapply(ACm$A,function(x)which(Am$Dosage==x))]-Cm$x[sapply(ACm$C,function(x)which(Cm$Dosage==x))]
BC=BCm$x-mu+Am$x[sapply(BCm$B,function(x)which(Bm$Dosage==x))]-Bm$x[sapply(BCm$C,function(x)which(Cm$Dosage==x))]
RMSE=0.01958648
###############
##  Graph    ##
###############B
plot(fullmod) #337 outlier
plot(fullmodel) #337 outlier
sort(fullmod$residuals)#337 -4.181322e-01,393 1.252207e-01
sort(fullmodel$residuals)# 337 -2.639963e-01, 393 6.138018e-02
interaction.plot(x.factor     = drugdata$A,
                 trace.factor = drugdata$B, 
                 response     = drugdata$yc, 
                 fun = mean,
                 type="b",
                 col=c("red","orange","yellow","green","blue","purple","pink","black"),  ### Colors for levels of trace var.
                 pch=15,             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

interaction.plot(x.factor     = drugdata$B,
                 trace.factor = drugdata$C, 
                 response     = drugdata$yc, 
                 fun = mean,
                 type="b",
                 col=c("red","orange","yellow","green","blue","purple","pink","black"),  ### Colors for levels of trace var.
                 pch=15,             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

interaction.plot(x.factor     = drugdata$A,
                 trace.factor = drugdata$C, 
                 response     = drugdata$yc, 
                 fun = mean,
                 type="b",
                 col=c("red","orange","yellow","green","blue","purple","pink","black"),  ### Colors for levels of trace var.
                 pch=15,             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

interaction.plot(x.factor     = drugdata$A,
                 trace.factor = drugdata$B, 
                 response     = drugdata$yn, 
                 fun = mean,
                 type="b",
                 col=c("red","orange","yellow","green","blue","purple","pink","black"),  ### Colors for levels of trace var.
                 pch=15,             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

interaction.plot(x.factor     = drugdata$B,
                 trace.factor = drugdata$C, 
                 response     = drugdata$yn, 
                 fun = mean,
                 type="b",
                 col=c("red","orange","yellow","green","blue","purple","pink","black"),  ### Colors for levels of trace var.
                 pch=15,             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

interaction.plot(x.factor     = drugdata$A,
                 trace.factor = drugdata$C, 
                 response     = drugdata$yn, 
                 fun = mean,
                 type="b",
                 col=c("red","orange","yellow","green","blue","purple","pink","black"),  ### Colors for levels of trace var.
                 pch=15,             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")


###############
##    Q3     ##
###############
drugdata$ynfit=predict(fullmodel,data=drugdata)
subdata=drugdata[drugdata$ynfit>=0.5,]
subdata=subdata[,-4:-7]
p=predict(fullmod,subdata)
frommin(p)#0.2211654
sort(p) #296
#A=10, B=10, C=100 minimize yc

###############
##    Q4     ##
###############
drugdata_coded
new=drugdata_coded[drugdata_coded$A==4:6 &drugdata$B==4:6&drugdata$C==4:6,]
View(new)
