library(readxl) 
library(openxlsx) 
library(spdep) 
library(rgdal) 
library(raster) 
library(corrplot)  
library(DescTools) 
library(nortest) 
library(car) 
library(spatialreg)
library(sf)
library(nortest) 
library(DescTools) 
library(lmtest)
library(tidyverse)
library(spgwr)
library(readr)
analisisGWR <- read_csv("F:/fahmi/Thesis/Analisis/PanselaThesis-20230124T142948Z-001/PL/PL_Raster/hasil/analisisGWR.csv")
View(analisisGWR)
analisis=analisisGWR
head(analisis)
dim(analisis)
peta<-readOGR(dsn="F:/fahmi/Thesis/Analisis/PanselaThesis-20230124T142948Z-001", layer = "AnalisisGWR")
analisis[!complete.cases(analisis),]
summary(analisis)
boxplot(analisis$y,main="Sebaran Peubah Y") 
boxplot(analisis$x1,main="Sebaran Peubah X1") 
boxplot(analisis$x2,main="Sebaran Peubah X2") 
boxplot(analisis$x3,main="Sebaran Peubah X3") 
plot(analisis$x1, analisis$y,
     xlab="X1", 
     ylab="Y",
     main="Scatter Plot of X1 and Y",
     pch=20, col="orange", cex=2)
reg.klasik = lm(y~x1, data=analisis)
lines.lm(reg.klasik, col=2, add=T)

plot(analisis$x2, analisis$y,
     xlab="X2", 
     ylab="Y",
     main="Scatter Plot of X2 and Y",
     pch=20, col="orange", cex=2)
reg.klasik = lm(y~x2, data=analisis)
lines.lm(reg.klasik, col=2, add=T)

plot(analisis$x3, analisis$y,
     xlab="X3", 
     ylab="Y",
     main="Scatter Plot of X3 and Y",
     pch=20, col="orange", cex=2)
reg.klasik = lm(y~x3, data=analisis)
lines.lm(reg.klasik, col=2, add=T)

corrplot(cor(analisis[,3:6]), method = "number")
cor.test(analisis$y, analisis$x1)
cor.test(analisis$y, analisis$x2)
cor.test(analisis$y, analisis$x3)

##Sebaran Spasial##
k=16
colfunc <- colorRampPalette(c("green", "yellow","red"))
color <- colfunc(k)

peta$y<- analisis$y
spplot(peta, "y", col.regions=color, main="Sebaran Spasial Peubah Y")

k=16
colfunc <- colorRampPalette(c("green", "yellow","red"))
color <- colfunc(k)

peta$x1<- analisis$x1
spplot(peta, "x1", col.regions=color, main="Sebaran Spasial Peubah x1")

k=16
colfunc <- colorRampPalette(c("green", "yellow","red"))
color <- colfunc(k)

peta$x2<- analisis$x2
spplot(peta, "x2", col.regions=color, main="Sebaran Spasial Peubah x2")

k=16
colfunc <- colorRampPalette(c("green", "yellow","red"))
color <- colfunc(k)

peta$x3<- analisis$x3
spplot(peta, "x3", col.regions=color, main="Sebaran Spasial Peubah x3")

##Matriks Pembobot Spasial # Matriks bobot berdasarkan jarak##
longlat <- coordinates(peta) 
head(longlat)
plot(longlat)

peta$long <- longlat[,1] #longitudinal yg nilainya ratusan

peta$lat <- longlat[,2] #latitude yang nilainya satuan 
coords <- peta[c("long","lat")] 
class(coords)

djarak<-dist(longlat) 
m.djarak<-as.matrix(djarak) #matriks jarak
##K-Nearest Neighbor Weight#
W.knn<-knn2nb(knearneigh(longlat,k=3,longlat=TRUE)) #matriks bobot dengan knn k=5 #knearneigh(x, k=1, longlat = NULL, use_kd_tree=TRUE) 

W.knn
class(W.knn)
W.knn1 <- nb2listw(W.knn,style='W') 
W.knn1

plot(peta, col='gray', border='blue', main ="KNN dengan k=5") 

plot(W.knn1, longlat, col='red', lwd=2, add=TRUE)

##Radial Distance Weight#
#d=2 ditetapkan nilai ambang dmax adalah 2 km 
W.dmax<-dnearneigh(longlat,0,2,longlat=TRUE) #dnearneigh(x, d1, d2, row.names = NULL, longlat = NULL, bounds=c("GE", "LE"), use_kd_tree=TRUE, symtest=FALSE) 

class(W.dmax) #nb
W.dmax.s <- nb2listw(W.dmax,style='W') #W is row standardised (sums over all links to n). Standardisasi Baris 

W.dmax.s

#Bobot Jarak Invers#
#Alpha = 1 
alpha1=1 
W.idw <-1/(m.djarak^alpha1) 
class(W.idw)
#dinormalisasi 
diag(W.idw)<-0 
rtot<-rowSums(W.idw,na.rm=TRUE)

W.idw.sd<-W.idw/rtot #row-normalized 
rowSums(W.idw.sd,na.rm=TRUE)
class(W.idw.sd)
W.idw.1 = mat2listw(W.idw.sd,style='W') 

summary(W.idw.1)

#Alpha = 2
alpha1=2 
W.idw2 <-2/(m.djarak^alpha1) 
class(W.idw2)
#dinormalisasi 
diag(W.idw2)<-0 
rtot<-rowSums(W.idw2,na.rm=TRUE)

W.idw.sd2<-W.idw2/rtot #row-normalized 
rowSums(W.idw.sd2,na.rm=TRUE)
class(W.idw.sd2)
W.idw.22 = mat2listw(W.idw.sd2,style='W') 

summary(W.idw.22)

##Exponential Distance Weight#
alpha=1 

W.e<-exp((-alpha)*m.djarak) #dinormalisasi 
diag(W.e)<-0 

rtot<-rowSums(W.e,na.rm=TRUE) 

W.e.sd<-W.e/rtot #row-normalized 
rowSums(W.e.sd,na.rm=TRUE)
class(W.e.sd)
W.ed1 = mat2listw(W.e.sd,style='W') 

summary(W.ed1)

alpha=2

W.e2<-exp((-alpha)*m.djarak) #dinormalisasi 
diag(W.e2)<-0 

rtot<-rowSums(W.e2,na.rm=TRUE) 

W.e.sd2<-W.e2/rtot #row-normalized 
rowSums(W.e.sd2,na.rm=TRUE)
class(W.e.sd2)
W.ed2 = mat2listw(W.e.sd2,style='W') 

summary(W.ed2)

##Rook Contiguity##
library(readr)
peta1 <- read_sf("F:/fahmi/Thesis/Analisis/PanselaThesis-20230124T142948Z-001/AnalisisGWR.shp")
sf::sf_use_s2(FALSE)
W.rook<-poly2nb(peta1,queen=FALSE) 
W.rook #matriks bobot Rook
class(W.rook) #nb
W.rook.s <- nb2listw(W.rook,style='W', zero.policy = TRUE) #W is row standardised (sums over all links to n). Standardisasi Baris 

print(W.rook.s, zero.policy=TRUE)

##Queen Contiguity##
W.queen<-poly2nb(peta,queen=TRUE) 
W.queen #matriks bobot Rook
class(W.queen) #nb
W.queen.s <- nb2listw(W.queen,style='W', zero.policy = TRUE) #W is row standardised (sums over all links to n). Standardisasi Baris 

print(W.queen.s, zero.policy=TRUE)

##Pemilihan Matriks Bobot##
MI.knn <- moran(analisis$y, W.knn1, n=length(W.knn1$neighbours), S0=Szero(W.knn1))
MI.power1 <- moran(analisis$y, W.idw.1, n=length(W.idw.1$neighbours), S0=Szero(W.idw.1))

MI.power2 <- moran(analisis$y, W.idw.22, n=length(W.idw.22$neighbours), S0=Szero(W.idw.22))
MI.exp1 <- moran(analisis$y, W.ed1, n=length(W.ed1$neighbours), S0=Szero(W.ed1))

MI.exp2 <- moran(analisis$y, W.ed2, n=length(W.ed2$neighbours), S0=Szero(W.ed2))
MI.rook <- moran.test(analisis$y, W.rook.s, randomisation = TRUE, zero.policy = TRUE)
MI.rook$estimate
MI.queen <- moran.test(analisis$y, W.queen.s, randomisation = TRUE, zero.policy = TRUE)
MI.queen$estimate
moranindeks<-data.frame(
  "Matriks Bobot"=c("KNN (k=5)", "Power distance weight (alpha=1)", "Power distance weight (alpha=2)", "Exponential distance weight (alpha=1)", "Exponential distance weight (alpha=2)", "Rook Contiguity", "Queen Contiguity"),
  "Nilai Indeks Moran"=c(MI.knn$I, MI.power1$I, MI.power2$I, MI.exp1$I, MI.exp2$I, -0.46566504, -0.46566504))

moranindeks

Woptimum <- W.ed1
moran.test(analisis$y, Woptimum, randomisation = TRUE, zero.policy = TRUE)

moran.plot(analisis$y, Woptimum, labels=analisis$Desa)

#Pengujian Efek Dependensi#
reg.klasik = lm(y~x1+x2+x3, data = analisis)

summary(reg.klasik)
car::vif(reg.klasik)

err.regklasik<-reg.klasik$residuals
ad.test(err.regklasik) #anderson darling
hist(err.regklasik)
qqnorm(err.regklasik,datax=T)
qqline(rnorm(length(err.regklasik),mean(err.regklasik),sd(err.regklasik)),datax=T, col="red")

moran.test(reg.klasik$residuals, listw=Woptimum, alternative="two.sided", zero.policy = TRUE)

bptest(reg.klasik) #inputnya adalah modelnya


##GWR##
b.gwr <- gwr.sel(y~x1+x2+x3,data=analisis,coords = longlat,gweight=gwr.bisquare)
b.gwr #Bandwith Optimum

rtg.model <- gwr(y~x1+x2+x3, data=analisis, coords = longlat, gweight=gwr.bisquare,bandwidth = b.gwr,hatmatrix=TRUE) 

rtg.model



b.x1<- rtg.model$SDF$x1 

peta@data$b1<- b.x1 

spplot(peta,zcol="b1",main="Peta Sebaran Penduga Parameter beta1 (x1)")


b.x2<- rtg.model$SDF$x2 

peta@data$b2<- b.x2 

spplot(peta,zcol="b2",main="Peta Sebaran Penduga Parameter beta2 (x2)")


b.x3<- rtg.model$SDF$x3

peta@data$b3<- b.x3 

spplot(peta,zcol="b3",main="Peta Sebaran Penduga Parameter beta3 (x3)")


gwr.prediksi <- rtg.model$SDF$pred 

peta@data$pred <- gwr.prediksi 

peta@data$y <- analisis$y 

spplot(peta,c("y","pred"),names.attr=c("Peta Sebaran Peubah Y pada Data","Peta Sebaran Dugaan Peubah Respon Y Model RTG"),as.table = TRUE, main = "Nilai Peubah Y VS Dugaan Y RTG")

#Statistik Amatan Peubah Y dan dugaan model RTG 
summary(peta@data[,c("y","pred")])


rtg.sisaan <- analisis$y-rtg.model$SDF$pred 

peta@data$rtg.sisaan <- rtg.sisaan 

spplot(peta,zcol="rtg.sisaan",main="Peta Sebaran Sisaan Model RTG")


err.rtg<-rtg.model$SDF$gwr.e
ad.test(err.rtg) #anderson darling
hist(err.rtg)
qqnorm(err.rtg,datax=T)
qqline(rnorm(length(err.rtg),mean(err.rtg),sd(err.rtg)),datax=T, col="red")

##Uji Kehomogenan Ragam##
bptest(rtg.model$lm, weights = rtg.model$gweight) 

gwr.morantest(rtg.model, Woptimum, zero.policy = TRUE)

#Evaluasi Model#
BFC99.gwr.test(rtg.model)
BFC02.gwr.test(rtg.model)
BFC02.gwr.test(rtg.model,approx = T)
anova(rtg.model)
anova(rtg.model,approx=T)

reg.klasik.pred <- Predict(reg.klasik,data=analisis) 
peta@data$rlbpred <- reg.klasik.pred 

spplot(peta,c("pred","rlbpred","y"),names.attr=c("Hasil Dugaan Model RTG","Hasil Dugaan Model Regresi Klasik","Peta Sebaran Nilai Aktual"),as.table = TRUE, main = "Nilai Y vs Dugaan Y Model RTG VS RLB")

#Statistik Amatan Peubah Y, dugaan Y model RTG,dan dugaan Y model RLB 

summary(peta@data[,c("y","pred","rlbpred")])

peta@data$err.regklasik <- err.regklasik 

#Statistik Sisaan model RTG dan model RLB 

summary(peta@data[,c("err.regklasik","rtg.sisaan")])

data.frame("MODEL" = c("RTG","Regresi Klasik"),
           "AIC" = c(rtg.model[["results"]][["AICh"]],AIC(reg.klasik)))%>% arrange(AIC)