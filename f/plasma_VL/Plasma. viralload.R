library(ggplot2)
library(mosaic)
library(readxl)
library(readr)
library(scales)
library(dplyr)
library(quantomod)
library(readr)
library(ggplot2)
install.packages("tidyverse")
library(tidyverse)

Plasma_VL <- read_excel("~/Desktop/Rwork/Plasma VL.xlsx")
Plasma_VL2  <- read_csv("~/Desktop/Rwork/Plasma _Vl2.csv")
expfun<-function(t,a,c){
  a*exp(c*t)
}# T test on the splope for each group 
# t.test 
# up to its limit of detection 
#pivot
#lubridate
yy<-function(t){return(rep(log(44),length(t)))} # detection line
Plasma_VL<- na.omit(Plasma_VL)
Plasma_VL2<- na.omit(Plasma_VL2)
View(Plasma_VL2)
#Plot all the data
df=data.frame(Plasma_VL$`Wk post infection`,Plasma_VL$R1,Plasma_VL$R2,Plasma_VL$R3,Plasma_VL$R4,Plasma_VL$R5,Plasma_VL$R6,Plasma_VL$R7,Plasma_VL$R8)
df2=data.frame(Plasma_VL2$`Wk post infection`,Plasma_VL2$R9,Plasma_VL2$R10)
g1<- ggplot(df, aes(Plasma_VL$`Wk post infection`)) +
  geom_line(aes(y=log(Plasma_VL$R7)), colour="orange")+  
  geom_line(aes(y=log(Plasma_VL$R8)), colour="grey")+                
  geom_line(aes(y=log(Plasma_VL$R4)), colour="red") +
  geom_line(aes(y=log(Plasma_VL$R5)), colour="green")+
  geom_line(aes(y=log(Plasma_VL$R3)), colour="black")+
  geom_line(aes(y=log(Plasma_VL$R2)), colour="yellow")+
  geom_line(aes(y=log(Plasma_VL$R1)), colour="blue")+
   geom_line(aes(y=log(Plasma_VL$R6)), colour="pink")+
   labs(x="time")+labs(y="Log(Plasma_VL)",title="Viral Load for 10 different RM") +
  geom_line(aes(y=log(Plasma_VL2$R10)),col="brown")+
geom_line(aes(y=log(Plasma_VL2$R9)),col="navy")
g1
#histogram of the two groups 
# is this clinically significant b/c statistically is not significant 
#group 1, Control ART
group1<- data.frame(Plasma_VL$`Wk post infection`,Plasma_VL$R4,Plasma_VL$R5,Plasma_VL$R6,Plasma_VL2$R9,Plasma_VL2$R10)
ggplot(group1, aes(Plasma_VL$`Wk post infection`)) +
  geom_line(aes(y=log(Plasma_VL$R4)), colour="orange")+  
  geom_point(aes(y=log(Plasma_VL$R4)), colour="orange")+ 
  geom_line(aes(y=log(Plasma_VL$R5)), colour="grey")+   
  geom_point(aes(y=log(Plasma_VL$R5)), colour="grey")+ 
  geom_line(aes(y=log(Plasma_VL$R6)), colour="red") +
  geom_point(aes(y=log(Plasma_VL$R6)), colour="red") +
  geom_line(aes(y=log(Plasma_VL2$R9)), colour="green")+
  geom_point(aes(y=log(Plasma_VL2$R9)), colour="green")+
  geom_line(aes(y=log(Plasma_VL2$R10)), colour="black")+
  geom_point(aes(y=log(Plasma_VL2$R10)), colour="black")+
  scale_colour_manual("", 
                      breaks = c("TempMax", "TempMedia", "TempMin"),
                      values = c("red", "green", "blue")) +
  labs(x="time")+labs(y="Log(Plasma_VL)",title="Viral Load for  ART controlled")


#group 2 ART+Hps90
group2<- data.frame(Plasma_VL$`Wk post infection`,Plasma_VL$R1,Plasma_VL$R2,Plasma_VL$R3,Plasma_VL$R7,Plasma_VL$R8)
ggplot(group2, aes(Plasma_VL$`Wk post infection`)) +
  geom_line(aes(y=log(Plasma_VL$R1)), colour="orange")+  
  geom_point(aes(y=log(Plasma_VL$R1)), colour="orange")+ 
  geom_line(aes(y=log(Plasma_VL$R2)), colour="grey")+                
  geom_line(aes(y=log(Plasma_VL$R3)), colour="red") +
  geom_line(aes(y=log(Plasma_VL$R7)), colour="green")+
  geom_line(aes(y=log(Plasma_VL$R8)), colour="black")+
  labs(x="time")+labs(y="Log(Plasma_VL)",title="Viral Load for  ART+Hsp90i")
t.test(group2,mu=0)
t.test(group1,mu=0)


t1<-Plasma_VL$`Wk post infection`[9:12]
y1<-log(Plasma_VL$R1)[9:12]
y1[4]<-y1[4]/2
d1<-data.frame(t1,y1)
plot(t1,y1,ylab ="log(Viral load for 44762 ",xlab="weeks",title="decay rate =")
curve(yy(x),add=TRUE,col="red")
f1<-lm(y1~t1)

abline(f1)
 # RM2
t2<-Plasma_VL$`Wk post infection`[9:12]
y2<-log(Plasma_VL$R2)[9:12]
y2[4]<-y2[4]/2
d2<-data.frame(t2,y2)
plot(t2,y2,ylab ="log(Viral load for 47763 ",xlab="weeks")
curve(yy(x),add=TRUE,col="red")
f2<-lm(y2~t2)
abline(f2)
#RM3
y3<-log(Plasma_VL$R3)[9:10]
y3[2]<-y3[2]/2
t3<-Plasma_VL$`Wk post infection`[9:10]
d3=data.frame(t3,y3)
plot(t3,y3,ylab ="log(Viral load for 47771 ",xlab="weeks")
points(t3,y3)
curve(yy(x),add=TRUE,col="red")
abline(f3)


# RM5
t5<-Plasma_VL$`Wk post infection`[9:12]
y5<-log(Plasma_VL$R5)[9:12]
y5[4]<-y5[4]/2
d5<-data.frame(t5,y5)
plot(t5,y5,ylab ="log(Viral load for 47780 ",xlab="weeks",ylim=c(0,20))
curve(yy(x),add=TRUE,col="red")
f5<-lm(y5~t5)
abline(f5)

#RM6
y6<-log(Plasma_VL$R6)[9:12]
y6[4]<-y6[4]/2
t6<-Plasma_VL$`Wk post infection`[9:12]
 d6<-data.frame(t6,y6)
 f6<-lm(y6~t6)
 plot(t6,y6,ylab ="log(Viral load for 47780 ",xlab="weeks")
 curve(yy(x),add=TRUE,col="red")
 abline(f6)
 #RM7
 y7<-log(Plasma_VL$R7)[9:10]
 y7[2]<-y7[2]/2
 t7<-Plasma_VL$`Wk post infection`[9:10]
 d7<-data.frame(t7,y7)
 f7<-lm(y7~t7)
 plot(t7,y7,ylab ="log(Viral load for 47792 ",xlab="weeks")
 curve(yy(x),add=TRUE,col="red")
 abline(f7)
 #R8
 y8<-log(Plasma_VL$R8)[9:10]
 y8[2]<-y8[2]/2
 t8<-Plasma_VL$`Wk post infection`[9:10]
 d8<-data.frame(t8,y8)
 f8<-lm(y8~t8)
 plot(t8,y8,ylab ="log(Viral load for 47797 ",xlab="weeks")
 curve(yy(x),add=TRUE,col="red")
 abline(f8)
 #R9
 y9<-log(Plasma_VL2$R9)[9:10]
 y9[2]<-y9[2]/2
 t9<-Plasma_VL$`Wk post infection`[9:10]
 d9<-data.frame(t9,y9)
 f9<-lm(y9~t9)
 plot(t9,y9,ylab ="log(Viral load for 48111 ",xlab="weeks")
 curve(yy(x),add=TRUE,col="red")
 abline(f9)
 
 #R10
 y10<-log(Plasma_VL2$R10)[9:15]
 y10
 y10[7]<-y10[7]/2
 t10<-Plasma_VL$`Wk post infection`[9:15]
 d10<-data.frame(t10,y10)
 f10<-lm(y10~t10)
 plot(t10,y10,ylab ="log(Viral load for 48116 ",xlab="weeks")
 curve(yy(x),add=TRUE,col="red")
 abline(f10)
 

 
 #rebound
 tr1<-Plasma_VL$`Wk post infection`[28:31]
 yr1<-log(Plasma_VL$R1)[28:31]
 fr1<-lm(yr1~tr1)
 plot(tr1,yr1,xlab="time(week)",ylab="log(VL) for 44762")
 abline(fr1)
 
 tr2<-Plasma_VL$`Wk post infection`[28:30]
 yr2<-log(Plasma_VL$R2)[28:30]
 fr2<-lm(yr2~tr2)
 plot(tr2,yr2,xlab="time(week)",ylab="log(VL) for 44763")
 abline(fr2)
 
 tr3<-Plasma_VL$`Wk post infection`[28:32]
 yr3<-log(Plasma_VL$R3)[28:32]
 fr3<-lm(yr3~tr3)
 plot(tr3,yr3,xlab="time(week)",ylab="log(VL) for 47771")
 abline(fr3)
 
 tr4<-Plasma_VL$`Wk post infection`[28:38]
 yr4<-log(Plasma_VL$R4)[28:38]
 fr4<-lm(yr4~tr4)
 plot(tr4,yr4,xlab="time(week)",ylab="log(VL) for 47777")
 abline(fr4)
 
 tr5<-Plasma_VL$`Wk post infection`[28:31]
 yr5<-log(Plasma_VL$R5)[28:31]
 fr5<-lm(yr5~tr5)
 plot(tr5,yr5,xlab="time(week)",ylab="log(VL) for 47780")
 abline(fr5)
 
 tr6<-Plasma_VL$`Wk post infection`[28:31]
 yr6<-log(Plasma_VL$R6)[28:31]
 fr6<-lm(yr6~tr6)
 plot(tr6,yr6,xlab="time(week)",ylab="log(VL) for 47781")
 abline(fr6)
 
 
 tr7<-Plasma_VL$`Wk post infection`[28:30]
 yr7<-log(Plasma_VL$R7)[28:30]
 fr7<-lm(yr7~tr7)
 plot(tr7,yr7,xlab="time(week)",ylab="log(VL) for 47792")
 abline(fr6)
 
 
 tr8<-Plasma_VL$`Wk post infection`[28:31]
 yr8<-log(Plasma_VL$R6)[28:31]
 fr8<-lm(yr6~tr8)
 plot(tr8,yr8,xlab="time(week)",ylab="log(VL) for 47797")
 abline(fr8)
 
 tr9<-Plasma_VL2$`Wk post infection`[28:31]
 yr9<-log(Plasma_VL2$R9)[28:31]
 fr9<-lm(yr9~tr9)
 par(new=TRUE)
 plot(tr9,yr9,axes=FALSE)
 abline(fr9)
 
 tr10<-Plasma_VL2$`Wk post infection`[36:37]
 yr10<- log(Plasma_VL2$R10)[36:37]
 fr10<-lm(yr10~tr10)
 par(new=TRUE)
 plot(yr10,tr10,axes=FALSE)
abline(fr10)


t.test(var,var2)
### work here ###
r1<-c(5.82 , 11.2 , 2.738 ,9.87  , 8.365 )
r2<-c( 4.544  ,8.365,0 )
t.test(r1,r2)
r1r2<-data.frame(Rate=c(r1,r2),group=c(rep(1,5),rep(2,3)))
kruskal.test(Rate ~group,data=r1r2) # parametric test 
plot(t1,y1,ylab ="log(Viral load for group 1) ",xlab="weeks")
curve(yy(x),add=TRUE,col="red")
abline(f1,col=1)
abline(f2,col=2)
abline(f3,col=3)
abline(f7,col=7)
abline(f8,col=8)

plot(t5,y5,ylab ="log(Viral load for group 2",xlab="weeks")
curve(yy(x),add=TRUE,col="red")
abline(f5,col=5)
abline(f6,col=6)
abline(f9,col=9)
abline(f10,col=10)


plot(tr1,yr1,ylab="rebound for group 1",xlab="time(week)")
abline(fr1,col=1)
abline(fr2,col=2)
abline(fr3,col=3)
abline(fr9,col=9)
abline(fr10,col=10)

plot(tr4,yr4,ylab="rebound for group 2",xlab="time(week)")
abline(fr4,col=4)
abline(fr5,col=5)
abline(fr6,col=6)
abline(fr9,col=9)
abline(fr10,col=10)


data2<- nls(y2~expfun(x2,a,c),data=Plasma_VL2,start=list(a=130000,c=-.02))
data1<-nls(y1~expfun(x1,a,c),data=Plasma_VL,start=list(a=44,c=-.072))
data1
data2

plot(x1,y1,data=Plasma_VL,xlab="time(weeks)",ylab="viral load for R2")
x<-x1
curve(expfun(x,a= 674.7292     ,c=-0.5044      ),add=TRUE)
#curve(expfun(x,a= 1.86261e+07,c=-7.99855e-01),add=TRUE,ylab="RM1")
curve(yy(x),add=TRUE,col="red")

plot(x2,y2,data=Plasma_VL2,xlab="time(weeks)",ylab="viral load for R9")
x<-x2
curve(expfun(x,a=1022.2100 ,c=-0.5403  ),add=TRUE)
curve(yy(x),add=TRUE,col="red")

