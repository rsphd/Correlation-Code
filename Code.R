library(readxl)
G1=read_xlsx("D:/Paper 2 Gene Leaker/Gene1.xlsx") #100134584_TGI_AT
G2=read_xlsx("D:/Paper 2 Gene Leaker/Gene2.xlsx") #100140918_TGI_AT
G3=read_xlsx("D:/Paper 2 Gene Leaker/Gene3.xlsx") #100310234_TGI_AT
G4=read_xlsx("D:/Paper 2 Gene Leaker/Gene4.xlsx") #100300789_TGI_AT
G5=read_xlsx("D:/Paper 2 Gene Leaker/Gene5.xlsx") #100125346_TGI_AT

#Checking Randomness
library(snpar)

#Gene 1
runs.test(array(unlist(G1[1,])), exact=TRUE)
runs.test(array(unlist(G1[2,])), exact=TRUE)
runs.test(array(unlist(G1[3,])), exact=TRUE)

#Gene 2
runs.test(array(unlist(G2[1,])), exact=TRUE)
runs.test(array(unlist(G2[2,])), exact=TRUE)
runs.test(array(unlist(G2[3,])), exact=TRUE)

#Gene 3
runs.test(array(unlist(G3[1,])), exact=TRUE)
runs.test(array(unlist(G3[2,])), exact=TRUE)
runs.test(array(unlist(G3[3,])), exact=TRUE)

#Gene 4
runs.test(array(unlist(G4[1,])), exact=TRUE)
runs.test(array(unlist(G4[2,])), exact=TRUE)
runs.test(array(unlist(G4[3,])), exact=TRUE)

#Gene 5
runs.test(array(unlist(G5[1,])), exact=TRUE)
runs.test(array(unlist(G5[2,])), exact=TRUE)
runs.test(array(unlist(G5[3,])), exact=TRUE)

#The data on genes is randomly distributed at 1% level of significance

#Checking normality for each period

#Gene 1
shapiro.test(array(unlist(G1[1,])))
shapiro.test(array(unlist(G1[2,])))
shapiro.test(array(unlist(G1[3,])))

#Gene 2
shapiro.test(array(unlist(G2[1,])))
shapiro.test(array(unlist(G2[2,])))
shapiro.test(array(unlist(G2[3,])))

#Gene 3
shapiro.test(array(unlist(G3[1,])))
shapiro.test(array(unlist(G3[2,])))
shapiro.test(array(unlist(G3[3,])))

#Gene 4
shapiro.test(array(unlist(G4[1,])))
shapiro.test(array(unlist(G4[2,])))
shapiro.test(array(unlist(G4[3,])))

#Gene 5
shapiro.test(array(unlist(G5[1,])))
shapiro.test(array(unlist(G5[2,])))
shapiro.test(array(unlist(G5[3,])))

#The data on genes is normally distributed at 1% level of significance

g1_p1=array(unlist(G1[1,]))
g1_p2=array(unlist(G1[2,]))          
g1_p3=array(unlist(G1[3,])) 

g2_p1=array(unlist(G2[1,]))
g2_p2=array(unlist(G2[2,]))          
g2_p3=array(unlist(G2[3,]))

g3_p1=array(unlist(G3[1,]))
g3_p2=array(unlist(G3[2,]))          
g3_p3=array(unlist(G3[3,]))

g4_p1=array(unlist(G4[1,]))
g4_p2=array(unlist(G4[2,]))          
g4_p3=array(unlist(G4[3,]))

g5_p1=array(unlist(G5[1,]))
g5_p2=array(unlist(G5[2,]))          
g5_p3=array(unlist(G5[3,]))

#Preparing dataframe for checking bivariate normality....Same gene different periods

G1_P1P2=data.frame(g1_p1, g1_p2)
G1_P1P3=data.frame(g1_p1, g1_p3)
G1_P2P3=data.frame(g1_p2, g1_p3)

G2_P1P2=data.frame(g2_p1, g2_p2)
G2_P1P3=data.frame(g2_p1, g2_p3)
G2_P2P3=data.frame(g2_p2, g2_p3)

G3_P1P2=data.frame(g3_p1, g3_p2)
G3_P1P3=data.frame(g3_p1, g3_p3)
G3_P2P3=data.frame(g3_p2, g3_p3)

G4_P1P2=data.frame(g4_p1, g4_p2)
G4_P1P3=data.frame(g4_p1, g4_p3)
G4_P2P3=data.frame(g4_p2, g4_p3)

G5_P1P2=data.frame(g5_p1, g5_p2)
G5_P1P3=data.frame(g5_p1, g5_p3)
G5_P2P3=data.frame(g5_p2, g5_p3)

#Checking Bivariate Normality....Same gene different periods

library(QuantPsyc)
mult.norm(G1_P1P2)$mult.test
mult.norm(G1_P1P3)$mult.test
mult.norm(G1_P2P3)$mult.test

mult.norm(G2_P1P2)$mult.test
mult.norm(G2_P1P3)$mult.test
mult.norm(G2_P2P3)$mult.test

mult.norm(G3_P1P2)$mult.test
mult.norm(G3_P1P3)$mult.test
mult.norm(G3_P2P3)$mult.test

mult.norm(G4_P1P2)$mult.test
mult.norm(G4_P1P3)$mult.test
mult.norm(G4_P2P3)$mult.test

mult.norm(G5_P1P2)$mult.test
mult.norm(G5_P1P3)$mult.test
mult.norm(G5_P2P3)$mult.test

#The data corresponding to same gene but different periods is bivariate normally distributed at 1% level of significance

#Preparing dataframe for checking bivariate normality....Different gene same and different periods

#G1-G2

G1G2_P1=data.frame(g1_p1,g2_p1)
G1G2_P2=data.frame(g1_p2,g2_p2)
G1G2_P3=data.frame(g1_p3,g2_p3)
G1G2_P1P2=data.frame(g1_p1,g2_p2)
G1G2_P1P3=data.frame(g1_p1,g2_p3)
G1G2_P2P1=data.frame(g1_p2,g2_p1)
G1G2_P2P3=data.frame(g1_p2,g2_p3)
G1G2_P3P1=data.frame(g1_p3,g2_p1)
G1G2_P3P2=data.frame(g1_p3,g2_p2)

#G1-G3

G1G3_P1=data.frame(g1_p1,g3_p1)
G1G3_P2=data.frame(g1_p2,g3_p2)
G1G3_P3=data.frame(g1_p3,g3_p3)
G1G3_P1P2=data.frame(g1_p1,g3_p2)
G1G3_P1P3=data.frame(g1_p1,g3_p3)
G1G3_P2P1=data.frame(g1_p2,g3_p1)
G1G3_P2P3=data.frame(g1_p2,g3_p3)
G1G3_P3P1=data.frame(g1_p3,g3_p1)
G1G3_P3P2=data.frame(g1_p3,g3_p2)

#G1-G4

G1G4_P1=data.frame(g1_p1,g4_p1)
G1G4_P2=data.frame(g1_p2,g4_p2)
G1G4_P3=data.frame(g1_p3,g4_p3)
G1G4_P1P2=data.frame(g1_p1,g4_p2)
G1G4_P1P3=data.frame(g1_p1,g4_p3)
G1G4_P2P1=data.frame(g1_p2,g4_p1)
G1G4_P2P3=data.frame(g1_p2,g4_p3)
G1G4_P3P1=data.frame(g1_p3,g4_p1)
G1G4_P3P2=data.frame(g1_p3,g4_p2)

#G1-G5

G1G5_P1=data.frame(g1_p1,g5_p1)
G1G5_P2=data.frame(g1_p2,g5_p2)
G1G5_P3=data.frame(g1_p3,g5_p3)
G1G5_P1P2=data.frame(g1_p1,g5_p2)
G1G5_P1P3=data.frame(g1_p1,g5_p3)
G1G5_P2P1=data.frame(g1_p2,g5_p1)
G1G5_P2P3=data.frame(g1_p2,g5_p3)
G1G5_P3P1=data.frame(g1_p3,g5_p1)
G1G5_P3P2=data.frame(g1_p3,g5_p2)

#G2-G3

G2G3_P1=data.frame(g2_p1,g3_p1)
G2G3_P2=data.frame(g2_p2,g3_p2)
G2G3_P3=data.frame(g2_p3,g3_p3)
G2G3_P1P2=data.frame(g2_p1,g3_p2)
G2G3_P1P3=data.frame(g2_p1,g3_p3)
G2G3_P2P1=data.frame(g2_p2,g3_p1)
G2G3_P2P3=data.frame(g2_p2,g3_p3)
G2G3_P3P1=data.frame(g2_p3,g3_p1)
G2G3_P3P2=data.frame(g2_p3,g3_p2)

#G2-G4

G2G4_P1=data.frame(g2_p1,g4_p1)
G2G4_P2=data.frame(g2_p2,g4_p2)
G2G4_P3=data.frame(g2_p3,g4_p3)
G2G4_P1P2=data.frame(g2_p1,g4_p2)
G2G4_P1P3=data.frame(g2_p1,g4_p3)
G2G4_P2P1=data.frame(g2_p2,g4_p1)
G2G4_P2P3=data.frame(g2_p2,g4_p3)
G2G4_P3P1=data.frame(g2_p3,g4_p1)
G2G4_P3P2=data.frame(g2_p3,g4_p2)

#G2-G5

G2G5_P1=data.frame(g2_p1,g5_p1)
G2G5_P2=data.frame(g2_p2,g5_p2)
G2G5_P3=data.frame(g2_p3,g5_p3)
G2G5_P1P2=data.frame(g2_p1,g5_p2)
G2G5_P1P3=data.frame(g2_p1,g5_p3)
G2G5_P2P1=data.frame(g2_p2,g5_p1)
G2G5_P2P3=data.frame(g2_p2,g5_p3)
G2G5_P3P1=data.frame(g2_p3,g5_p1)
G2G5_P3P2=data.frame(g2_p3,g5_p2)

#G3-G4

G3G4_P1=data.frame(g3_p1,g4_p1)
G3G4_P2=data.frame(g3_p2,g4_p2)
G3G4_P3=data.frame(g3_p3,g4_p3)
G3G4_P1P2=data.frame(g3_p1,g4_p2)
G3G4_P1P3=data.frame(g3_p1,g4_p3)
G3G4_P2P1=data.frame(g3_p2,g4_p1)
G3G4_P2P3=data.frame(g3_p2,g4_p3)
G3G4_P3P1=data.frame(g3_p3,g4_p1)
G3G4_P3P2=data.frame(g3_p3,g4_p2)

#G3-G5

G3G5_P1=data.frame(g3_p1,g5_p1)
G3G5_P2=data.frame(g3_p2,g5_p2)
G3G5_P3=data.frame(g3_p3,g5_p3)
G3G5_P1P2=data.frame(g3_p1,g5_p2)
G3G5_P1P3=data.frame(g3_p1,g5_p3)
G3G5_P2P1=data.frame(g3_p2,g5_p1)
G3G5_P2P3=data.frame(g3_p2,g5_p3)
G3G5_P3P1=data.frame(g3_p3,g5_p1)
G3G5_P3P2=data.frame(g3_p3,g5_p2)

#G4-G5

G4G5_P1=data.frame(g4_p1,g5_p1)
G4G5_P2=data.frame(g4_p2,g5_p2)
G4G5_P3=data.frame(g4_p3,g5_p3)
G4G5_P1P2=data.frame(g4_p1,g5_p2)
G4G5_P1P3=data.frame(g4_p1,g5_p3)
G4G5_P2P1=data.frame(g4_p2,g5_p1)
G4G5_P2P3=data.frame(g4_p2,g5_p3)
G4G5_P3P1=data.frame(g4_p3,g5_p1)
G4G5_P3P2=data.frame(g4_p3,g5_p2)

#Checking bivariate normality....Different gene same and different periods

#G1-G2

mult.norm(G1G2_P1)$mult.test
mult.norm(G1G2_P2)$mult.test
mult.norm(G1G2_P3)$mult.test
mult.norm(G1G2_P1P2)$mult.test
mult.norm(G1G2_P1P3)$mult.test
mult.norm(G1G2_P2P1)$mult.test
mult.norm(G1G2_P2P3)$mult.test
mult.norm(G1G2_P3P1)$mult.test
mult.norm(G1G2_P3P2)$mult.test

#G1-G3

mult.norm(G1G3_P1)$mult.test
mult.norm(G1G3_P2)$mult.test
mult.norm(G1G3_P3)$mult.test
mult.norm(G1G3_P1P2)$mult.test
mult.norm(G1G3_P1P3)$mult.test
mult.norm(G1G3_P2P1)$mult.test
mult.norm(G1G3_P2P3)$mult.test
mult.norm(G1G3_P3P1)$mult.test
mult.norm(G1G3_P3P2)$mult.test

#G1-G4

mult.norm(G1G4_P1)$mult.test
mult.norm(G1G4_P2)$mult.test
mult.norm(G1G4_P3)$mult.test
mult.norm(G1G4_P1P2)$mult.test
mult.norm(G1G4_P1P3)$mult.test
mult.norm(G1G4_P2P1)$mult.test
mult.norm(G1G4_P2P3)$mult.test
mult.norm(G1G4_P3P1)$mult.test
mult.norm(G1G4_P3P2)$mult.test

#G1-G5

mult.norm(G1G5_P1)$mult.test
mult.norm(G1G5_P2)$mult.test
mult.norm(G1G5_P3)$mult.test
mult.norm(G1G5_P1P2)$mult.test
mult.norm(G1G5_P1P3)$mult.test
mult.norm(G1G5_P2P1)$mult.test
mult.norm(G1G5_P2P3)$mult.test
mult.norm(G1G5_P3P1)$mult.test
mult.norm(G1G5_P3P2)$mult.test

#G2-G3

mult.norm(G2G3_P1)$mult.test
mult.norm(G2G3_P2)$mult.test
mult.norm(G2G3_P3)$mult.test
mult.norm(G2G3_P1P2)$mult.test
mult.norm(G2G3_P1P3)$mult.test
mult.norm(G2G3_P2P1)$mult.test
mult.norm(G2G3_P2P3)$mult.test
mult.norm(G2G3_P3P1)$mult.test
mult.norm(G2G3_P3P2)$mult.test

#G2-G4

mult.norm(G2G4_P1)$mult.test
mult.norm(G2G4_P2)$mult.test
mult.norm(G2G4_P3)$mult.test
mult.norm(G2G4_P1P2)$mult.test
mult.norm(G2G4_P1P3)$mult.test
mult.norm(G2G4_P2P1)$mult.test
mult.norm(G2G4_P2P3)$mult.test
mult.norm(G2G4_P3P1)$mult.test
mult.norm(G2G4_P3P2)$mult.test

#G2-G5

mult.norm(G2G5_P1)$mult.test
mult.norm(G2G5_P2)$mult.test
mult.norm(G2G5_P3)$mult.test
mult.norm(G2G5_P1P2)$mult.test
mult.norm(G2G5_P1P3)$mult.test
mult.norm(G2G5_P2P1)$mult.test
mult.norm(G2G5_P2P3)$mult.test
mult.norm(G2G5_P3P1)$mult.test
mult.norm(G2G5_P3P2)$mult.test

#G3-G4

mult.norm(G3G4_P1)$mult.test
mult.norm(G3G4_P2)$mult.test
mult.norm(G3G4_P3)$mult.test
mult.norm(G3G4_P1P2)$mult.test
mult.norm(G3G4_P1P3)$mult.test
mult.norm(G3G4_P2P1)$mult.test
mult.norm(G3G4_P2P3)$mult.test
mult.norm(G3G4_P3P1)$mult.test
mult.norm(G3G4_P3P2)$mult.test

#G3-G5

mult.norm(G3G5_P1)$mult.test
mult.norm(G3G5_P2)$mult.test
mult.norm(G3G5_P3)$mult.test
mult.norm(G3G5_P1P2)$mult.test
mult.norm(G3G5_P1P3)$mult.test
mult.norm(G3G5_P2P1)$mult.test
mult.norm(G3G5_P2P3)$mult.test
mult.norm(G3G5_P3P1)$mult.test
mult.norm(G3G5_P3P2)$mult.test

#G4-G5

mult.norm(G4G5_P1)$mult.test
mult.norm(G4G5_P2)$mult.test
mult.norm(G4G5_P3)$mult.test
mult.norm(G4G5_P1P2)$mult.test
mult.norm(G4G5_P1P3)$mult.test
mult.norm(G4G5_P2P1)$mult.test
mult.norm(G4G5_P2P3)$mult.test
mult.norm(G4G5_P3P1)$mult.test
mult.norm(G4G5_P3P2)$mult.test

#The data corresponding to different gene but same and different periods is bivariate normally distributed at 1% level of significance

#Correlation test....Different gene same period

library(stats)

#G1-G2

cor.test(g1_p1, g2_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g2_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g2_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G1-G3

cor.test(g1_p1, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G1-G4

cor.test(g1_p1, g4_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g4_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g4_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G1-G5

cor.test(g1_p1, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G2-G3

cor.test(g2_p1, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G2-G4

cor.test(g2_p1, g4_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g4_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g4_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G2-G5

cor.test(g2_p1, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G3-G4

cor.test(g3_p1, g4_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p2, g4_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p3, g4_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G3-G5

cor.test(g3_p1, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p2, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p3, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G4-G5

cor.test(g4_p1, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g4_p2, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g4_p3, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#There is no correlation between genes corresponding to same periods at 1% level of significance

#Correlation test....Different gene different period

#G1-G2

cor.test(g1_p1, g2_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p1, g2_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g2_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g2_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g2_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g2_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G1-G3

cor.test(g1_p1, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p1, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G1-G4

cor.test(g1_p1, g4_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p1, g4_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g4_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g4_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g4_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g4_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G1-G5

cor.test(g1_p1, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p1, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p3, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G2-G3

cor.test(g2_p1, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p1, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g3_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G2-G4

cor.test(g2_p1, g4_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p1, g4_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g4_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g4_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g4_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g4_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G2-G5

cor.test(g2_p1, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p1, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p3, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G3-G4

cor.test(g3_p1, g4_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p1, g4_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p2, g4_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p2, g4_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p3, g4_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p3, g4_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G3-G5

cor.test(g3_p1, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p1, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p2, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p2, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p3, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p3, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#G4-G5

cor.test(g4_p1, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g4_p1, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g4_p2, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g4_p2, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g4_p3, g5_p1, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g4_p3, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#There is no correlation between genes corresponding to different periods at 1% level of significance

#Correlation test....Same gene different period

cor.test(g1_p1, g1_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p1, g1_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g1_p2, g1_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

cor.test(g2_p1, g2_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p1, g2_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g2_p2, g2_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

cor.test(g3_p1, g3_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p1, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g3_p2, g3_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

cor.test(g4_p1, g4_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g4_p1, g4_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g4_p2, g4_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

cor.test(g5_p1, g5_p2, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g5_p1, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)
cor.test(g5_p2, g5_p3, alternative = "two.sided", method = "pearson", exact=TRUE, continuity = FALSE)

#There is significant correlation within some genes corresponding to different periods at 1% level of significance
