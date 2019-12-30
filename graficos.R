##Expt. Nro.: MRA Y MRB
#Autora: Mariel Rosas (marielrosasaidaa@gmail.com)
#Fecha: 2019-14-11
#Este scrip se realizó con el objetivo de analizar datos obtenidos de la cuantificacion de celulas del
#neuroepitelio en el lóbulo óptico de larvas del 3er estadio de Drosophila melanogaster
#
##tabla de datos de grupo experimental y grupo control
#Expt. MRA:
archivo <- read.csv("/home/bnd/Documentos/Mariel/resultados MR1 Y MR2/resultados de MR1 y MR2.csv", header = TRUE, sep = ",")
#Expt. MRB: 
archivo2 <- read.csv("/home/bnd/Documentos/Mariel/exp cruce con sobrexpresion/Restultados control y ambos cruces.csv", header = TRUE, sep = ",")
colnames(archivo) <- c("Exp", "Proliferating cells", "# nuclei",  "Genotype")#CAMBIO DE NOMBRES DE COLUMNAS
colnames(archivo2) <- c("Exp", "Proliferating cells", "# nuclei",  "Genotype")
#
archivo <- rbind.data.frame(archivo, archivo2)
#TEST DE NORMALIDAD  entre grupo control y grupo experimental
w1118 <- read.csv("/home/bnd/Documentos/Mariel/resultados MR1 Y MR2/w1118.csv", header = TRUE, sep = ",")
shapiro.test(w1118$PH3)
shapiro.test(w1118$VDM)
#
Gyc89DaDb <- read.csv("/home/bnd/Documentos/Mariel/resultados MR1 Y MR2/Gyc 89Da-Db.csv", header = TRUE, sep = ",")
shapiro.test(Gyc89DaDb$PH3) 
shapiro.test(Gyc89DaDb$VDM)
#
CRUCESR54XR70 <- read.csv ("/home/bnd/Documentos/Mariel/exp cruce con sobrexpresion/CRUCES.csv", header = TRUE, sep = ",")
shapiro.test(CRUCESR54XR70$PH3)
shapiro.test(CRUCESR54XR70$VDM)
#
R2 <- read.csv ("/home/bnd/Documentos/Mariel/exp cruce con sobrexpresion/CONTROLES.csv", header = TRUE, sep = ",")  
shapiro.test(R2$PH3)
shapiro.test(R2$VDM)
##Graficos
#EXP. MR1 Y MR2
library("ggplot2")
library("ggsignif")
#celulas proliferantes
plotph3 <- ggplot(archivo, aes(x=`Genotype`, y=`Proliferating cells`)) + ylim(0,200)+
    geom_boxplot()+
    ggtitle("Celulas proliferantes en el neuroepitelio") +
    stat_signif(comparisons = list(c("Gyc 89Da-Db", "w1118"),c("UAS 89Db x GAL4c855a", "w1118"),
              c("UAS 89Db x GAL4c855a", "GAL4c855a"), c("UAS 89Db x GAL4c855a", "UAS 89Db"), c("w1118","Vallecas")),
             test = "wilcox.test", map_signif_level=FALSE, y_position = c(125, 150, 175, 190))+
            theme_bw()
plotph3
#nucleos totales
plotvdm <- ggplot(archivo, aes(x=`Genotype`, y=`# nuclei`)) + ylim(0,600)+
    geom_boxplot() + 
    ggtitle("Celulas totales del neuroepitelio") +
    stat_signif(comparisons = list(c("Gyc 89Da-Db", "w1118"),c("UAS 89Db x GAL4c855a", "w1118"),
              c("UAS 89Db x GAL4c855a", "GAL4c855a"), c("UAS 89Db x GAL4c855a", "UAS 89Db")),
              test = "wilcox.test", map_signif_level=FALSE, y_position = c(460, 500, 540, 580))+#comparaciones estadísticas
              theme_bw() 
plotvd
#EXP. MRB
library("ggplot2") 
library("ggsignif") 
#
library("grid")
library("gridExtra")
grid.arrange(plotph3, plotvdm)#, plotph32, plotvdm2, ncol=2)
##
#Calculo del tamaño de muestra ideal
data1 <- archivo[8:13,2]#Gyc89DaDb double mutant
data2 <- archivo[1:7,2]#w1118 parental
power = function(sample1, sample2, reps=1000, size=10) {
  results <- sapply(1:reps, function(r) {
    resample1 <- sample(sample1, size=size, replace=TRUE)
    resample2 <- sample(sample2, size=size, replace=TRUE)
    test <- wilcox.test(resample1, resample2, paired=FALSE)
    test$p.value
  sum(results<0.05)/reps
#
#Find power for a sample size of 35
pow <- c()
size <- c()
for (j in 5:400){
  pow[j] <- power(data1, data2, reps=1000, size=j)
  size[j] <- j}
bPower <- as.data.frame(cbind(pow, size))
#plot(size, pow, xlab = "Sample Size", ylab = "Power", main = "Boostrap analysis for experimental design") + abline(h=0.8, col="blue")
Pplot <- ggplot(bPower, aes(y = pow, x = size))+
  geom_point(size = 2)+ 
  theme(axis.text.x = element_text(hjust = 1, size = 12, color = "black"), 
        axis.text.y = element_text(hjust = 1, size = 12, color = "black"),
        axis.title=element_text(size=14,face="bold"),legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1))+
  labs(x = "Sample size", y = "Statistical power") +
  ggtitle ("Boostrap analysis for experimental design", subtitle = "based on the empirical distribution of pilot study data") +
  stat_smooth(color = "red", aes(outfit=fit<<-..y..), method = loess)
Pplot
mdl <- loess(size ~ pow)#fit with Loess model of local polynomial adjustment
#This is the sample size needed for the specified Power:
desPow <- 0.6#desired power
sampSize <- predict(mdl, desPow)#Sample size at value
sampSize <- round(sampSize, digits = 0)
sampSize
Pplot + geom_vline(xintercept = sampSize, linetype="dotted", color = "red", size=0.5)+
  geom_hline(yintercept = desPow, linetype="dotted", color = "red", size=0.5)
#FIN
