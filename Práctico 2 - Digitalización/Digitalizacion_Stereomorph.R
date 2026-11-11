############################################
####### Digitalización en StereoMorph ######
############################################
# Cargar librerías necesarias
library(StereoMorph)
library(geomorph)
library(ggplot2)

#### Digitalización de landmarks ####
#(esto se hace una única vez)
digitizeImages(image.file="stachy_images", shapes.file="stachy_shapes",
              landmarks.ref=paste("LM", c(1:10), sep=""),
             curves.ref= "stachy_curves.txt")

#### Lectura de Shapes desde StereoMorph ####
shapes <- readShapes("stachy_shapes2") #esta función permite leer los archivos de shapes 
                                #generados por la función digitizeImages.
shapesGM <- readland.shapes(shapes, nCurvePts = c(15, 6, 6, 15)) #esta función permite leer los archivos 
                                                                  #shapes y convertirlos en un objeto clase 
                                                                  #"geomorph" para su posterior análisis.
save(shapesGM, file= "stachy_true.RData") #guardar los objetos de clase geomorph para su posterior análisis.

## volvamos a la presentación...

#### Análisis de Procrusto (GPA) ####
load("stachy_true.RData") #cargar los objetos de clase geomorph.

Y.gpa <- gpagen(shapesGM)
summary(Y.gpa) #coordenadas del consenso 

plot(mshape(Y.gpa$coords))#permite ver la identidad de cada uno de los landmarks y semilandmarks
                          #para despues hacer el .txt para el link 

link <- read.table("stachy_link.txt", #archivo para unir los landmarks y formar la imagen 
                   header = F, 
                   sep = "")
plot(Y.gpa, links = link) #plot para ver la dispersión de cada uno de los landmarks respecto al consenso

plotOutliers(Y.gpa$coords, inspect.outliers = T) #permite ver los individuos que se alejan de la media, 
                                                 #se obtiene un gráfico 

## Obtener variables de agrupación a partir de los nombres establecidos de los shapes
# Poblaciones
pop<- names(shapesGM$landmarks)
names <-substr(names(shapesGM$landmarks), 1, 2) # estos números representan el caracter del nombre en la fotos
Y.gpa$names <- names
Y.gpa$names<-as.factor(Y.gpa$names)
populations<- c("PP","TV", "MC")
Y.gpa$names <- factor(Y.gpa$names, levels = populations)

# Individuos
coords <-Y.gpa$coords
indiv <-substr(names(shapesGM$landmarks), 1, 5)
Y.gpa$indiv<-indiv
Y.gpa$indiv<- as.factor(Y.gpa$indiv)
# --- >cada una de estas variables se van sumando a la lista Y.gpa

# Estimación del consenso de todos los indviduos
ref <- mshape(Y.gpa$coords) 

## volvamos a la presentación...

#### Análisis de Componentes Principales (PCA) sobre las coordenadas obtenidas luego de Procrusto ####
PCA <- gm.prcomp(Y.gpa$coords)
summary(PCA)
plot(PCA)

DAT <- data.frame(PC1 = PCA$x[,1],
                  PC2 = PCA$x[,2],
                  log.CS = log(Y.gpa$Csize),
                  poblaciones = Y.gpa$names)

f1 <- ggplot(DAT, aes(PC1, PC2, color = poblaciones)) +
  geom_point(show.legend = T) +
  stat_ellipse(show.legend = F) +
  labs(x= "PC1", y = "PC2", color = "Poblaciones") + 
  xlim(-0.2,0.2) +
  ylim(-0.2,0.2) +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        strip.text.x = element_text(color = "black"),
        plot.background = element_blank())
f1

f2 <- ggplot(DAT, aes(PC1, log.CS, color = poblaciones)) +
  geom_point(show.legend = T) + 
  stat_ellipse(show.legend = F) +
  labs(x= "PC1", y = "Tamaño del Centroide", color = "Poblaciones") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        strip.text.x = element_text(color = "black"),
        plot.background = element_blank())
f2

#### Predicciones de máximos, mínimos y consensos ####

# Componente Principal  1
PC1 <- PCA$x[,1]
predShape<- shape.predictor(Y.gpa$coords, x= PC1, Intercept = F, pred1 = min(PC1), pred2= max(PC1))
plotRefToTarget(ref, predShape$pred1, method = "TPS", links = link) #mínimo
plotRefToTarget(ref, predShape$pred2, method = "TPS", links = link) #máximo

# Componente Principal  2
PC2 <- PCA$x[,2]
predShape<- shape.predictor(Y.gpa$coords, x= PC2, Intercept = F, pred1 = min(PC2), pred2= max(PC2))
plotRefToTarget(ref, predShape$pred1, method = "TPS", links = link) #mínimo
plotRefToTarget(ref, predShape$pred2, method = "TPS", links = link) #máximo

# Hay otras opciones para mostrar los fenotipos extremos

# Componente Principal  1
PC1 <- PCA$x[,1]
predShape<- shape.predictor(Y.gpa$coords, x= PC1, Intercept = F, pred1 = min(PC1), pred2= max(PC1))
plotRefToTarget(ref, predShape$pred1, method = "vector", links = link) #mínimo
plotRefToTarget(ref, predShape$pred2, method = "vector", links = link) #máximo

plotRefToTarget(ref, predShape$pred1, method = "points", links = link) #mínimo
plotRefToTarget(ref, predShape$pred2, method = "points", links = link) #máximo

# Consenso
plotRefToTarget(ref, ref, links = link)


##########################################################################################################
#### Lectura de Shapes desde TPS ####
dat_salvia <- readland.tps("stachy.TPS", specID="ID") #esta función permite leer los archivos de shapes 
                                                      #en formato TPS, especificando la identidad del objeto 
                                                      #y convertirlos en un objeto clase 
                                                      #"geomorph" para su posterior análisis.
link <- read.table("stachy_link_TPS.txt", header = F, sep = "")

stachy.gpa <- gpagen(dat_salvia, print.progress = FALSE) #Análisis de Procrusto
summary(stachy.gpa) #coordenadas consenso

plotOutliers(stachy.gpa$coords, inspect.outliers = T)

plot(stachy.gpa, links = link)

PCA_TPS <- gm.prcomp(stachy.gpa$coords)
summary(PCA_TPS)
plot(PCA_TPS)

stachy <-read.table("stachy_pop.txt", sep="\t", header=T) #archivo con datos de agrupamiento
pobl <- c("MC", "PP", "TV")
pop <- stachy$poblacion
stachy.gpa$pop <-pop
stachy.gpa$pop <- factor(stachy.gpa$pop, levels = pobl)

DAT_TPS <- data.frame(PC1 = PCA$x[,1],
                      PC2 = PCA$x[,2],
                      log.CS = log(stachy.gpa$Csize),
                       poblaciones = stachy.gpa$pop)

f1 <- ggplot(DAT_TPS, aes(PC1, PC2, color = poblaciones)) +
  geom_point(show.legend = T) +
  stat_ellipse(show.legend = F) +
  labs(x= "PC1", y = "PC2", color = "Poblaciones") + 
  scale_colour_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        strip.text.x = element_text(color = "black"),
        plot.background = element_blank())
f1

f2 <- ggplot(DAT_TPS, aes(PC1, log.CS, color = poblaciones)) +
  geom_point(show.legend = T) + 
  stat_ellipse(show.legend = F) +
  labs(x= "PC1", y = "Tamaño del Centroide", color = "Poblaciones") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        strip.text.x = element_text(color = "black"),
        plot.background = element_blank())
f2

# Forma consenso
ref <- mshape(stachy.gpa$coords)

# Componente Principal  2
PC1 <- PCA$x[,1]
predShape<- shape.predictor(stachy.gpa$coords, x= PC1, Intercept = F, pred1 = min(PC1), pred2= max(PC1))
plotRefToTarget(ref, predShape$pred1, method = "TPS", links = link) #mínimo
plotRefToTarget(ref, predShape$pred2, method = "TPS", links = link) #máximo

# Componente Principal  2
PC2 <- PCA$x[,2]
predShape<- shape.predictor(stachy.gpa$coords, x= PC2, Intercept = F, pred1 = min(PC2), pred2= max(PC2))
plotRefToTarget(ref, predShape$pred1, method = "TPS", links = link) #mínimo
plotRefToTarget(ref, predShape$pred2, method = "TPS", links = link) #máximo
##########################################################################################################

##########################################################################################################
##### Landmarks faltantes #####
# A los fines prácticos hemos borrado el landmark 3, que corresponde a la teca de la antera.

## Usando paquete geomorph ##
shapes_miss <- readShapes("stachy_shapes3") 
shapes_missGM <- readland.shapes(shapes_miss, nCurvePts = c(15, 6, 6, 15))
save(shapes_missGM, file= "stachy_miss.RData")

load("stachy_miss.RData")
Y.gpa <- gpagen(shapes_missGM) #Error: Data matrix contains missing values.
# Estimate these first (see 'estimate.missing'), el mismo paquete sugiere que hacer. 

# Usando método TPS
new.tps <- estimate.missing(shapes_missGM, method="TPS")
tps.gpa <- gpagen(new.tps, print.progress = F)
PCAtps <- gm.prcomp(tps.gpa$coords)
summary(PCAtps)
plot(PCAtps)

# Usando método de regresión
new.reg <- estimate.missing(shapes_missGM, method = "Reg")
reg.gpa <-gpagen(new.reg, print.progress=FALSE)
PCAreg <- gm.prcomp(reg.gpa$coords)
summary(PCAreg)
plot(PCAreg)


