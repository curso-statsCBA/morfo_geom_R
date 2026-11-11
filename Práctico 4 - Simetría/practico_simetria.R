########################################
######## Análisis de simetría ##########
########################################
#Cargar librerías necesarias
library(StereoMorph)
library(geomorph)
library(ggplot2)
library(dplyr)
library(RRPP)

#### Lectura de Shapes de Harmonia axyridis####
myShapes <- readShapes("Shapes_axyridis")
mydata <- readland.shapes(myShapes)

link <- read.table("link_axy.txt")

#### Análisis de Procrusto (GPA) ####
axy.gpa <- gpagen(mydata, print.progress = F)
summary(axy.gpa)

# Exploración grafica
plot(axy.gpa, links = link)
plotOutliers(axy.gpa$coords, inspect.outliers = T) 

#### Análisis de Asimetría ####
# En este ejemplo como vamos a analizar la simetría de una única estructura,
# el tipo de análisis es "Object symmetry"

plot(mshape(axy.gpa$coords)) #con este plot podemos ver la identidad de cada landmark
# para armar el objeto que nos permite decirle al geomorph que vamos a trabajar en simetría.

izde <- cbind(c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 23, 27), 
              c(21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 24, 28))
# los landmarks 22, 25, 26, 1 no se incluyen debido a que pertenecen
#a la línea media del individuo, o línea de simetría.

axy.gpa$land.pairs <- izde # sirve para Object symmetry
# se agrega a la lista del análisis de procrusto la matriz
# que contiene números correspondientes a pares de landmarks emparejados
# a ambos lados de la línea de simetría.

dimnames(axy.gpa$coords)[[3]] # esta es una manera diferente de obtener los nombres
                              # de los shapes, a través del doble corchete

# Ahora podemos obtener las diferentes variables de agrupación
# DRO-AgB_050X [sitio de muestreo] [manejo] [individuo] 
axy.gpa$ind <- substr(dimnames(axy.gpa$coords)[[3]], 9, 11) # individuo
axy.gpa$man <- substr(dimnames(axy.gpa$coords)[[3]], 5, 6) # manejo
axy.gpa$rep <- rep(c(1,2), 75) # este vector indica a qué replica pertenece cada objeto
                               # y tiene que ser igual al número total de individuos.

#### Análisis de asimetría direccional y fluctuante ####
axy.sym <- bilat.symmetry(axy.gpa$coords, ind= axy.gpa$ind, replicate = axy.gpa$rep, 
                      object.sym = TRUE, land.pairs = axy.gpa$land.pairs, 
                      print.progress = FALSE, iter = 9999)
# cuando trabajamos con "object symmetry", object.sym = TRUE
# cuando trabajamos con "matching symmetry", object.sym = FALSE y se agrega un nuevo argumento
# "side", que indica qué objeto pertenece a qué grupo de simetría.

## Procrustes ANOVA
summary(axy.sym) 
## ind = variación entre individuos, 
# indican diferencias morfológicas entre los inividuos de la muestra.
## side = mide asimetría direccional,
# permite ver si hay una tendencia de que un lado sea diferente del otro
## ind:side = mide ASIMETRÍA FLUCTUANTE, 
# permite ver si hay variación individual en el grado de asimetría.
## ind:slide:replicate = permite ver si la diferencia que vemos es prpia de los individuos
# o es un error de medición.

plot(axy.sym) # ver presentación.

###############################################################
# Preparación de datos para análisis estadísticos posteriores #
###############################################################
# Como por cada individuo hay dos réplicas,
# debemos elegir una de ellas para poder realizar análisis posteriores.
ret <-seq(from= 1, to= 150, by = 2) # en este caso nos quedamos con los individuos impares.

# Para poder hacer los análisis, creamos un objeto nuevo con los reultados que obtuvimos
# de los análisis de Procrusto y simetría anteriores.
gdf.asym <- geomorph.data.frame(sym = axy.sym$symm.shape, # Componente simétrico de la forma
                                asym = axy.sym$asymm.shape, # Componente asimétrico de la forma
                                FA = axy.sym$FA.component, # Componente asimetría fluctuante
                                manejo = as.factor(axy.gpa$man[ret]),
                                size = axy.gpa$Csize[ret],
                                ind = axy.gpa$ind[ret])

## Realizamos un Procrustes ANOVA para responder preguntas biológicas, 
# en este caso vamos a trabajar únicamente con el componente simétrico, 
# porque nos interesa saber cuánto influye el tamaño de las vaquitas y el manejo
# en la simetría de los individuos 
fit1 <- procD.lm(sym ~ size + manejo, data = gdf.asym, print.progress = FALSE,
                 iter = 9999, # número de iteraciones para el test estadístico
                 RRPP = TRUE,
                 turbo = FALSE)
summary(fit1)

 
#armamos un vector que nos permite separar ambos manejos para visualizar mejor en el grafico
colores1 <- as.character(factor(gdf.asym$manejo, labels = c("orange2", "cyan3")))
allom.plot <- plotAllometry(fit1, size = gdf.asym$size,  logsz = TRUE, method = "RegScore", 
                            pch = 19, col = colores1, xlab = "Log(CS)")

#### Predicciones de máximos, mínimos y consensos ####
preds <- shape.predictor(fit1$GM$fitted, x= allom.plot$RegScore, Intercept = FALSE, 
                         predmin = min(allom.plot$RegScore), 
                         predmax = max(allom.plot$RegScore)) 
M <- mshape(axy.sym$symm.shape)
plot(M, links = link)
plotRefToTarget(M, preds$predmin, mag=5, links = link) 
plotRefToTarget(M, preds$predmax, mag=5, links = link)


##################################################################################################################
########################################
#### A  practicar!!! ####
########################################
# Para esta parte práctica vamos a usar datos que ya están pre-cargados en el paquete geomorph
# Tienen que cargar el set de datos "mosquito" y luego,
# ejecutar el análisis de Procrustes y simetría correspondiente.

data(mosquito)

######################
#### Missing data ####
######################
# En especímenes simétricos, las ubicaciones a un lado del eje de simetría pueden utilizarse para estimar
# las ubicaciones de los puntos de referencia ausentes en el otro lado del eje.
myShapes_m <- readShapes("Shapes_axyridis_miss") # algunos shapes no tiene el landmark 5 
mydata_m <- readland.shapes(myShapes_m)
link <- read.table("link_axy.txt")
axy_m.gpa <- gpagen(mydata_m, print.progress = F)
# Error en gpagen(mydata_m, print.progress = F): 
# Data matrix contains missing values. Estimate these first (see 'estimate.missing').

# Identificar especímenes con al menos un NA
individuos_NA <- which(
  sapply(mydata_m$landmarks,
         function(x) any(is.na(x)))
)
individuos_NA

# Usando método de reconstrucción por espejo
# El landmark homólogo al 5 es el número 18, entonces tenemos que relfejar ese landmark para que 
# lo podamos utilizar como reemplazo del landmark faltante.
# Sabemos la identidad de los individuos que tienen landmarks faltantes, entonces vamos a repararlos de a uno:

landmarks_eje <- c(22, 25, 26, 1) #Definimos los landmarks que están sobre el eje de simetría para todos los ejemplares.

#### REPARACIÓN DEL ESPÉCIMEN DRO-AgC_553X [5] ####
# Aislamos al espécimen 5 desde la lista
axy553 <- mydata_m$landmarks[[5]]  
# Calculamos el centro local de simetría de ESTE individuo
# calculamos para cada uno, porque fueron digitalizados de manera individual.
centro_553 <- mean(
  axy553[landmarks_eje, 1], 
  na.rm = TRUE
)
# Creamos su versión espejo usando su propio centro
espejo553 <- axy553
espejo553[, 1] <- 2 * centro_553 - axy553[, 1]

axy553_reparado <- axy553 
axy553_reparado[5, ] <- espejo553[18, ] # Reemplazamos el landmark 5 dañado, espejando el landmark 18

mydata_m$landmarks[[5]] <- axy553_reparado # Devolvemos el ejemplar ya arreglado a la lista general 

#### REPARACIÓN DEL ESPÉCIMEN RSI-CvB_045X [23] ####
axy045 <- mydata_m$landmarks[[23]]
centro_045 <- mean(
  axy045[landmarks_eje, 1],
  na.rm = TRUE
)
espejo045 <- axy045
espejo045[, 1] <- 2 * centro_045 -
  axy045[, 1]
axy045_reparado <- axy045
axy045_reparado[5, ] <- espejo045[18, ]
mydata_m$landmarks[[23]] <- axy045_reparado


# Ahora vamso a comparar como se ven los landmarks originales y los reparados,
#para ver si hay diferencias significativas y ver si nuestra metodologíá es válida.
original <- mydata$landmarks[[5]] # este es el especímen original sin el landmark faltante
reparado <- mydata_m$landmarks[[5]]

par(mfrow = c(1, 2))

plotRefToTarget(original, original, links = link)
mtext("Original")

plotRefToTarget(original, reparado, links = link)
mtext("Original vs. reconstruido")

par(mfrow = c(1, 1))

# Individuo 5
lm5_original <- original[5, ]# Coordenada real del landmark 5
lm5_reparado <- reparado[5, ]# Coordenada reconstruida
error_5 <- sqrt(sum( # Error euclidiano
  (lm5_original - lm5_reparado)^2
))

# Individuo 23
error_23 <- sqrt(sum(
  (mydata$landmarks[[23]][5, ] -
     mydata_m$landmarks[[23]][5, ])^2
))

data.frame(
  Individuo = c(5, 23),
  Error = c(error_5, error_23)
)
