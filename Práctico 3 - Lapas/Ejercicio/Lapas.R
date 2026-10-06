library(geomorph)
library(StereoMorph)
library(ggplot2)
library(ggrepel)
library(MASS)
library(randomForest)
library(mclust)
library(vegan)

# ============================================================
# 0. PROCESADO DE FOTOS
# ============================================================

# digitizeImages(
#   image.file = 'images',
#   shapes.file = 'shapes',
#   landmarks.ref = paste("LM", c(1:2), sep = ""),
#   curves.ref = "curves.txt"
# )

# ============================================================
# 1. CARGAR SHAPES
# ============================================================

shapes <- readShapes("Práctico 3 - Lapas/Ejercicio/shapes")

shapesGM <- readland.shapes(
  shapes,
  nCurvePts = c(15, 15)
)

link <- read.table(
  "Práctico 3 - Lapas/Ejercicio/link.txt",
  header = FALSE,
  sep = ""
)


# ============================================================
# 2. GENERALIZED PROCRUSTES ANALYSIS (GPA)
# ============================================================

Y.gpa <- gpagen(shapesGM)


# ============================================================
# 3. IDENTIFICAR GRUPOS E INDIVIDUOS
# ============================================================

# Nombres:
# L01, L02... = Litoral
# S01, S02... = Sublitoral
# B01, B02... = Basurero

grupo <- substr(
  names(shapesGM$landmarks), 1, 1)

indiv <- substr(
  names(shapesGM$landmarks), 2, 3)

# Convertir códigos a nombres descriptivos
grupo <- factor(
  grupo,
  levels = c("L", "S", "B"),
  labels = c("Litoral", "Sublitoral", "Basurero")
)

# Comprobar cantidad de individuos por grupo
table(grupo)


# Guardar información en el objeto GPA
Y.gpa$grupo <- grupo
Y.gpa$indiv <- indiv


# ============================================================
# 4. FORMA MEDIA
# ============================================================

ref <- mshape(Y.gpa$coords)


# ============================================================
# 5. PRINCIPAL COMPONENTS ANALYSIS (PCA)
# ============================================================

PCA <- gm.prcomp(Y.gpa$coords)

summary(PCA)


# ============================================================
# 6. CONSENSUS / REFERENCE SHAPE
# ============================================================

plotRefToTarget(
  ref,
  ref,
  links = link,
  label = TRUE
)


# ============================================================
# 7. DATA FRAME PARA GRAFICOS
# ============================================================

DAT <- data.frame(
  PC1 = PCA$x[, 1],
  PC2 = PCA$x[, 2],
  PC3 = PCA$x[, 3],
  log.CS = log(Y.gpa$Csize),
  grupo = grupo,
  indiv = names(shapesGM$landmarks)
)


# ============================================================
# 8. PCA: LITORAL VS SUBLITORAL + BASURERO
# ============================================================

ggplot(
  DAT,
  aes(
    x = PC1,
    y = PC2,
    color = grupo
  )
) +  
  # Elipses solamente para los grupos conocidos
  stat_ellipse(
    data = subset(DAT, grupo != "Basurero"),
    aes(group = grupo),
    level = 0.95,
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  geom_text_repel(
    aes(label = indiv),
    size = 3,
    show.legend = FALSE
  ) +
  theme_classic()


# ============================================================
# 9. PROCRUSTES ANOVA
#    LITORAL VS SUBLITORAL
# ============================================================

# Seleccionar solamente los grupos conocidos
conocidos <- grupo %in% c(
  "Litoral",
  "Sublitoral"
)

# Coordenadas Procrustes de los individuos conocidos
Y.conocidos <- Y.gpa$coords[, , conocidos]

# Factor de grupo
grupo.conocidos <- droplevels(
  grupo[conocidos]
)

# Modelo Procrustes
modelo <- procD.lm(
  Y.conocidos ~ grupo.conocidos,
  iter = 999
)

# Resultados
summary(modelo)


# ============================================================
# 10. SILUETAS PARA FORMAS MAX Y MIN DE CADA EJE
# ============================================================

# ----------------------------
# PC1
# ----------------------------

PC1 <- PCA$x[, 1]

M <- mshape(Y.gpa$coords)

preds <- shape.predictor(
  Y.gpa$coords,
  x = PC1,
  Intercept = FALSE,
  pred1 = min(PC1),
  pred2 = max(PC1),
  pred3 = mean(PC1)
)

plotRefToTarget(
  M,
  preds$pred1,
  links = link,
  mag = 1,
  method = "points"
)

mtext("PC1 - Min.")

plotRefToTarget(
  M,
  preds$pred2,
  links = link,
  mag = 1,
  method = "points"
)

mtext("PC1 - Max.")

plotRefToTarget(
  M,
  preds$pred3,
  links = link,
  mag = 1,
  method = "points"
)

mtext("PC1 - Mean")


# ----------------------------
# PC2
# ----------------------------

PC2 <- PCA$x[, 2]

M <- mshape(Y.gpa$coords)

preds <- shape.predictor(
  Y.gpa$coords,
  x = PC2,
  Intercept = FALSE,
  pred1 = min(PC2),
  pred2 = max(PC2),
  pred3 = mean(PC2)
)

plotRefToTarget(
  M,
  preds$pred1,
  links = link,
  mag = 1,
  method = "points"
)

mtext("PC2 - Min.")

plotRefToTarget(
  M,
  preds$pred2,
  links = link,
  mag = 1,
  method = "points"
)

mtext("PC2 - Max.")

plotRefToTarget(
  M,
  preds$pred3,
  links = link,
  mag = 1,
  method = "points"
)

mtext("PC2 - Mean")


# ============================================================
# 11. Analisis Discriminante
# ============================================================

### A. Separar individuos con grupo conocido

DAT.conocidos <- subset(DAT, grupo != "Basurero")

# Eliminar niveles no utilizados
DAT.conocidos$grupo <- droplevels(DAT.conocidos$grupo)

# Cantidad de individuos por grupo
table(DAT.conocidos$grupo)


### B. Ajustar LDA con PC1 y PC2 

LDA <- lda(
  grupo ~ PC1 + PC2,
  data = DAT.conocidos
)

LDA


### C. Validación cruzada leave-one-out

LDA.cv <- lda(
  grupo ~ PC1 + PC2,
  data = DAT.conocidos,
  CV = TRUE
)

# Matriz de confusión
table(
  Real = DAT.conocidos$grupo,
  Predicho = LDA.cv$class
)

# Accuracy
accuracy.LDA <- mean(
  LDA.cv$class == DAT.conocidos$grupo
)

accuracy.LDA


### D. Identificar individuos mal clasificados (Si los hubiera)

DAT.conocidos$LDA_pred <- LDA.cv$class

DAT.conocidos[
  DAT.conocidos$grupo != DAT.conocidos$LDA_pred,
  c("indiv", "grupo", "LDA_pred")
]


### E. Separar ejemplares del Basurero 

DAT.basurero <- subset(DAT, grupo == "Basurero")


### F. Clasificar Basurero con el LDA 

pred.B <- predict(
  LDA,
  newdata = DAT.basurero[, c("PC1", "PC2")]
)


### G. Resultados de clasificación

resultado.LDA <- data.frame(
  indiv = DAT.basurero$indiv,
  prediccion = pred.B$class,
  prob.Litoral = pred.B$posterior[, "Litoral"],
  prob.Sublitoral = pred.B$posterior[, "Sublitoral"]
)

resultado.LDA

### 8. Scores de la función discriminante

LDA.scores <- predict(LDA)$x

DAT.LDA <- data.frame(
  indiv = DAT.conocidos$indiv,
  grupo = DAT.conocidos$grupo,
  LD1 = LDA.scores[, 1]
)

DAT.LDA

ggplot(DAT.LDA, aes(x = LD1, y = 0, color = grupo)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = indiv),
    size = 3,
    show.legend = FALSE
  ) +
  theme_classic() +
  labs(
    x = "Función discriminante 1 (LD1)",
    y = NULL
  )


### =========================================================
### 12. RANDOM FOREST
### =========================================================

### A. Ajustar Random Forest ---------------------------------

set.seed(123)

RF <- randomForest(
  grupo ~ PC1 + PC2,
  data = DAT.conocidos,
  ntree = 1000,
  importance = TRUE
)


### B. Evaluación del modelo ---------------------------------

# Matriz de confusión
RF$confusion

# Accuracy OOB
accuracy.RF <- 1 - RF$err.rate[nrow(RF$err.rate), "OOB"]

accuracy.RF


### C. Importancia de las variables --------------------------

importance(RF)

varImpPlot(RF)


### D. Clasificar ejemplares del Basurero --------------------

clasificacion.RF <- predict(
  RF,
  newdata = DAT.basurero[, c("PC1", "PC2")]
)

# Probabilidades de pertenencia
prob.RF <- predict(
  RF,
  newdata = DAT.basurero[, c("PC1", "PC2")],
  type = "prob"
)


### E. Resultados de clasificación ---------------------------

resultado.RF <- data.frame(
  indiv = DAT.basurero$indiv,
  prediccion = clasificacion.RF,
  prob.Litoral = prob.RF[, "Litoral"],
  prob.Sublitoral = prob.RF[, "Sublitoral"]
)

resultado.RF


### =========================================================
### PCA FINAL - GRUPOS CONOCIDOS + BASURERO
### =========================================================

ggplot() +
  
  # Elipses de los grupos conocidos
  stat_ellipse(
    data = subset(DAT, grupo != "Basurero"),
    aes(x = PC1, y = PC2, color = grupo),
    level = 0.95,
    linewidth = 0.8
  ) +
  
  # Individuos conocidos
  geom_point(
    data = subset(DAT, grupo != "Basurero"),
    aes(x = PC1, y = PC2, color = grupo),
    size = 3
  ) +
  
  # Ejemplares del Basurero
  geom_point(
    data = DAT.basurero,
    aes(x = PC1, y = PC2),
    shape = 17,
    size = 3.5,
    color = "black"
  ) +
  
  # Nombres de los ejemplares del Basurero
  geom_text_repel(
    data = DAT.basurero,
    aes(x = PC1, y = PC2, label = indiv),
    size = 3,
    color = "black"
  ) +
  
  # Nombres de los ejemplares conocidos
  geom_text_repel(
    data = subset(DAT, grupo != "Basurero"),
    aes(x = PC1, y = PC2, label = indiv, color = grupo),
    size = 2.5,
    show.legend = FALSE
  ) +
  
  theme_classic() +
  
  labs(
    x = paste0("PC1 (", round(summary(PCA)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(PCA)$importance[2,2] * 100, 1), "%)"),
    color = "Grupo"
  )


### =========================================================
### 13. GAUSSIAN ADMIXTURE ANALYSIS
### =========================================================

### A. Buscar mejor modelo

variables <- c("PC1", "PC2")
BIC <- mclustBIC(DAT[, variables])
plot(BIC)
summary(BIC)
BIC

### Best model
mod1 <- Mclust(DAT[, variables], x = BIC)
summary(mod1, parameters = TRUE)

### second model
mod2 <- Mclust(DAT[, variables], G = 2, modelNames = "EEE")
summary(mod2, parameters = TRUE)

plot(mod1, what = "classification")

plot(mod2, what = "classification")

### B. comparación entre la clasificación a priori y con GMA
table(DAT$grupo, mod1$classification)

colores <- as.character(DAT$grupo)
colores[colores == "Litoral"] <- "#e35e06"
colores[colores == "Basurero"] <- "black"
colores[colores == "Sublitoral"] <- "#0258af"

plot(DAT[, c("PC1", "PC2")], type = "n")
ordihull(
  DAT[, c("PC1", "PC2")], groups = mod1$classification,
  draw = "polygon", col = "grey", lwd = 0.1
)
points(DAT[, c("PC1", "PC2")], col = colores, pch = 19)
