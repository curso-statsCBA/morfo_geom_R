#Práctico 8 - Introducción a Métodos Comparativos Filogenéticos (y no filogenéticos)

rm(list=ls(all=TRUE))
library(StereoMorph)
library(geomorph)
library(ape)
library(phytools)
library(magrittr)
library(morphospace)
library(mvMORPH)
library(Morpho)
library(geiger)
library(readxl)
library(ggplot2)
library(viridis)
library(Hmisc)
library(svglite)

#Manipulación de filogenias
#Ultamétrico y no ultramétrico

#load("Beak_full_shape_coordinates_sps_means.RData")# contains shapes_a object
#load("Beak_non_allometric_shape_coordinates_sps_means.RData")# contains shapes_na object
#load("Beak_Centroid_Size_sps_means.RData") 

#Obtenemos las coordenadas crudas
Shapesfurn<-readShapes(file='Shapes')

Shapes_Beak_Furnariids_with_sliding_semilandmarks<-readland.shapes(Shapesfurn,nCurvePts=c(10,10,13))
#Shapes_Neurocranium_Furnariids_without_sliding_semilandmarks<-readland.shapes(Shapesfurn,nCurvePts=NULL)

myGPA <- gpagen(Shapes_Beak_Furnariids_with_sliding_semilandmarks) # GPA perfomed with minimized bending energy
plot(myGPA)

#Generar el arbol recortado con las 85 especies
Species<-factor(c("Geositta_cunicularia","Geositta_punensis","Sclerurus_caudacutus","Sclerurus_rufigularis","Sclerurus_scansor","Dendrocincla_anabatina","Dendrocincla_homochroa","Dendrocincla_merula","Dendrocolaptes_certhia","Dendrocolaptes_platyrostris","Drymornis_bridgesii","Glyphorynchus_spirurus","Glyphorynchus_spirurus","Lepidocolaptes_affinis","Lepidocolaptes_angustirostris","Lepidocolaptes_angustirostris","Lepidocolaptes_falcinellus","Lepidocolaptes_fuscicapillus","Sittasomus_griseicapillus","Xiphocolaptes_major","Xiphorhynchus_flavigaster","Xiphorhynchus_guttatus","Xiphorhynchus_obsoletus","Xiphorhynchus_spixii","Berlepschia_rikeri","Berlepschia_rikeri","Xenops_minutus","Xenops_minutus","Xenops_minutus","Xenops_rutilans","Pygarrhichas_albogularis","Pygarrhichas_albogularis","Pygarrhichas_albogularis","Ochetorhynchus_phoenicurus","Premnornis_guttuliger","Pseudocolaptes_boissonneautii","Pseudocolaptes_lawrencii","Tarphonomus_certhioides","Cinclodes_fuscus","Cinclodes_fuscus","Cinclodes_nigrofumosus","Cinclodes_oustaleti","Cinclodes_patagonicus","Cinclodes_patagonicus","Cinclodes_patagonicus","Furnarius_figulus","Furnarius_leucopus","Furnarius_leucopus","Furnarius_rufus","Furnarius_rufus","Furnarius_rufus","Furnarius_rufus","Furnarius_torridus","Limnornis_curvirostris","Phleocryptes_melanops","Phleocryptes_melanops","Phleocryptes_melanops","Upucerthia_dumetaria","Upucerthia_dumetaria","Anumbius_annumbi","Anumbius_annumbi","Anumbius_annumbi","Aphrastura_spinicauda","Aphrastura_spinicauda","Asthenes_baeri","Asthenes_dorbignyi","Asthenes_modesta","Asthenes_pyrrholeuca","Certhiaxis_cinnamomeus","Certhiaxis_cinnamomeus","Coryphistera_alaudina","Cranioleuca_demissa","Cranioleuca_gutturata","Leptasthenura_aegithaloides","Leptasthenura_aegithaloides","Leptasthenura_platensis","Leptasthenura_striata","Limnoctites_rectirostris","Phacellodomus_ruber","Phacellodomus_ruber","Phacellodomus_striaticollis","Phacellodomus_striaticollis","Pseudoseisura_gutturalis","Pseudoseisura_lophotes","Pseudoseisura_lophotes","Pseudoseisura_lophotes","Schoeniophylax_phryganophilus","Synallaxis_cinerascens","Synallaxis_frontalis","Synallaxis_ruficapilla","Synallaxis_rutilans","Synallaxis_scutata","Synallaxis_spixi","Margarornis_squamiger","Premnoplex_brunnescens","Anabacerthia_lichtensteini","Anabacerthia_ruficaudata","Anabazenops_dorsalis","Ancistrops_strigilatus","Automolus_infuscatus","Automolus_melanopezus","Automolus_ochrolaemus","Automolus_ochrolaemus","Automolus_paraensis","Clibanornis_rubiginosus","Dendroma_rufa","Heliobletus_contaminatus","Philydor_erythrocercum","Philydor_pyrrhodes","Syndactyla_rufosuperciliata","Syndactyla_rufosuperciliata","Thripadectes_holostictus","Thripadectes_rufobrunneus"))

#Cargar la topología de Harvey et al 2020 y excluir las especies presentes en los datos de morfometria geometrica
tree0<-read.tree("Harvey_et_al_2020.tre")
tree_l <- drop.tip(tree0, tip = tree0$tip.label[!tree0$tip.label%in%levels(Species)])
plot(tree_l,cex=0.4)

#Exportamos la filogenia recortada para usarla en las prácticas siguientes
save(tree_l,file="tree_85.RData")

#Obtenemos las coordenadas de Procrustes consenso de cada especie y hacemos que el orden de los nombres siga el de la filogenia
consensus_l<-expected_shapes(myGPA$coords, Species)#expected_shapes es una función del paquete morphospace
consensus_l<-consensus_l[,,dimnames(consensus_l)[[3]]%in%tree_l$tip.label]
shapes_a <- consensus_l[,,tree_l$tip.label]

#Importante: comprobar si los valores consenso de las especies en la matriz de coordenadas de Procrustes tienen el mismo orden que en la filogenia (debe dar TRUE) para evitar resultados espurios en distintos análisis filogenéticos.
identical(dimnames(shapes_a)[[3]], tree_l$tip.label)

all(tree_l$tip.label %in% dimnames(consensus_l)[[3]])

#Creamos el objeto R con las coordenadas para usarlo en futuras sesiones
save(shapes_a,file="shapes_a_85_PICO.RData")
#Podemos exportar en .csv también
write.csv(shapes_a,"shapes_a_PICO.csv")

#Extraemos el Tamaño del Centroide (CS) obtenido previamente en el GPA
sizes_l<-tapply(X=myGPA$Csize, INDEX = Species, FUN=mean)
sizes_l<-sizes_l[names(sizes_l)%in%tree_l$tip.label]
#Reordenar los nombres de las especies en función de la filogenia
sizes_l <- c(sizes_l[tree_l$tip.label])

#Chequear que esta ok el orden de los nombres
identical(names(sizes_l), tree_l$tip.label)
#y tambien chequeamos que no se hayan introducido NAs por coerción
all(tree_l$tip.label %in% names(sizes_l))

#Obtenemos los residuos sin alometría de las coordenadas de Procrustes (Esto puede hacerse con una regresión lineal no-filogenética si no tienen filogenia)
shapes_na<-lm(two.d.array(shapes_a) ~ log(sizes_l)) %>% detrend_shapes(tree = tree_l) %>% arrayspecs(p = nrow(shapes_a), k = 2)

#Señal filogenética
#K-mult con geomorph
kmult<-physignal(shapes_a, tree_l, iter = 999, seed = NULL, print.progress = FALSE)
plot(kmult)
summary(kmult)

kmult1Zscoresna<-physignal.z(shapes_a, tree_l, iter = 999, seed = NULL, print.progress = FALSE)
kmult1Zscoresna

#Test de alometria filogenético
SkullAllometry <- geomorph.data.frame(shapes_a, sizes_l, phy = tree_l)
SkullAllometry.Dpgls <- procD.pgls(shapes_a ~ log(sizes_l), phy = tree_l, data = SkullAllometry, 
                                   iter = 999)
summary(SkullAllometry.Dpgls)
plot(SkullAllometry.Dpgls)

#Test de alometría no filogenético
SkullProclm <- geomorph.data.frame(shapes_a,sizes_l)
skullproclm<-procD.lm(shapes_a ~ log(sizes_l), data = SkullProclm,iter = 999)
summary(skullproclm)
plot(skullproclm)

#Análisis de Componentes Principales clásico
#Con alometria
fullshapecoords<-two.d.array(shapes_a)
pca<-prcomp(two.d.array(shapes_a))
plot(pca$x,col=Species)# colores sin nombres
summary(pca)
todos_los_pcs<-data.frame(pca$x)
todos_los_pcs <- todos_los_pcs[rownames(todos_los_pcs)%in%tree_l$tip.label,][tree_l$tip.label,]
PC1<-setNames(todos_los_pcs$PC1, rownames(todos_los_pcs))
PC1a5<-todos_los_pcs[,1:5]
rownames(pca$x) <- dimnames(shapes_a)[[3]]

#Sin colores y con nombres
plot(pca$x[, 1], pca$x[, 2],
     col = "black",
     xlab = "PC1",
     ylab = "PC2",
     pch = 19)

text(pca$x[, 1], pca$x[, 2],
     labels = rownames(pca$x),
     pos = 3,
     cex = 0.6)

#Con colores y nombres
plot(pca$x[, 1], pca$x[, 2],
     col = Species,
     pch = 19,
     xlab = "PC1",
     ylab = "PC2")

text(pca$x[, 1], pca$x[, 2],
     labels = rownames(pca$x),
     pos = 3,
     cex = 0.6)

#PCA Sin alometria
Non_Allometric_Shapecoords<-two.d.array(shapes_na)
pca_na<-prcomp(two.d.array(shapes_na))
plot(pca_na$x,col=Species)
summary(pca_na)
todos_los_pcs_na<-data.frame(pca_na$x)
todos_los_pcs_na <- todos_los_pcs_na[rownames(todos_los_pcs_na)%in%tree_l$tip.label,][tree_l$tip.label,]
PC1_na<-setNames(todos_los_pcs_na$PC1, rownames(todos_los_pcs_na))
PC1a5_na<-todos_los_pcs_na[,1:5]

#Con colores y nombres
plot(pca_na$x[, 1], pca_na$x[, 2],
     col = Species,
     pch = 19,
     xlab = "PC1",
     ylab = "PC2")

text(pca_na$x[, 1], pca_na$x[, 2],
     labels = rownames(pca_na$x),
     pos = 3,
     cex = 0.6)

#Colorear por clados o variable de intéres
#Variable ecológica
data <- read_excel("Furnariidae_107_sps_Ecological_Variables.xlsx", 
                   sheet = "Type_of_Nest")
data <- as.data.frame(data)
rownames(data) <- data[,1]
data <- data[,-1]
data <- data[rownames(data) %in% tree_l$tip.label, ][tree_l$tip.label, ]

#Creamos una función para plotear el PCA con distintas variables

plot_pca_grupos <- function(pca, grupo, titulo = "",
                            paleta = NULL) {
  
  grupo <- factor(grupo)
  
  if (is.null(paleta)) {
    colores <- rainbow(nlevels(grupo))
  } else {
    colores <- paleta[levels(grupo)]
  }
  
  col_puntos <- colores[as.character(grupo)]
  
  plot(pca$x[, 1], pca$x[, 2],
       col = col_puntos,
       pch = 19,
       xlab = "PC1",
       ylab = "PC2",
       main = titulo)
  
  text(pca$x[, 1], pca$x[, 2],
       labels = rownames(pca$x),
       pos = 3,
       cex = 0.6)
  
  legend("topright",
         legend = levels(grupo),
         col = colores[levels(grupo)],
         pch = 19,
         cex = 0.7)
}

#Usamos tipo de nido
Type_of_Nest<-as.factor(data$Type_of_Nest)
names(Type_of_Nest)<-rownames(data)

colores_nido <- c(
  "0" = "#0072B2",  # Cavidad
  "1" = "#009E73"   # Domo
)


plot_pca_grupos(
  pca,
  Type_of_Nest,
  "PCA por tipo de nido",
  paleta = colores_nido
)

#Filomorfoespacios (filogenético). Se puede hacer con phytools o con geomorph
library(phytools)
#Creo un sinmap con la variable nido usando el vector factor "Type_of_Nest" que genere antes
treesimmNest<-make.simmap(tree_l,Type_of_Nest,model="ER",nsim = 1)#esto genera una sola reconstruccion de estados ancestrales estocastica. Lo ideal de cara a una publicación es generar 99 y correr OUM y BMM sobre las 99 
col1 <- setNames(c("#009E73", "#CC79A7"), c("0", "1"))
plot(treesimmNest,fsize=0.5, colors = col1)

#Con nombres de especies
phylomorphospace(treesimmNest,todos_los_pcs,label=c("horizontal"),bty="l",fsize=0.6,node.by.map=TRUE,node.size=c(0,1.2),xlab="PC 1",
                 ylab="PC 2",colors = col1)

#Sin nombres de especies y con mapeo de indice kinetico (semi-cuantitativa, categorias del 0 al 7) y formas de los tips en funcion del tipo de nido
library(viridis)

# KI as a factor with seven nivels 0:7
KI <- as.factor(data$Kinetic_Index)
levels(KI) <- as.character(0:7)
names(KI) <- rownames(data)

# Simmap
set.seed(123)
treesimmki <- make.simmap(tree_l, KI, model = "ER", nsim = 1)

# Paleta viridis con 8 colores
cols <- setNames(viridis(8), as.character(0:7))

# Plot simmap with discret viridis palette
plot(treesimmki, fsize = 0.4, colors = cols)
add.simmap.legend(colors = cols, prompt = FALSE,
                  x = 0, y = 0.9*max(nodeHeights(treesimmki)),
                  fsize = 0.7, vertical = FALSE)

phylomorphospace(treesimmki,PC1a5,label=c("horizontal"),bty="l",ftype="off",node.by.map=TRUE,node.size=c(0,0),xlab="PC 1",
                 ylab="PC 2", colors = cols)

#Nido
#Type_of_Nest<-as.factor(data$Type_of_Nest)
#names(Type_of_Nest)<-rownames(data)

# Elegir forma para los tips (cuadrado nro 21 y circulo nro 22)
pch_vals <- setNames(c(21, 22), c("0", "1"))

for (state in levels(Type_of_Nest)) {
  tips_state <- names(Type_of_Nest)[Type_of_Nest== state]
  coords <- todos_los_pcs[tips_state, ]
  
  # color según Tipo de Nido (usamos los colores ya definidos para el primer filomorfoespacio)
  cols_tip <- cols[as.character(KI[tips_state])]
  
  points(coords,
         pch = pch_vals[state],
         bg = cols_tip,
         col = "black",
         cex = 2.3,
         lwd = 1)
}

#Podemos hacerlo para habitat (categorica con tres categorias)
data2 <- read_excel("Furnariidae_107_sps_Ecological_Variables.xlsx", 
                   sheet = "Primary_Habitat")
data2 <- as.data.frame(data2)
rownames(data2) <- data2[,1]
data2 <- data2[,-1]
data2 <- data2[rownames(data2) %in% tree_l$tip.label, ][tree_l$tip.label, ]
Primary_Habitat<-as.factor(data2$Primary_Habitat)
names(Primary_Habitat)<-rownames(data2)
treesimmHabitat<-make.simmap(tree_l,Primary_Habitat,model="ER",nsim = 1) 
cols3 <- setNames(c("#E69F00","green","deeppink"), c("0", "1", "2"))
plot(treesimmHabitat,fsize=0.4, colors = cols3)

phylomorphospace(treesimmki,PC1a5,label=c("horizontal"),bty="l",ftype="off",node.by.map=TRUE,node.size=c(0,0),xlab="PC 1",
                 ylab="PC 2", colors = cols)

#Tres formas distintas para las categorias de Habitat Primario
pch_vals <- setNames(c(21, 22,25), c("0", "1","2"))

#Habitat primario
for (state in levels(Primary_Habitat)) {
  tips_state <- names(Primary_Habitat)[Primary_Habitat== state]
  coords <- todos_los_pcs[tips_state, ]
  
# color según Habitat Primario (usamos los colores ya definidos)
  cols_tip <- cols[as.character(KI[tips_state])]
  
  points(coords,
         pch = pch_vals[state],
         bg = cols_tip,
         col = "black",
         cex = 2.3,
         lwd = 1)
}


#Fenograma (filogenético)
library(phytools)
phenogram(treesimmNest,PC1,spread.labels=TRUE, spread.cost=c(1,0.7),link=0.2,fsize=0.5,colors = col1)
phenogram(treesimmki,PC1,spread.labels=TRUE, spread.cost=c(1,0.7),link=0.2,fsize=0.5,colors = cols)
phenogram(treesimmHabitat,PC1,spread.labels=TRUE, spread.cost=c(1,0.7),link=0.2,fsize=0.5,colors = cols3)

#Regresiones lineales Filogenéticas y No filogenéticas: MANOVAs y MANCOVAs
library(geomorph)
library(mvMORPH)
library(phytools)
library(viridis)
library(readxl)

#Tipo de nido MANCOVA filogenético
SkullNAC_Nest <- geomorph.data.frame(shapes_a, data$Type_of_Nest,data$Kinetic_Index, phy = tree_l)
SkullNAC.Dpgls_Nest <- procD.pgls(shapes_a ~ data$Type_of_Nest + data$Type_of_Nest:data$Kinetic_Index
                                  , phy = tree_l, data = SkullNAC_Nest, iter = 999)
summary(SkullNAC.Dpgls_Nest)

#Tipo de nido MANCOVA no filogenético
SkullProclm_Nest <- geomorph.data.frame(shapes_a, data$Type_of_Nest,data$Kinetic_Index)
skullproclm_Nest<-procD.lm(shapes_a ~ data$Type_of_Nest + data$Type_of_Nest:data$Kinetic_Index, data = SkullProclm,iter = 999)
summary(skullproclm_Nest)

#Graficar la relacion entre tres variables
install.packages("Hmisc")
library(Hmisc)
library(ggplot2)
ggplot(data = data.frame(Type_of_Nest = data$Type_of_Nest,
                         PC1 = PC1,
                         Kinetic_Index = data$Kinetic_Index),
       aes(x = factor(Type_of_Nest), y = PC1, color = Kinetic_Index)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 18) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) +
  labs(x = "Type of Nest", y = "Shape (PC 1)", color = "Kinetic Index") +
  theme_minimal()

#Graficar variación asociada a la variable correlacionada
#Supplementary Figures S5, S6, S8 and S9
mean.shape <- mshape(shapes_a)
coef_vector <- SkullNAC.Dpgls_Nest$pgls.coefficients[2, ]
coef_matrix_raw <- as.matrix(t(coef_vector))
coef_matrix <- arrayspecs(coef_matrix_raw, p = 33, k = 2)
coef_matrix <- coef_matrix[, , 1]
dim(coef_matrix)
[1] 33  2
scaling_factor <- 1
target.shape <- mean.shape + coef_matrix * scaling_factor
plotRefToTarget(M1 = mean.shape,
                M2 = target.shape,
                method = "vector")

range_var <- range(data$Type_of_Nest)
delta_shape <- coef_matrix * diff(range_var)
dist_procrustes <- sqrt(sum(delta_shape^2))
dist_procrustes


#############E-PGLS###########################################################

#Grafico Santi

####Variables ecológicas multivariadas########################################
#PLS filogenetico y PLS no filogenetico con variables ecológicas multivariadas
library(readxl)
data3 <- read_excel("Furnariidae_107_sps_Ecological_Variables.xlsx", 
                    sheet = "Dietary_Items")
data3 <- as.data.frame(data3)
rownames(data3) <- data3[,1]
data3 <- data3[,-1]
data3 <- data3[rownames(data3)%in%tree_l$tip.label,][tree_l$tip.label,]


# Eliminamos la columnad del indice kinetico porque el PLS es entre la forma y la variable ecológica solamente
eco_vars <- data3[, colnames(data3) != "Kinetic_Index"]

pls_phylo <- phylo.integration(A = shapes_a, 
                               A2 = eco_vars, 
                               phy = tree_l, 
                               iter = 999)

summary(pls_phylo)

scores_morpho <- pls_phylo$XScores[, 1]

scores_diet  <- pls_phylo$YScores[, 1]


library(ggplot2)
install.packages("svglite")
library(svglite)

# Crear el objeto para plotear
KI <- data3$Kinetic_Index
df_plot <- data.frame(scores_morpho, scores_diet, KI)

# Scatterplot of PLS scores (morphology vs diet), colored by kinetic index
p <- ggplot(df_plot, aes(x = scores_morpho, y = scores_diet, color = KI)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, aes(group = 1), color = "black") + # Add global linear regression line
  scale_color_viridis_c() + # Continuous color scale for kinesis
  theme_minimal(base_size = 14) +
  labs(x = "PLS x-score", y = "PLS y-score", color = "Kinesis")

p

# dev.off()

# Grabar en formato SVG
ggsave("PLS_forma_items_alimentarios.svg", 
       plot = p,
       device = "svg", width = 6, height = 5, units = "in")

######################
#Análisis de Componentes Principales filogenéticos. Se puede hacer con phytools o con geomorph.
#PCA filogenético con allometria:
library(phytools)
phylPCA<-phyl.pca(tree_l,fullshapecoords)
plot(phylPCA$S,main="",col=Species)

Phy_PCA<-data.frame(phylPCA$S)
Phy_PCA<- Phy_PCA[rownames(Phy_PCA)%in%tree_l$tip.label,][tree_l$tip.label,]
Phy_PC1<-setNames(Phy_PCA$PC1, rownames(Phy_PCA))
Phy_PC1aPC5<-Phy_PCA[,1:5]

write.csv(Phy_PCA,"Phylogenetic_PCA.csv")

#PCA filogenético sin allometria:
Phy_PCA_N_A<-phyl.pca(tree_l,Non_Allometric_Shapecoords)
par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(Phy_PCA_N_A$S,main="",col=Species)
Phy_PCA_N_A_PCs<-data.frame(Phy_PCA_N_A$S)
Phy_PCA_N_A_PCs<- Phy_PCA_N_A_PCs[rownames(Phy_PCA_N_A_PCs)%in%tree_l$tip.label,][tree_l$tip.label,]
Phy_PCA_N_A_PC1<-setNames(Phy_PCA_N_A_PCs$PC1, rownames(Phy_PCA_N_A_PCs))
Phy_PCA_N_A_PC1aPC5<-Phy_PCA_N_A_PCs[,1:5]

write.csv(Phy_PCA_N_A_PCs,"PCA_Phy_N_A_PCs.csv")

#####################################################################################################
#Test de Hipótesis de modularidad evolutiva##########################################################
#Para este análisis cargamos set de datos del cráneo completo de furnáridos para 83 especies
#de Stefanini et al. 2025
load("Skull_full_shape_coordinates_83_sps_means.RData")
load("tree_83_sps.RData")

#Hypothesis 1 Two modules: Beak (rostrum + nares) and neurocranium
land.gpsH1<-c("A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A")
H1a<- phylo.modularity(shapes_a, partition.gp=land.gpsH1, phy=tree_l,CI = TRUE, iter=999)

#Hypothesis 2 Two modules: Nares, (rostrum + neurocranium)
land.gpsH2<-c("B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","A","A","A","A","A","A","A","A","A","A","A")
H2a<- phylo.modularity(shapes_a, partition.gp=land.gpsH2, phy=tree_l,CI = TRUE, iter=999)

#Hypothesis 3 Two modules: Beak, (nares + neurocranium)
land.gpsH3<-c("A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B")
H3 <- phylo.modularity(shapes_a, partition.gp=land.gpsH3, phy=tree_l,CI = TRUE, iter=999)

#Hypothesis 4 Three modules: Nares, beak and neurocranium
land.gpsH4<-c("A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","C","C","C","C","C","C","C","C","C","C","C")
H4 <- phylo.modularity(shapes_a, partition.gp=land.gpsH4, phy=tree_l,CI = TRUE, iter=999)

#Compare
testcr <- compare.CR(H1, H2, H3, H4, CR.null = TRUE)
summary(testcr)
