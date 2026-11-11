#Práctico 10 - Simulacion de modelos macroevolutivos
# BM, OU1, OUM y Early Burst (EB)
# Trayectorias fenotipicas + disparidad entre linajes contemporaneos
# ============================================================

rm(list=ls(all=TRUE))
library(ape)
library(ggplot2)

set.seed(20261020)

# -------------------------
# 1. Filogenia y parametros
# -------------------------
#Simulamos 300 especies
n_tips <- 300
tree <- ape::rcoal(n_tips)
tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
Tmax <- max(ape::node.depth.edgelength(tree))
root_node <- n_tips + 1L

# Resolucion de la simulacion y numero de replicas para disparidad
dt <- 0.01
n_replicates <- 100

# Parametros separados para cada modelo
sigma_BM <- 0.30

theta_OU1 <- 0
alpha_OU1 <- 8
sigma_OU1 <- 0.30

theta1_OUM <- -0.6
theta2_OUM <- 0.6
alpha_OUM <- 1.5
sigma_OUM <- 0.20

sigma_EB <- 0.45
beta_EB <- -6  # desaceleracion marcada: divergencia temprana y trayectorias amesetadas

root_value <- 0

# -------------------------
# 2. Definir dos regimenes para OU de dos óptimos. 
# -------------------------
# Para que el patron sea claro asignamos theta1 y theta2 a los dos
# grandes clados que parten de la raiz. Asi ambos regimenes actuan desde
# el comienzo de la evolución del clado.
root_children <- tree$edge[tree$edge[, 1] == root_node, 2]

if (length(root_children) != 2L) {
  stop("Se espera una raiz bifurcante con dos clados descendientes")
}

clade_node <- root_children[2]
clade_tips <- ape::extract.clade(tree, clade_node)$tip.label

message(
  "OUM: regimen theta1 = ", theta1_OUM, " en un clado y theta2 = ",
  theta2_OUM, " en el otro. El clado theta2 contiene ",
  length(clade_tips), " especies."
)

# Devuelve TRUE si node es el nodo del clado o desciende de el.
is_descendant_or_self <- function(node, ancestor, edge_matrix) {
  if (node == ancestor) return(TRUE)
  current <- node
  repeat {
    idx <- which(edge_matrix[, 2] == current)
    if (length(idx) == 0) return(FALSE)
    current <- edge_matrix[idx[1], 1]
    if (current == ancestor) return(TRUE)
  }
}

# -------------------------
# 3. Simular un caracter sobre todas las ramas
# -------------------------
simulate_model <- function(model, tree, dt, root_value,
                           sigma_BM, sigma_OU1, alpha_OU1, theta_OU1,
                           sigma_OUM, alpha_OUM, theta1_OUM, theta2_OUM,
                           sigma_EB, beta_EB, clade_node) {
  ntip <- length(tree$tip.label)
  n_total <- ntip + tree$Nnode
  edge_mat <- tree$edge
  depth <- ape::node.depth.edgelength(tree)
  node_values <- rep(NA_real_, n_total)
  node_values[ntip + 1L] <- root_value

  # Procesar primero las ramas cuyos nodos parentales estan mas cerca de la raiz
  ord <- order(depth[edge_mat[, 1]])
  traj <- vector("list", length(ord))

  for (ii in seq_along(ord)) {
    e <- ord[ii]
    parent <- edge_mat[e, 1]
    child <- edge_mat[e, 2]
    L <- tree$edge.length[e]
    nstep <- max(1L, ceiling(L / dt))
    h <- L / nstep
    times <- seq(depth[parent], depth[parent] + L, length.out = nstep + 1L)
    x <- numeric(nstep + 1L)
    x[1] <- node_values[parent]

    # Regimen OUM: rama bajo el nodo del clado usa theta2.
    in_clade <- is_descendant_or_self(child, clade_node, edge_mat)
    theta_branch <- if (in_clade) theta2_OUM else theta1_OUM

    for (j in seq_len(nstep)) {
      z <- rnorm(1)
      tmid <- (times[j] + times[j + 1L]) / 2

      if (model == "BM") {
        x[j + 1L] <- x[j] + sigma_BM * sqrt(h) * z

      } else if (model == "OU1") {
        decay <- exp(-alpha_OU1 * h)
        x[j + 1L] <- theta_OU1 + (x[j] - theta_OU1) * decay +
          sigma_OU1 * sqrt((1 - exp(-2 * alpha_OU1 * h)) /
                             (2 * alpha_OU1)) * z

      } else if (model == "OUM") {
        decay <- exp(-alpha_OUM * h)
        x[j + 1L] <- theta_branch + (x[j] - theta_branch) * decay +
          sigma_OUM * sqrt((1 - exp(-2 * alpha_OUM * h)) /
                             (2 * alpha_OUM)) * z

      } else if (model == "EB") {
        # Early Burst: proceso tipo BM con tasa de evolucion dependiente del tiempo.
        # Varianza instantanea: sigma^2(t) = sigma_EB^2 * exp(beta_EB * t).
        # Con beta_EB < 0, los incrementos se hacen gradualmente menores.
        # No hay atraccion hacia un optimo: EB no es un modelo OU.
        sigma_t <- sigma_EB * exp(beta_EB * tmid / 2)
        x[j + 1L] <- x[j] + sigma_t * sqrt(h) * z
      }
    }

    node_values[child] <- x[length(x)]
    traj[[ii]] <- data.frame(
      model = model, edge = e, parent = parent, child = child,
      time = times, trait = x
    )
  }

  list(
    trajectory = do.call(rbind, traj),
    tips = data.frame(
      species = tree$tip.label,
      trait = node_values[seq_len(ntip)],
      model = model
    ),
    nodes = node_values
  )
}

simulate_one <- function(model) {
  simulate_model(
    model = model, tree = tree, dt = dt, root_value = root_value,
    sigma_BM = sigma_BM,
    sigma_OU1 = sigma_OU1, alpha_OU1 = alpha_OU1, theta_OU1 = theta_OU1,
    sigma_OUM = sigma_OUM, alpha_OUM = alpha_OUM,
    theta1_OUM = theta1_OUM, theta2_OUM = theta2_OUM,
    sigma_EB = sigma_EB, beta_EB = beta_EB,
    clade_node = clade_node
  )
}

models <- c("BM", "OU1", "OUM", "EB")

# -------------------------
# 4. Figura 1: una simulacion de trayectorias por modelo
# -------------------------
set.seed(20261020)
simulations <- lapply(models, simulate_one)
names(simulations) <- models
trajectory <- do.call(rbind, lapply(simulations, `[[`, "trajectory"))

p_trajectory <- ggplot(trajectory,
                       aes(x = time, y = trait, group = edge)) +
  geom_line(linewidth = 0.35, alpha = 0.7, colour = "#2878B5") +
  facet_wrap(~model, ncol = 2, scales = "free_y") +
  labs(
    x = "Tiempo desde el origen (unidades de tiempo)",
    y = "Valor del caracter",
    title = "Trayectorias fenotipicas simuladas sobre una filogenia",
    subtitle = paste0("n = ", n_tips,
                      " especies; misma topologia; semilla = 20261020")
  ) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
print(p_trajectory)
ggsave("figura_1_trayectorias_corregida.png", p_trajectory,
       width = 10, height = 7, dpi = 300)

# -------------------------
# 5. Figura 2: disparidad media a traves del tiempo
# -------------------------
# En cada instante, obtener exactamente un estado por linaje contemporaneo.
# Evitamos los instantes exactos de especiacion para no duplicar estados nodales.
times_grid <- seq(0.01, Tmax - 0.01, length.out = 100)
node_depth <- ape::node.depth.edgelength(tree)
internal_nodes <- (n_tips + 1L):(n_tips + tree$Nnode)
node_times <- node_depth[internal_nodes]
times_grid <- times_grid[
  vapply(times_grid, function(tt) all(abs(node_times - tt) > 1e-7),
         logical(1))
]

get_lineage_values <- function(traj_df, t) {
  vals <- numeric(0)
  for (e in unique(traj_df$edge)) {
    d <- traj_df[traj_df$edge == e, ]
    t0 <- min(d$time)
    t1 <- max(d$time)
    # Solo ramas que contienen t en su interior, no en los extremos.
    if (t > t0 + 1e-10 && t < t1 - 1e-10) {
      vals <- c(vals, approx(d$time, d$trait, xout = t)$y)
    }
  }
  vals
}

set.seed(20261021)
disp_list <- vector("list", length(models) * n_replicates)
counter <- 1L

for (m in models) {
  message("Disparidad: simulando ", m, " (", n_replicates, " replicas)")
  for (r in seq_len(n_replicates)) {
    sim_r <- simulate_one(m)
    d_r <- sim_r$trajectory
    for (tt in times_grid) {
      vals <- get_lineage_values(d_r, tt)
      if (length(vals) >= 2L) {
        disp_list[[counter]] <- data.frame(
          model = m, replicate = r, time = tt,
          disparity = stats::var(vals),
          n_lineages = length(vals)
        )
        counter <- counter + 1L
      }
    }
  }
}

disp_reps <- do.call(rbind, disp_list[seq_len(counter - 1L)])
disp_reps$model <- factor(disp_reps$model, levels = models)

summary_groups <- split(
  disp_reps,
  interaction(disp_reps$model, disp_reps$time, drop = TRUE)
)
disp_summary <- do.call(rbind, lapply(summary_groups, function(d) {
  data.frame(
    model = as.character(d$model[1]),
    time = d$time[1],
    mean = mean(d$disparity, na.rm = TRUE),
    lower = unname(quantile(d$disparity, 0.025, na.rm = TRUE)),
    upper = unname(quantile(d$disparity, 0.975, na.rm = TRUE))
  )
}))
disp_summary$model <- factor(disp_summary$model, levels = models)
disp_summary <- disp_summary[order(disp_summary$model, disp_summary$time), ]

p_disparity <- ggplot(disp_summary, aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#7FB3D5", alpha = 0.35) +
  geom_line(colour = "#21618C", linewidth = 1) +
  facet_wrap(~model, ncol = 2, scales = "free_y") +
  labs(
    x = "Tiempo desde el origen (unidades de tiempo)",
    y = "Varianza entre linajes contemporaneos",
    title = "Disparidad fenotipica a traves del tiempo",
    subtitle = paste0(n_replicates,
                      " replicas por modelo; media e intervalo central del 95%")
  ) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
print(p_disparity)
ggsave("figura_2_disparidad_corregida.png", p_disparity,
       width = 10, height = 7, dpi = 300)

# -------------------------
# 6. Exportar tablas
# -------------------------
write.csv(trajectory, "trayectorias_una_simulacion.csv", row.names = FALSE)
write.csv(disp_reps, "disparidad_todas_las_replicas.csv", row.names = FALSE)
write.csv(disp_summary, "disparidad_resumen.csv", row.names = FALSE)

# -------------------------
# 7. Ideas de ejercicios
# -------------------------
# 1. Cambiar n_tips (p.ej., 24, 96, 300).
# 2. En OU1, cambiar alpha_OU1 y observar la fuerza de atraccion al optimo.
# 3. En OUM, cambiar theta1_OUM y theta2_OUM y comparar los dos clados.
# 4. Cambiar alpha_OUM y sigma_OUM para ver la relacion atraccion/ruido.
# 5. En EB, comparar beta_EB = -6, -3, 0 y +1.
#    beta = 0 equivale a BM con tasa constante; beta < 0 desacelera la evolucion.
# 6. Aumentar n_replicates para estabilizar el resumen de disparidad.
#
# Interpretacion:
# BM: varianza suele aumentar con el tiempo.
# OU1: la varianza tiende a estabilizarse bajo un optimo compartido.
# OUM: la separacion depende de los dos optimos y de la asignacion de regimenes.
# EB: con beta muy negativo se espera mucha divergencia temprana y trayectorias
#     progresivamente mas horizontales. No hay atraccion hacia un optimo; los linajes
#     conservan las diferencias acumuladas al principio.
#
# Nota: son simulaciones didacticas. La disparidad es la varianza entre
# linajes contemporaneos y no una metrica DTT clasica de disparidad de subclados.
