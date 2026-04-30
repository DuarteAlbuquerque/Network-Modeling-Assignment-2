# =============================================================================
# Assignment 2 – Task 1: Geometrically Weighted Terms in ERGMs
# Glasgow School Data (Wave 1)
# =============================================================================

library(ergm)
library(network)
library(sna)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(42)

DATA <- "../Glasgow/"

# =============================================================================
# Load data
# =============================================================================
f1_mat <- as.matrix(read.csv(paste0(DATA, "f1.csv"), header = FALSE))
net1   <- as.network(f1_mat, directed = TRUE)

demographic <- read.csv(paste0(DATA, "demographic.csv"))
alcohol     <- read.csv(paste0(DATA, "alcohol.csv"))

net1 %v% "gender"  <- demographic$gender
net1 %v% "age"     <- demographic$age
net1 %v% "alcohol" <- alcohol$Wave.1

cat("=== Glasgow Wave 1 Network ===\n")
cat("Nodes           :", network.size(net1), "\n")
cat("Edges           :", network.edgecount(net1), "\n")
cat("Density         :", round(network.density(net1), 4), "\n")

outdeg <- degree(net1, cmode = "outdegree")
cat("Out-degree range:", min(outdeg), "-", max(outdeg), "\n")
cat("Mean out-degree :", round(mean(outdeg), 3), "\n\n")


# =============================================================================
# Part (1) – Change statistic of gwodegree
# =============================================================================
#
# (1.i)  Analytical derivation
# ---------------------------------------------------------------------------
# The GWD network statistic is:
#
#   s_GWD(x, alpha) = e^alpha * sum_{k=1}^{n} {1 - (1 - e^{-alpha})^k} * D_k^out(x)
#
# Adding a directed edge to a sender node that currently has out-degree k:
#   • D_k^out   decreases by 1   (node leaves degree class k)
#   • D_{k+1}^out increases by 1  (node enters degree class k+1)
#
# Change statistic:
#   delta_GWD(k, alpha)
#     = s_GWD(x + edge) - s_GWD(x)
#     = e^alpha * [ {1-(1-e^{-alpha})^{k+1}} - {1-(1-e^{-alpha})^k} ]
#     = e^alpha * [ (1-e^{-alpha})^k - (1-e^{-alpha})^{k+1} ]
#     = e^alpha * (1-e^{-alpha})^k * [1 - (1-e^{-alpha})]
#     = e^alpha * (1-e^{-alpha})^k * e^{-alpha}
#
#   => delta_GWD(k, alpha) = (1 - e^{-alpha})^k
#
# ---------------------------------------------------------------------------

delta_GWD <- function(k, alpha) (1 - exp(-alpha))^k

# (1.ii) Plot the change statistic
# ---------------------------------------------------------------------------
alpha_vals  <- c(-1, -0.5, 0, 0.1, 0.2, 0.3, 0.5, log(2), 1, 2, 3)
alpha_labs  <- c("-1", "-0.5", "0", "0.1", "0.2", "0.3", "0.5", "ln2", "1", "2", "3")
k_range     <- 0:max(outdeg)

df_plot <- do.call(rbind, lapply(seq_along(alpha_vals), function(i) {
  data.frame(
    k      = k_range,
    delta  = delta_GWD(k_range, alpha_vals[i]),
    alpha  = factor(alpha_labs[i], levels = alpha_labs),
    group  = ifelse(alpha_vals[i] < 0, "negative", "non-negative")
  )
}))

# Figure 1 – negative decay values (one figure each)
for (a_name in alpha_labs[alpha_vals < 0]) {
  df_sub <- df_plot %>% filter(alpha == a_name)
  p <- ggplot(df_sub, aes(x = k, y = delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(colour = "steelblue", linewidth = 0.9) +
    geom_point(colour = "steelblue", size = 2) +
    labs(
      title    = bquote("Change statistic," ~ alpha == .(a_name)),
      subtitle = "Negative decay value",
      x        = "Current out-degree k",
      y        = expression(delta[GWD](k, alpha))
    ) +
    theme_bw(base_size = 13)
  fname <- paste0("task1_fig_neg_alpha", gsub("\\.", "p", gsub("-", "m", a_name)), ".pdf")
  ggsave(fname, p, width = 6, height = 4)
  cat("Saved", fname, "\n")
}

# Figure 2 – all non-negative decay values in one figure
df_nonneg <- df_plot %>% filter(group == "non-negative")
cols      <- RColorBrewer::brewer.pal(max(3, nlevels(droplevels(df_nonneg$alpha))), "Set1")
# fallback if RColorBrewer not available
tryCatch({
  library(RColorBrewer)
}, error = function(e) NULL)

p_nonneg <- ggplot(df_nonneg, aes(x = k, y = delta, colour = alpha, linetype = alpha)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_colour_brewer(palette = "Set1", name = expression(alpha)) +
  scale_linetype_discrete(name = expression(alpha)) +
  labs(
    title    = "Change statistic - non-negative decay values",
    x        = "Current out-degree k",
    y        = expression(delta[GWD](k, alpha))
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

ggsave("task1_fig_nonneg_alpha.pdf", p_nonneg, width = 7, height = 5)
cat("Saved task1_fig_nonneg_alpha.pdf\n\n")


# =============================================================================
# Part (2) – Effect of gwodegree coefficient on out-degree distribution
# =============================================================================

# (2.i) M1: edges-only reference model
cat("=== M1: edges-only model ===\n")
M1 <- ergm(net1 ~ edges, control = control.ergm(seed = 42))
summary(M1)

# (2.ii) M2: fix gwodegree coefficient at -2 and +2, decay = 0.3
#   Re-estimate edges so expected density matches observed density.
#   Using offset() fixes the coefficient of gwodegree while ergm estimates edges.

cat("\n=== M2a: gwodegree coef = -2, decay = 0.3 ===\n")
M2a <- ergm(
  net1 ~ edges + offset(gwodegree(0.3, fixed = TRUE)),
  offset.coef  = c(-2),
  control      = control.ergm(seed = 42)
)
summary(M2a)

cat("\n=== M2b: gwodegree coef = +2, decay = 0.3 ===\n")
M2b <- ergm(
  net1 ~ edges + offset(gwodegree(0.3, fixed = TRUE)),
  offset.coef  = c(2),
  control      = control.ergm(seed = 42)
)
summary(M2b)

# (2.iii) Simulate 100 networks from each model
NSIM <- 100

cat("\nSimulating", NSIM, "networks from M1 ...\n")
sim_M1  <- simulate(M1,  nsim = NSIM, output = "network")

cat("Simulating", NSIM, "networks from M2a ...\n")
sim_M2a <- simulate(M2a, nsim = NSIM, output = "network")

cat("Simulating", NSIM, "networks from M2b ...\n")
sim_M2b <- simulate(M2b, nsim = NSIM, output = "network")

# Compute node count by out-degree for each simulated network
outdeg_counts <- function(sim_list, label) {
  k_max <- max(outdeg)          # use observed max as a common ceiling
  do.call(rbind, lapply(seq_len(NSIM), function(i) {
    degs <- degree(sim_list[[i]], cmode = "outdegree")
    tab  <- tabulate(degs + 1L, nbins = k_max + 1L)   # +1: degree 0 -> index 1
    data.frame(outdegree = 0:k_max, count = tab, sim = i, model = label)
  }))
}

df_sim <- bind_rows(
  outdeg_counts(sim_M1,  "M1: edges only"),
  outdeg_counts(sim_M2a, "M2a: theta_GWD = -2"),
  outdeg_counts(sim_M2b, "M2b: theta_GWD = +2")
) %>%
  mutate(model = factor(model,
                        levels = c("M1: edges only",
                                   "M2a: theta_GWD = -2",
                                   "M2b: theta_GWD = +2")))

# Observed node counts by out-degree
obs_tab <- data.frame(
  outdegree = 0:max(outdeg),
  count     = tabulate(outdeg + 1L, nbins = max(outdeg) + 1L)
)

# Violin plot
p_violin <- ggplot(df_sim, aes(x = factor(outdegree), y = count, fill = model)) +
  geom_violin(scale = "width", alpha = 0.65, colour = NA) +
  geom_point(
    data      = obs_tab,
    aes(x = factor(outdegree), y = count),
    inherit.aes = FALSE,
    shape = 23, size = 3, fill = "black", colour = "white"
  ) +
  facet_wrap(~ model, ncol = 1) +
  scale_fill_manual(
    values = c("M1: edges only"       = "#377EB8",
               "M2a: theta_GWD = -2"  = "#E41A1C",
               "M2b: theta_GWD = +2"  = "#4DAF4A")
  ) +
  labs(
    title    = "Out-degree distribution: 100 simulated networks vs. observed",
    subtitle = "Black diamonds = observed node counts",
    x        = "Out-degree",
    y        = "Number of nodes",
    fill     = "Model"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold")
  )

ggsave("task1_outdeg_violin.pdf", p_violin, width = 9, height = 11)
cat("\nSaved task1_outdeg_violin.pdf\n")

# Summary statistics for interpretation
cat("\n=== Summary: mean out-degree in simulations ===\n")
df_sim %>%
  group_by(model, sim) %>%
  summarise(mean_outdeg = sum(outdegree * count) / sum(count), .groups = "drop") %>%
  group_by(model) %>%
  summarise(
    mean = round(mean(mean_outdeg), 3),
    sd   = round(sd(mean_outdeg), 4),
    .groups = "drop"
  ) %>%
  print()

cat("\nObserved mean out-degree:", round(mean(outdeg), 3), "\n")

cat("\n=== Task 1 complete. ===\n")

