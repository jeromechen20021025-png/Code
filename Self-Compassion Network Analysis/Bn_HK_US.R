# Uncomment the next 3 lines if you need to update R. This will open an installer that will prompt you to confirm the update:
# install.packages("installr")
# library(installr)
# updateR()

# Run this to check R Version:
# R.version.string # R should be 4.4.2 for bnlearn to work. 
# If it isn't, I suggest uninstalling older versions and installing this one

# Load necessary packages:
if (!require(pacman)) install.packages("pacman")
if (!require(writexl)) install.packages("writexl")
if (!require(BNSL)) install.packages("BNSL")
if (!require(philentropy)) install.packages("philentropy")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(pacman)
pacman::p_load(networktools, bootnet, IsingFit, qgraph, mgm, NetworkComparisonTest, dplyr, readxl, huge, bnlearn, Rgraphviz, entropy, gRain, BNSL, philentropy, gt, webshot2)
BiocManager::install("Rgraphviz", ask = FALSE, force = TRUE)

# Load necessary libraries
library(qgraph)
library(bootnet)
library(dplyr)
library(readxl)
library(networktools)
library(bnlearn)
library(writexl)
library(Rgraphviz)

# Load the dataset
rm(list = ls()) # Clear environment
gc()

set.seed(123) # For reproducibility

data <- read_excel("internal.xlsx")
x <- data

# Function to estimate and plot the network
estimate_and_plot_network <- function(data, group_name, color_vector) {
  network <- estimateNetwork(data, default = "EBICglasso", 
                             tuning = 0.5, lambda.min.ratio = 0.01) # Adjust parameters
  
  return(network)
  print(network)
  
  pdf(file = paste0(group_name, ".pdf"))
  qgraph(network$graph, layout = "spring", color = color_vector)
  dev.off()
}

# Function to estimate and plot Bayesian network
estimate_and_plot_bn <- function(data, group_name) {
  bn <- hc(data)  # Using hill-climbing algorithm for Bayesian network estimation
  pdf(file = paste0(group_name, "_bn.pdf"))
  graphviz.plot(bn, main = paste("Bayesian Network for", group_name), layout = "dot")
  dev.off()
  return(bn)
}

### US
xcurr <- x %>% filter(region == 5) %>% select(CS_K:SCSY_o) 
color_vector_us <- c(rep("Blue", 12), rep("Red", 12), rep("Green", 12), rep("Yellow", 12), rep("Red", 4))
us <- estimate_and_plot_network(xcurr, "us", color_vector_us)

# Apply goldbricker function in the R-package networktools on the weighted adjacency matrix
us_gb <- goldbricker(xcurr)

# Check if suggested reductions are character
if (!is.null(names(us_gb$suggested_reductions)) && all(sapply(us_gb$suggested_reductions, is.character))) {
  xcurr_reduced <- net_reduce(data = xcurr, badpairs = us_gb, method = "best_goldbricker")
} else {
  warning("No suggested reductions found or they are not in character format.")
  xcurr_reduced <- xcurr
}

us_reduced <- estimate_and_plot_network(xcurr_reduced, "us_reduced", color_vector_us)

cat("Network has", us$nPerson, "cases.")
W <- us$graph
W # Weighted adjacency matrix
centralityPlot(us, include = "Strength", scale = "z-scores", orderBy = "Strength")
centralityPlot(us, include = "ExpectedInfluence", scale = "z-scores", orderBy = "ExpectedInfluence")
centralityPlot(us, include = "All", scale = "z-scores", orderBy = "Strength")
boot <- bootnet(xcurr, default = "EBICglasso", nBoots = 1000, nCores = 8, type = "case", statistics = c("strength", "closeness", "betweenness", "ExpectedInfluence"))
corStability(boot)
plot(boot, statistics = "strength")

# New graphs
boot1 <- bootnet(us, nBoots = 1000, nCores = 8, statistics = c("edge", "strength", "expectedInfluence"))
plot(boot1, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
plot(boot1, "strength")
plot(boot1, "expectedInfluence")

### Hong Kong
xcurr <- x %>% filter(region == 3) %>% select(CS_K:SCSY_o) 
color_vector_hk <- c(rep("Blue", 12), rep("Red", 12), rep("Green", 12), rep("Yellow", 12), rep("Red", 4))
hk <- estimate_and_plot_network(xcurr, "hk", color_vector_hk)

# Apply goldbricker function in the R-package networktools on the weighted adjacency matrix
hk_gb <- goldbricker(xcurr)

# Check if suggested reductions are character
if (!is.null(names(hk_gb$suggested_reductions)) && all(sapply(hk_gb$suggested_reductions, is.character))) {
  xcurr_reduced <- net_reduce(data = xcurr, badpairs = hk_gb, method = "best_goldbricker")
} else {
  warning("No suggested reductions found or they are not in character format.")
  xcurr_reduced <- xcurr
}

hk_reduced <- estimate_and_plot_network(xcurr_reduced, "hk_reduced", color_vector_hk)

cat("Network has", hk$nPerson, "cases.")
W <- hk$graph
W # Weighted adjacency matrix
centralityPlot(hk, include = "Strength", scale = "z-scores", orderBy = "Strength")
centralityPlot(hk, include = "ExpectedInfluence", scale = "z-scores", orderBy = "ExpectedInfluence")
centralityPlot(hk, include = "All", scale = "z-scores", orderBy = "Strength")
boot <- bootnet(xcurr, default = "EBICglasso", nBoots = 1000, nCores = 8, type = "case", statistics = c("strength", "closeness", "betweenness", "ExpectedInfluence"))
corStability(boot)
plot(boot, statistics = "strength")

# New graphs
boot1 <- bootnet(hk, nBoots = 1000, nCores = 8, statistics = c("edge", "strength", "expectedInfluence"))
plot(boot1, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
plot(boot1, "strength")
plot(boot1, "expectedInfluence")

### Comparison between High_ability and Non_high_ability
centralityPlot(list(US = us, HK = hk), include = c("Strength", "ExpectedInfluence"), scale = "z-scores", orderBy = "Strength")

## Tests
difference_value <- function(NCT) {
  diff_edges <- NCT$einv.pvals  # Use the entire edge invariance p-values table without filtering
  
  for (i in 1:nrow(diff_edges)) {
    var_1 <- as.character(diff_edges[i, 1])
    var_2 <- as.character(diff_edges[i, 2])
    
    # Debugging output to check for NAs and the variables being checked
    cat("Checking edge:", var_1, "and", var_2, "\n")
    
    # Check if var_1 and var_2 are not NA
    if (!is.na(var_1) && !is.na(var_2)) {
      # Check if the variables exist in the network matrices
      if (var_1 %in% rownames(NCT$nw1) && var_2 %in% colnames(NCT$nw1)) {
        value_net_1 <- NCT$nw1[var_1, var_2]
        value_net_2 <- NCT$nw2[var_1, var_2]
        
        abs_difference <- abs(value_net_1 - value_net_2)
        p_value <- diff_edges$`p-value`[i]
        
        cat("Test Edge", i, "\n----\n")
        cat(var_1, "and", var_2)
        cat("\nNetwork 1:", value_net_1, "\nNetwork 2:", value_net_2)
        cat("\nAbsolute difference:", abs_difference, "with p-value =", p_value, "\n----\n")
      } else {
        cat("Test Edge", i, "\n----\n")
        cat(var_1, "and", var_2)
        cat("\nEdge not found in one or both networks.")
        cat("\n----\n")
      }
    } else {
      cat("Test Edge", i, "\n----\n")
      cat(var_1, "and", var_2)
      cat("\nOne or both variables are NA.")
      cat("\n----\n")
    }
  }
}

niter <- 1000
nct_tests <- NCT(us, hk, it = niter, test.edges = TRUE, test.centrality = TRUE, p.adjust.methods = "BH", progressbar = TRUE)
#summary(nct_tests)
plot(nct_tests, what = "strength")
plot(nct_tests, what = "network")
plot(nct_tests, what = "edge")
plot(nct_tests, what = "centrality")
#difference_value(nct_tests)

### Bayesian Network Comparison

# Combine data with the group indicator for bnlearn comparison
data_bn <- x %>% select(CS_K:SCSY_o, region)

# Function to calculate centrality measures for Bayesian networks
calculate_centrality_measures <- function(bn, data, criterion) {
  nodes <- nodes(bn)
  
  # Calculate degree centrality
  arc_strength <- arc.strength(bn, data = data, criterion = criterion)
  degree <- sapply(nodes, function(node) {
    sum(arc_strength$from == node | arc_strength$to == node)
  })
  
  # Calculate betweenness centrality
  bn_arcs <- arcs(bn)
  graph <- igraph::graph_from_data_frame(bn_arcs, directed = TRUE)
  betweenness <- igraph::betweenness(graph, v = nodes, directed = TRUE)
  
  # Calculate closeness centrality
  closeness <- igraph::closeness(graph, v = nodes, mode = "all")
  
  centrality_measures <- data.frame(Degree = degree, Betweenness = betweenness, Closeness = closeness)
  rownames(centrality_measures) <- nodes
  return(centrality_measures)
}

# Function to plot centrality measures
plot_centrality_measures <- function(centrality_measures, group_name) {
  par(mfrow = c(1, 3))
  
  barplot(centrality_measures$Degree, main = paste("Degree Centrality for", group_name), las = 2)
  barplot(centrality_measures$Betweenness, main = paste("Betweenness Centrality for", group_name), las = 2)
  barplot(centrality_measures$Closeness, main = paste("Closeness Centrality for", group_name), las = 2)
  
  par(mfrow = c(1, 1))
}

# us Bayesian Network
data_bn_us <- data_bn %>% filter(region == 5) %>% select(-region)
us_bn <- estimate_and_plot_bn(data_bn_us, "us")

# Apply goldbricker function in the R-package networktools on the weighted adjacency matrix
us_bn_gb <- goldbricker(data_bn_us)

# Check if suggested reductions are character
if (!is.null(names(us_bn_gb$suggested_reductions)) && all(sapply(us_bn_gb$suggested_reductions, is.character))) {
  data_bn_us_reduced <- net_reduce(data = data_bn_us, badpairs = us_bn_gb, method = "best_goldbricker")
} else {
  warning("No suggested reductions found or they are not in character format.")
  data_bn_us_reduced <- data_bn_us
}

us_bn_reduced <- estimate_and_plot_bn(data_bn_us_reduced, "us_reduced")

# Calculate and plot centrality measures
us_centrality <- calculate_centrality_measures(us_bn, data_bn_us, "bic-g")
plot_centrality_measures(us_centrality, "US Group")
us_centrality_reduced <- calculate_centrality_measures(us_bn_reduced, data_bn_us_reduced, "bic-g")
plot_centrality_measures(us_centrality_reduced, "US Group Reduced")

# Bootstrapping and plotting results
boot_us_bn <- boot.strength(data_bn_us, algorithm = "hc", R = 1000)
boot_avg_us_bn <- averaged.network(boot_us_bn)
pdf(file = "us_bn_bootstrapped.pdf")
graphviz.plot(boot_avg_us_bn, main = "Bootstrapped Bayesian Network for US Group", layout = "dot")
dev.off()

# Obtain edge strength 



create_strength_table <- function(strength_data) {
  strength_data <- strength_data[order(-abs(strength_data$strength)), ]
  n_rows <- ceiling(nrow(strength_data) / 3)
  
  make_column <- function(start_idx, end_idx) {
    if(start_idx > nrow(strength_data)) return(NULL)
    end_idx <- min(end_idx, nrow(strength_data))
    subset <- strength_data[start_idx:end_idx, ]
    list(
      number = start_idx:end_idx,
      direction = paste(subset$from, "\u2192", subset$to),
      weight = subset$strength
    )
  }
  
  cols <- list(
    make_column(1, n_rows),
    make_column(n_rows + 1, 2 * n_rows),
    make_column(2 * n_rows + 1, 3 * n_rows)
  )
  
  result <- data.frame(
    check.names = FALSE,
    " " = cols[[1]]$number,
    "Association direction" = cols[[1]]$direction,
    "Weight" = cols[[1]]$weight
  )
  
  if(!is.null(cols[[2]])) {
    result[["  "]] = cols[[2]]$number
    result[["Association direction "]] = cols[[2]]$direction
    result[["Weight "]] = cols[[2]]$weight
  }
  
  if(!is.null(cols[[3]])) {
    result[["   "]] = cols[[3]]$number
    result[["Association direction  "]] = cols[[3]]$direction
    result[["Weight  "]] = cols[[3]]$weight
  }
  
  gt(result) %>%
    tab_header(title = "Estimated edge weights in Bayesian network") %>%
    fmt_number(columns = contains("Weight"), decimals = 2) %>%
    tab_style(
      style = cell_borders(sides = "all", color = "gray", weight = px(1)),
      locations = cells_body()
    ) %>%
    tab_options(
      table.border.top.style = "hidden",
      table.border.bottom.style = "hidden"
    )
}

strength_us_bn <- arc.strength(us_bn, data_bn_us)
strength_us_bn_table <- create_strength_table(strength_us_bn) # Equivalent to content in Table 4 of Sheida Moradi et al.
gtsave(strength_us_bn_table, filename = "strength_us_bn_table.png")

# plot DAG with edge thickness as strength (Fig.4 in Sheida Moradi et al.)
pdf(file = "strength_us_bn.pdf")
strength.plot(us_bn,strength_us_bn, main = paste("Bayesian Network for US (Edge thickness represent strength)"), layout = "dot")
dev.off()

# Hong Kong Bayesian Network
data_bn_hk <- data_bn %>% filter(region == 3) %>% select(-region)
hk_bn <- estimate_and_plot_bn(data_bn_hk, "hk")

# Apply goldbricker function in the R-package networktools on the weighted adjacency matrix
hk_bn_gb <- goldbricker(data_bn_hk)

# Check if suggested reductions are character
if (!is.null(names(hk_bn_gb$suggested_reductions)) && all(sapply(hk_bn_gb$suggested_reductions, is.character))) {
  data_bn_hk_reduced <- net_reduce(data = data_bn_hk, badpairs = hk_bn_gb, method = "best_goldbricker")
} else {
  warning("No suggested reductions found or they are not in character format.")
  data_bn_hk_reduced <- data_bn_hk
}

hk_bn_reduced <- estimate_and_plot_bn(data_bn_hk_reduced, "hk_reduced")

# Calculate and plot centrality measures
hk_centrality <- calculate_centrality_measures(hk_bn, data_bn_hk, "bic-g")
plot_centrality_measures(hk_centrality, "HK Group")
hk_centrality_reduced <- calculate_centrality_measures(hk_bn_reduced, data_bn_hk_reduced, "bic-g")
plot_centrality_measures(hk_centrality_reduced, "HK Group Reduced")

# Bootstrapping and plotting results
boot_hk_bn <- boot.strength(data_bn_hk, algorithm = "hc", R = 1000)
boot_avg_hk_bn <- averaged.network(boot_hk_bn)
pdf(file = "hk_bn_bootstrapped.pdf")
graphviz.plot(boot_avg_hk_bn, main = "Bootstrapped Bayesian Network for HK Group", layout = "dot")
dev.off()

# Obtain edge strength 
strength_hk_bn <- arc.strength(hk_bn, data_bn_hk)
strength_hk_bn_table <- create_strength_table(strength_hk_bn) # Equivalent to content in Table 4 of Sheida Moradi et al.
gtsave(strength_hk_bn_table, filename = "strength_hk_bn_table.png")
# Plot DAG with edge thickness as strength (Fig.4 in Sheida Moradi et al.)
pdf(file = "strength_hk_bn.pdf")
strength.plot(hk_bn, strength_hk_bn, main = paste("Bayesian Network for HK (Edge thickness represent strength)"), layout = "dot")
dev.off()

# Checking if the learned networks are indeed DAG:
# US
acyclic(us_bn)
directed(us_bn)
acyclic(boot_avg_us_bn)
directed(boot_avg_us_bn) # There are no undirected edges.
undirected.arcs(boot_avg_us_bn) # There are no undirected edges.
print("US BN and its bootsrapped average are DAGs!")

# HK
directed(hk_bn)
acyclic(hk_bn)
directed(boot_avg_hk_bn) # There are no undirected edges.
acyclic(boot_avg_hk_bn) # There are no undirected edges.
print("HK BN and its bootsrapped average are DAGs!")

# Network Comparison Test
NCT(data_bn_us, data_bn_hk) # Notice this is implemented using the data used for each BN. NCT function can either take data sets or objects estimated using estimateNetwork from bootnet. It doesn't work on bn objects.

# Alternatives: Compare DAG network figures or function native to bnlearn
graphviz.compare(us_bn, hk_bn, 
                 main = c("", "Comparison from US to HK (pictured)"), 
                 diff = "from-first", 
                 sub = c("", "Edges missing or with different direction in US but present in HK are in red. \nEdges present in US but not present in HK are in dashed blue. \nEdges present in both are black."))

graphviz.compare(hk_bn, us_bn,
                 main = c("", "Comparison from HK to US (pictured)"), 
                 diff = "from-first", 
                 sub = c("", "Edges missing or with different direction in HK but present in US are in red. \nEdges present in HK but not present in US are in dashed blue. \nEdges present in both are black."))

all.equal(us_bn, hk_bn) # Checks if two networks are identical
hamming(us_bn, hk_bn) # Hamming distance between two networks
shd(us_bn, hk_bn) # Structural Hamming distance (SHD) between two networks
bnlearn::compare(target = us_bn, current = hk_bn) # Display differences in networks
# We can also list which edges are in each category:
bnlearn::compare(target = us_bn, current = hk_bn, arcs = TRUE)

################################################################################
##################### New Additions Below ######################################
################################################################################
us_bn_fit <- bn.fit(us_bn,data_bn_us)
hk_bn_fit <- bn.fit(hk_bn,data_bn_hk)

# Functions 
run_shd <- function(bn1, bn2) {
  return(shd(bn1, bn2))
}

run_hamming <- function(bn1, bn2) {
  return(hamming(bn1, bn2))
}

run_graph_edit_distance <- function(bn1, bn2) {
  bn1_igraph <- as.igraph(bn1)
  bn2_igraph <- as.igraph(bn2)
  return(igraph::graph.edit.distance(bn1_igraph, bn2_igraph))
}

run_arc_strength <- function(bn, data) {
  return(arc.strength(bn, data = data))
}

run_centrality_measures <- function(igraph_obj) {
  list(
    degree = igraph::degree(igraph_obj),
    betweenness = igraph::betweenness(igraph_obj),
    closeness = igraph::closeness(igraph_obj)
  )
}

run_cpquery <- function(bn, event, evidence, n) {
  return(cpquery(bn, event = event, evidence = evidence, n = n))
}

run_bootstrap <- function(data, R, algorithm, cluster = 8) {
  bootstrap_strength <- boot.strength(data, R = R, algorithm = algorithm)
  return(bootstrap_strength)
}

run_ave_boot_bn <- function(bootstrap_strength, threshold) {
  return(averaged.network(bootstrap_strength, threshold = threshold))
}

compare_arcs <- function(bn1, bn2) {
  # Use bnlearn's compare() function
  comparison <- bnlearn::compare(bn1, bn2)
  
  # Extract metrics
  tp <- comparison$tp  # True positives: arcs present in both
  fp <- comparison$fp  # False positives: arcs in bn2 but not in bn1
  fn <- comparison$fn  # False negatives: arcs in bn1 but not in bn2
  
  # Return results as a list
  list(TP = tp, FP = fp, FN = fn)
}

compare_cpdag <- function(bn1, bn2) {
  cpdag1 <- cpdag(bn1)
  cpdag2 <- cpdag(bn2)
  identical(cpdag1, cpdag2)  # Returns TRUE if CPDAGs are equivalent
}

compare_parameters <- function(fit1, fit2, node) {
  params1 <- coef(fit1[[node]])
  params2 <- coef(fit2[[node]])
  list(Node = node, Parameters1 = params1, Parameters2 = params2)
}

# Calculate Structural Hamming Distance (SHD) between US and HK Bayesian Networks
shd_value <- run_shd(us_bn, hk_bn)
cat("The Structural Hamming Distance between the US and Hong Kong Bayesian Networks is:", shd_value, "\n")

# Calculate Hamming Distance between US and HK Bayesian Networks
hamming_value <- run_hamming(us_bn, hk_bn)
cat("The Hamming Distance between the US and Hong Kong Bayesian Networks is:", hamming_value, "\n")

# Parameter Comparison: Arc Strength Comparison between US and HK Bayesian Networks
us_arc_strength <- run_arc_strength(us_bn, data = data_bn_us)
hk_arc_strength <- run_arc_strength(hk_bn, data = data_bn_hk)
cat("Arc Strengths for the US Bayesian Network:\n")
print(us_arc_strength)
cat("Arc Strengths for the Hong Kong Bayesian Network:\n")
print(hk_arc_strength)

# Likelihood and Score Comparison between US and HK Bayesian Networks
us_score <- score(us_bn, data = data_bn_us, type = "bic-g")
hk_score <- score(hk_bn, data = data_bn_hk, type = "bic-g")
cat("The BIC Score for the US Bayesian Network is:", us_score, "\n")
cat("The BIC Score for the Hong Kong Bayesian Network is:", hk_score, "\n")

# Simulation-Based Validation: Generate Simulated Data for Bayesian Network Validation
simulated_obs <- 10^6 # simulate 1 million observations
#US
simulated_data_us <- rbn(us_bn_fit, n = simulated_obs)
simulated_bn_us <- hc(simulated_data_us)
# Calculate Structural Hamming Distance between the original and simulated US Bayesian Networks
shd_simulated_us <- shd(us_bn, simulated_bn_us)
cat("The Structural Hamming Distance between the original and simulated US Bayesian Networks is:", shd_simulated_us, "\n")
#HK
simulated_data_hk <- rbn(hk_bn_fit, n = simulated_obs)
simulated_bn_hk <- hc(simulated_data_hk)
# Calculate Structural Hamming Distance between the original and simulated US Bayesian Networks
shd_simulated_hk <- shd(hk_bn, simulated_bn_hk)
cat("The Structural Hamming Distance between the original and simulated HK Bayesian Networks is:", shd_simulated_us, "\n")
# Calculate Structural Hamming Distance  between simulated BNs
shd_simulated_us_hk <- shd(simulated_bn_us, simulated_bn_hk)
cat("The Structural Hamming Distance between the simulated US and HK Bayesian Networks is:", shd_simulated_us_hk, "\n")
# Plot BNs from simulated dta with arc strength:
pdf(file = "simulated_strength_us_bn.pdf")
strength.plot(simulated_bn_us, strength = run_arc_strength(simulated_bn_us, data = simulated_data_us), main = paste("Bayesian Network for Simulated US (Edge thickness represent strength)"), layout = "dot")
dev.off()
pdf(file = "simulated_strength_us_bn.pdf")
strength.plot(simulated_bn_hk, strength = run_arc_strength(simulated_bn_hk, data = simulated_data_hk), main = paste("Bayesian Network for Simulated HK (Edge thickness represent strength)"), layout = "dot")
dev.off()

# Probabilistic Measures: Conditional Probability Comparison
# Query the conditional probability for a specific event in the US Bayesian Network
# The query below is the probability, given the fitted US BN, that CS_CM is higher than 17 (its 3rd quartile), 
# conditional on PSA_total being higher than its median value of 56. This is obtained from simulating 10k samples from the fitted BN
cp_us <- cpquery(us_bn_fit, event = (CS_CM >= 17), evidence = (PSA_total >= 56), n = 10000)
cat("The Conditional Probability of CS_CM being above 17 given PSA_total above 56 in the US Bayesian Network is:", cp_us, "\n")

# Query the conditional probability for a specific event in the HK Bayesian Network
cp_hk <- cpquery(hk_bn_fit, event = (CS_CM >= 17), evidence = (PSA_total >= 56), n = 10000)
cat("The Conditional Probability of CS_CM being above 17 given PSA_total is above 56 in the HK Bayesian Network is:", cp_hk, "\n")

# Probabilistic Measures: simulated values of one value, conditional on values of another variable
# US
cpdist_us <- cpdist(us_bn_fit, nodes = "CS_CM", evidence = (PSA_total >= 56), n = 10000)
hist(cpdist_us)
# HK
cpdist_hk <- cpdist(hk_bn_fit, nodes = "CS_CM", evidence = (PSA_total >= 56), n = 10000)
hist(cpdist_hk)

# Node-Level Comparison: Centrality Analysis
# Convert Bayesian networks to igraph objects for centrality calculation
us_igraph <- as.igraph(us_bn)
hk_igraph <- as.igraph(hk_bn)

# Calculate centrality measures for the US Bayesian Network
degree_us <- igraph::degree(us_igraph)
betweenness_us <- igraph::betweenness(us_igraph)
closeness_us <- igraph::closeness(us_igraph)
# Print centrality measures for the US Bayesian Network
cat("Degree Centrality for the US Bayesian Network:\n")
print(degree_us)
cat("Betweenness Centrality for the US Bayesian Network:\n")
print(betweenness_us)
cat("Closeness Centrality for the US Bayesian Network:\n")
print(closeness_us)

# Calculate centrality measures for the Hong Kong Bayesian Network
degree_hk <- igraph::degree(hk_igraph)
betweenness_hk <- igraph::betweenness(hk_igraph)
closeness_hk <- igraph::closeness(hk_igraph)
# Print centrality measures for the Hong Kong Bayesian Network
cat("Degree Centrality for the Hong Kong Bayesian Network:\n")
print(degree_hk)
cat("Betweenness Centrality for the Hong Kong Bayesian Network:\n")
print(betweenness_hk)
cat("Closeness Centrality for the Hong Kong Bayesian Network:\n")
print(closeness_hk)

# Bootstrapping and Stability: Bootstrap with Confidence Thresholds
threshold <- 0.85  # Threshold for arc inclusion (e.g., 85% strength)
# Bootstrap the Bayesian network structure for the US Bayesian Network
bootstrap_us <- run_bootstrap(data_bn_us, R = 1000, algorithm = "hc")
# Print bootstrap results for the US Bayesian Network
cat("Bootstrap Strength Results for the US Bayesian Network:\n")
print(bootstrap_us)
# Plot the averaged network with bootstrap confidence
us_avg_network <- run_ave_boot_bn(bootstrap_us, threshold = threshold)
cat("Plotting Averaged Network for the US Bayesian Network with Bootstrap Confidence:\n")
pdf(file = "us_avg_bootstrap_network.pdf")
graphviz.plot(us_avg_network, main = "Averaged Bayesian Network for US with Bootstrap Confidence")
dev.off()
# Bootstrap the Bayesian network structure for the HK Bayesian Network
bootstrap_hk <- run_bootstrap(data_bn_hk, R = 1000, algorithm = "hc")
# Print bootstrap results for the HK Bayesian Network
cat("Bootstrap Strength Results for the HK Bayesian Network:\n")
print(bootstrap_hk)
# Plot the averaged network with bootstrap confidence
hk_avg_network <- run_ave_boot_bn(bootstrap_hk, threshold = threshold)
cat("Plotting Averaged Network for the HK Bayesian Network with Bootstrap Confidence:\n")
pdf(file = "hk_avg_bootstrap_network.pdf")
graphviz.plot(hk_avg_network, main = "Averaged Bayesian Network for HK with Bootstrap Confidence")
dev.off()

# Visual and Graphical Comparisons: Graph Overlay
# Compare the US and Hong Kong Bayesian Networks by overlaying them
graphviz.compare(us_avg_network, hk_avg_network, 
                 main = c("", "Comparison from US to HK Bootstrapped Averaged Networks"),
                 diff = "from-first", 
                 sub = c("", "Edges missing or with different direction in US but present in HK are in red. \nEdges present in US but not present in HK are in dashed blue. \nEdges present in both are black."))

graphviz.compare(hk_avg_network, us_avg_network,
                 main = c("", "Comparison from HK to US Bootstrapped Averaged Networks"),
                 diff = "from-first", 
                 sub = c("", "Edges missing or with different direction in HK but present in US are in red. \nEdges present in HK but not present in US are in dashed blue. \nEdges present in both are black."))

# Structural comparisons
# Compare if the networks are structurally identical
all.equal(us_avg_network, hk_avg_network)

#Performance-based comparisons
#Cross-validation
# Perform 10-fold cross-validation for predictive performance
cv_us <- bn.cv(data_bn_us, us_bn, method = "k-fold", k = 10)
cv_hk <- bn.cv(data_bn_hk, hk_bn, method = "k-fold", k = 10)

print("Cross-validation scores:")
print(cv_us)
print(cv_hk)

#Log-Likelihood and Scoring Functions

# Log-likelihood comparison
loglik_us <- logLik(us_bn_fit, data_bn_us)
loglik_hk <- logLik(hk_bn_fit, data_bn_hk)
aic_us <- AIC(us_bn_fit, data_bn_us)
bic_us <- BIC(us_bn_fit, data_bn_us)
aic_hk <- AIC(hk_bn_fit, data_bn_hk)
bic_hk <- BIC(hk_bn_fit, data_bn_hk)

print(paste0("US> Log-likelihood: ", round(loglik_us,2), "; AIC: ", round(aic_us,2), "; BIC: ", round(bic_us,2),"."))
print(paste0("HK> Log-likelihood: ", round(loglik_hk,2), "; AIC: ", round(aic_hk,2), "; BIC: ", round(bic_hk,2),"."))

# Comparison with Reference Network
# Find Consensus Network (taking only arcs that show up in both)
# Extract arcs from both networks
us_arcs <- arcs(us_bn)
hk_arcs <- arcs(hk_bn)
colnames(us_arcs) <- colnames(hk_arcs) <- c("from", "to")
us_arcs <- us_arcs[order(us_arcs[, 1], us_arcs[, 2]), ]
hk_arcs <- hk_arcs[order(hk_arcs[, 1], hk_arcs[, 2]), ]
common_arcs <- merge(as.data.frame(us_arcs), as.data.frame(hk_arcs))
nodes <- unique(c(common_arcs$from, common_arcs$to))
consensus_bn <- empty.graph(nodes)
arcs(consensus_bn) <- as.matrix(common_arcs)
is_dag <- directed(consensus_bn) && acyclic(consensus_bn)
print(is_dag)
graphviz.plot(consensus_bn, main = "Consensus Network")

#Robustness and Sensitivity Analysis
# Average edge strength
avg_strength_us <- mean(bootstrap_us$strength)
avg_strength_hk <- mean(bootstrap_hk$strength)

print(paste("Average edge strength for us_bn:", round(avg_strength_us,4)))
print(paste("Average edge strength for hk_bn:", round(avg_strength_hk,4)))

#Graphical Comparisons
# Adjust outer margins to create space for the overarching title
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2), oma = c(0, 0, 5, 0))  # Add space at the top via oma
# Plot the two networks
strength.plot(us_bn, us_arc_strength, main = "US")  # Plot the first Bayesian network
strength.plot(hk_bn, hk_arc_strength, main = "HK")  # Plot the second Bayesian network
# Add a single overarching title
mtext("Bayesian Networks with Arc Strength for US and HK", side = 3, line = 1, outer = TRUE, cex = 1.5)

#Mutual Information Between Networks
# Compute mutual information matrices
mi_us <- mi_matrix(data_bn_us)
mi_hk <- mi_matrix(data_bn_hk)
# Compare the matrices
difference_mi <- mi_us - mi_hk
print(difference_mi)

#Path-Based Similarity
# Example: Check if a causal path from CS_K to SCSY_m exists in both networks
path_in_us <- path.exists(us_bn, from="CS_K", to="SCSY_cm", direct = FALSE)
path_in_hk <- path.exists(hk_bn, from="CS_K", to="SCSY_cm", direct = FALSE)

# Compare paths
print(path_in_us == path_in_hk)

#Conditional Independence Tests
# Test independence of SCSY_cm and SCSY_sk, conditional on CS_I
ci_test_us <- ci.test("SCSY_cm", "SCSY_sk", "CS_I", data = data_bn_us, test = "mi-g")
ci_test_hk <- ci.test("SCSY_cm", "SCSY_sk", "CS_I", data = data_bn_hk, test = "mi-g")

print(ci_test_us)
print(ci_test_hk)

#Predictive Performance on Held-Out Data
# Split data into train and test sets
train_idx <- sample(seq_len(nrow(x)), size = 0.7 * nrow(x))
train_data <- x[train_idx, ]
test_data <- x[-train_idx, ]

data_bn_us_train <- train_data %>% filter(region == 5) %>% select(CS_K:SCSY_o)
data_bn_us_test <- test_data %>% filter(region == 5) %>% select(CS_K:SCSY_o)
data_bn_hk_train <- train_data %>% filter(region == 3) %>% select(CS_K:SCSY_o)
data_bn_hk_test <- test_data %>% filter(region == 3) %>% select(CS_K:SCSY_o)

# Fit networks to training data
us_bn_fit <- bn.fit(us_bn, data_bn_us_train)
hk_bn_fit <- bn.fit(hk_bn, data_bn_hk_train)

# Predict and compare log-likelihood on test data
loglik_us <- logLik(us_bn_fit, data_bn_us_test)
loglik_hk <- logLik(hk_bn_fit, data_bn_hk_test)

print(paste("Log-likelihood (us_bn):", loglik_us))
print(paste("Log-likelihood (hk_bn):", loglik_hk))

# True Positives, False Positives, and False Negatives
arc_comparison <- compare_arcs(us_bn, hk_bn)
cat("Arc Comparison Results for BNs:\n")
print(arc_comparison)
arc_comparison <- compare_arcs(simulated_bn_us, simulated_bn_hk)
cat("Arc Comparison Results for BNs with simulated data:\n")
print(arc_comparison)
arc_comparison <- compare_arcs(boot_avg_us_bn,boot_avg_hk_bn)
cat("Arc Comparison Results for averaged bootstrapped BNs:\n")
print(arc_comparison)

# Equivalence Class Analysis (CPDAG)
cpdag_comparison <- compare_cpdag(us_bn, hk_bn)
cat("Are the equivalence classes (CPDAGs) of the two BNs identical?:", cpdag_comparison, "\n")
cpdag_comparison <- compare_cpdag(boot_avg_us_bn, boot_avg_hk_bn)
cat("Are the equivalence classes (CPDAGs) of the two averaged bootstrapped BNs identical?:", cpdag_comparison, "\n")

# Parameter-Level Comparisons
# Example: Compare parameters of a specific node
node_to_compare <- "PSA_total"  # Replace with another node name
parameter_comparison <- compare_parameters(us_bn_fit, hk_bn_fit, node_to_compare) 
cat("Parameter Comparison for Node:", node_to_compare, "\n") # The parameters here show the upstream connections (nodes leading to node_to_compare)
cat("Respective parameters for us_bn_fit and hk_bn_fit respectively:", "\n")
print(parameter_comparison)
#These retrieve the numbers that gie the strength in the side by side graphs. A higher value indicates the value in question has a higher causal strength to the one being evaluated

## Average Bootstrapped Networks' Density
# Convert bnlearn networks to igraph objects
us_igraph <- as.igraph(us_avg_network)
hk_igraph <- as.igraph(hk_avg_network)

# Calculate graph density using edge_density
us_density <- edge_density(us_igraph)
hk_density <- edge_density(hk_igraph)

# Print the densities
cat("Graph Density for US Averaged Network:", us_density, "\n")
cat("Graph Density for HK Averaged Network:", hk_density, "\n")