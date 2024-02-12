# Define bnlearn functions  

###################################################
## FUNCTION 1: fit.the.model

# The Bayesian network model is constructed using the `fit.the.model()` function, 
# which employs a constrained learning algorithm based on conditional independence
# testing, specifically the semi-interleaved HITON-PC method. In this process, the 
# algorithm identifies "parent-child" relationships within the network, where nodes
# represent either phenotypic traits or genetic variants, and edges signify associations
# between these nodes. Notably, traits can have parent nodes that are either genetic 
# variants or other traits, but they can only serve as parents to other traits, adhering
# to the constraint that genetic variants do not act on traits (since genotypes are constant 
# across an individual's lifespan).

# To enforce this constraint, a "blacklist" is created using the `tiers2blacklist()`
# function, preventing arcs from being directed towards specific nodes. This restriction
# aims to guide the learning process by ensuring that known causal relationships are 
# inferred in the correct direction from genetic variants to traits, while also allowing
# for customization to blacklist other traits if needed (e.g., you can customize this
# to force/restrict a trait->trait relationship).

# After the networks are learned, the nodes are categorized into subsets for 
# visualization, and the network structures are determined by maximizing the Bayesian
# Information Criteria (BIC). This approach facilitates the construction of Bayesian 
# networks that capture probabilistic relationships between traits and genetic variants, 
# with the learned structures reflecting potential causal associations. 

# Further details on the methodology can be found at
# http://www.bnlearn.com/research/genetics14/.

# data: The data set to be used
# alpha: The significance level for conditional independence tests 

fit.the.model <- function(data, alpha) { # Start fit.the.model
  
  # Creates an empty list (cpc) with a length equal to the length of the traits vector
  cpc <- vector(length(traits), mode = "list")
  # Assign names to the elements of the cpc list based on the values in the traits vector
  names(cpc) <- traits

  # Identify the parents of each trait (may be variants or other traits)
  for (t in seq_along(traits)) { # Start for 
   cpc[[t]] = learn.nbr(data[, c(traits, genes)], node = traits[t], debug = FALSE,
                   method = "si.hiton.pc", test = "cor", alpha = alpha)
  } # End for 

  # Merge the relevant variables to use for learning
  nodes <- unique(c(traits, unlist(cpc)))
  
  # Set up blacklist which specifies that variants cannot depend on traits
  blacklist <- tiers2blacklist(list(nodes[!(nodes %in% traits)], traits))
   
  # Build the Bayesian network
  bn <- hc(data[, nodes], blacklist = blacklist)

  # Return Bayesian network 
  return(bn)

} # End fit.the.model

###################################################

###################################################

## FUNCTION 2: xval.the.model

# The `xval.the.model()` function performs model training with n-fold 
# cross-validation (in the case of this example, 5-fold cross-validation).
 
# During this process, the data set is divided into multiple partitions, 
# with each partition serving as the test set while the remaining data is used
# for training. This process is repeated iteratively to ensure that all 
# data points are included in the test set at least once.

# During each fold of cross-validation, the following steps are performed:

# (1) Data Splitting: The dataset is divided into a training set (dtraining) 
# and a test set (dtest).
# (2) Model Fitting: A Bayesian network model is fitted to the training data 
# using the `fit.the.model()` function. This model captures the probabilistic 
# relationships between phenotypic traits and genetic variants.
# (3) Prediction: The model is used to predict the values of phenotypic traits
# on the test set. These predictions are stored in the prediction matrix.
# (4) Posterior Estimation: Posterior estimates are computed for each trait based
# on the test data. These estimates are stored in the post matrix.
# (5) Correlation Computation: The correlations between the predicted and observed 
# values for each trait are calculated, both before (predcor) and after (postcor) 
# cross-validation. These correlations provide a measure of the model's predictive performance.

# Finally, the function returns various results, including the predicted values, 
# posterior estimates, observed values, and correlation coefficients for each trait.
# Additionally, it provides the learned models for each fold of cross-validation. 
# This process allows for the assessment of how well the Bayesian network model 
# generalizes to unseen data and provides insights into its predictive capabilities.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  } # End redraw.label.ggnet
                                                                                       
# data: The data set to be used 
# k_crossval: The number of folds for cross-validation 
# cluster: The number of clusters to be used for parallel processing
# alpha: The significance level for conditional independence tests 
# ridge: A logical value indicating whether to use a custom threshold for network visualization

# Create a function to train the model with n fold cross validation, 
# 100-100/n% train and 100/n% test 
# Hold out each partition iteratively so that no data are totally held out
xval.the.model = function(data, k_crossval = 10, cluster, alpha, ridge) { # Start xval.the.model
  
  # Note: In cases of small sample sizes, k may need to be lowered
  
  # Store number of observations 
  n <- nrow(data)
  
  # Create a numeric vector named predcor initialized with a length
  # equal to the number of elements in the traits vector 
  # This vector will be used to store correlation values between
  # predicted and observed values for each trait
  predcor <- numeric(length(traits))
  names(predcor) <- traits
  
  # Repeat, creating a postcor vector that will be used to store 
  # correlation values for posterior estimates
  postcor <- numeric(length(traits))
  names(postcor) <- traits
  
  # Shuffle the data to get unbiased splits
  kcv <- split(sample(n), seq_len(k_crossval))
  # Store the length of each test set
  kcv.length <- sapply(kcv, length)
  
  # Train network models and support the optional use of ridge regression for model re-fitting under certain conditions
  predicted <- parLapply(kcv, cl = cluster, function(test) { # Start predicted
    
    # Create a matrix to store the predicted values
    pred <- matrix(0, nrow = length(test), ncol = length(traits))
    colnames(pred) <- traits
    # Create a matrix to store posterior estimates
    post = matrix(0, nrow = length(test), ncol = length(traits))
    colnames(post) = traits
    
    # Print a message for each cross-validation fold
    cat("* beginning cross-validation fold.\n")
    
    # Split training and test data
    dtraining <- data[-test, ]
    dtest <- data[test, ]
    
    # Fit the model on the training data (using function created in chunk above)
    model <- fit.the.model(dtraining, alpha = alpha)
    fitted <- bn.fit(model, dtraining[, nodes(model)])
    # Optional re-fit with ridge regression
    if (ridge) { # Start if
      
      # Load penalized library
      library(penalized)
      
      # Run ridge 
      for (no in nodes(fitted)) { # Start for
        
        node.parents <- parents(fitted, no)
        
        if (length(node.parents) < 3)
          next
        
        opt.lambda <- optL2(response = dtraining[, no],
                            penalized = dtraining[, node.parents],
                            model = "linear", trace = FALSE,
                            minlambda2 = 10e-5, maxlambda = 500)$lambda
        fitted[[no]] <- penalized(response = dtraining[, no],
                                  penalized = dtraining[, node.parents],
                                  model = "linear", trace = FALSE,
                                  lambda1 = 0, lambda2 = opt.lambda)
        
      } # End for
      
    } # End if -> then 
    
    # Subset the test data
    dtest <- dtest[, nodes(model)]
    
    # Print message reporting the number of nodes in the Bayesian network model
    cat("  > model has", length(nodes(model)), "nodes.\n")
    
    # Predict each trait in turn, given all the parents
    # For each trait, use the fitted Bayesian network model (fitted) to make 
    # predictions (predict) based on the test data (dtest) for the nodes
    # in the model
    # The predictions are stored in the pred matrix, where each column corresponds to a trait
    for (t in traits)
      pred[, t] <- predict(fitted, node = t, data = dtest[, nodes(model)])
    
    # Iterate over the rows of the test data (dtest)
    # For each row, calculate the posterior estimates for all traits (traits) 
    # in the Bayesian network
    # Note that this uses the cpdist function to estimate the conditional
    # probability distribution of each trait given the evidence (values 
    # of other variables) from the test data
    # The resulting posterior estimates are stored in the post matrix, 
    # with each column representing a trait
    for (i in seq(nrow(dtest)))
      post[i, traits] = colMeans(cpdist(fitted, nodes = traits,
                                        evidence = as.list(dtest[i, names(dtest) %in% genes]),
                                        method = "lw", n = 1000))
    
    # Return a list with the pred and post results
    return(list(model = fitted, pred = pred, post = post))
    
  }) # End predicted
  
  
  # Summarize results from cross-validation and calculate and print statistics 
  # on the model's performance 
  posterior <- do.call(rbind, lapply(predicted, `[[`, "post"))
  causal <- do.call(rbind, lapply(predicted, `[[`, "pred"))
  cat("* overall cross-validated correlations:\n")
  
  # Iterate over each trait in the traits vector to calculate 
  # and print performance statistics for each trait separately
  for (t in traits) {
    predcor[t] = cor(causal[, t], data[unlist(kcv), t])
    cat("  > PREDCOR(", t, "):", predcor[t], "\n")
    postcor[t] = cor(posterior[, t], data[unlist(kcv), t])
    cat("  > POSTCOR(", t, "):", postcor[t], "\n")
  } # End for
  
  # Construct and return a list that summarizes various aspects of the 
  # cross-validation results and model performance
  return(list(predicted = causal, posterior = posterior,
              observed = data[unlist(kcv), t], predcor = predcor, postcor = postcor,
              models = lapply(predicted, `[[`, "model")))
  
} # End xval.the.model              

###################################################

###################################################

## FUNCTION 3: run_plot_graph

# This code defines an R function called `run_plot_graph` that, calling the functions 
# above, learns and visualizes the Bayesian networks in the data, including customizing 
# node labels and appearance, and then returns results related to the network analysis. 
# Note that the following function uses parallel computing which involves breaking down
# a complex task into smaller subtasks that can be executed simultaneously, 
# or in parallel, by multiple processors or computers to speed up computationally
# intense calculations.

# data: The dataset to be used for plotting the networks
# k_crossval: The number of folds for cross-validation
# k_iterations: The number of times the entire process (which includes k-fold cross-validation) is repeated.
# alpha: The significance level for conditional independence tests
# use.custom.threshold: A logical value indicating whether to use a custom threshold for network visualization
# custom.threshold: The custom threshold value to be used if use.custom.threshold is set to TRUE
# ncluster: The number of clusters for parallel computation
# ...: Additional arguments to be passed to other functions

# Define function to plot the networks
run_plot_graph <- function(data = df,
                           k_crossval = 10,
                           k_iterations = 5,
                           alpha = 0.01,
                           use.custom.threshold = FALSE,
                           custom.threshold = 0.90,
                           ncluster = 8,
                           custom_labels = custom_labels,
                           ...) { # Start run_plot_graph
  
  #### CLUSTER
  # Create a cluster (cl) for parallel computation 
  cl <- makeCluster(ncluster)
  
  # Load library(bnlearn) and library(lme4) in the cluster 
  invisible(clusterEvalQ(cl, library(bnlearn)))
  invisible(clusterEvalQ(cl, library(lme4)))
  
  # Export the traits, genes, and fit.the.model objects to the cluster 
  # to make them available in the parallel environment
  clusterExport(cl = cl, c("traits", "genes", "fit.the.model"))
  
  # Initialize the pr001 variable as a vector of length k with the mode set to "list"
  pr001 <- vector(k_iterations, mode = "list")
  
  # Perform cross-validation (xval.the.model) k times and store the results in the pr001 list
  for (i in seq_along(pr001)) { # Start for 
    pr001[[i]] <- xval.the.model(data, 
                                 cluster = cl, 
                                 alpha = alpha, 
                                 ridge = FALSE,
                                 k_crossval = k_crossval) 
  } # End for 
  
  # Stop the cluster
  stopCluster(cl)
  #### CLUSTER
  
  # Calculate the average prediction correlation (pred.summary) and average 
  # posterior correlation (post.summary) based on cross-validation results
  pred.summary <- sapply(pr001, `[[`, "predcor")
  print(rowMeans(pred.summary))
  post.summary <- sapply(pr001, `[[`, "postcor")
  print(rowMeans(post.summary))
  
  # Average the network structures
  # Initialize a list (arclist) to store the network structures
  arclist <- list()
  # Extract the models from pr001 and append their arcs to the arclist
  for (i in seq_along(pr001)) { # Start for
    # Extract the models
    run <- pr001[[i]]$models
    for (j in seq_along(run))
      arclist[[length(arclist) + 1]] = arcs(run[[j]])
  } # End for
  # Compute the arc strengths
  # The nodes (genes and traits) present in the network structures are 
  # identified and stored in the nodes variable
  nodes <- unique(unlist(arclist))
  
  # Blacklist
  # Note: default boot.strength() or custom.strength() runs with 
  # "cpdag = TRUE" which means reversible arcs can have positive 
  # strength in both directions (so blacklist won't work as we expect)
  # Use "cpdag = FALSE" to truly make the blacklist work as expected 
  strength <- custom.strength(arclist, nodes = nodes, cpdag = FALSE)
  # Estimate the threshold and average the networks
  averaged <- averaged.network(strength)
  
  # Select relevant nodes (those with a degree greater than 0) 
  # and store in relevant.nodes
  relevant.nodes <- nodes(averaged)[sapply(nodes, degree, object = averaged) > 0]
  # Create a subgraph (averaged2) based on the relevant nodes
  averaged2 <- subgraph(averaged, relevant.nodes)
  # Extract the arc strengths (strength2) between relevant nodes 
  strength2 <- strength[(strength$from %in% relevant.nodes) & (strength$to %in% relevant.nodes), ]
  
  # Print information about the threshold value and minimum strength greater than the threshold 
  cat("threshold: ", attr(strength2, "threshold"), "\n")
  
  # Customize node labels 
  # Store the current node names from averaged
  current_nodes <- nodes(averaged)
  # Create a custom labels dictionary with all nodes from averaged2
  c.labels <- custom_labels2[current_nodes]
  # Relabel nodes in averaged2 using the custom labels
  nodes(averaged2) <- c.labels[nodes(averaged2)] 
  
  # Update the 'from' and 'to' columns in the strength object with the custom labels
  strength2$from <- c.labels[strength2$from]
  strength2$to <- c.labels[strength2$to]
  
  # Update the attribute names in the strength and strength2 objects
  attr(strength2, "nodes") <- c.labels[attr(strength2, "nodes")] 
  attr(strength, "nodes") <- c.labels[attr(strength, "nodes")]  # Newly added to correct bug, 12-13-2023
  
  # Need to also updated results$averaged and results$traits
  nodes(averaged) <- c.labels[nodes(averaged)]
  traits <- c.labels[traits]
  
  # strength2 retains nodes that are not in relevant.nodes as attributes 
  # Since the plotting code below uses the attributes of strength2 and averaged2,
  # we need to update the attributes to reflect only relevant nodes
  # Extract unique nodes from both "from" and "to" columns
  nodes_from_strength <- unique(c(strength$from, strength$to))
  nodes_from_strength2 <- unique(c(strength2$from, strength2$to))
  # Create a list of nodes to keep (nodes that appear in either data frame)
  nodes_to_keep <- intersect(nodes_from_strength, nodes_from_strength2)
  # Update the attributes of both data frames to only include nodes_to_keep
  attr(strength2, "nodes") <- attr(strength2, "nodes")[nodes_to_keep]
  
  # Print the minimum strength observed that is > threshold
  t <- attr(strength2, "threshold")
  v <- strength2$strength
  cat("min strength > threshold: ", min(v[v > t]), "\n")
  cat("strength: ", sort(strength2$strength))
  
  # Plot the network 
  # The appearance of nodes and edges is customized based on the data
  if (use.custom.threshold) {
    cat("Using custom threshold of ", custom.threshold, "\n")
    gR = strength.plot(
      averaged2, 
      strength2, 
      shape = "rectangle",
      layout = "fdp",
      render = F,
      threshold = custom.threshold
    )
  } else {
    gR = strength.plot(
      averaged2,
      strength2,
      shape = "rectangle",
      layout = "fdp",
      render = F
    )
  }
  
  # Customize plot
  nodeRenderInfo(gR)$fill = "lightblue"
  nodeRenderInfo(gR)$col = "darkblue"
  nodeRenderInfo(gR)$fill[traits] = "limegreen"
  nodeRenderInfo(gR)$col[traits] = "darkgreen"
  a <- arcs(bnlearn::subgraph(averaged, traits))
  a <- as.character(interaction(a[, "from"], a[, "to"], sep = "~"))
  edgeRenderInfo(gR)$col = "grey"
  edgeRenderInfo(gR)$col[a] = "darkgreen"
  
  # Render plot (gR) using Rgraphviz::renderGraph
  renderGraph(gR)
  
  # Create a list named results that contains the averaged2 network, 
  # strength2, averaged network, traits variables, and threshold value
  results <-
    list(
      averaged2 = averaged2,
      strength2 = strength2,
      averaged = averaged,
      traits = traits,
      threshold = t
    )
  
  # Return results
  return(results)
  
} # End run_plot_graph

###################################################

###################################################

## FUNCTION 4: match.arcs.and.directions

# This code defines an R function called `match.arcs.and.directions` which is 
# designed to match arcs in a network with their corresponding direction 
# coefficients in a strengths data frame, providing flexibility to either 
# keep or discard additional information based on the keep argument.


# arcs: Arcs within a Bayesian network representing relationships 
# between nodes (genes or traits)
# nodes: Nodes within a Bayesian network, which can represent 
# either genes or traits
# strengths: Strength coefficients associated with the arcs in 
# the Bayesian network, indicating the strength and direction
# of relationships
# keep: A logical value indicating whether to keep additional
# information when matching arcs and directions in the 
# Bayesian network

# Define function to match arcs and directions 
match.arcs.and.directions <- function (arcs, 
                                       nodes,
                                       strengths, 
                                       keep = FALSE) { # Start match.arcs.and.directions
  if (nrow(strengths) < nrow(arcs))
    stop("insufficient number of strength coefficients.")
  a_hash <- interaction(arcs[, "from"], arcs[, "to"])
  s_hash <- interaction(strengths[, "from"], strengths[, "to"])
  if (keep) {
    s <- strengths[match(a_hash, s_hash), , drop = FALSE]
    coef = s$direction
  } else {
    s <- strengths[match(a_hash, s_hash), "direction"]
    from <- strengths[match(a_hash, s_hash), "from"]
    to <- strengths[match(a_hash, s_hash), "to"]
    names(s) <- paste0(from, "~", to)
    coef <- s
  }
  if (any(is.na(coef))) {
    missing = apply(arcs[is.na(coef), , drop = FALSE], 1,
                    function(x) {
                      paste(" (", x[1], ", ", x[2], ")", sep = "")
                    })
    stop(
      "the following arcs do not have a corresponding direction coefficients:",
      missing,
      "."
    )
  }
  
  # Return 
  return(s)
} # End match.arcs.and.directions

###################################################

###################################################

## FUNCTION 5: redraw.graph.labels

# This code defines the `redraw.graph.labels` function, which enhances the graph 
# created above by adding labels to its edges. These labels make the graph more
# informative by conveying information about the strength and direction of
# relationships represented by each edge in the graph. 

# Strength: a measure of confidence of that edge while fixing the rest of
# the network structure and is defined as the empirical frequency a specific 
# edge is observed over a set of networks learned from iterations 
# (i.e., the number of times the edge was present out of the total number of 
# iterations). 

# Direction: the probability of the edge’s direction conditional on the edge’s
# presence within the network (i.e., the number of times the edge traveled in 
# a specific direction out of the total number of iterations
# in which it was present). 


# averaged2: Bayesian network learned through functions above (relevant nodes only)
# strength2: Bayesian network strengths learned through the functions above
# averaged: Bayesian network learned through functions above (all nodes)
# traits: Phenotypes of interest
# custom.threshold: Custom threshold to customize line thickness/type

# Define function to redraw graph labels 
redraw.graph.labels <-
  function(averaged2,
           strength2,
           averaged,
           traits,
           custom.threshold) { # Start redraw.graph.labels
    gR <- strength.plot(
      averaged2,
      strength2,
      shape = "rectangle",
      layout = "fdp",
      render = F,
      threshold = custom.threshold
    )
    x <- averaged2
    str <-
      match.arcs.and.directions(
        arcs = x$arcs,
        nodes = names(x$nodes),
        strengths = strength2
      )
    str2 <- bnlearn:::match.arcs.and.strengths(
      arcs = x$arcs,
      nodes = names(x$nodes),
      strengths = strength2
    )
    
    # Customize edge labels 
    labels2 <- paste0(round(str2, 2), ":", round(str, 2))
    names(labels2) <- names(str)
    
    # Begin plot
    gR <- Rgraphviz::layoutGraph(gR,
                                 edgeAttrs = list(label = labels2),
                                 layoutType = "fdp")
    
    # Customize plot
    nodeRenderInfo(gR)$fill[genes] = "grey"
    nodeRenderInfo(gR)$fill = "lightblue"
    nodeRenderInfo(gR)$col = "darkblue"
    nodeRenderInfo(gR)$fill[traits] = "limegreen"
    nodeRenderInfo(gR)$col[traits] = "darkgreen"
    a = arcs(bnlearn::subgraph(averaged2, traits))
    a = as.character(interaction(a[, "from"], a[, "to"], sep = "~"))
    edgeRenderInfo(gR)$col = "grey"
    edgeRenderInfo(gR)$col[a] = "darkgreen"
    Rgraphviz::renderGraph(gR)
  } # End redraw.graph.labels

###################################################

###################################################

## FUNCTION 6: redraw.label.ggnet

# The following function, `redraw.label.ggnet`, is used for an additional 
# redrawing of the Bayesian network but in a ggnetwork/ggplot2 framework, 
# enhancing it with labels and color-coding for interpretation in the context
# of Bayesian analysis. Of note, this function calls in `ggrepel()` so 
# that node and edge labels do not overlap, further supporting the user
# in creating publication quality figures. 


# averaged2: Bayesian network learned through functions above (relevant nodes only)
# strength2: Bayesian network strengths learned through the functions above
# averaged: Bayesian network learned through functions above (all nodes)
# traits: Phenotype of interest
# ...: Additional arguments to be passed to other functions

# Define function to redraw the graph in ggnetwork/ggplot2 framework
redraw.label.ggnet <- function(averaged2, 
                               strength2,
                               averaged,
                               traits,
                               df_vertex_table=df_vertex_table,
                               ...) { # Start redraw.label.ggnet
  
  
  x <- averaged2
  
  # Relies on bnlearn code to return a directed graph
  # Get the directions and strengths
  vec_dir <- match.arcs.and.directions(
    arcs = x$arcs,
    nodes = names(x$nodes),
    strengths = strength2
  )
  vec_str <- bnlearn:::match.arcs.and.strengths(
    arcs = x$arcs,
    nodes = names(x$nodes),
    strengths = strength2
  )
  
  # Create data.frame with node/directions/strengths
  # Join this to the data.frame created by ggnetwork
  df_dir_str <- data.frame(do.call(rbind, str_split(names(vec_dir), "~")),
                           vec_dir,
                           vec_str,
                           stringsAsFactors = F)
  names(df_dir_str)
  
  # Convert NEL graph to iGraph to ggnetwork data.frame
  # This step 'shuffles' the edge orders so vec_dir and vec_str are not in same order
  # as df_ggnetwork
  graph_igraph <- igraph::igraph.from.graphNEL(as.graphNEL(averaged2))
  df_ggnetwork <- ggnetwork(graph_igraph, ...)
  df_ggnetwork$vertex.names <- as.character(df_ggnetwork$name)
  df_ggnetwork$xend <- as.numeric(df_ggnetwork$xend)
  df_ggnetwork$yend <- as.numeric(df_ggnetwork$yend)
  df_ggnetwork$x <- as.numeric(df_ggnetwork$x)
  df_ggnetwork$y <- as.numeric(df_ggnetwork$y)
  
  # Create a data.frame of node coordinates as this is needed to determine the
  # destination node/end vertex
  df_node_coord <- df_ggnetwork %>% 
    filter(is.na(weight)) %>% 
    select(vertex.names, x, y) %>%
    transmute(vertex.names.end = vertex.names, xend = x, yend = y)
  
  # Join back to form end vertex
  # This step relies on fuzzyjoin with L2 norm < 0.05 since ggnetwork shifts the node 
  # positions around randomly
  df_ggnetwork <- df_ggnetwork %>%
    fuzzyjoin::distance_left_join(df_node_coord,
                                  by = c("xend", "yend"),
                                  max_dist = 0.05)
  
  # Join the direction/probability table based on source node and end node
  df_ggnetwork <- df_ggnetwork %>% left_join(df_dir_str,
                                             by = c("vertex.names" = "X1", "vertex.names.end" = "X2"))
  
  # Join node type table
  df_ggnetwork <- df_ggnetwork %>% left_join(df_vertex_table, by = "vertex.names")
  
  # Create a column for node colors based on the node type
  # Note that in this example, we only have two node types (symptom/SNP)
  df_ggnetwork <- df_ggnetwork %>%
    mutate(node_color = ifelse(grepl("^rs", vertex.names), "tan3", "steelblue"))
  
  # Fix the data.frame column names, fuzzyjoin appended .x and .y to column names
  names(df_ggnetwork)[names(df_ggnetwork) == "xend.x"] <- "xend"
  names(df_ggnetwork)[names(df_ggnetwork) == "yend.x"] <- "yend"
  
  # Plot with ggplot and ggnetwork geoms
  ggplot(df_ggnetwork, aes(
    x = x,
    y = y,
    xend = xend,
    yend = yend)) +
    geom_edges(aes(linetype = (vec_str > 0.9)),
               color = "grey50",
               arrow = arrow(length = unit(8, "pt"), type = "open")) +
    geom_nodes(size = 10,
               aes(shape = type, color = node_color, fill = node_color),
               stroke = 0.5,   # Set the size of the border
               color = "black") +  # Set the color of the border
    geom_edgelabel_repel(aes(label = paste0(
      round(vec_str, 2), ":", round(vec_dir, 2)))) +
    geom_nodelabel_repel(
      aes(label = vertex.names),
      force = 1,
      fontface = "bold",
      box.padding = unit(1.2, "lines"),
      color = "black") +
    scale_color_manual(values = c("tan3" = "tan3", "steelblue" = "steelblue")) +
    scale_fill_manual(values = c("tan3" = "tan3", "steelblue" = "steelblue")) +
    scale_linetype_manual(
      "Node Signif.",
      values = c("FALSE" = "dotted", "TRUE" = "solid"),
      label = c("<= 90 Strength", "> 0.90 Strength")) +
    scale_shape_manual(values = c(21:23)) + theme_blank(legend.position = "none")
  
} # End redraw.label.ggnet

###################################################

