find_nn_in_chunks <- function(query.emb, ref.emb, k = 10, chunk_size = 10000, n_cores = 10) {
  # Split query embeddings into chunks
  chunks = split(1:nrow(query.emb), 
                 ceiling(seq_along(1:nrow(query.emb))/chunk_size))
  
  # Process chunks in parallel
  nn.chunks = parallel::mclapply(chunks, function(idx) {
    query.emb.chunk = query.emb[idx, ]
    RANN::nn2(ref.emb, query.emb.chunk, k = k)
  }, mc.cores = n_cores)
  
  # Combine results
  list(
    nn.idx = do.call(rbind, lapply(nn.chunks, function(x) x$nn.idx)),
    nn.dists = do.call(rbind, lapply(nn.chunks, function(x) x$nn.dists))
  )
}

# Function to get predicted labels from nearest neighbors
get_pred_labels <- function(nn_idx, ref_labels) {
  query_labels = matrix(nrow = nrow(nn_idx), ncol = ncol(nn_idx))
  
  # Get labels for k nearest neighbors
  for (i in 1:nrow(nn_idx)) {
    neighbor_labels = ref_labels[nn_idx[i,]]
    query_labels[i,] = neighbor_labels
  }
  
  # Get most frequent label and probability
  pred.probs = apply(query_labels, 1, function(x) {
    sort(table(x), decreasing = TRUE)[1]
  })
  pred.labels = apply(query_labels, 1, function(x) {
    names(sort(table(x), decreasing = TRUE))[1]
  })
  
  list(
    labels = pred.labels,
    probs = pred.probs
  )
}


transfer_labels <- function(ref, ref.reduction, query, query.reduction, dims, ref.label.names, k = 10, n_cores = 10) {
  # Get embeddings
  ref.emb = Embeddings(ref, reduction = ref.reduction)[, dims]
  query.emb = Embeddings(query, reduction = query.reduction)[, dims]
  
  # Find nearest neighbors
  nn = find_nn_in_chunks(query.emb, ref.emb, k = k, n_cores = n_cores)
  
  # Get reference labels
  results = lapply(ref.label.names, function(label.name) {
    ref_labels = FetchData(ref, vars = c(label.name), cells = rownames(ref.emb))[[1]]
    ref_labels = as.character(ref_labels)
    
    # Get predicted labels
    pred = get_pred_labels(nn$nn.idx, ref_labels)
    
    # Prepare results
    results = data.frame(
      row.names = rownames(query.emb),
      pred_label = pred$labels,
      pred_prob = pred$probs / k
    )
    colnames(results) = paste0(label.name, c("_pred", "_prob"))
    return(results)
  }) %>% do.call(cbind, .)
  
  return(results)
}
