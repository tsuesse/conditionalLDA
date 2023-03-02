
#***************************************************************************
# obtain_sample_SRS ----
#***************************************************************************
#' @description This function takes random samples of grouping blocks (blocks).
#' @param data dataset with observations, dataframe (maybe also matrix with conversion in function)
#' @param block name of a column in data in which identifiers of aggregation blocks are stored
#' @param ns number of blocks (aggregated) to be sampled
#' @param equal.n.var equal number of groups based on this variable
#' @export


obtain_sample_SRS <- function(data, block, ns, equal.n.var=NULL) {
  
  blocks_units<-as.integer(data[,block])
  duplicated_blocks<-duplicated(blocks_units)
  
  blocks<- unique(blocks_units)
  
  if(is.null(equal.n.var)){
    
    s_blocks <- sort(sample(blocks, ns))  
    s <- (1:a)[is.element(blocks_units, s_blocks)] # sample of blocks
    return(tibble::lst(s, s_blocks))
    
    }else{
  
  data_blocks<-data.frame(blocks=blocks,equal.n.var=data[,equal.n.var][!duplicated_blocks])
  
    a <- nrow(data)
 
    stratified <- data_blocks %>%
    group_by(equal.n.var) %>%
    sample_n(size=ns)
  
  s_blocks<-  sort(stratified$blocks)
  s <- (1:a)[is.element(blocks_units, s_blocks)] # sample of blocks
    
    # table(stratified[,"equal.n.var"])
  ind<-duplicated(data[s,block])
  
  return(tibble::lst(s, s_blocks, ns_per_group=table(data[s,equal.n.var][!ind])))
  
  } 
 
 
  
  
  # check
  # all.equal(unique(data[s,"field"]),s_blocks)
  # data[s,"croptype"]
  # sum(!duplicated_blocks)
  # length(blocks)
  # ind<-duplicated(data[s,"field"])
  # table(data[s,"croptype"][!ind])

  
} # end function


#***************************************************************************
# obtain_sample_cluster_adj ----
#***************************************************************************
#' @description This function takes random samples of grouping blocks (blocks) using a spatial cluster approach.
#' @param data dataset with observations, dataframe (maybe also matrix with conversion in function)
#' @param block name of a column in data in which identifiers of aggregation blocks are stored
#' @param ns number of blocks (aggregated) to be sampled
#' @param coords geopgraphical coordinates
#' @param distance maximimum distance
#' @param prop_clust proportion of sample that form initial clusters
#' @param equal.n.var equal number of groups based on this variable
#' @export

obtain_sample_cluster_adj <- function(data, block, ns, coords, distance = 1500,
                                      prop_clust = 0.1,equal.n.var=NULL, dist_mat_blocks=NULL) {
# obtain training data with cluster sampling, prop_clust is automatically increased by 0.01
# if there are to little blocks within given distance of cluster centers
prop_clust <- prop_clust - 0.01
hilf_s_cluster <- NULL
while (is.null(hilf_s_cluster$s)) {
  prop_clust <- prop_clust + 0.01
  hilf_s_cluster <- tryCatch(
    {
      obtain_sample_cluster(
        data = data, block = block, ns = ns,
        coords = coords, distance = distance,
        prop_clust = prop_clust,equal.n.var=equal.n.var,dist_mat_blocks=dist_mat_blocks)
    },
    error = function(cond) {
      message <- cat("Number of blocks close to random fields was less than required.\n")
      message <- cat("prop_clust will be increased by 0.01 to ", prop_clust + 0.01, "\n")
      return(NULL)
    }
  )
}

return(hilf_s_cluster)
}

#***************************************************************************
# obtain_sample_cluster ----
#***************************************************************************
# This function draws a cluster sample, i.e. a given proportion of the desired
# sample size is sampled randomly, and further observations are added based on
# there distance to the clustersamples. Blocks within a given distance are included.

#' @param data dataset with observations, dataframe (maybe also matrix with conversion in function)
#' @param aggr name of a column in data in which identifiers of aggregation blocks are stored
#' @param n number of blocks (aggregated) to be sampled
#' @param coords names of columns in data in which coords for each observation are stored
#' @param distance distance around central cluster blocks within which further blocks are incorporated into sample

obtain_sample_cluster <- function(data, block, ns, coords, distance = 1500,
                                  prop_clust = 0.1, equal.n.var=NULL,dist_mat_blocks=NULL) {
  
  
  blocks_units<-as.integer(data[,block])
  duplicated_blocks<-duplicated(blocks_units)
  
  blocks<- unique(blocks_units)
  
  # block to block distance
  #  dist_mat_blocks
  
  n_clust <- ceiling(ns * prop_clust) # number of clusters (central cluster blocks)
  n_within <- ns - n_clust # remaining number of blocks
  
  # again: is sorting necessary?
  
  if(is.null(equal.n.var)){
    s_blocks <- sort(sample(blocks, n_clust))  
  }else{
    
    data_blocks<-data.frame(blocks=blocks,equal.n.var=data[,equal.n.var][!duplicated_blocks])
    
    # n_units <- nrow(data)
    
    stratified <- data_blocks %>%
      group_by(equal.n.var) %>%
      sample_n(size=n_clust)
    
    s_blocks<-  sort(stratified$blocks)
    
    # table(stratified[,"equal.n.var"])
  } 
  
  
  
  DMsr<- dist_mat_blocks[s_blocks,-s_blocks]
  
  neighbors <- which(DMsr < distance, arr.ind = T)
  neighbors <- neighbors[neighbors[, 1] != neighbors[, 2], ]
  
  # obtain unsampled blocks
  blocks_not_sampled <- setdiff(blocks,s_blocks)
  # obtain blocks that are within certain distance
  blocks_not_sampled_close <- unique(blocks_not_sampled[neighbors[, 2]])
  
  # all.equal(sort(c(blocks_not_sampled,s_blocks)),blocks)
  
  # all.equal(data_blocks[,1],1:nrow(data_blocks))
  
  
  
  if (length(blocks_not_sampled_close) < n_within) {
    stop("Number of blocks close to random fields is less than required.\n", 
         "Increase prop_clust or n or distance!")
  }
  
  if(is.null(equal.n.var)){
    
    blocks_within <- try(sample(blocks_not_sampled_close, n_within))
  
  }else{
    
    data_blocks_not_sampled<- data_blocks[blocks_not_sampled_close,]
    
    # n_units <- nrow(data)
    
    stratified <- try(data_blocks_not_sampled %>%
      group_by(equal.n.var) %>%
      sample_n(size=n_within))
    
    blocks_within<-  try(sort(stratified$blocks))
    
    # table(stratified[,"equal.n.var"])
  } 
  
  
  
  if (!inherits(blocks_within, "try-error") && length(intersect(blocks_within, s_blocks)) == 0) {
    cat("Successful Cluster Sampling\n")
    
    s_blocks <- c(s_blocks, blocks_within)
    a <- nrow(data) 
    s <- (1:a)[is.element(blocks_units, s_blocks)]
    
    if(is.null(equal.n.var)){
      return(tibble::lst(s, s_blocks, n_clust, n_within, prop_clust))
    }else{
      
      ind<-duplicated(data[s,block])
      ns_per_group<-table(data[s,equal.n.var][!ind])
      
      return(tibble::lst(s, s_blocks, n_clust, n_within, prop_clust,
                         ns_per_group=ns_per_group)) 
    
    }
  
  
  }else{
    
    cat("Issues with sampling, not enough blocks within distance\n")
    return(tibble::lst(s=NULL))  
    
  }
} # end


#***************************************************************************
# lda_predmaj ----
#***************************************************************************
# This function predicts classes based on given lda fit and assigns predictions 
# to all pixels of a field depending on a majority vote within the field
#*********************************
#' @param object lda fit
#' @param newdata test data to be predicted
#' @param fac identifier for fields/blocks
#' @export

lda_predmaj <- function(object, newdata, fac = NULL) {
  
  
  
  majority <- function(x) {
    levels(x)[nnet::which.is.max(table(x))]
  }
  
  majority_filter <- function(x, fac) {
    for (lev in levels(fac)) {
      x[fac == lev] <- majority(x[fac == lev])
    }
    x
  }
  
  pred <- predict(object, newdata = newdata, type = "response")
  if (!is.null(fac)){ 
    pred <- majority_filter(pred$class, newdata[, fac]) 
  }
  
  return(pred)
}

#***************************************************************************
# calc_BIC ----
#***************************************************************************
# This function calculates BIC based on the loglikelihoof of the model, number of
# training observations and number of estimated parameters.
#*****************************
# function that calculates BIC
calc_BIC <- function(LLmodel, ns_units, number_pars) {
  BIC <- (-2 * LLmodel + log(ns_units) * number_pars)
}
