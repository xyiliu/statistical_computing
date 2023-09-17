#' @name DistributionBasedFeatureExtractor
#' @title Extract Features from Distribution of Lengths
#' @description
#' 4 feature extraction methods that facilitate classifier learning
#' from genomic variants are included in this function.
#' @param values A dataframe of lengths to be discretized with 2 columns, sample_id and lengths.
#' @param labels A dataframe of labels with 2 columns, sample_id and labels. Only needed for supervised scenarios.
#' @param breakpoint_type  A string of type of breakpoints to generate.
#' Available options are: equal, quantile, supervised, clustering. Defaults to 'quantile'.
#' @param n_bins An integer or "auto" of number of features to generate.
#' An integer explicitly selects given number of features. "auto" is only for supervised.
#' @param prefix A string of the name
#' @param bw A float number of bandwidth for kde. Used only in supervised breakpoints. Defaults to 0.5.
#' @param log_scale A bool whether log scale should be used on the x axis to calculate the bins.
#' Used only in supervised breakpoints. Defaults to True.
#' @param resolution An integer of resolution of kde. Used only in supervised breakpoints. Defaults to 200.
#' @param value_range Range of values. NULL is allowed.
#' @return  A dataframe with all the features, combined with labels for supervised scenarios.
#' @import mclust
#' @import dplyr
#' @importFrom Rcpp sourceCpp
#' @examples
#' \dontrun{
#' DistributionBasedFeatureExtractor(values=df, labels=y)
#' }
#' @exportPattern "^[[:alpha:]]+"
#'
#'
#' @export
DistributionBasedFeatureExtractor <- function(values,
                                              labels,
                                              breakpoint_type="quantile",
                                              n_bins=4,
                                              prefix="dbfe",
                                              bw=0.5,
                                              resolution=200,
                                              log_scale=TRUE,
                                              value_range=NULL)
{


  #  if (length(bins) != 0){
  if ((n_bins == 1) && (breakpoint_type != "supervised")){
    bins = c(0, Inf)
  } else {
    check_values(values)
    if (!is.null(value_range)){
      values = filter_values_to_range(values, value_range)
    }
    if (breakpoint_type == "equal"){
      bins = generate_equal_bins(values, n_bins, log_scale)
    }
    else if (breakpoint_type == "quantile"){
      bins = generate_quantile_bins(values, n_bins)
    }
    else if (breakpoint_type == "cluster"){
      bins = generate_clustering_bins(values, n_bins)
    }
    else if (breakpoint_type == "supervised"){
      res = generate_supervised_bins(values, labels, max_bins="auto", bw=bw, resolution=resolution, log_scale=T)
      bins = res$bins
    }
    else{
      stop("'Unsupported transformation type! Only 'quantile', 'clustering', and 'supervised' ")
    }
  }
  #  }
  result_tab = create_sample_features_from_bins(values, bins, prefix="dbfe")
  out_df = merge(result_tab, labels, by=0, all.x=TRUE)
  out_df = within(out_df, rm("sample_id"))
  return(out_df)
}







generate_equal_bins <- function(lengths, n_bins=4, log_scale=TRUE){

  lengths_min = min(lengths$var_len)
  lengths_max = max(lengths$var_len)
  if(log_scale==TRUE){
    breakpoints = logspace(log10(lengths_min), log10(lengths_max), n = n_bins + 1)[2:n_bins]
  } else{
    breakpoints = linspace(lengths_min, lengths_max, n = n_bins + 1)[2:n_bins]
  }
  result = breakpoints_to_bins(breakpoints)

  return(result)
}

generate_quantile_bins <- function(lengths, n_bins = 4){
  quantiles = quantile(lengths$var_len, probs=seq(0, 1, 1/n_bins))
  # the 0% quantile is the smallest value, and the 100% quantile is the largest value
  last_quantile = length(quantiles)

  breakpoints = round(quantiles[-c(1,last_quantile)])

  result = breakpoints_to_bins(breakpoints)

  return(result)
}

generate_clustering_bin <- function(lengths, no_cluster){

  library(mclust)
  model = Mclust(lengths$var_len, no_cluster, model = "V") ##Gaussian mixture fitted by EM

  lengths$Cluster = model$classification
  group_min <- lengths %>%
    group_by(Cluster) %>%
    summarise(Min = min(var_len)
    )
  breakpoints <- group_min$Min[-1]
  result = breakpoints_to_bins(breakpoints)
  return(list(result = result, model = model))
}

generate_supervised_bins <- function(lengths, classes, max_bins="auto", bw=0.5, resolution=200, log_scale=T, only_peaks=FALSE) {
  density = calculate_density(lengths, classes, bw, resolution, log_scale)
  d1_mean = aggregate(Density ~ Length, data = density[density$Class==1,], mean)
  d0_mean = aggregate(Density ~ Length, data = density[density$Class==0,], mean)

  density_diff = d1_mean - d0_mean
  density_diff[,"Length"] = d1_mean[,"Length"]

  breakpoints = calculate_breakpoints_and_peaks_from_density_diff(density_diff)$breakpoints
  peaks = calculate_breakpoints_and_peaks_from_density_diff(density_diff)$peaks
  bins = breakpoints_to_bins(breakpoints)
  n_bins = calculate_n_bins(max_bins, density_diff, peaks, bins)
  return(list(bins=bins, density=density))
}

calculate_n_bins = function(max_bins, density_diff, peaks, bins)
{
  if ((max_bins == 'auto') | (max_bins == 'sd'))
  {
    n_bins = sum((abs(peaks) >= sd(density_diff$Density))*1)
  }
  else if (max_bins == round(max_bins)) #check if the max_bins is an integer
  {
    n_bins = max_bins
  }
  else if (is.double(max_bins)) # if not an int, check if it's a double
  {
    n_bins = sum((abs(peaks) >= max_bins * max(abs(density_diff$Density)))*1)
  }
  else if (max_bins == "all")
  {
    n_bins == len(bins)
  }
  else
  {
    stop("Unsupported \"max_bins\" value, but this should have been detected earlier")
  }
  return(n_bins)
}


calculate_breakpoints_and_peaks_from_density_diff = function(density_diff)
{
  breakpoints = c()
  peaks = c()

  curr_peak = 0
  for (i in 1: (dim(density_diff)[1] - 1)) {
    if (abs(density_diff$Density[i]) > abs(curr_peak)) {
      curr_peak = density_diff$Density[i]
    }
    if ( density_diff$Density[i] * density_diff$Density[i+1] < 0 ) {
      breakpoints = c(breakpoints, (density_diff$Length[i]+density_diff$Length[i+1])/2)
      peaks = c(peaks, curr_peak)
      curr_peak = 0
    }
  }

  peaks = c(peaks, curr_peak)
  return (list(breakpoints = breakpoints, peaks=peaks))
}


calculate_density = function(lengths, classes, bw=0.5, resolution=200, log_scale=TRUE)
{

  density = data.frame()
  if (log_scale==TRUE){
    x_grid = seq(log(min(lengths[,2])), log(max(lengths[,2])), length.out=resolution)
  }
  else{
    x_grid = seq(min(lengths[,2]), max(lengths[,2]), length.out=resolution)
  }

  density = rbind(density, calculate_density_in_class(x_grid, lengths, classes, 1, bw, log_scale))
  density = rbind(density, calculate_density_in_class(x_grid, lengths, classes, 0, bw, log_scale))

  return (density)
}


calculate_density_in_class <- function(x, values, classes, c, bw, log_scale){
  sample = rownames(classes[which(classes$label==c),])
  values_c = values[values$sample_id %in% c(sample),"var_len"]

  if (log_scale==TRUE){
    values_c = log(values_c)
  }

  kde_skl <- density(values_c, bw=bw, kernel="gaussian")
  density_func = approxfun(kde_skl$x, kde_skl$y)
  density <- density_func(x)

  density_len <- length(density)
  if (log_scale==TRUE){
    result = data.frame("Length" = exp(x), "Density" = density, "Class" = rep(c, density_len))
  } else {
    result = data.frame("Length" = x, "Density" = density, "Class" = rep(c, density_len))
  }
  return (result)
}

create_sample_features_from_bins <- function(lengths, bins, prefix="dbfe"){
  num_bins = length(bins)-1
  # add a column to df indicating what the class are they in
  generated_bins_df = lengths %>%
    group_by(sample_id) %>%
    mutate(bin = cut(var_len, breaks = bins, labels = c(1:num_bins)))

  feature_names = c()
  for (i in 1:(length(bins)-1)){
    feature_name = paste(prefix, "_", bins[i], "_", bins[i+1])
    feature_names = c(feature_names, feature_name)
  }
  # group_by sample ids and count lengths in each bins
  each_bin_count = generated_bins_df %>% distinct() %>% count(bin)
  ## get the uniq sample names
  sample_id_uniq = unique(lengths$sample_id)
  values = c()
  for (j in sample_id_uniq){
    # each row is a sample, we first create a vector with all 0
    #represent the counts of variant lengths in each bins for jth sample
    each_sample_value = rep(0, num_bins)
    # get the jth sample's bin assignment
    each_sample_df = each_bin_count[each_bin_count$sample_id==j, ]
    each_sample_bin = each_sample_df$bin # get the bin assignment
    each_sample_counts = each_sample_df$n # get the counts
    ## will figure out how to prevent a nested loop later but this is all I could do
    ##at this point
    ## default is 4 bins
    for (i in 1:length(each_sample_bin)){
      # updating the bin counts if that bin's is contating at least one counts
      each_sample_value[as.integer(each_sample_bin[i])] = each_sample_counts[i]
    }
    values = rbind(values, each_sample_value)
  }
  #values = rbind(values, each_sample_value)
  #df <- data.frame("feature_names" = feature_names, "values"=values)
  rownames(values) = sample_id_uniq
  colnames(values) = feature_names
  return (values)
}


breakpoints_to_bins <- function(breakpoints){
  bins = c(0, breakpoints, Inf) # set the first number of the bins to be zero
  return(bins)
}


is_float_vec <- function(range_vec){
  return (is.numeric(range_vec) && !all(range_vec %% 1 == 0, na.rm = TRUE))
}

is_int_vec <- function(range_vec){
  return (mean(range_vec - floor(range_vec)) == 0)
}

filter_values_to_range <- function (values, value_range){
  if (value_range[1] >= value_range[2]){
    stop("Minimum of desired range must be smaller than maximum.")
  }

  if (!is.vector(value_range) | (length(value_range)!=2)){
    stop("Incorrect Range or type! Expect a vector of length 2")
  }

  if (is_float_vec(value_range)){
    min_quantile = quantile(values$var_len, value_range[1])
    max_quantile = quantile(values$var_len, value_range[2])
    values = values[values$var_len >= min_quantile, ]
    values = values[values$var_len <= max_quantile, ]
  }
  else if (is_int_vec(value_range)){
    values = values[values$var_len >= value_range[1], ]
    values = values[values$var_len <= value_range[2], ]
  }
  else{
    stop("Unexpected type of range. Expected a pair of ints or floats")}
  return (values)
}


check_values <- function(values){
  if (any(values$var_len) < 0){
    stop("'All values must be greater than 0!'")
  }
}


