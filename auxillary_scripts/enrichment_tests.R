#' Test if elements of a query vector are enriched in a test vector using Fisher's exact test
#'
#' @param test A character vector. Should only contain unique elements. 
#' @param query A character vector. Should only contain unique elements. 
#' @param universe A vector containing all the elements that test and query could contain. Should only contain unique elements. 
#' @param alternative The alternative hypothesis for the test. Default is "greater".
#' @param return_overlap A logical value indicating whether to return the overlap between test and query as a character with elements separated by ";".
#' @return A data.frame with the results of Fisher's exact test
#' @export
fisher_test_vectors = function(test, query, universe, alternative = "greater", return_overlap = F){
  
  # Check if there are any duplicated elements in test, query or universe
  if(anyDuplicated(test)){stop("test cannot contain duplicated elements")}
  if(anyDuplicated(query)){stop("query cannot contain duplicated elements")}
  if(anyDuplicated(universe)){stop("universe cannot contain duplicated elements")}
  
  # Perform Fisher tests
  ft_result = fisher.test(
    x = factor(universe %in% test, levels = c("TRUE", "FALSE")), 
    y = factor(universe %in% query, levels = c("TRUE", "FALSE")), 
    alternative = alternative)
  
  # Calculate the size of the overlaps of the test and the universe with the query
  test_overlap_size = sum(test %in% query)
  universe_overlap_size = sum(universe %in% query)
  
  # Calcualte the proportion of the test vector which overlaps the query
  test_overlap_proportion = test_overlap_size/length(test)
  
  # Calculate the proportion of the universe which overlaps the query
  universe_overlap_proportion = universe_overlap_size/length(universe)
  
  # Calculate the enrichment of the query in the test relative to the universe
  relative_enrichment = test_overlap_proportion/universe_overlap_proportion
  
  # Create a data.frame with the results
  result_df = data.frame(test_size = length(test), query_size = length(query), universe_size = length(universe),
    test_overlap_size, universe_overlap_size, test_overlap_proportion, universe_overlap_proportion, relative_enrichment, p_value = ft_result$p.value)
  
  # Add the intersection of test and query if requested
  if(return_overlap){
    test_query_overlap = paste(intersect(test, query), collapse = "; ")
     if(length(test_query_overlap) == 0){test_query_overlap = NA}
    result_df$test_query_overlap = test_query_overlap
  }
  
  # Return results data.frame
  return(result_df)
}

#' Test if a list of query vectors is enriched in a test vector using Fisher's exact test
#'
#' @param test A character vector. Should only contain unique elements. 
#' @param query_list A list of character vectors to test for enrichment in test. Each vector should only contain unique elements. Names of list are included in the results.
#' @param universe A vector containing all the elements that test and query_list could contain. Should only contain unique elements. 
#' @param alternative The alternative hypothesis for the test. Default is "greater".
#' @param p_adjust_method Method to use to adjust p-values. Default is "fdr".
#' @param return_overlap A logical value indicating whether to return the overlap between test and query as a character with elements separated by ";". 
#' @return A data.frame with the resuls of Fisher's exact test
#' @export
fisher_test_apply = function(test, query_list, universe, alternative = "greater", p_adjust_method = "fdr", return_overlap = F, ncores = 1){
  
  # Create cluster if ncores greater than 1
  if(ncores > 1){
    cl = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl, ncores)
    `%dopar%` = foreach::`%dopar%`
    on.exit(parallel::stopCluster(cl))
  } else {
    `%dopar%` = foreach::`%do%`
  }
  
  # Create a list with results for each vector in query_list
  results_list = foreach::foreach(query = query_list) %dopar% {
    
    enrichmentTests::fisher_test_vectors(test = test, query = query, universe = universe, 
      alternative = alternative, return_overlap = return_overlap)
    
  }
  
  # Name results_list with names of query_list
  names(results_list) = names(query_list)
    
  # Combine results_list into a single data.frame
  results_df = dplyr::bind_rows(results_list, .id = "query_name")
  results_df$q_value = p.adjust(results_df$p_value, method = p_adjust_method)
  
  # Swap q_value and test_query_overlap columns if return_overlap is set to TRUE
  if(return_overlap){
    results_df = dplyr::relocate(results_df, q_value, .before = test_query_overlap)
  }
  
  # Return results data.frame
  return(results_df)
}

# Create a function which will perform a Chi-squared test for the overlap of CpG sites in a set of query regions 
# with a set of subject regions compared to a set of control regions
query_subject_cpg_overlap_test = function(query_regions, subject_regions, control_regions, background_cpgs, alternative = "greater"){
  
  # Check that alternative is one of the permitted values
  match.arg(arg = alternative, choices = c("two.sided", "greater", "less"), several.ok = F)
  
  # Find CpG sites overlapping the query, subject and control regions
  query_regions_cpgs = subsetByOverlaps(background_cpgs, query_regions, ignore.strand = T)$name
  subject_regions_cpgs = subsetByOverlaps(background_cpgs, subject_regions, ignore.strand = T)$name
  control_regions_cpgs = subsetByOverlaps(background_cpgs, control_regions, ignore.strand = T)$name
  
  # Find the number of CpGs in common between query_regions and subject_regions and control_regions and subject_regions
  query_subject_overlap = intersect(query_regions_cpgs, subject_regions_cpgs)
  control_subject_overlap = intersect(control_regions_cpgs, subject_regions_cpgs)
  
  # Perform a chi-square test
  chi_test_result = prop.test(
    x = c(length(query_subject_overlap), length(control_subject_overlap)), 
    n = c(length(query_regions_cpgs), length(control_regions_cpgs)), 
    correct = F, alternative = alternative
  )
  
  results_df = data.frame(broom::tidy(chi_test_result))
  results_df$query_subject_overlap = length(query_subject_overlap)
  results_df$control_subject_overlap = length(control_subject_overlap)
  
  return(results_df)
  
}

# Create a function which will iterate through a list of GRanges and test for enrichment among a set of query_regions
query_subject_cpg_overlap_test_apply = function(query_regions, control_regions, background_cpgs, gr_list, alternative = "greater"){
  system.time({overlap_test_results = lapply(gr_list, function(x) 
   query_subject_cpg_overlap_test(query_regions = query_regions, subject_regions = x, 
     control_regions = control_regions, background_cpgs = background_cpgs, alternative = alternative))})
  overlap_test_results = dplyr::bind_rows(overlap_test_results, .id = "tf_name")
  overlap_test_results$q_value = p.adjust(overlap_test_results$p.value, method = "fdr")
  overlap_test_results = dplyr::mutate(overlap_test_results, enrichment = estimate1 / estimate2)
}

#' Pairwise intersection lengths
#'
#' Find intersection lengths of all pairs of elements in a list of vectors
#'
#' @param lst A list of vectors
#' @param proportion If proportion is set to TRUE, intersection is giving as a proportion of the 
#' total number of elements from the column group. Default value is FALSE. 
#' @return A matrix with the intersection lengths for all pairs of elements in a list. If list has names, row names and column names are the names of list elements.
#' @export
intersect_lengths_all_pairwise = function(lst, proportion = F){
# Returns the length of the intersection for all pairwise combination of elements in a list as a matrix
  intersection_matrix = matrix(NA, length(lst), length(lst), dimnames = list(names(lst), names(lst)))
  for (r in 1:nrow(intersection_matrix)){
    for (c in 1:ncol(intersection_matrix)){
      if(c <= r){intersection_matrix[r, c] = length(intersect(lst[[c]], lst[[r]]))}
    }
    intersection_matrix[upper.tri(intersection_matrix)] = t(intersection_matrix)[upper.tri(t(intersection_matrix))]
  }
  
  if(proportion){
    unique_elements = lengths(lapply(lst, unique))
    intersection_matrix = sweep(intersection_matrix, 2, unique_elements, FUN = "/")
  }
  
  return(intersection_matrix)
}