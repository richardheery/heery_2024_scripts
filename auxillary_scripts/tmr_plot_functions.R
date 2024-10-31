# Create a function which will add regions upstream and downstream to transcripts
expand_transcripts = function (grl, expand_upstream = 0, expand_downstream = 0){
    
  # Get transcripts from GRL
  transcripts = unlist(reduce(grl))
    
  # Create a temporary GRanges to store start of upstream expansions
  tss_temp = resize(transcripts, width = 1, fix = "start")
  
  # For expansion in expand_upstream
  upstream_sections = foreach(expansion = c(expand_upstream[1], diff(expand_upstream)), label = expand_upstream) %do% {
    
    # Create an upstream section for the relevant expansion
    upstream_section = promoters(tss_temp, upstream = expansion, downstream = 0)
    upstream_section$region = paste0("TSS-", label)
    
    # Update TSS temp
    tss_temp =  adjust_gr(tss_temp, upstream = expansion)
    
    # Return The upstream section
    upstream_section$transcript_id = names(upstream_section)
    upstream_section
    
  }
  
  # Convert downstream_sections into a GRanges
  upstream_sections = unlist(GRangesList(upstream_sections))
  
  # Create a temporary GRanges to store start of downstream expansions
  tes_temp = resize(transcripts, width = 1, fix = "end")
  
  # For expansion in expand_downstream
  downstream_sections = foreach(expansion = c(expand_downstream[1], diff(expand_downstream)), label = expand_downstream) %do% {
    
    # Create an downstream section for the relevant expansion
    downstream_section = plyranges::shift_downstream(promoters(tes_temp, upstream = 0, downstream = expansion), 1)
    downstream_section$region = paste0("TES+", label)
    
    # Update TES temp
    tes_temp = resize(adjust_gr(gr = tes_temp, upstream = 0, downstream = expansion), width = 1, fix = "end")
    
    # Return The downstream section
    downstream_section$transcript_id = names(downstream_section)
    downstream_section
    
  }
  
  # Convert downstream_sections into a GRanges
  downstream_sections = unlist(GRangesList(downstream_sections))
  
  # Create a data.frame from grl
  grl_df = data.frame(unlist(grl))
  complete_df = dplyr::bind_rows(grl_df, data.frame(downstream_sections), 
      data.frame(upstream_sections))
  complete_df = dplyr::arrange(complete_df, transcript_id, start)
  complete_gr = makeGRangesFromDataFrame(complete_df, keep.extra.columns = TRUE, 
      seqinfo = seqinfo(grl))
  complete_gr$region = ifelse(!is.na(complete_gr$rank), 
    paste(complete_gr$region, complete_gr$rank, sep = "_"), complete_gr$region)
  complete_gr$rank = NULL
  complete_grl = GRangesList(split(complete_gr, complete_gr$transcript_id))
  return(unlist(complete_grl[names(grl)]))
}

plot_tmr_regions = function(tmrs, transcript_regions_gr, regions_filter, normalize = T, title = "Distribution of TMRs"){
  
  # Get title size of each class of region
  region_sizes = sapply(split(transcript_regions_gr, transcript_regions_gr$region), function(x)
    sum(width(reduce(x, ignore.strand = T))))
  
  # Find overlaps between tmrs and transcript_regions_gr
  overlaps_df = data.frame(findOverlaps(tmrs, transcript_regions_gr))
  
  # Add names of transcripts, direction and region to overlap results
  overlaps_df$transcript1 = tmrs$ID[overlaps_df$queryHits]
  overlaps_df$direction = tmrs$direction[overlaps_df$queryHits]
  overlaps_df$transcript2 = transcript_regions_gr$transcript_id[overlaps_df$subjectHits]
  overlaps_df$region = transcript_regions_gr$region[overlaps_df$subjectHits]
  
  # Filter for overlaps for TMRs and their associated transcript
  overlaps_df = dplyr::filter(overlaps_df, transcript1 == transcript2)
  
  # Count number of TMRs per region
  overlaps_summary = summarize(group_by(overlaps_df, direction, region), count = n())
  overlaps_summary = tidyr::complete(tibble(overlaps_summary), direction, region, fill = list(count = 0))
  
  # Normalize count by coverage of region
  overlaps_summary$normalized_count = overlaps_summary$count/region_sizes[overlaps_summary$region]*1e6
  overlaps_summary = data.frame(overlaps_summary)
  
  # Filter for desired regions
  tmr_overlaps_count_filtered = dplyr::filter(overlaps_summary, region %in% regions_filter)
  regions_filter = gsub("intron_", "Intron ", regions_filter)
  regions_filter = gsub("exon_", "Exon ", regions_filter)
  regions_filter = gsub("500", "250", regions_filter)
  regions_filter = gsub("1000", "750", regions_filter)
  regions_filter = gsub("1500", "1250", regions_filter)
  regions_filter = gsub("2000", "1750", regions_filter)
  
  # Convert region into a factor
  tmr_overlaps_count_filtered$region = gsub("intron_", "Intron ", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = gsub("exon_", "Exon ", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = gsub("500", "250", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = gsub("1000", "750", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = gsub("1500", "1250", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = gsub("2000", "1750", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = factor(tmr_overlaps_count_filtered$region, levels = regions_filter)
  
  # Update count to normalized count if normalize is TRUE
  ylab = "TMR Count"
  if(normalize){
    tmr_overlaps_count_filtered$count = tmr_overlaps_count_filtered$normalized_count
    ylab = "Number of TMRs per MB"
  }
  
  # Add a column indicating if regions are upstream, in the transcribed region or downstream
  tmr_overlaps_count_filtered$region_class = "Transcribed Region"
  tmr_overlaps_count_filtered$region_class[grep("TSS", tmr_overlaps_count_filtered$region)] = "Upstream of TSS"
  tmr_overlaps_count_filtered$region_class[grep("TES", tmr_overlaps_count_filtered$region)] = "Downstream of TES"
  tmr_overlaps_count_filtered$region_class = factor(tmr_overlaps_count_filtered$region_class, 
    levels = c("Upstream of TSS", "Transcribed Region", "Downstream of TES"))
  
  # Create plot of number of TMRs per region
  regions_plot_normalized = ggplot(tmr_overlaps_count_filtered, aes(x = region, y = count, fill = direction)) +
    geom_col(position = "dodge") + scale_fill_manual(values = c(colour_list$purple_and_gold_light), drop = F)
  regions_plot_normalized = customize_ggplot_theme(regions_plot_normalized, 
    title = NULL, 
    xlab = NULL, ylab = ylab,
    fill_colors = c(colour_list$purple_and_gold_light), x_labels_angle = 45) +
    scale_fill_manual(values = c(colour_list$purple_and_gold_light), drop = F)
  regions_plot_normalized + facet_grid(~region_class, scales = "free_x", space = "free_x") + 
    theme(panel.spacing = unit(0,'lines'), strip.background = element_blank(),
      strip.text = element_text(size = 16))
  
}

# Make a function which will bin relative TMRs
bin_relative_tmrs = function(tmrs, width, transcripts_subset = NULL){
  
  # If transcripts_subset provided, subet for TMRs associated with these transcripts
  if(!is.null(transcripts_subset)){
      tmrs = tmrs[tmrs$transcript_id %in% transcripts_subset]
  }
  
  # Convert the clusters to relative ranges
  relative_ranges_tmrs = 
    methodicalFinal::relative_ranges(tmrs, reference_positions = GRanges(tmrs$tss_location))

  # Bin the relative ranges into 500 bp windows
  binned_ranges_tmrs = bin_relative_ranges(relative_ranges = relative_ranges_tmrs, bin_start = -width, 
    bin_end = width, bin_step = 500, category = tmrs$direction)
  
  # Convert to long format with separate rows for negative and positive TMRs
  binned_ranges_tmrs = reshape2::melt(binned_ranges_tmrs, id.vars = "bin_center")
  
  return(binned_ranges_tmrs)
  
}

