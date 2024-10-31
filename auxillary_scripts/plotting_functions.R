library(ggplot2)

#' Customize the theme of a ggplot object
#'
#' @param plot A ggplot object
#' @param base_theme The base theme to use for the plot. Default is theme_bw().
#' @param title Title for the plot
#' @param plot_title_size Size of the plot title in pts. Default is 24.
#' @param xlab The x-axis title
#' @param ylab The y-axis title
#' @param axis_title_size Size of the axes titles in pts. Default is 20.
#' @param axis_text_size Size of the axes text in pts. Default is 18.
#' @param x_labels Optional labels to use for the x-axis
#' @param format_x_labels Logical value indicating whether to format the X axis labels by replacing "_" with a space and using title case. Default is FALSE.
#' @param x_labels_angle Number of degrees to rotate x-axis labels. Default is 0. 
#' @param show_legend A logical value indicating whether to show the legend or not. Default is TRUE. 
#' @param fill_title The title for fill in the legend
#' @param color_title The title for color in the legend
#' @param colors An optional vector of colors to use for the color mapping. 
#' @param fill_colors An optional vector of colors to use for the fill. 
#' @param fill_labels An optional vector of labels for the fill legend
#' @param legend_title_size Size of the legend titles in pts. Default is 20.
#' @param legend_text_size Size of the legend text in pts. Default is 18.
#' @param legend_key_size Size of the legend key in cm. Default is 1.
#' @param format_fill_labels Logical value indicating whether to format the legend labels by replacing "_" with a space and using title case. Default is FALSE.
#' @param scale_x A scale to use for the x-axis. If none provided, one is chosen by default.
#' @param scale_y A scale to use for the y-axis. If none provided, one is chosen by default. 
#' @param facet A variable to use for faceting. Default is not to facet. 
#' @param facet_nrow Number of rows to use for faceting. 
#' @param facet_ncol Number of rows to use for faceting. 
#' @param facet_scales scale for facet_wrap(). Default is "free". 
#' @param strip_text_size Size of the facet labels. Default is 20.
#' @param facet_labels Labels to use for panels when facetting
#' @return A ggplot object
#' @export
customize_ggplot_theme = function(plot, base_theme = theme_bw(), title = NULL, plot_title_size = 24,
  xlab = NULL, ylab = NULL, axis_title_size = 20, axis_text_size = 18, x_labels = NULL, format_x_labels = F, x_labels_angle = 0,  
  show_legend = T, fill_title = NULL, color_title = NULL, fill_colors = NULL, colors = NULL, fill_labels = waiver(), 
  legend_title_size = 20, legend_text_size = 18, legend_key_size = 1, format_fill_labels = F, 
  scale_x = NULL, scale_y = NULL, 
  facet = NULL, facet_nrow = NULL, facet_ncol = NULL, facet_scales = "free", strip_text_size = 20, facet_labels = NULL){
  
  plot = plot +
    base_theme +
    theme(plot.title = element_text(hjust = 0.5, size = plot_title_size), 
	    axis.title = element_text(size = axis_title_size), axis.text = element_text(size = axis_text_size), 
  	  legend.title = element_text(size = legend_title_size), legend.text = element_text(size = legend_text_size), 
      legend.key.size = unit(legend_key_size, "cm"), 
      strip.text = element_text(size = strip_text_size)) +
    labs(title = title, x = xlab, y = ylab, fill = fill_title,  color = color_title)
  
  # Remove legend if specified
  if(!show_legend){
    plot = plot + theme(legend.position = "None")
  }
  
  # Rotate x-axis lables if specified
  if(x_labels_angle != 0){
    plot = plot + theme(axis.text.x = element_text(angle = x_labels_angle, hjust = 1))
  }
  
  # Use x_labels if provided
  if(!is.null(x_labels)){
    plot = plot + scale_x_discrete(labels = x_labels)
  }
  
  # Format x axis labels if specified
  if(format_x_labels){
    xlabels = ggplot_build(plot)$layout$panel_params[[1]]$x$get_labels()
    new_xlabels = stringr::str_to_title(gsub("_", " ", xlabels))
    plot = plot + scale_x_discrete(labels = new_xlabels)
  }
  
  # Find the type of geom the plot uses
  plot_geom = class(plot$layers[[1]]$geom)
  
  # Change fill colours and labels if specified
  if(!is.null(fill_colors)){
    plot = plot + scale_fill_manual(values = fill_colors, labels = fill_labels)
  } else {
    plot = plot + scale_fill_discrete(labels = fill_labels)
  }
  
  # Change colors if specified
  if(!is.null(colors)){
    plot = plot + scale_color_manual(values = colors)
  }
  
  # If scale_x and scale_y not provided. they are inferred
  if(is.null(scale_x)){
    if("GeomBar" %in% plot_geom & !"GeomCol" %in% plot_geom){
      scale_x = scale_x_continuous(expand = c(0,0), labels = scales::comma) 
    }
  }
  
  if(is.null(scale_y)){
    if("GeomBar" %in% plot_geom){
      scale_y = scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::comma)
    }
  }
  
  plot = plot + scale_x + scale_y
  
  if(!is.null(facet)){
    if(!is.null(facet_labels)){
      original_facet_labels = sort(unique(plot$data[[facet]]))
      labeller = as_labeller(setNames(facet_labels, original_facet_labels))
    } else {
      labeller = "label_value"
    }
    plot = plot + facet_wrap(facet, nrow = facet_nrow, ncol = facet_ncol, labeller = labeller, scales = facet_scales)
  }
  
  return(plot)
  
}

#' Encode p-values with symbols
#'
#' @param p_values A vector of p-values
#' @param cutpoints A vector of cutpoints for significance. Default is c(0, 0.001, 0.01, 0.05, 1)
#' @param symbol
#' @return A vector of symbols representing p-values
#' @export
sig_sym = function(p_values, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbol = "*"){
  
  symbols = sapply(3:0, function(x) paste(rep(symbol, x), collapse = ""))
  
  significance_symbols = as.character(symnum(p_values, corr = FALSE, na = FALSE, cutpoints = cutpoints, symbols = symbols))
  return(significance_symbols)
}

#' Create a heatmap from a matrix without clustering
#'
#' @param mat A matrix or object which can ne coerced to a matrix.
#' @param colors A vector if colors to use for heatmap. Default is colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Reds"))(100). 
#' Can also just supply a single colour and it will be interpolated using colorRampPalette().
#' @param reverse_colors A logical value indicating whether to reverse the supplied or created colors e.g. if all values are negative. Default is FALSE.
#' @param breaks Vector of breaks associated with colors. Default is NA.
#' @param mask_diagonal Logical value indicating whether main diagonal should be masked. Default is T.
#' @param decimal_places Number of decimal places to display. Default is 2.
#' @param title Title of the plot. Default is no title. 
#' @param row_labels Labels to use for rows. Default is to use row names of mat.
#' @param col_labels Labels to use for columns. Default is to use column names of mat.
#' @param angle_col Angle to rotate column labels. Default is 270. 
#' @param title_size Fontsize for the title. Default is 18.
#' @param label_size Fontsize for the labels. Default is 15.  
#' @param number_size Fontsize for numbers. Default is 18.
#' @param filename Optional filename to save plot with. 
#' @param file_dimensions A vector of length two with the width and height of the file in inches. 
#' Default is a width of 6 and height of 6.  
#' @param return_ggplot A logical value indicating whether the heatmap should be returned as a ggplot object
#' @return A pheatmap object or ggplot object dericed from it depending on 
#' @export
heatmap_without_clustering = function(mat, colors = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Reds"))(100),
  reverse_colors = F, breaks = NA, mask_diagonal = T, decimal_places = 2, title = "", row_labels = NULL, col_labels = NULL, 
  angle_col = 270, title_size = 18, label_size = 15, number_size = 18, filename = NA, file_dimensions = c(6, 6), return_ggplot = F){
  
  # Interpolate colors if just a single colour provided
  if(length(colors) == 1){
    colors = colorRampPalette(c("white", colors))(100)
  }
  
  # Reverse colurs if specified
  if(reverse_colors){colors = rev(colors)}
  
  # Ensure mat is a matrix
  mat = as.matrix(mat)
  
  # Round values to specified number of places and convert to a character
  char_mat = matrix(as.character(round(mat, decimal_places)), ncol = ncol(mat))
  
  # Mask the main diagonal if specified
  if(mask_diagonal){
    diag(mat) = NA
    diag(char_mat) = ""
  } 
  
  # Use row names of mat as row labels if none provided
  if(is.null(row_labels)){
    row_labels = row.names(mat)
  }
  
  # Use column names of mat as row labels if none provided
  if(is.null(col_labels)){
    col_labels = colnames(mat)
  }
  
  # Create heatmap
  heatmap = pheatmap::pheatmap(mat = mat, cluster_rows = F, cluster_cols = F, 
    color = colors, breaks = breaks, na_col = "white", legend = F, angle_col = angle_col, 
    display_numbers = char_mat, number_color = "black", border_color = "black", 
    main = title, labels_row = row_labels, labels_col = col_labels, 
    fontsize = title_size, fontsize_row = label_size, fontsize_col = label_size, fontsize_number = number_size, 
    filename = filename, width = file_dimensions[1], height = file_dimensions[2])
  
  # Reset graphics device if saving plot
  if(.Device == "pdf"){dev.off()}
  
  # Return heatmap
  if(return_ggplot){
    return(ggplotify::as.ggplot(heatmap$gtable))
  } else {
    return(heatmap)
  }
  
}

colour_list = list(
  ggplot_red_blue = c("#F8766D", "#00BFC4"),
  two_blues = c("#1F8AC0", "#104C91"),
  red_navy_pair = c("#D01C1FFF", "#4B878BFF"),
  gold_green_red = c("#DDCC77", "#117733", "#882255"),
  three_greens = c("#c4e6c3", "#4da284", "#1d4f60"),
  three_blues = c("#d1eeea", "#68abb8", "#2a5674"),
  three_reds = c("#FCBBA1", "#FB6A4A", "#A50F15"),
  three_purples = c("#BCBDDC", "#807DBA", "#54278F"),
  rainbox_six = c("#9E0142", "#FDAE61", "#FEE08B", "#ABDDA4", "#3288BD", "#5E4FA2"),
  nine_greens = RColorBrewer::brewer.pal(9, "YlGn"),
  purple_and_gold_dark = c("#7B5C90", "#bfab25"),
  purple_and_gold_light = c("#A28CB1", "#D2C465")
)

#' Save plots to a PDF file
#'
#' @param plotlist A list of ggplot objects
#' @param titles A vector of titles for the plots
#' @return A ggplot object
#' @export
title_ggplots = function(plotlist, titles){
  
  if(length(plotlist) != length(titles)){stop("plotlist and titles must be the same length")}
    
  lapply(seq_along(plotlist), function(x) plotlist[[x]] + ggtitle(titles[[x]]))
}

#' Save a list of plots to a PDF file
#'
#' @param plotlist A list of ggplot objects
#' @param nrows Number of rows in each sheet
#' @param ncols Number of columns in each sheet
#' @param width Width of each page. Default is 16.
#' @param height Heigth of each page. Default is 9.
#' @param filename Name of output file
#' @return Nothing
#' @export
pdf_save = function(plotlist, nrows = 1, ncols = 1, filename, width = 16, height = 9){
  on.exit(expr = dev.off())
  pdf(filename, width = width, height = height, title = basename(filename))
  for(x in seq(1, length(plotlist)/prod(nrows, ncols))) {
    print(ggpubr::ggarrange(plotlist = plotlist[x:(x+prod(nrows, ncols)-1)], nrow = nrows, ncol = ncols))
  }
  dev.off()
  on.exit()
}
