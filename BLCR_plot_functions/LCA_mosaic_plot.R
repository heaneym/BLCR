library(ggplot2)
library(ggmosaic)

#This function creates a list of ggplots, each list item being a mosaic plot corresponding to an 
#item extracted from LCA.

#Parameters are:

# itemprob - the item probability parameters (theta) from LCA formatted as a list of matrices of length m,
# with each matrix having dimension G x K_j, (G = #groups, K_j = #levels for item j)

# classprob - class probability parameters/mixing proportions (pi) from LCA formatted as a vector of length G.

# group_labels - the group labels indicated on the x-axis. This is a vector of length G. If = NULL, default = letters A,B,C,D,...

# item_titles - vector of length M with custom titles for each item. only displayed if show_title = TRUE
# If no names are assigned to the theta list, the titles displayed are Item 1, Item 2, etc.

# show_title - logical indicating whether the title is displayed. default is true, and titles are the names of the items taken from theta.

# show_y_axis_numbers - logical specifying whether to display a y-axis scale.

# custom_colours - vector of length max(K) specifying a custom colour palette from 1:K_j for each item j. 


plot_mosaic_gg_LCA <- function(
    itemprob,
    classprob,
    group_labels = NULL,
    item_titles = NULL,
    show_title = TRUE,
    show_y_axis_numbers = FALSE,
    custom_colours = NULL
){
  G <- length(classprob)
  M <- length(itemprob)
  K <- unlist(lapply(itemprob, ncol))
  xmax <- cumsum(classprob)
  xmin <- c(0, head(xmax, -1))
  mosaic_plot_list <- list()
  
  if (is.null(group_labels)) {
    group_labels <- LETTERS[1:G]
  }
  if (length(group_labels) != G) {
    stop("Length of group_labels must match number of groups (G).")
  }
  
  if (is.null(item_titles)) {
    if (is.null(names(itemprob))) {
      item_titles <- paste("Item", 1:M)
    } else {
      item_titles <- names(itemprob)
    }
  } else {
    if (length(item_titles) != M) {
      stop("Length of item_titles must match number of items (M).")
    }
  }

  max_K <- max(K)
  if (is.null(custom_colours)) {
    extended_colours <- c(
      "#0072B2",  
      "#E69F00",  
      "#D55E00",  
      "#009E73",  
      "#F0E442",  
      "#56B4E9",  
      "#CC79A7",  
      "#999999"   
    )
    if (max_K > length(extended_colours)) {
      colour_func <- colorRampPalette(extended_colours)
      extended_colours <- colour_func(max_K)
    }
    
    cb_palette <- setNames(
      extended_colours[1:max_K],
      paste0("Response ", 1:max_K)
    )
  } else {
    if (length(custom_colors) < max_K) {
      stop(paste0("custom_colors must have at least ", max_K, " colors."))
    }
    cb_palette <- setNames(
      custom_colors[1:max_K],
      paste0("Response ", 1:max_K)
    )
  }
  
  for (j in 1:M){
    plot_data <- data.frame()
    for (g in 1:G) {
      y_cum <- c(0, cumsum(itemprob[[j]][g, ]))
      for (c in 1:K[j]) {
        plot_data <- rbind(plot_data, data.frame(
          group = group_labels[g],
          category = paste0("Response ", c),
          xmin = xmin[g],
          xmax = xmax[g],
          ymin = y_cum[c],
          ymax = y_cum[c + 1]
        ))
      }
    }
    
    p <- ggplot(plot_data) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = category), color = "black") +
      scale_fill_manual(
        values = cb_palette,
        labels = 1:max(K[j]),
        guide  = "legend"     
      ) +
      scale_x_continuous(
        breaks = (xmin + xmax) / 2,
        labels = group_labels
      ) +
      labs(x = "", y = "", fill = "Response") +
      theme_minimal() +                     
      theme(
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        axis.text.y = if (show_y_axis_numbers) element_text() else element_blank(),
        axis.ticks.y = if (show_y_axis_numbers) element_line() else element_blank(),
        panel.grid = element_blank()        
      )
    
    if(show_title) {
      p <- p + labs(title = item_titles[j]) + 
        theme(plot.title = element_text(hjust = 0.5))
    }
    
    mosaic_plot_list[[j]] <- p
  }
  
  return(mosaic_plot_list)
}

