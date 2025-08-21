library(ggplot2)
library(ggmosaic)


plot_mosaic_gg_LCA <- function(
    itemprob, 
    classprob, 
    group_labels = NULL, 
    show_title = TRUE, 
    show_y_axis_numbers = FALSE
){
  G <- length(classprob)
  M <- length(itemprob)
  K <- unlist(lapply(itemprob, ncol))
  xmax <- cumsum(classprob)
  xmin <- c(0, head(xmax, -1))
  mosaic_plot_list <- list()
  
  # Default group labels: A, B, C, ...
  if (is.null(group_labels)) {
    group_labels <- LETTERS[1:G]
  }
  if (length(group_labels) != G) {
    stop("Length of group_labels must match number of groups (G).")
  }
  
  # Colorblind-friendly palette
  cb_palette <- c(
    "Response 1" = "#0072B2",
    "Response 2" = "#E69F00",
    "Response 3" = "#D55E00"
  )
  
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
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category), color = "black") +
      scale_fill_manual(values = cb_palette, labels = c(1, 2, 3)) +
      scale_x_continuous(
        breaks = (xmin + xmax) / 2,
        labels = group_labels
      ) +
      labs(x = "", y = "", fill = "Response") +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        axis.text.y = if (show_y_axis_numbers) element_text() else element_blank(),
        axis.ticks.y = if (show_y_axis_numbers) element_line() else element_blank()
      )
    if(show_title) {
      p <- p + labs(title = names(itemprob)[j]) + theme(plot.title = element_text(hjust = 0.5))
    }
    mosaic_plot_list[[j]] <- p
  }
  return(mosaic_plot_list)
}


save_mosaic_plots <- function(plot_list, itemprob, dir = ".", base_filename = "mosaic_plot", height = 4, width = 6, dpi = 300) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  item_names <- names(itemprob)
  if (is.null(item_names)) {
    item_names <- paste0("item", seq_along(plot_list))
  }
  for (i in seq_along(plot_list)) {
    safe_item_name <- gsub("[^A-Za-z0-9_]", "_", item_names[i])
    filename <- paste0(base_filename, "_", safe_item_name, ".png")
    filepath <- file.path(dir, filename)
    ggsave(filename = filepath, plot = plot_list[[i]], width = width, height = height, dpi = dpi)
    message("Saved: ", filepath)
  }
}

