#Function for saving a list of ggplots at once - convenient for a list of mosaic plots from LCA.

save_plot_list <- function(plot_list, 
                           plot_names = NULL,
                           dir = ".", 
                           base_filename = "plot",
                           format = "png",
                           height = 4, 
                           width = 6, 
                           dpi = 300) {
  
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  valid_formats <- c("png", "pdf", "svg", "jpg", "jpeg", "tiff", "bmp", "eps", "ps")
  format <- tolower(format)
  if (!format %in% valid_formats) {
    stop("Invalid format. Must be one of: ", paste(valid_formats, collapse = ", "))
  }
  if (is.null(plot_names)) {
    plot_names <- names(plot_list)
  }
  if (is.null(plot_names)) {
    plot_names <- paste0("plot", seq_along(plot_list))
  }
  for (i in seq_along(plot_list)) {
    safe_name <- gsub("[^A-Za-z0-9_]", "_", plot_names[i])
    filename <- paste0(base_filename, "_", safe_name, ".", format)
    filepath <- file.path(dir, filename)
    
    ggsave(filename = filepath, 
           plot = plot_list[[i]], 
           width = width, 
           height = height, 
           dpi = dpi,
           device = format)
    
    message("Saved: ", filepath)
  }
  
  invisible(filepath)
}