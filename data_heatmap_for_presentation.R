library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
heatmap_cols <- viridis(3, option = 'D', direction = 1)
png("presentation_example_heatmap.png", width = 8, height = 2, units = "in", res = 600)
pheatmap(Y_CSHQ[which(rowSums(Y_CSHQ) == 48)[c(1,2,3,4,10)],], 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         main = "All 5 participants have identical CSHQ total = 48",
         color = heatmap_cols, 
         breaks = c(0.5, 1.5, 2.5, 3.5), 
         legend_breaks = c(1, 2, 3),
         cellwidth = 14,
         cellheight = 18,
         legend_labels = c("1 (Rarely)", "2 (Sometimes)", "3 (Usually)"),
         fontsize = 6,
         fontsize_col = 10,
         fontsize_row = 10,
         labels_row = 1:5)
dev.off()

pheatmap(Y_CSHQ_reduced, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         color = heatmap_cols,
         breaks = c(0.5, 1.5, 2.5, 3.5), 
         legend_breaks = c(1, 2, 3))
pheatmap(Y_CSHQ_reduced[order(CSHQ_4group_minVI_cluster$cl, rowSums(Y_CSHQ_reduced)),], 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         color = heatmap_cols,
         breaks = c(0.5, 1.5, 2.5, 3.5), 
         legend_breaks = c(1, 2, 3),
         gaps_row = c(61,103,136))


#Trying the same thing in ggplot


df_long <- as.data.frame(Y_CSHQ_reduced) |>
  mutate(.row = row_number()) |>
  pivot_longer(-.row, names_to = ".col", values_to = "val") |>
  mutate(
    .col = factor(.col, levels = colnames(Y_CSHQ_reduced)),
    .row = factor(.row, levels = rev(seq_len(nrow(Y_CSHQ_reduced))))
  )

p1_discrete <- ggplot(mutate(df_long, val = factor(val, levels = c(1,2,3))),
                      aes(x = .col, y = .row, fill = val)) +
  geom_tile(color = NA, show.legend = FALSE) +
  scale_fill_manual(values = heatmap_cols, drop = FALSE, labels = c("1 (Rarely)", "2 (Sometimes)", "3 (Usually)")) +
  coord_cartesian() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )
p1_discrete

ggsave(filename = 'full_data_heatmap_unsorted.png', plot = p1_discrete, dpi = 600)







row_order <- order(CSHQ_4group_minVI_cluster$cl, -rowSums(Y_CSHQ_reduced))          
Y_ord <- Y_CSHQ_reduced[row_order,, drop = FALSE]                                   
gaps <- c(61, 103, 136)                                                              

idx <- c(0, gaps, nrow(Y_ord))                                                       
blocks <- Map(function(a,b) seq(a+1, b), head(idx, -1), tail(idx, -1))               
rev_rows <- unlist(rev(blocks))                                                      
Y_rev <- Y_ord[rev_rows, , drop = FALSE]                                             


df_rev <- as.data.frame(Y_rev) |>
  mutate(.row = row_number()) |>
  pivot_longer(-.row, names_to = ".col", values_to = "val") |>
  mutate(
    .col = factor(.col, levels = colnames(Y_CSHQ_reduced)),
    .row = factor(.row, levels = seq_len(nrow(Y_rev))),
    val  = factor(val, levels = c(1,2,3))
  )                                                                                   


p_base <- ggplot(df_rev, aes(.col, .row, fill = val)) +
  geom_tile() +
  scale_fill_manual(values = heatmap_cols, drop = FALSE, labels = c("1 (Rarely)", "2 (Sometimes)", "3 (Usually)")) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )                                                                                  


boundaries <- c(0, 12, 45, 87, 148)                                                 
bar_y_num  <- boundaries + 0.5                                                      


row_levels <- levels(df_rev$.row)
n_rows <- length(row_levels)


bar_y_fac <- factor(pmin(as.integer(boundaries + 1), n_rows), levels = seq_len(n_rows))  


mid_rows <- (head(boundaries, -1) + tail(boundaries, -1)) / 2
label_y_fac <- factor(as.integer(round(mid_rows + 0.5)), levels = seq_len(n_rows))      

group_labels <- c("Group D (Elevated Difficulties)", "Group C (Short/Light)", "Group B (Anxious)", "Group A (Comparatively Healthy)")                            


x_min <- 0.5         
tick_len <- 2.0       
label_offset <- 0.40  

p_rev <- p_base +
  coord_cartesian(clip = "off") +                                                         
  theme(plot.margin = margin(5.5, 5.5, 5.5, 65)) +                                       
  geom_segment(
    data = data.frame(y = bar_y_fac),
    aes(x = x_min - tick_len, xend = x_min, y = y, yend = y),
    inherit.aes = FALSE, linewidth = 1.6, color = "black", lineend = "butt"
  ) +
  geom_text(
    data = data.frame(y = label_y_fac, lab = group_labels),
    aes(x = x_min - tick_len - label_offset, y = y, label = lab),
    inherit.aes = FALSE, hjust = 1, vjust = 0.5, size = 4.2, color = "black"
  )


p_rev                                        


ggsave(filename = 'full_data_heatmap_sorted.png', plot = p_rev, dpi = 600)





































row_order <- order(CSHQ_4group_minVI_cluster$cl, -rowSums(Y_CSHQ_reduced))          
Y_ord <- Y_CSHQ_reduced[row_order,, drop = FALSE]                                   
gaps <- c(61, 103, 136)                                                              

idx <- c(0, gaps, nrow(Y_ord))                                                       
blocks <- Map(function(a,b) seq(a+1, b), head(idx, -1), tail(idx, -1))               
rev_rows <- unlist(rev(blocks))                                                      
Y_rev <- Y_ord[rev_rows, , drop = FALSE]                                             

df_rev <- as.data.frame(Y_rev) |>
  mutate(.row = row_number()) |>
  pivot_longer(-.row, names_to = ".col", values_to = "val") |>
  mutate(
    .col = factor(.col, levels = colnames(Y_CSHQ_reduced)),
    .row = factor(.row, levels = seq_len(nrow(Y_rev))),
    val  = factor(val, levels = c(1,2,3))
  )                                                                                   

p_base <- ggplot(df_rev, aes(.col, .row, fill = val)) +
  geom_tile() +
  scale_fill_manual(values = heatmap_cols, drop = FALSE, labels = c("1 (Rarely)", "2 (Sometimes)", "3 (Usually)")) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )                                                                                  

boundaries <- c(0, 12, 45, 87, 148)                                                 
bar_y_num  <- boundaries + 0.5                                                      

row_levels <- levels(df_rev$.row)
n_rows <- length(row_levels)

bar_y_fac <- factor(pmin(as.integer(boundaries + 1), n_rows), levels = seq_len(n_rows))  

mid_rows <- (head(boundaries, -1) + tail(boundaries, -1)) / 2
label_y_fac <- factor(as.integer(round(mid_rows + 0.5)), levels = seq_len(n_rows))      

group_labels <- c("Group D\n(Elevated)", "Group C\n(Short/Light)", "Group B\n(Anxious)", "Group A\n(Healthy)")                            

x_min <- 0.5         
tick_len <- 2.0       
label_offset <- 0.50  # Increased slightly for better spacing

p_rev <- p_base +
  coord_cartesian(clip = "off") +                                                         
  theme(plot.margin = margin(5.5, 5.5, 5.5, 85)) +  # Increased left margin from 65 to 85                                     
  geom_segment(
    data = data.frame(y = bar_y_fac),
    aes(x = x_min - tick_len, xend = x_min, y = y, yend = y),
    inherit.aes = FALSE, linewidth = 1.6, color = "black", lineend = "butt"
  ) +
  geom_text(
    data = data.frame(y = label_y_fac, lab = group_labels),
    aes(x = x_min - tick_len - label_offset, y = y, label = lab),
    inherit.aes = FALSE, hjust = 1, vjust = 0.5, size = 3.8, color = "black"  # Reduced size from 4.2 to 3.5
  )

p_rev                                        

ggsave(filename = 'full_data_heatmap_sorted.png', plot = p_rev, 
       width = 12, height = 8, dpi = 600, units = "in")

