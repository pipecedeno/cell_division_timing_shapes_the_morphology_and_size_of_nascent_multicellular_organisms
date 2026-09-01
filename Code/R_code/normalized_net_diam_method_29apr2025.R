
library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(stringi)
# library(glue)
library(purrr)
library(stringr)
library(ggpubr)
library(tidyverse)
library(cowplot)
library(igraph)
library(gridExtra)
library(grid)


# Function to grow an exponential tree
grow_exponential_tree <- function(network) {
  # Get the node and edge list and rename nodes to create a new edge list
  nodes <- V(network)
  edges <- as_edgelist(network)
  
  highest_node <- max(nodes)
  num_nodes <- length(nodes)
  
  # Create new nodes and add them to the network
  new_nodes <- (num_nodes+1):(num_nodes+num_nodes)
  network <- add_vertices(network, length(new_nodes))
  
  # Make the connection of the old branch to the new branch
  lowest_node <- min(new_nodes)
  network <- add_edges(network, c(1, lowest_node))
  
  # Only process renamed edges if there are any edges
  if (nrow(edges) > 0) {
    renamed_edges <- matrix(NA, nrow = nrow(edges), ncol = 2)
    for (i in 1:nrow(edges)) {
      renamed_edges[i, 1] <- highest_node + edges[i, 1]
      renamed_edges[i, 2] <- highest_node + edges[i, 2]
    }
    
    # Add the renamed edges to the network
    network <- add_edges(network, t(renamed_edges))
  }
  
  return(network)
}


# Plotting networks to test function ####
# Create initial network with 1 node
snowflake <- make_empty_graph(n = 1, directed = FALSE)

# Set up plotting area
par(mfrow = c(2, 3))

# First plot - initial single node
plot(snowflake, 
     vertex.size = 10,
     vertex.label.cex = 0.8,
     vertex.color = "lightblue",
     vertex.label.color = "black",
     edge.color = "gray",
     main = "Initial")

# Second plot - after first growth
snowflake <- grow_exponential_tree(snowflake)
plot(snowflake, 
     vertex.size = 10,
     vertex.label.cex = 0.8,
     vertex.color = "lightblue",
     vertex.label.color = "black",
     edge.color = "gray",
     main = "First Growth")

# Third plot - after second growth
snowflake <- grow_exponential_tree(snowflake)
plot(snowflake, 
     vertex.size = 10,
     vertex.label.cex = 0.8,
     vertex.color = "lightblue",
     vertex.label.color = "black",
     edge.color = "gray",
     main = "Second Growth")

# Fourth plot - after third growth
snowflake <- grow_exponential_tree(snowflake)
plot(snowflake, 
     vertex.size = 10,
     vertex.label.cex = 0.8,
     vertex.color = "lightblue",
     vertex.label.color = "black",
     edge.color = "gray",
     main = "Third Growth")

snowflake <- grow_exponential_tree(snowflake)
plot(snowflake, 
     vertex.size = 10,
     vertex.label.cex = 0.8,
     vertex.color = "lightblue",
     vertex.label.color = "black",
     edge.color = "gray",
     main = "Fourth Growth")

snowflake <- grow_exponential_tree(snowflake)
plot(snowflake, 
     vertex.size = 10,
     vertex.label.cex = 0.8,
     vertex.color = "lightblue",
     vertex.label.color = "black",
     edge.color = "gray",
     main = "Fifth Growth")


# Simulating growing network ####
# Create initial network with 2 nodes and one edge
snowflake <- make_empty_graph(n = 2, directed = FALSE)
snowflake <- add_edges(snowflake, c(1, 2))

num_cells <- c(vcount(snowflake))
diameter <- c(diameter(snowflake))

num_gen <- 8

for (i in 1:num_gen) {
  snowflake <- grow_exponential_tree(snowflake)
  
  diameter <- c(diameter, diameter(snowflake))
  num_cells <- c(num_cells, vcount(snowflake))
}

num_cells <- as.numeric(num_cells)
diameter <- as.numeric(diameter)


# Making normalized network metric ####
# Fit a linear model to log-transformed data
p <- coef(lm(diameter ~ log10(num_cells)))

# Create plot 1: Basic relationship
p1=ggplot(data = data.frame(num_cells, diameter), aes(x = num_cells, y = diameter))+
  geom_point(col="cornflowerblue")+
  geom_line(col="cornflowerblue")+
  labs(x = "Number of cells", y = "Network diameter")+
  theme_classic(base_size = 10)+
  NULL
p1

# Create plot 2: Log scale with fitted line
p2 <- ggplot(data = data.frame(num_cells, diameter), aes(x = num_cells, y = diameter))+
  geom_point(col="cornflowerblue")+
  geom_line(col="cornflowerblue")+
  geom_line(aes(y = p[1] + p[2] * log10(num_cells)), linetype = "dashed", color = "black")+
  scale_x_log10()+
  labs(x = expression(log[10](Number~of~cells)), y = "Network diameter")+
  theme_classic(base_size = 10)+
  NULL
p2

# Create plot 3: Normalized network diameter
p3 <- ggplot(data = data.frame(num_cells, norm_diam = diameter / (p[2] * log10(num_cells) + p[1])), 
             aes(x = num_cells, y = norm_diam))+
  geom_line(col="cornflowerblue")+
  ylim(0, 2)+
  labs(x = "Number of cells", y = "Normalized network diameter") +
  theme_classic(base_size = 10)+
  NULL
p3

supp_fig_4=plot_grid(p1, p2, p3, labels=c('A', 'B', 'C'), ncol=3, align='hv')
supp_fig_4

# ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig4_normalized_net_diam_method_29apr2025.png',
#        plot=supp_fig_4, dpi='retina', width=9, height=3)


# Adding strains data ####


strain_diam_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/fig_4_network_growth_with_fragmentation/mean_diameter_values_3oct2025.csv", header=TRUE)
colnames(strain_diam_df)

exponential_df=data.frame(strain="Exponential",
                          num_nodes=num_cells,
                          mean_diameter=diameter)

temp_strain_df=rbind(data.frame(strain=c("Grande", "Petite w/o Delay", "Petite"),
                                num_nodes=c(2, 2, 2),
                                mean_diameter=c(1,1,1)),
                     strain_diam_df[strain_diam_df$num_nodes<=512 & strain_diam_df$num_nodes>1, c("strain", "num_nodes", "mean_diameter")])

network_df=rbind(exponential_df, temp_strain_df)
network_df$strain=factor(network_df$strain, levels=c("Exponential", "Grande", "Petite w/o Delay", "Petite"))

p1_strain=ggplot(network_df, aes(x = num_nodes, y = mean_diameter, color=strain))+
  geom_line()+
  geom_point(data=network_df[network_df$strain=='Exponential',],
             aes(x=num_nodes, y=mean_diameter))+
  scale_color_manual(values=c("cornflowerblue", "#AE93BEFF", "#B4DAE5FF", "#F0D77BFF"))+
  labs(x = "Number of cells", y = "Network diameter")+
  theme_classic(base_size = 10)+
  guides(col='none')+
  NULL
p1_strain

p2_strain <- ggplot(network_df, aes(x = num_nodes, y = mean_diameter, color=strain))+
  geom_line()+
  geom_point(data=network_df[network_df$strain=='Exponential',],
             aes(x=num_nodes, y=mean_diameter))+
  geom_line(data=network_df[network_df$strain=='Exponential',],
            aes(y = p[1] + p[2] * log10(num_nodes)), linetype = "dashed", color = "black")+
  scale_x_log10()+
  labs(x = "Number of cells", y = "Network diameter")+
  scale_color_manual(values=c("cornflowerblue", "#AE93BEFF", "#B4DAE5FF", "#F0D77BFF"))+
  theme_classic(base_size = 10)+
  guides(col='none')+
  NULL
p2_strain

network_df$norm_diameter=network_df$mean_diameter/ (p[2] * log10(network_df$num_nodes) + p[1])

p3_strain <- ggplot(network_df, aes(x = num_nodes, y = norm_diameter, color=strain))+
  geom_line()+
  ylim(0.5, 1.5)+
  labs(x = "Number of cells", y = "Normalized network diameter") +
  scale_color_manual(values=c("cornflowerblue", "#AE93BEFF", "#B4DAE5FF", "#F0D77BFF"))+
  theme_classic(base_size = 10)+
  guides(col='none')+
  NULL
p3_strain


# temporary plot to extract the legend from
temp_plot <- ggplot(network_df, aes(x = num_nodes, y = norm_diameter, color=strain))+
  geom_line()+
  ylim(0.5, 1.5)+
  labs(x = "Number of cells", y = "Normalized network diameter") +
  scale_color_manual(values=c("cornflowerblue", "#AE93BEFF", "#B4DAE5FF", "#F0D77BFF"))+
  theme_classic(base_size = 10)+
  theme(legend.position = "bottom",           # Position at bottom
        legend.direction = "horizontal")+      # Make it horizontal
  labs(color='Strain')+
  NULL

# Extract the legend
shared_legend <- ggpubr::get_legend(temp_plot)


# Create the plot grid without legend
plots_grid <- plot_grid(p1_strain, p2_strain, p3_strain, ncol = 3)

# Combine the plots with the shared legend
final_plot <- plot_grid(plots_grid, shared_legend, ncol = 1, rel_heights = c(1, 0.1))
final_plot

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig_5_normalized_net_diam_method_28aug2026.png',
       plot=final_plot, dpi='retina', width=9, height=3)

