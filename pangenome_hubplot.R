
# 1. Install and Load necessary packages
# (Run install.packages if needed)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

# 2. Load your data
# Assuming your data frame in R is already loaded and named 'sv_cover'
# If you need to load it from the RDS file again:
sv <- readRDS("/Users/yuhenghuang/Documents/Postdoc_UCI/Pangenome_TE/Pangenome_TE_Results/seq_sv_big_on_sv_cover.rds")
sv_cover <- sv
#test <- readRDS("/Users/yuhenghuang/Documents/Postdoc_UCI/Pangenome_TE/Results/GCA_003402055.1_5_5_full.rds")
#sv_cover <- sv |> 
#  slice(1:50000)
  #  filter(cover == 1)
#write_rds(sv_cover, "/Users/yuhenghuang/Documents/Postdoc_UCI/Pangenome_TE/Results/sv_cover.rds")

#sv_cover <- rename(sv_cover, C1 = C1.x, C8 = C8.x, V8 = V8.x, p1 = p1.x, p8 = p8.x)
# 3. Define Column Names based on your head() output
col_sv1 <- "V1"
col_sv2 <- "V8"
col_p1 <- "p1"
col_p8 <- "p8"
col_sv1_attribute <- "C1" # For potential coloring info later
col_sv2_attribute <- "C8" # For potential coloring info later

# Define the threshold
p_threshold <- 0.85
min_sequences_for_node <- 3 # New: Minimum sequences for a node to be plotted

# --- 1. Identify 'Mutual Nesting' Edges (Strong Reciprocal Links) ---
# These edges define the groups/components
#mutual_edges_df <- sv_cover %>%
#  filter(.data[[col_p1]] > p_threshold & .data[[col_p8]] > p_threshold) %>%
#  select(from = .data[[col_sv1]], to = .data[[col_sv2]]) %>%
#  filter(from != to) # Ignore self-comparisons if they exist

mutual_edges_df <- sv_cover |> 
  filter(p1 > p_threshold & p8 > p_threshold) |> 
  select(from = V1, to = V8)  |> 
  filter(from != to)
print(paste("Found", nrow(mutual_edges_df), "mutual nesting links (p1 >", p_threshold, "& p8 >", p_threshold, ")"))


# --- 2. Identify All Unique Sequences ---
# We need all sequences present in the data, even if they don't form mutual pairs
all_sv_ids_in_data <- unique(c(sv_cover[[col_sv1]], sv_cover[[col_sv2]]))
node_list_all_seqs <- tibble(name = all_sv_ids_in_data)

# --- 3. Build Graph Based ONLY on Mutual Nesting & Find Components ---
# This graph helps group sequences based on the strong condition
mutual_graph <- tbl_graph(nodes = node_list_all_seqs, edges = mutual_edges_df, directed = FALSE)

# Assign a component ID to each sequence based on the mutual graph structure
# Sequences not part of any mutual pair will form components of size 1
node_component_mapping <- mutual_graph %>%
  activate(nodes) %>%
  mutate(component_id = group_components()) %>%
  as_tibble() %>%
  select(name, component_id) # Map: Original Sequence Name -> Component ID

# --- 4. Summarize Components to Define Final Nodes ---
# Get attribute information (like C1/C8) for potential coloring
attributes_1 <- sv_cover %>% select(name = .data[[col_sv1]], node_attribute = .data[[col_sv1_attribute]])
attributes_2 <- sv_cover %>% select(name = .data[[col_sv2]], node_attribute = .data[[col_sv2_attribute]])
all_attributes <- bind_rows(attributes_1, attributes_2) %>% distinct(name, .keep_all = TRUE)

# Combine component mapping with attributes
component_nodes_data <- node_component_mapping %>%
  left_join(all_attributes, by = "name") %>%
  mutate(node_attribute = factor(replace_na(node_attribute, "Unknown")))

# Create the summary table for the final graph nodes
component_summary_all <- component_nodes_data %>%
  group_by(component_id) %>%
  summarise(
    n_sequences = n(), # Node size aesthetic
    # For Color: Use component ID as a factor. Could also use most frequent attribute etc.
    sequences_in_component = list(name) # Keep track of original sequences
  ) %>%
  ungroup()

 # Filter Components (Nodes) by Size (before TE family assignment)
component_summary_size_filtered <- component_summary_all %>%
  filter(n_sequences >= min_sequences_for_node)

print(paste("Components after size filtering (n_sequences >", min_sequences_for_node, "):", nrow(component_summary_size_filtered)))

# --- Check TE Family Consistency Within Components ---
# Load TE information
sv_te <- readRDS("/Users/yuhenghuang/Documents/Postdoc_UCI/Pangenome_TE/Pangenome_TE_Results/seq_sv_big_on_set_cover.rds")

# Print structure of sv_te to understand the data
print("\n=== TE Data Structure ===")
print(paste("Columns in sv_te:", paste(names(sv_te), collapse = ", ")))
print(paste("Number of rows in sv_te:", nrow(sv_te)))
print("First few rows of sv_te:")
print(head(sv_te, 3))

# Get TE family information for sequences
# User specified: V8 in sv_te contains TE family information
# V1 in sv_te should contain sequence IDs (matching sequences in sv_cover)
# Create a mapping of sequence ID to TE family
te_mapping <- sv_te %>% 
  select(sequence_id = V1, te_family = V8) %>%
  filter(!is.na(sequence_id), !is.na(te_family)) %>%
  distinct(sequence_id, .keep_all = TRUE)

print(paste("\nNumber of sequences with TE information:", nrow(te_mapping)))
print(paste("Number of unique TE families:", n_distinct(te_mapping$te_family)))

# --- Load and Process TE Taxonomy Data ---
te_taxonomy_file <- "/Users/yuhenghuang/Documents/Postdoc_UCI/names_D_mel_TEs"
te_taxonomy_lines <- readLines(te_taxonomy_file)

# Parse the lines to extract family and superfamily
# Format: >"family_id"#"superfamily_id"
# Remove the leading ">"
te_taxonomy_clean <- str_remove(te_taxonomy_lines, "^>")

# Split by "#"
te_taxonomy_df <- tibble(raw_line = te_taxonomy_clean) %>%
  separate(raw_line, into = c("te_family", "te_superfamily"), sep = "#", extra = "merge", fill = "right") %>%
  # Ensure distinct mapping
  distinct(te_family, .keep_all = TRUE)

print("\n=== TE Taxonomy Data ===")
print(head(te_taxonomy_df))

# Join TE information with component mapping
component_with_te <- node_component_mapping %>%
  left_join(te_mapping, by = c("name" = "sequence_id")) %>%
  # Add superfamily info
  left_join(te_taxonomy_df, by = "te_family")

# Assign TE families to components based on dominant family (>85% threshold)
# Now also including superfamily info
te_family_assignment <- component_with_te %>%
  filter(!is.na(te_family)) %>%  # Only check sequences with TE information
  group_by(component_id, te_family, te_superfamily) %>%
  summarise(n_in_family = n(), .groups = "drop") %>%
  group_by(component_id) %>%
  mutate(
    total_sequences = sum(n_in_family),
    family_proportion = n_in_family / total_sequences
  ) %>%
  # Find the dominant family (highest proportion)
  slice_max(family_proportion, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Assign TE family if dominant family has >95% of sequences, otherwise mark for removal
  mutate(
    assigned_te_family = if_else(family_proportion >= 0.95, te_family, NA_character_),
    assigned_te_superfamily = if_else(family_proportion >= 0.95, te_superfamily, NA_character_),
    has_dominant_family = family_proportion >= 0.95
  ) %>%
  select(component_id, assigned_te_family, assigned_te_superfamily, family_proportion, total_sequences, has_dominant_family)

# Summary statistics
total_components_with_te <- nrow(te_family_assignment)
components_with_dominant_family <- sum(te_family_assignment$has_dominant_family, na.rm = TRUE)
components_to_remove <- total_components_with_te - components_with_dominant_family

print("\n=== TE Family Assignment Analysis ===")
print(paste("Total components with TE information:", total_components_with_te))
print(paste("Components with dominant TE family (>95%):", components_with_dominant_family))
print(paste("Components to be removed (no dominant family):", components_to_remove))
print(paste("Percentage with dominant family:", round(100 * components_with_dominant_family / total_components_with_te, 2), "%"))

# Show examples of components being removed
if (components_to_remove > 0) {
  print("\nExamples of components being removed (no dominant TE family >85%):")
  print(te_family_assignment %>% 
    filter(!has_dominant_family) %>% 
    arrange(desc(total_sequences)) %>% 
    head(10))
}

# --- Integrate TE Family Assignment and Filter Components ---
# Join TE family assignment to component summary
component_summary_with_te <- component_summary_size_filtered %>%
  left_join(te_family_assignment %>% select(component_id, assigned_te_family, assigned_te_superfamily, has_dominant_family), 
            by = "component_id") %>%
  # Remove components without dominant TE family
  filter(has_dominant_family == TRUE | is.na(has_dominant_family)) %>%
  # Remove components that don't have TE information at all (optional - you may want to keep them)
  # For now, we'll keep components without TE info but mark them
  mutate(
    te_family = factor(replace_na(assigned_te_family, "Unknown")),
    te_superfamily = factor(replace_na(assigned_te_superfamily, "Unknown")),
    component_id_factor = factor(component_id)
  ) %>%
  select(-has_dominant_family, -assigned_te_family, -assigned_te_superfamily)

# Final filtered component summary (after size and TE family filtering)
filtered_component_summary <- component_summary_with_te

# Get the IDs of the components that passed all filters
valid_component_ids <- filtered_component_summary$component_id

print("\n=== Final Component Summary (after size and TE family filtering) ===")
print(paste("Retained", nrow(filtered_component_summary), "nodes after all filtering"))
print(paste("Number of unique TE families:", n_distinct(filtered_component_summary$te_family)))
print(paste("Number of unique TE superfamilies:", n_distinct(filtered_component_summary$te_superfamily)))
print("TE family distribution:")
print(table(filtered_component_summary$te_family))
print("TE superfamily distribution:")
print(table(filtered_component_summary$te_superfamily))
print("\nFiltered Component Summary (Nodes for Final Plot):")
print(head(filtered_component_summary))

# --- 5. Identify 'One-Way' Similarity Edges Between Original Sequences ---
oneway_edges_df <- sv_cover %>%
  filter(
    (.data[[col_p1]] > p_threshold & .data[[col_p8]] <= p_threshold) | # p1 high, p8 low
      (.data[[col_p1]] <= p_threshold & .data[[col_p8]] > p_threshold)   # p1 low, p8 high
  ) %>%
  # Define direction: if p1>T, edge goes from sv1->sv2; if p8>T, edge goes from sv2->sv1
  mutate(
    from_seq = if_else(.data[[col_p1]] > p_threshold, .data[[col_sv1]], .data[[col_sv2]]),
    to_seq = if_else(.data[[col_p1]] > p_threshold, .data[[col_sv2]], .data[[col_sv1]])
  ) %>%
  select(from_seq, to_seq) %>%
  filter(from_seq != to_seq) # Ignore self-loops

print(paste("Found", nrow(oneway_edges_df), "potential one-way links (p1 > T or p8 > T, but not both)."))

# --- 6. Map One-Way Sequence Edges to Component IDs ---
seq_to_comp_lookup <- setNames(node_component_mapping$component_id, node_component_mapping$name)

component_edges_all_df <- oneway_edges_df %>% # Renamed to indicate it's 'all' before filtering
  mutate(
    from_comp = seq_to_comp_lookup[from_seq],
    to_comp = seq_to_comp_lookup[to_seq]
  ) %>%
  filter(!is.na(from_comp) & !is.na(to_comp)) %>% # Ensure sequences were found in mapping
  filter(from_comp != to_comp) %>%
  distinct(from_comp, to_comp)

print(paste("Mapped to", nrow(component_edges_all_df), "directed edges between component nodes."))
print("Component Edges (before TE family filtering):")
print(head(component_edges_all_df))
# Note: Edges will be filtered after TE family filtering (see below)

# Graph will be created after TE family filtering (see below)


# Re-filter edges based on final valid component IDs (after TE family filtering)
filtered_component_edges_df <- component_edges_all_df %>%
  filter(from_comp %in% valid_component_ids & to_comp %in% valid_component_ids) %>%
  select(from = from_comp, to = to_comp) # Rename for graph creation

print("\n=== Final Edge Summary ===")
print(paste("Edges after all filtering:", nrow(filtered_component_edges_df)))
print("Component Edges (Edges for Final Plot):")
print(head(filtered_component_edges_df))

# --- 7. Create the Final Aggregate Graph ---
# Check if we have nodes before creating graph
if (nrow(filtered_component_summary) == 0) {
  stop("No nodes remaining after filtering. Cannot create graph.")
}

# Ensure all edge endpoints exist in the node set
valid_node_ids <- filtered_component_summary$component_id

# Debug: Check edge endpoints
if (nrow(filtered_component_edges_df) > 0) {
  edges_from_valid <- filtered_component_edges_df$from %in% valid_node_ids
  edges_to_valid <- filtered_component_edges_df$to %in% valid_node_ids
  print(paste("Edges with valid 'from' nodes:", sum(edges_from_valid), "out of", nrow(filtered_component_edges_df)))
  print(paste("Edges with valid 'to' nodes:", sum(edges_to_valid), "out of", nrow(filtered_component_edges_df)))
  
  # Show problematic edges if any
  if (any(!edges_from_valid) || any(!edges_to_valid)) {
    print("Problematic edges (referencing non-existent nodes):")
    print(filtered_component_edges_df[!(edges_from_valid & edges_to_valid), ])
  }
}

# Filter edges to only include those where both endpoints exist
filtered_component_edges_df <- filtered_component_edges_df %>%
  filter(from %in% valid_node_ids & to %in% valid_node_ids)

print(paste("Final edges after ensuring node matches:", nrow(filtered_component_edges_df)))
print(paste("Number of nodes:", nrow(filtered_component_summary)))

# Ensure component_id types match between nodes and edges
# Get the original type of component_id
component_id_type <- class(filtered_component_summary$component_id)[1]

if (nrow(filtered_component_edges_df) > 0) {
  # Convert edge endpoints to match component_id type
  if (component_id_type == "factor") {
    filtered_component_edges_df$from <- factor(filtered_component_edges_df$from, levels = levels(filtered_component_summary$component_id))
    filtered_component_edges_df$to <- factor(filtered_component_edges_df$to, levels = levels(filtered_component_summary$component_id))
  } else if (component_id_type == "numeric" || component_id_type == "integer") {
    filtered_component_edges_df$from <- as.numeric(filtered_component_edges_df$from)
    filtered_component_edges_df$to <- as.numeric(filtered_component_edges_df$to)
  } else {
    filtered_component_edges_df$from <- as.character(filtered_component_edges_df$from)
    filtered_component_edges_df$to <- as.character(filtered_component_edges_df$to)
  }
  
  # Final check - remove any edges that still don't match
  filtered_component_edges_df <- filtered_component_edges_df %>%
    filter(from %in% valid_node_ids & to %in% valid_node_ids) %>%
    filter(!is.na(from) & !is.na(to))
}

# Create graph - ensure all nodes in edges are in nodes dataframe
if (nrow(filtered_component_edges_df) == 0) {
  print("Warning: No edges remaining. Creating graph with nodes only.")
  final_aggregate_graph_filtered <- tbl_graph(
    nodes = filtered_component_summary,
    edges = NULL,
    directed = TRUE,
    node_key = "component_id"
  )
} else {
  # Get all unique node IDs from edges
  edge_node_ids <- unique(c(filtered_component_edges_df$from, filtered_component_edges_df$to))
  
  # Check if all edge nodes are in the nodes dataframe
  missing_nodes <- edge_node_ids[!edge_node_ids %in% valid_node_ids]
  if (length(missing_nodes) > 0) {
    print(paste("Warning: Found", length(missing_nodes), "nodes in edges that are not in nodes dataframe. Removing those edges."))
    filtered_component_edges_df <- filtered_component_edges_df %>%
      filter(from %in% valid_node_ids & to %in% valid_node_ids)
  }
  
  # Double check - ensure no edges reference missing nodes
  if (nrow(filtered_component_edges_df) > 0) {
    all_edge_nodes <- unique(c(filtered_component_edges_df$from, filtered_component_edges_df$to))
    if (!all(all_edge_nodes %in% valid_node_ids)) {
      stop("Error: Some edges still reference nodes not in the nodes dataframe. This should not happen.")
    }
  }
  
  # Ensure nodes dataframe has component_id as first column or as a clear identifier
  # Make sure component_id is clean (no NAs, proper type)
  filtered_component_summary <- filtered_component_summary %>%
    filter(!is.na(component_id)) %>%
    distinct(component_id, .keep_all = TRUE)
  
  # Re-filter edges one more time to match cleaned nodes
  valid_node_ids <- filtered_component_summary$component_id
  filtered_component_edges_df <- filtered_component_edges_df %>%
    filter(from %in% valid_node_ids & to %in% valid_node_ids) %>%
    filter(!is.na(from) & !is.na(to))
  
  print(paste("Final node count:", nrow(filtered_component_summary)))
  print(paste("Final edge count:", nrow(filtered_component_edges_df)))
  
  # Create graph - ensure nodes are properly set up
  if (nrow(filtered_component_edges_df) > 0) {
    # Get all nodes from edges
    all_nodes_in_edges <- unique(c(filtered_component_edges_df$from, filtered_component_edges_df$to))
    print(paste("Unique nodes in edges:", length(all_nodes_in_edges)))
    print(paste("Nodes in nodes dataframe:", length(valid_node_ids)))
    
    # Ensure all edge nodes are in nodes
    if (!all(all_nodes_in_edges %in% valid_node_ids)) {
      stop("Critical error: Not all edge nodes are in nodes dataframe")
    }
  }
  
  # Create the graph - include ALL nodes, even isolated ones
  # The issue is that tbl_graph with node_key can have problems when edges don't reference all nodes
  # Solution: Create graph with all nodes explicitly specified as vertices
  
  if (nrow(filtered_component_edges_df) > 0) {
    # Prepare vertices dataframe with ALL nodes
    # graph_from_data_frame requires 'name' column for vertices
    vertices_df <- filtered_component_summary %>%
      mutate(name = as.character(component_id)) %>%
      select(name, everything())
    
    # Prepare edges dataframe - use character names that match vertices
    edges_df <- filtered_component_edges_df %>%
      mutate(
        from = as.character(from),
        to = as.character(to)
      ) %>%
      filter(!is.na(from) & !is.na(to))
    
    # Create graph using igraph - this will include ALL vertices we specify
    # even if they're not in edges
    g <- graph_from_data_frame(
      edges_df,
      directed = TRUE,
      vertices = vertices_df  # This ensures ALL nodes are included
    )
    
    # Convert to tbl_graph
    final_aggregate_graph_filtered <- as_tbl_graph(g) %>%
      activate(nodes) %>%
      # Ensure component_id is the correct type (it should be preserved from vertices_df)
      mutate(component_id = as.character(component_id))
      
  } else {
    # No edges - create graph from nodes only
    final_aggregate_graph_filtered <- tbl_graph(
      nodes = filtered_component_summary,
      edges = NULL,
      directed = TRUE,
      node_key = "component_id"
    )
  }
}

print("Final Filtered Aggregate Graph Summary:")
print(final_aggregate_graph_filtered)


# --- 8. Plot the Final Aggregate Graph ---
set.seed(99) # Seed for layout reproducibility

# --- Generate Colors for TE Families ---
num_te_families <- n_distinct(filtered_component_summary$te_family)
family_colors <- character(0) # Initialize
if (num_te_families > 0) {
  te_family_levels <- levels(filtered_component_summary$te_family)
  # Generate different colors for different families (no limit)
  family_colors <- setNames(
    scales::hue_pal()(num_te_families),
    te_family_levels
  )
}

# --- Generate Colors for TE Superfamilies ---
# Define User-Specified Palettes
dna_palette <- c("#56b4e9", "#0072b2", "#5d95c7", "#a2d2ff")
line_palette <- c("#009e73", "#b3c074", "#b0d146", "#9adab2ff")
ltr_palette <- c("#d55e00", "#e69f00",  "#cc79a7")
unknown_color <- "#999999"

te_superfamily_levels <- levels(filtered_component_summary$te_superfamily)
superfamily_colors <- character(0) # Initialize

if (length(te_superfamily_levels) > 0) {
  superfamily_colors <- character(length(te_superfamily_levels))
  names(superfamily_colors) <- te_superfamily_levels
  
  # Indices for cycling through palettes
  dna_idx <- 1
  line_idx <- 1
  ltr_idx <- 1
  
  for (sf in te_superfamily_levels) {
    # Ensure safe string for checking
    sf_upper <- toupper(sf)
    
    if (is.na(sf) || sf == "Unknown") {
      superfamily_colors[sf] <- unknown_color
    } else if (grepl("DNA", sf_upper) || grepl("RC", sf_upper) || grepl("HELITRON", sf_upper)) {
      # DNA Group
      superfamily_colors[sf] <- dna_palette[((dna_idx - 1) %% length(dna_palette)) + 1]
      dna_idx <- dna_idx + 1
    } else if (grepl("LINE", sf_upper)) {
      # LINE Group
      superfamily_colors[sf] <- line_palette[((line_idx - 1) %% length(line_palette)) + 1]
      line_idx <- line_idx + 1
    } else if (grepl("LTR", sf_upper)) {
      # LTR Group
      superfamily_colors[sf] <- ltr_palette[((ltr_idx - 1) %% length(ltr_palette)) + 1]
      ltr_idx <- ltr_idx + 1
    } else {
      # Other/Unknown
      superfamily_colors[sf] <- unknown_color
    }
  }
}

# Create custom layout that prevents overlap of large nodes while keeping them centered
# Use 'kk' layout to maintain the original centered appearance
base_layout <- create_layout(final_aggregate_graph_filtered, layout = 'kk')

# Calculate layout coordinate scale
x_range <- diff(range(base_layout$x))
y_range <- diff(range(base_layout$y))
layout_scale <- max(x_range, y_range, 1)

# Calculate relative node sizes (proportional to sqrt of n_sequences for radius)
max_n_sequences <- max(base_layout$n_sequences)
node_size_ratio <- sqrt(base_layout$n_sequences / max_n_sequences)

# Convert node sizes to layout coordinate units
# Use a smaller radius to keep nodes closer together, but ensure no overlap
max_node_radius_in_layout <- layout_scale * 0.08
node_radii_layout <- node_size_ratio * max_node_radius_in_layout

# Find the largest node and ensure it doesn't overlap with other large nodes
node_order <- order(base_layout$n_sequences, decreasing = TRUE)
if (length(node_order) >= 2) {
  idx_largest <- node_order[1]  # Largest node
  
  # Check the largest node against other large nodes (top 5 largest)
  n_large_nodes <- min(5, length(node_order))
  
  for (j in 2:n_large_nodes) {
    idx_other <- node_order[j]
    
    # Calculate current distance
    dx <- base_layout$x[idx_other] - base_layout$x[idx_largest]
    dy <- base_layout$y[idx_other] - base_layout$y[idx_largest]
    current_dist <- sqrt(dx^2 + dy^2)
    
      # Minimum distance needed (sum of radii + minimal padding to prevent overlap)
      min_dist <- (node_radii_layout[idx_largest] + node_radii_layout[idx_other]) * 1.05  # 5% padding - just enough to prevent overlap
    
    # If overlapping or too close, push them apart
    if (current_dist < min_dist) {
      # If nodes are at same position, use horizontal separation to keep them centered
      if (current_dist < 1e-10) {
        dx <- 1
        dy <- 0
      } else {
        # Normalize direction
        dx <- dx / current_dist
        dy <- dy / current_dist
      }
      
      # Calculate movement needed
      move_distance <- (min_dist - current_dist) / 2
      
      # Move nodes apart symmetrically to keep them centered
      # Move the other node more to preserve the largest node's central position
      base_layout$x[idx_other] <- base_layout$x[idx_other] + dx * move_distance * 1.5
      base_layout$y[idx_other] <- base_layout$y[idx_other] + dy * move_distance * 1.5
      base_layout$x[idx_largest] <- base_layout$x[idx_largest] - dx * move_distance * 0.5
      base_layout$y[idx_largest] <- base_layout$y[idx_largest] - dy * move_distance * 0.5
    }
  }
}

# --- 9. Plot and Save (Two Plots: Family and Superfamily) ---

# Function to create plot
create_te_plot <- function(graph, layout, color_var, color_palette, legend_title) {
  ggraph(graph, layout = layout) +
    geom_edge_link(
      arrow = arrow(length = unit(2.5, 'mm')),
      start_cap = circle(4, 'mm'),
      end_cap = circle(4, 'mm'),
      alpha = 0.4, color = '#808080'
    ) +
    geom_node_point(
      aes(size = n_sequences, color = !!sym(color_var)),
      alpha = 0.9
    ) +
    scale_size_area(max_size = 20, name = "Sequences per Group") +
    # Only apply scale_color_manual if there are colors to apply
    (if(length(color_palette) > 0) scale_color_manual(values = color_palette, name = legend_title) else NULL) +
    theme_graph(base_family = 'sans', background = "white") +
    theme(
      legend.key = element_rect(fill = "white", colour = NA),
      legend.title = element_text(size = 12, face = "bold"), 
      legend.text = element_text(size = 10),
      legend.position = "right"
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 5), title = legend_title, order = 1),
      size = guide_legend(title = "Sequences per Group", order = 2)
    )
}

# 1. Plot colored by TE Family
plot_family <- create_te_plot(final_aggregate_graph_filtered, base_layout, "te_family", family_colors, "TE Family")

# Display and Save Family Plot
if (nrow(filtered_component_summary) > 0) {
  print("Saving TE Family plot...")
  ggsave(paste0("/Users/yuhenghuang/Documents/Postdoc_UCI/Pangenome_TE/Pangenome_TE_Results/TE_100bp_aggregate_filtered_gt", min_sequences_for_node, "_family_95percent.pdf"),
         plot = plot_family, width = 15, height = 10, dpi = 500)
} else {
  print(paste("No nodes to plot after filtering for n_sequences >", min_sequences_for_node))
}

# 2. Plot colored by TE Superfamily
plot_superfamily <- create_te_plot(final_aggregate_graph_filtered, base_layout, "te_superfamily", superfamily_colors, "TE Superfamily")

# Display and Save Superfamily Plot
if (nrow(filtered_component_summary) > 0) {
  print("Saving TE Superfamily plot...")
  ggsave(paste0("/Users/yuhenghuang/Documents/Postdoc_UCI/Pangenome_TE/Pangenome_TE_Results/TE_100bp_aggregate_filtered_gt", min_sequences_for_node, "_superfamily_95percent.pdf"),
         plot = plot_superfamily, width = 12, height = 10, dpi = 500)
}
