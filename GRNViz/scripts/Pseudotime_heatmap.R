#!/usr/bin/env Rscript

# Load required libraries with error handling
required_packages <- c("ggplot2", "pheatmap", "plotly", "jsonlite")

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("[ERROR] Required package '%s' is not installed.\n", pkg), file = stderr())
        cat(sprintf("[ERROR] Please install it using: install.packages('%s')\n", pkg), file = stderr())
        quit(status = 1)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Logging functions
log_error <- function(message) {
    cat(sprintf("[ERROR] %s\n", message), file = stderr())
}

log_info <- function(message) {
    cat(sprintf("[INFO] %s\n", message))
}

args <- commandArgs(trailingOnly = TRUE)

# Validate command line arguments
if (length(args) < 6) {
    log_error(sprintf("Insufficient arguments. Expected 6 arguments, got %d", length(args)))
    log_error("Usage: Rscript Pseudotime_heatmap.R <expression_file> <trajectory_file> <sif_file> <output_file> <top_num> <sampled_num>")
    quit(status = 1)
}

# Function to detect file separator
detect_separator <- function(file_path) {
    tryCatch({
        first_line <- readLines(file_path, n=1)  # Read first line only
        if (grepl(",", first_line)) {
            return(",")
        } else if (grepl("\t", first_line)) {
            return("\t")
        } else if (grepl(" ", first_line)) {
            return(" ")
        } else {
            return(NULL)  # Return NULL if detection fails
        }
    }, error = function(e) {
        log_error(sprintf("Failed to detect separator in file: %s", file_path))
        log_error(sprintf("Error: %s", conditionMessage(e)))
        return(NULL)
    })
}

# Validate input files exist
input_expression_file <- args[1]
input_trajectory_file <- args[2]
sif_file <- args[3]
output_file_path <- args[4]

for (f in c(input_expression_file, input_trajectory_file, sif_file)) {
    if (!file.exists(f)) {
        log_error(sprintf("Input file not found: %s", f))
        quit(status = 1)
    }
}

# Validate numeric arguments
top_num <- suppressWarnings(as.numeric(args[5]))
sampled_num <- suppressWarnings(as.numeric(args[6]))

if (is.na(top_num) || top_num <= 0) {
    log_error(sprintf("Invalid top_num value: '%s'. Must be a positive number.", args[5]))
    quit(status = 1)
}

if (is.na(sampled_num) || sampled_num <= 0) {
    log_error(sprintf("Invalid sampled_num value: '%s'. Must be a positive number.", args[6]))
    quit(status = 1)
}

log_info(sprintf("Processing expression file: %s", input_expression_file))
log_info(sprintf("Processing trajectory file: %s", input_trajectory_file))
log_info(sprintf("Processing SIF file: %s", sif_file))
log_info(sprintf("Output file: %s", output_file_path))
log_info(sprintf("Top genes: %d", top_num))
log_info(sprintf("Sampled cells: %d", sampled_num))

# Detect separator from SIF file
sep_detected <- detect_separator(sif_file)

# Load data based on detected separator
if (!is.null(sep_detected)) {
    tryCatch({
        tenet_sif <- read.table(sif_file, header=FALSE, sep=sep_detected, stringsAsFactors=FALSE)
        log_info(sprintf("Successfully loaded SIF file with %d rows and %d columns", nrow(tenet_sif), ncol(tenet_sif)))
    }, error = function(e) {
        log_error(sprintf("Failed to read SIF file: %s", sif_file))
        log_error(sprintf("Error: %s", conditionMessage(e)))
        quit(status = 1)
    })
} else {
    log_error("Cannot detect separator. Please check the file format.")
    quit(status = 1)
}

# Function to replace empty list ([]) with empty string ("") in JSON file for title values
replace_empty_titles_in_file <- function(file_path) {
  tryCatch({
    # Read file
    json_content <- readLines(file_path, warn = FALSE)

    # Concatenate file content into a single string
    json_content <- paste(json_content, collapse = "\n")

    # Replace empty list with empty string (only for title values)
    json_content <- gsub('"title":\\[\\]', '"title":""', json_content)

    # Also replace empty lists in axis titles if present
    json_content <- gsub('"zaxis":\\{"title":\\[\\]\\}', '"zaxis":{"title":""}', json_content)
    json_content <- gsub('"xaxis":\\{"title":\\[\\]\\}', '"xaxis":{"title":""}', json_content)
    json_content <- gsub('"yaxis":\\{"title":\\[\\]\\}', '"yaxis":{"title":""}', json_content)

    # Write modified content back to file
    writeLines(json_content, file_path)
    log_info(sprintf("Successfully cleaned JSON file: %s", file_path))
  }, error = function(e) {
    log_error(sprintf("Failed to clean JSON file: %s", file_path))
    log_error(sprintf("Error: %s", conditionMessage(e)))
  })
}

pseudotime_heatmap <- function(matrix, selected_gene, gene_list,
                               pseudotime, span=0.5, use_pseudotime_origin = T, use_z_score = T,
                               order_pseudo_score=F, p_min=-0.7, p_max=0.7, p_legend=T,
                               max_min = T, out_result = F, target_average = F, target, pseudo_decrease=T) {

  tryCatch({
    #==================== Subset & order by pseudotime =================

    # Check and adjust data length
    n_samples <- length(pseudotime)
    log_info(sprintf("Processing %d samples with %d selected genes", n_samples, length(selected_gene)))

    # Validate selected genes exist in matrix
    missing_genes <- selected_gene[!selected_gene %in% colnames(matrix)]
    if (length(missing_genes) > 0) {
      log_info(sprintf("Warning: %d genes not found in expression matrix, skipping them", length(missing_genes)))
      selected_gene <- selected_gene[selected_gene %in% colnames(matrix)]
    }

    if (length(selected_gene) == 0) {
      log_error("No valid genes found in expression matrix")
      quit(status = 1)
    }

    sub_matrix <- matrix[, selected_gene, drop=FALSE]
    sub_matrix <- cbind(data.frame(pseudotime=pseudotime), sub_matrix)

    # Remove NA and Inf values
    valid_rows <- which(!is.infinite(sub_matrix$pseudotime) & !is.na(sub_matrix$pseudotime))
    if (length(valid_rows) == 0) {
      log_error("No valid pseudotime values found (all NA or Inf)")
      quit(status = 1)
    }
    sub_matrix <- sub_matrix[valid_rows, ]
    log_info(sprintf("After removing invalid pseudotime values: %d rows remaining", nrow(sub_matrix)))

    # Sort by pseudotime
    sub_matrix <- sub_matrix[order(sub_matrix$pseudotime), ]

    # Ensure all columns have the same length
    n_rows <- nrow(sub_matrix)
    if (!use_pseudotime_origin) {
      sub_matrix$pseudotime <- 1:n_rows
    }

    #=================== Target average heatmap =================
    if (target_average == T) {
      average_matrix <- data.frame(pseudotime = sub_matrix$pseudotime)

      for (j in 1:length(selected_gene)) {
        current_gene <- selected_gene[j]

        # Construct target_list
        if (ncol(gene_list) >= 3) {
          target_list <- subset(gene_list, gene_list[,1] == gsub("[.]", "-", current_gene))
          if (nrow(target_list) > 0) {
            target_list <- gsub("-", ".", target_list[,3])
          } else {
            target_list <- current_gene
          }
        } else {
          target_list <- current_gene
        }

        # Calculate average
        if (length(target_list) == 1) {
          temp_values <- sub_matrix[[target_list]]
        } else {
          valid_targets <- target_list[target_list %in% colnames(sub_matrix)]
          if (length(valid_targets) == 0) {
            temp_values <- rep(NA, nrow(sub_matrix))
          } else {
            temp_values <- rowMeans(sub_matrix[, valid_targets, drop=FALSE], na.rm = TRUE)
          }
        }

        average_matrix[[paste0(current_gene, " (", length(target_list), ")")]] <- temp_values
      }

      sub_matrix <- average_matrix
    }

    #==================== Z-score transformation =================
    if(use_z_score == T) {
      for (i in 2:dim(sub_matrix)[2]) {
        col_sd <- sd(sub_matrix[, i], na.rm = TRUE)
        col_mean <- mean(sub_matrix[, i], na.rm = TRUE)
        if (is.na(col_sd) || col_sd == 0) {
          log_info(sprintf("Warning: Column %d has zero or NA standard deviation, skipping z-score", i))
        } else {
          sub_matrix[, i] <- (sub_matrix[, i] - col_mean) / col_sd
        }
      }
    }

    #==================== Normalization with loess =================
    normalized_data <- list()
    normalized_data[[1]] <- sub_matrix[, 1]  # Store pseudotime values

    # Perform normalization for each gene
    for (i in 2:ncol(sub_matrix)) {
      x <- sub_matrix[, 1]  # pseudotime
      y <- sub_matrix[, i]  # expression values

      # Remove NA, NaN, Inf values
      valid_idx <- which(is.finite(x) & is.finite(y))
      if (length(valid_idx) > 2) {  # Need at least 3 points for loess
        tryCatch({
          lo <- loess(y[valid_idx] ~ x[valid_idx], span = span)
          normalized_data[[i]] <- predict(lo, newdata = x)
        }, error = function(e) {
          # Use original data if loess fails
          log_info(sprintf("Warning: loess failed for column %d, using original data", i))
          normalized_data[[i]] <<- y
        })
      } else {
        # Use original data if insufficient valid data
        log_info(sprintf("Warning: Insufficient valid data for column %d (%d points), using original", i, length(valid_idx)))
        normalized_data[[i]] <- y
      }
    }

    # Convert list to matrix
    normalized_matrix <- do.call(cbind, normalized_data)
    colnames(normalized_matrix) <- colnames(sub_matrix)
    normalized_sub_matrix_sorted <- t(normalized_matrix)

    #==================== Pseudo score ordering =================
    if(order_pseudo_score == T) {
      pseudo_score <- data.frame()

      for (i in 2:dim(normalized_sub_matrix_sorted)[1]) {
        pseudo_score[i, 1] <- sum(normalized_sub_matrix_sorted[1, ] * normalized_sub_matrix_sorted[i, ], na.rm = TRUE)
      }

      normalized_sub_matrix_sorted_pseudo_ordered <- cbind(pseudo_score, normalized_sub_matrix_sorted)
      rownames(normalized_sub_matrix_sorted_pseudo_ordered) <- rownames(normalized_sub_matrix_sorted)
      temp <- normalized_sub_matrix_sorted_pseudo_ordered[-1, ]
      temp <- temp[order(temp[, 1], decreasing = pseudo_decrease), ]
      normalized_sub_matrix_sorted_pseudo_ordered <- rbind(normalized_sub_matrix_sorted_pseudo_ordered[1, ], temp)
      final_matrix <- normalized_sub_matrix_sorted_pseudo_ordered
      pheat_start <- 2
    } else {
      final_matrix <- normalized_sub_matrix_sorted
      pheat_start <- 1
    }

    #==================== Create pheatmap =================
    tryCatch({
      if(max_min == T) {
        s1 <- pheatmap(final_matrix[2:dim(final_matrix)[1], pheat_start:dim(final_matrix)[2]],
                       breaks = seq(p_min, p_max, length.out=100), show_colnames=F, cluster_rows = F, cluster_cols = F, legend = p_legend)
      } else {
        s1 <- pheatmap(final_matrix[2:dim(final_matrix)[1], pheat_start:dim(final_matrix)[2]],
                       show_colnames=F, cluster_rows = F, cluster_cols = F, legend = p_legend)
      }
      log_info("Successfully created pheatmap")
    }, error = function(e) {
      log_error(sprintf("Failed to create pheatmap: %s", conditionMessage(e)))
    })

    #==================== Convert to Plotly =================
    # Ensure 'z' is a matrix
    z_matrix <- as.matrix(final_matrix[2:dim(final_matrix)[1], pheat_start:dim(final_matrix)[2]])

    fig <- plot_ly(
      z = z_matrix,
      x = colnames(final_matrix)[pheat_start:dim(final_matrix)[2]],
      y = rownames(final_matrix)[2:dim(final_matrix)[1]],
      type = "heatmap",
      colorscale = "Viridis"
    )

    # Set proper title format
    layout <- list(
      title = "Heatmap",
      xaxis = list(
        title = "Pseudotime",
        titlefont = list(size = 14)
      ),
      yaxis = list(
        title = "Genes",
        titlefont = list(size = 14)
      ),
      colorbar = list(
        title = list(
          text = "Expression Level",
          font = list(size = 12)
        )
      )
    )

    fig <- fig %>% layout(layout)

    # Convert to JSON directly
    fig_json <- plotly_json(fig, jsonedit = FALSE, pretty = TRUE)

    # Save JSON file
    tryCatch({
      write(fig_json, output_file_path)
      log_info(sprintf("Successfully saved heatmap to: %s", output_file_path))
    }, error = function(e) {
      log_error(sprintf("Failed to save output file: %s", output_file_path))
      log_error(sprintf("Error: %s", conditionMessage(e)))
      quit(status = 1)
    })

  }, error = function(e) {
    log_error("Failed to create pseudotime heatmap")
    log_error(sprintf("Error: %s", conditionMessage(e)))
    quit(status = 1)
  })
}

count_degree <- function(data, decreasing=T, degree="out") {
    # Check data format (number of columns)
    num_cols <- ncol(data)

    if (num_cols == 3) {  # 3-column format (source, weight, target)
        if (degree == "out") {
            att <- as.data.frame(table(data[,1]))[
                order(as.data.frame(table(data[,1]),)[,2],
                decreasing = decreasing),]
        }
        if (degree == "in") {
            att <- as.data.frame(table(data[,3]))[
                order(as.data.frame(table(data[,3]),)[,2],
                decreasing = decreasing),]
        }
    } else if (num_cols == 2) {  # 2-column format (gene, weight)
        # Ignore degree parameter for 2-column format
        att <- data.frame(
            Var1 = data[,1],
            Freq = data[,2]
        )
        # Sort by weight (Freq)
        att <- att[order(att$Freq, decreasing=decreasing),]
    } else {
        log_error(sprintf("Input data must have 2 or 3 columns, got %d columns", num_cols))
        quit(status = 1)
    }

    return(att)
}

# Load input files with error handling
tryCatch({
    expression_data <- read.csv(input_expression_file)
    log_info(sprintf("Successfully loaded expression file with %d rows and %d columns", nrow(expression_data), ncol(expression_data)))
}, error = function(e) {
    log_error(sprintf("Failed to read expression file: %s", input_expression_file))
    log_error(sprintf("Error: %s", conditionMessage(e)))
    quit(status = 1)
})

tryCatch({
    trajectory <- read.table(input_trajectory_file, quote="\"", comment.char="")
    log_info(sprintf("Successfully loaded trajectory file with %d rows", nrow(trajectory)))
}, error = function(e) {
    log_error(sprintf("Failed to read trajectory file: %s", input_trajectory_file))
    log_error(sprintf("Error: %s", conditionMessage(e)))
    quit(status = 1)
})

# Calculate degree from SIF data
count_TENET <- count_degree(tenet_sif, degree = "out")
selected_gene <- count_TENET$Var1
log_info(sprintf("Found %d genes in SIF file", length(selected_gene)))

# Select only numeric columns
numeric_columns <- sapply(expression_data, is.numeric)
numeric_expression_data <- expression_data[, numeric_columns, drop=FALSE]

if (ncol(numeric_expression_data) == 0) {
    log_error("No numeric columns found in expression data")
    quit(status = 1)
}
log_info(sprintf("Selected %d numeric columns from expression data", ncol(numeric_expression_data)))

# Calculate top N genes
actual_top_num <- min(top_num, nrow(numeric_expression_data))
if (actual_top_num < top_num) {
    log_info(sprintf("Requested top %d genes, but only %d available", top_num, actual_top_num))
}
top_genes <- head(order(rowSums(numeric_expression_data, na.rm = TRUE), decreasing = TRUE), actual_top_num)

# Sample top N genes
sampled_expression_data <- numeric_expression_data[top_genes, , drop=FALSE]

# Sample top N cells from trajectory
actual_sampled_num <- min(sampled_num, nrow(trajectory))
if (actual_sampled_num < sampled_num) {
    log_info(sprintf("Requested %d samples, but only %d available", sampled_num, actual_sampled_num))
}
sampled_trajectory <- trajectory[1:actual_sampled_num, ]

# Select only genes that exist in expression data column names
valid_genes <- selected_gene[selected_gene %in% colnames(sampled_expression_data)]

if (length(valid_genes) == 0) {
    log_error("No matching genes found between SIF file and expression data")
    log_error("Please check that gene names in SIF file match column names in expression file")
    quit(status = 1)
}
log_info(sprintf("Found %d valid genes matching between SIF and expression data", length(valid_genes)))

# Call the main function
log_info("Starting pseudotime heatmap generation...")
pseudotime_heatmap(matrix = sampled_expression_data, gene_list = tenet_sif,
                   selected_gene = as.character(valid_genes), pseudotime = sampled_trajectory,
                   span=0.7, use_pseudotime_origin = T, use_z_score = T,
                   order_pseudo_score=T, max_min = T, p_min=-0.7, p_max=0.7, p_legend=T,
                   out_result =F, target_average = T, pseudo_decrease=T)

log_info("Pseudotime heatmap generation completed successfully")
