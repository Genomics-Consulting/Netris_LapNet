# PurIST ----
#' PurIST
#'
#' Purity Independent Subtyping of Tumors
#'
#' @param dat \code{vector} A single sample gene expression profile 
#'     with gene symbols in names
#' 
#' @return probability
#'
#' @examples
#' probabibity <- PurIST(dat)
#'
#' @author Heewon Seo, \email{Heewon.Seo@uhnresearch.ca}
#'
#' @export
#' @importFrom gtools inv.logit
PurIST_binary <- function(dat)
{
  biomarkers = c("GPR87", "REG4", 
                 "KRT6A", "ANXA10",
                 "BCAR3", "GATA6",
                 "PTGES", "CLDN18",
                 "ITGA3", "LGALS4",
                 "C16orf74", "DDC", 
                 "S100A2", "SLC40A1", 
                 "KRT5", "CLRN3")
  
  has_genes = is.element(biomarkers, names(dat))
  if (!all(has_genes)) {
    stop(sprintf("Please provide the following genes\n%s",
                 paste(biomarkers, collapse = ", ")))
  }
  
  coefficient <- c(1.994, 2.031, 1.618, 0.922, 1.059, 0.929, 2.505, 0.485)
  intercept <- -6.815
  
  if (dat["GPR87"] > dat["REG4"])     { penal_1 <- 1 } else { penal_1 <- 0 }
  if (dat["KRT6A"] > dat["ANXA10"])   { penal_2 <- 1 } else { penal_2 <- 0 }
  if (dat["BCAR3"] > dat["GATA6"])    { penal_3 <- 1 } else { penal_3 <- 0 }
  if (dat["PTGES"] > dat["CLDN18"])   { penal_4 <- 1 } else { penal_4 <- 0 }
  if (dat["ITGA3"] > dat["LGALS4"])   { penal_5 <- 1 } else { penal_5 <- 0 }
  if (dat["C16orf74"] > dat["DDC"])   { penal_6 <- 1 } else { penal_6 <- 0 }
  if (dat["S100A2"] > dat["SLC40A1"]) { penal_7 <- 1 } else { penal_7 <- 0 }
  if (dat["KRT5"] > dat["CLRN3"])     { penal_8 <- 1 } else { penal_8 <- 0 }
  
  penal <- c(penal_1, penal_2, penal_3, penal_4, penal_5, penal_6, penal_7, penal_8)
  
  return( gtools::inv.logit(intercept + sum(coefficient * penal)) )
}

# PurIST Classifier (Log Ratio Model)
# Based on: Rashid et al., Clin Cancer Res 2020, derived from Moffitt subtypes.
# This version correctly implements the weighted log2 ratio of the gene pairs
# and includes the necessary intercept term.

#' PurIST
#'
#' Purity Independent Subtyping of Tumors (PurIST) classifier.
#' Calculates the probability of a sample belonging to the Basal-like subtype
#' based on the weighted log2 ratio of 8 gene pairs.
#'
#' @param dat A numeric vector representing a single sample's gene expression
#'            profile, with HUGO gene symbols in the names (must be normalized, e.g., TPM).
#' @param epsilon A small pseudocount added to all expression values to prevent errors
#'                when taking the log of zero or dividing by zero. Standard value is 1.
#' @return A numeric value (0 to 1) representing the probability of being the Basal-like subtype.
#'
#' @export
#' @importFrom gtools inv.logit
PurIST_LogRatio <- function(dat, epsilon = 1) {
  
  # 1. Define Classifier Parameters (from PurIST Coefficients table)
  # The intercept and coefficients derived from the logistic regression model
  INTERCEPT_B0 <- -6.815 
  
  purist_model <- data.frame(
    GeneA = c("GPR87", "KRT6A", "BCAR3", "PTGES", "ITGA3", "C16orf74", "S100A2", "KRT5"),
    GeneB = c("REG4", "ANXA10", "GATA6", "CLDN18", "LGALS4", "DDC", "SLC40A1", "CLRN3"),
    # Beta coefficients provided in your PurIST Coefficients table
    Beta = c(1.994, 2.031, 1.618, 0.922, 1.059, 0.929, 2.505, 0.485)
  )
  
  # All 16 gene markers needed for the calculation
  biomarkers <- unique(c(purist_model$GeneA, purist_model$GeneB))
  
  # 2. Input Validation
  # Check if all 16 required genes are present in the input data
  has_genes <- is.element(biomarkers, names(dat))
  if (!all(has_genes)) {
    stop(sprintf("Missing required PurIST genes: %s", 
                 paste(biomarkers[!has_genes], collapse = ", ")))
  }
  
  # Check for NA/Inf values in the expression data
  if (any(is.na(dat[biomarkers])) || any(dat[biomarkers] < 0)) {
    stop("Expression data must be non-negative and contain no NA values for PurIST genes.")
  }
  
  # 3. Calculate Weighted Logit Score (S)
  # Start the logit score with the intercept term
  logit_score <- INTERCEPT_B0 
  
  # Loop through all 8 gene pairs
  for (i in 1:nrow(purist_model)) {
    
    gene_a <- purist_model$GeneA[i]
    gene_b <- purist_model$GeneB[i]
    beta <- purist_model$Beta[i]
    
    # Extract expression and apply pseudocount
    expr_A <- dat[gene_a]
    expr_B <- dat[gene_b]
    
    # Calculate the log2 ratio (core of the PurIST model)
    # log2((GeneA + epsilon) / (GeneB + epsilon))
    log2_ratio <- log2((expr_A + epsilon) / (expr_B + epsilon))
    
    # Add the weighted term to the logit score
    logit_score <- logit_score + (beta * log2_ratio)
  }
  
  # 4. Convert Logit Score (S) to Probability (PurIST Score)
  # Uses the inverse logit function: P = 1 / (1 + exp(-S))
  # Note: Requires the 'gtools' package (inv.logit is the same)
  # If gtools is not installed, you can use: 1 / (1 + exp(-logit_score))
  purist_probability <- 1 / (1 + exp(-logit_score))
  
  # The result is the probability of being Basal-like (Basal-like > 0.5, Classical <= 0.5)
  return(purist_probability)
}


#' PurIST_Batch_Predict
#'
#' Applies the PurIST classifier to an entire gene expression matrix.
#' Handles samples with invalid (negative) expression data by assigning NA.
#'
#' @param expression_matrix A numeric matrix where rows are genes and columns are samples.
#'                          Row names must be HUGO gene symbols. Values must be normalized (e.g., TPM).
#' @param epsilon A small pseudocount (default 1).
#' @return A data frame with columns: Sample, PurIST_Probability (Basal-like), and Subtype_Call.
#'
#' @export
PurIST_Batch_Predict <- function(expression_matrix, epsilon = 1) {
  
  if (is.null(rownames(expression_matrix))) {
    stop("The expression matrix must have gene symbols as row names.")
  }
  
  samples <- colnames(expression_matrix)
  probabilities <- numeric(length(samples))
  names(probabilities) <- samples
  
  for (s in samples) {
    # Extract data for the current sample (column)
    sample_dat <- expression_matrix[, s]
    names(sample_dat) <- rownames(expression_matrix)
    
    # Check for non-positive values in the *full* sample vector before calling the core function
    if (any(is.na(sample_dat)) || any(sample_dat < 0)) {
      warning(sprintf("Sample %s contains NA or negative expression values and will be skipped.", s))
      probabilities[s] <- NA
      next
    }
    
    # Apply the single-sample classifier
    tryCatch({
      probabilities[s] <- PurIST_LogRatio(sample_dat, epsilon = epsilon)
    }, error = function(e) {
      # Catch errors primarily caused by missing genes in the matrix row names
      warning(sprintf("Classification error for sample %s: %s", s, e$message))
      probabilities[s] <- NA
    })
  }
  
  # Create the final results data frame
  results <- data.frame(
    Sample = samples,
    PurIST_Probability = probabilities,
    Subtype_Call = ifelse(probabilities > 0.5, "Basal-like", "Classical"),
    row.names = NULL
  ) %>%
    # Handle NA probabilities for the subtype call
    mutate(Subtype_Call = ifelse(is.na(PurIST_Probability), "No-Call (Invalid Data)", Subtype_Call))
  
  return(results)
}

