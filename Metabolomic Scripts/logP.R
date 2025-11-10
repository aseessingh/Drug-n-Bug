# -----------------------------------------------
#            Calculate LogP from SMILES
# -----------------------------------------------

# --- Load necessary libraries ---
library(rcdk)
library(rJava)
library(readr)

# --- INPUT: Set the path for your CSV file here ---
input_file <- "input.csv"   # Path to your CSV file (must contain CHEMICAL_NAME and SMILES columns)
output_file <- "output_with_logP.csv"  # Path to output the result

# --- Read input data ---
df <- read_csv(input_file)

# Check if required columns exist
if(!all(c("CHEMICAL_NAME", "SMILES") %in% colnames(df))) {
  stop("The input CSV file must contain 'CHEMICAL_NAME' and 'SMILES' columns.")
}

# --- Function to calculate logP for a given SMILES string ---
calculate_logP <- function(smiles) {
  if (is.na(smiles)) {
    return(NA)  # Return NA if the SMILES string is NA
  }
  molecule <- parse.smiles(smiles)[[1]]  # Convert SMILES to molecule object
  if (is.null(molecule)) {
    return(NA)  # Return NA if the SMILES is invalid
  }
  logP <- get.desc(molecule, "logP")     # Calculate logP
  return(logP)
}

# --- Apply the function to the SMILES column ---
df$logP <- sapply(df$SMILES, calculate_logP)

# --- View the results (optional) ---
print(df)

# --- Save the results to a new CSV file ---
write_csv(df, output_file)

# Message to confirm the process is complete
message("LogP calculation complete. The results are saved in '", output_file, "'")
