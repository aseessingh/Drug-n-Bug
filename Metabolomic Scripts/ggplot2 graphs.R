library(ggplot2)
library(readr)
library(ggrepel)
library(dplyr)

create_custom_plot <- function(file_name, x, y, color, label, plot_title, color_title, output_file) {
  # Read the data from the file
  data <- as.data.frame(read_csv(file_name, show_col_types = FALSE))
  
  # Standardize the ontology column
  data$ontology <- tolower(trimws(data$ontology))
  
  # Calculate the top 0.1% cutoff values for x and y
  x_cutoff <- quantile(data[[x]], 0.999, na.rm = TRUE)
  y_cutoff <- quantile(data[[y]], 0.999, na.rm = TRUE)
  
  # Filter data for labeling based on the cutoffs
  label_data <- data %>%
    filter(!is.na(.data[[x]]) & !is.na(.data[[y]])) %>%
    filter(abs(.data[[x]]) > 0.1 & abs(.data[[y]]) > 0.1 | .data[[x]] >= x_cutoff | .data[[y]] >= y_cutoff)
  
  # Data below the cutoff for grey points
  below_cutoff_data <- data %>%
    filter(!is.na(.data[[x]]) & !is.na(.data[[y]])) %>%
    filter(abs(.data[[x]]) < 0.1 | abs(.data[[y]]) < 0.1)
  
  # Identify extreme points on the x and y axes
  extreme_x <- data %>%
    filter(.data[[x]] == min(.data[[x]], na.rm = TRUE) | .data[[x]] == max(.data[[x]], na.rm = TRUE))
  extreme_y <- data %>%
    filter(.data[[y]] == min(.data[[y]], na.rm = TRUE) | .data[[y]] == max(.data[[y]], na.rm = TRUE))
  
  # Combine label data with extreme points
  label_data <- bind_rows(label_data, extreme_x, extreme_y) %>%
    distinct()  # Remove duplicates
  
  plot <- ggplot() +
    # Plot the data below the cutoff in grey
    geom_point(data = below_cutoff_data, aes(x = .data[[x]], y = .data[[y]]), color = "grey", size = 3, alpha = 0.6) +
    # Plot all data points
    geom_point(data = data, aes(x = .data[[x]], y = .data[[y]], color = .data[[color]], shape = .data[["ontology"]]), size = 3, alpha = 0.6) +
    # Label the filtered data
    geom_text_repel(data = label_data, aes(x = .data[[x]], y = .data[[y]], label = .data[[label]]), 
                    size = 2.8, 
                    vjust = -1, 
                    hjust = 1, 
                    max.overlaps = Inf,  # Allow unlimited overlaps
                    nudge_x = 0.05,      # Nudge labels slightly to the right
                    nudge_y = 0.05) +    # Nudge labels slightly up
    scale_color_gradient(low = "blue", high = "red") +  # Use gradient for continuous values
    scale_shape_manual(values = c(16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,26)) +  # Adding more shapes
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +  # Horizontal axis at y = 0
    geom_vline(xintercept = 0, color = "black", linetype = "solid") +
    labs(title = plot_title,
         x = x,
         y = y,
         color = color_title,
         shape = "ontology") +
    theme_minimal(base_size = 11) +
    theme(
      legend.title = element_text(face = "bold"),  # Ensure legend titles are visible
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    scale_x_continuous() +  # Automatically set x-axis limits
    scale_y_continuous()    # Automatically set y-axis limits
  
  # Save the plot
  ggsave(output_file, plot, width = 10, height = 8)
}
lm_case1 <- create_custom_plot(file_name = "lm_case1.csv",x = "Beta_Dose",y="Beta_HbA1c",color = "Beta_Serum",label = "metabolite",plot_title = "Linear Regression Plot-All Group 3 Metformin Intake",color_title = "Beta_Serum",output_file = "lm_case1.png" )
lm_case2 <- create_custom_plot(file_name = "lm_case2.csv",x = "Beta_Dose",y="Beta_HbA1c",color = "Beta_Serum",label = "metabolite",plot_title = "Linear Regression Plot-All Group 3 Metformin Intake with Dose Data",color_title = "Beta_Serum",output_file = "lm_case2.png" )
lm_case3 <- create_custom_plot(file_name = "lm_case3.csv",x = "Beta_Dose",y="Beta_HbA1c",color = "Beta_Serum",label = "metabolite",plot_title = "Linear Regression Plot-All Group 3 Metformin Intake with Dose and Serum Data",color_title = "Beta_Serum",output_file = "lm_case3.png" )
lm_case4 <- create_custom_plot(file_name = "lm_case4.csv",x = "Beta_Dose",y="Beta_HbA1c",color = "Beta_Serum",label = "metabolite",plot_title = "Linear Regression Plot-All Group 3 Metformin Intake Pos vs Neg",color_title = "Beta_Serum",output_file = "lm_case4.png" )
metad_case1 <- create_custom_plot(file_name = "metad_case1.csv",x = "Corr_Dose",y="Corr_HbA1c",color = "Corr_Serum",label = "CHEMICAL_NAME",plot_title = "MetadeconfoundR-All Group 3 Metformin Intake",color_title = "Corr_Serum",output_file = "metad_case1.png" )
metad_case2 <- create_custom_plot(file_name = "metad_case1.csv",x = "Corr_Dose",y="Corr_HbA1c",color = "Corr_Serum",label = "CHEMICAL_NAME",plot_title = "MetadeconfoundR-All Group 3 Metformin Intake with Dose Data",color_title = "Corr_Serum",output_file = "metad_case2.png" )
metad_case3<- create_custom_plot(file_name = "metad_case3.csv",x = "Corr_Dose",y="Corr_HbA1c",color = "Corr_Serum",label = "CHEMICAL_NAME",plot_title = "MetadeconfoundR-Group 3 Metformin Intake with Dose and Serum Data",color_title = "Corr_Serum",output_file = "metad_case3.png" )
metad_case4 <- create_custom_plot(file_name = "metad_case4.csv",x = "Corr_Dose",y="Corr_HbA1c",color = "Corr_Serum",label = "CHEMICAL_NAME",plot_title = "MetadeconfoundR-All Group 3 Metformin Intake Pos vs Neg",color_title = "Corr_Serum",output_file = "metad_case4.png" )

#Plot Function 
create_custom_plot_log10 <- function(file_name, x, y, color, label, plot_title, color_title, output_file) {
  # Read the data from the file
  data <- as.data.frame(read_csv(file_name, show_col_types = FALSE))
  
  # Standardize the ontology column
  data$ontology <- tolower(trimws(data$ontology))
  
  # Calculate the top 0.1% cutoff values for x and y
  x_cutoff <- quantile(data[[x]], 0.999, na.rm = TRUE)
  y_cutoff <- quantile(data[[y]], 0.999, na.rm = TRUE)
  
  # Filter data for labeling based on the presence of x and y values or top 0.1% cutoffs
  label_data <- data %>%
    filter((.data[[x]] != 0 & .data[[y]] != 0) | (.data[[x]] >= x_cutoff | .data[[y]] >= y_cutoff))
  
  # Define a set of 27 shapes
  shapes <- c(16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 26)
  
  plot <- ggplot() +
    # Plot all data points
    geom_point(data = data, aes(x = .data[[x]], y = .data[[y]], color = .data[[color]], shape = .data[["ontology"]]), size = 3, alpha = 0.6) +
    # Label the filtered data
    geom_text_repel(data = label_data, aes(x = .data[[x]], y = .data[[y]], label = .data[[label]]), 
                    size = 2.8, 
                    vjust = -1, 
                    hjust = 1, 
                    max.overlaps = Inf,  # Allow unlimited overlaps
                    nudge_x = 0.05,      # Nudge labels slightly to the right
                    nudge_y = 0.05) +    # Nudge labels slightly up
    scale_color_gradient(low = "blue", high = "red") +  # Use gradient for continuous values
    scale_shape_manual(values = shapes) +  # Use the defined set of shapes
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +  # Horizontal axis at y = 0
    geom_vline(xintercept = 0, color = "black", linetype = "solid") +
    labs(title = plot_title,
         x = x,
         y = y,
         color = color_title,
         shape = "ontology") +
    theme_minimal(base_size = 11) +
    theme(
      legend.title = element_text(face = "bold"),  # Ensure legend titles are visible
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    scale_x_continuous() +  # Automatically set x-axis limits
    scale_y_continuous()    # Automatically set y-axis limits
  
  # Save the plot
  ggsave(output_file, plot, width = 10, height = 8)
}

lm_case3_log10 <- create_custom_plot_log10(file_name = "lm_case3.csv",x = "Beta_Dose",y="Beta_HbA1c",color = "Beta_Serum",label = "metabolite",plot_title = "Linear Regression Plot-All Group 3 Metformin Intake with Dose and Serum Data",color_title = "Beta_Serum",output_file = "lm_case3_log10.png" )
metad_case3_log10 <- create_custom_plot_log10(file_name = "metad_case3.csv",x = "Corr_Dose",y="Corr_HbA1c",color = "Corr_Serum",label = "CHEMICAL_NAME",plot_title = "MetadeconfoundR-Group 3 Metformin Intake with Dose and Serum Data",color_title = "Corr_Serum",output_file = "metad_case3_log10.png" )
lm_case3_AG_log10 <- create_custom_plot_log10(file_name = "lm_case3_AG.csv",x = "Beta_Dose",y="Beta_AG",color = "Beta_Serum",label = "metabolite",plot_title = "Linear Regression Plot-All Group 3 Metformin Intake with Dose and Serum Data",color_title = "Beta_Serum",output_file = "lm_case3_AG_log10.png" )
metad_case4_log10 <- create_custom_plot_log10(file_name = "metad_case4.csv",x = "Corr_Dose",y="Corr_HbA1c",color = "Corr_Serum",label = "CHEMICAL_NAME",plot_title = "MetadeconfoundR-All Group 3 Metformin Intake Pos vs Neg",color_title = "Corr_Serum",output_file = "metad_case4_log10.png" )
