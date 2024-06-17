# Load necessary libraries
library(ggplot2)
library(reshape2)

# Function to read files and create a combined data frame
read_and_combine_files <- function(file_list) {
  combined_data <- data.frame()
  
  for (file in file_list) {
    # Read the file
    data <- read.table(file, header = FALSE, sep = "\t")
    
    # Select the 8th column (assumed to be V8)
    data_subset <- data[, 8, drop = FALSE]
    
    # Add a column for the file name
    data_subset$file <- basename(file)
    
    # Combine data
    combined_data <- rbind(combined_data, data_subset)
  }
  
  return(combined_data)
}

# Get the list of text files in the current directory
file_list <- list.files(pattern = "\\.txt$")

# Read and combine data from all files
combined_data <- read_and_combine_files(file_list)

# Convert file column to factor
combined_data$file <- as.factor(combined_data$file)

# Create a swarm plot
ggplot(combined_data, aes(x = file, y = V8, color = file)) +
  geom_quasirandom() +
  labs(title = "Swarm Plot of Column 8 Values from Multiple Files",
       x = "File",
       y = "Column 8 Values") +
  theme_minimal()

# Save the plot to a file
ggsave("swarm_plot.png")
