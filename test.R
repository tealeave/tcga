# Load the required library
library(TCGAbiolinks)

# Download the subtype data into a variable
subtype_data <- TCGAbiolinks::PanCancerAtlas_subtypes()

# Print the column names to the console
print(colnames(subtype_data))