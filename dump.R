#### Step Two: Data Import
#$ [Insert Data Here]
DATA <- as.data.frame(Bloodsmear_data_full)
#$ [Optional: Data Cleaning]
DATA$Timepoint <- replace(DATA$Timepoint, DATA$Timepoint == 17, 7) #okay that worked

#Subsetting
DATA$Timepoint <- as.factor(DATA$Timepoint)
levels(DATA$Timepoint)


#Splitting the datasets by timepoint
Day2_Data = subset(DATA, DATA$Timepoint == 2)
Day2_x100s = Day2_Data[c(1:6, 11:13)]
Day2_x1000s = Day2_Data[1:10]

Day7_Data = subset(DATA, DATA$Timepoint == 7)
Day7_x100s = Day7_Data[c(1:6, 11:13)]
Day7_x1000s = Day7_Data[1:10]

Day14_Data = subset(DATA, DATA$Timepoint == 14)
Day14_x100s = Day14_Data[c(1:6, 11:13)]
Day14_x1000s = Day14_Data[1:10]

Day21_Data = subset(DATA, DATA$Timepoint == 21)
Day21_x100s = Day21_Data[c(1:6, 11:13)]
Day21_x1000s = Day21_Data[1:10]

#Then I need to split them further by temp..... wait. I should have split them by temp for the impacts of day, because day is still in the model....
















comp_clean <- function(acomp_object, identifier) {
  # Get the original variable name as a string
  original_name <- deparse(substitute(acomp_object))
  
  # Append the suffix to create a new name
  new_name <- paste0(original_name, "_cleaned")
  
  # Clean the acomp_object
  BDL_MNAR_REPLACE <- 0.001 # assign a value to replace BDL and MNARs
  print("Original acomp_object:")
  print(acomp_object)
  
  acomp_object <- zeroreplace(acomp_object, BDL_MNAR_REPLACE) # replace zeros
  print("After zeroreplace:")
  print(acomp_object)
  
  acomp_object[is.na(acomp_object)] <- BDL_MNAR_REPLACE # replace NAs
  print("After replacing NAs:")
  print(acomp_object)
  
  # Assign the cleaned object to the new name in the global environment
  assign(new_name, acomp_object, envir = .GlobalEnv)
  
  # Explore the composition again
  clean_Ident <- paste0(identifier, "_cleaned_comp_explore_graphs_") # this is to identify the output from your specific run or hypothesis
  pdf(file = paste(clean_Ident, ".pdf", sep = "_"))
  
  # PDF code
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  text(5, 8, "Composition after treating missing values")
  text(5, 7, identifier)
  
  plot(acomp_object)
  title("Untransformed Cleaned Composition", outer = TRUE, line = -1)
  plot(clr(acomp_object)) # clr transformation (central log)
  title("CLR-transformed Cleaned Composition", outer = TRUE, line = -1)
  plot(ilr(acomp_object)) # ilr transformation (isometric log) (check)
  title("ILR-transformed Cleaned Composition", outer = TRUE, line = -1)
  
  dev.off()
  print("Graphs Made!")
}

# Example usage:
# comp_clean(Spooky_acomp, "TestIdentifier")