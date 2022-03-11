packages <- c("MplusAutomation","MASS","stringr")
# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages],Ncpus=4)
# }

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Change
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)


filenames <-list.files(paste0(wd,"/Templates"))

for (j in filenames){
  createModels(paste0(wd,"/Templates/",j))
}

#Change number of starts etc
### Run all models within a directory and log output
### Skipping models with existing output files
### replaceOutfile ="modifiedDate"
### Checks whether there is an existing output file for a given input file
### If there is, it checks to see whether the date the input file was modified is newer than output
### If input is newer then the model is run, otherwise it is skipped: replaceOutfile = "modifiedDate"
runModels(wd,recursive=TRUE,showOutput = FALSE,replaceOutfile = "never",
          logFile = paste0(wd,"/ComparisonLog.txt"))