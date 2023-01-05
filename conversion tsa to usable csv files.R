library(tidyverse)
library(here)

## empty list of dataframes
df_list <- list()

## list of TSA files
files <- list.files(path = "tsa files", full.names = TRUE, pattern = ".TSA$")

#loop trough TSA file list
#warnings may be ignored
for (i in 1:length(files)) {
  
  #read each file
  d <- readLines(files[i],  encoding = "windows-1252")
    
    #convert to tibble
    d <- as_tibble(d) %>%
      slice(which.max(value == "#TRIAL BEGIN") : n()) %>%
      separate(col = value, into = c("var","value"))  %>%
      filter(var %in% c('year','study','interventionEvent','interventionTotal','controlEvent','controlTotal'))
    
    #select only relevant rows and convert to a usable data format
    d <- d %>%
      mutate(id = rep(1:(dim(d)[1]/6),each = 6)) %>%
      pivot_wider(names_from = var) %>%
      mutate(across(c(interventionEvent,interventionTotal,controlEvent,controlTotal),as.numeric)) %>%
      mutate(RCT=paste(study,year,sep="")) %>%
      select(RCT,interventionEvent,interventionTotal,controlEvent,controlTotal)
    
    #add to list
    df_list[[i]] <- d
  
}

#rename to original file names
names(df_list) <- files
names(df_list) <- gsub("tsa files/","",gsub(".TSA","",files))

#identify meta-analyses (MAs) with <3 studies because this leads to an error in the analyses.
small_MAs <- sapply(df_list, function(x) { 
  if (nrow(x) < 3) return(x)
})
small_MAs <- Filter(Negate(is.null), small_MAs)

#limit df_list to include only MAs with >2 studies for a given outcome
large_MAs <- sapply(df_list, function(x) { 
  if (nrow(x) > 2) return(x)
})
large_MAs <- Filter(Negate(is.null), large_MAs)

#write to seperate csv's
for(i in names(large_MAs)){
  write.csv(large_MAs[[i]], quote = FALSE, row.names = FALSE, file = paste0("csv files/",i,".csv"))
}

