library(tidyverse)
library(here)

#empty list of dataframes
df_list <- list()

#list of TSA files
files <- list.files(path = "tsa files", full.names=TRUE, pattern = ".TSA$")

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

#save list
write_rds(df_list, "data/df_list1-1080.data.rds")

#write to seperate csv's
for(i in names(df_list)){
  write.csv(df_list[[i]], quote = FALSE, row.names = FALSE, file = paste0("csv files/",i,".csv"))
}
