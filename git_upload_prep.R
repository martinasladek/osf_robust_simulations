for(i in list.dirs()[stringr::str_detect(list.dirs(), "simulation_exports")]){
  
 writeLines("", con = paste0(i, "/export_location.txt"))
  
}
