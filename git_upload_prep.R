for(i in list.dirs()[stringr::str_detect(list.dirs(), "simulation_exports")]){
  
 writeLines(
   paste0(
          "All files can be accessed at: https://sussex.box.com/s/wgsemizl7ulv4zt53zxi0tyr38n6myeb \n\nContact Martina at m.sladekova@sussex.ac.uk if you encounter problems accessing the files"
          ), 
   con = paste0(i, "/export_location.txt"))
  
}
