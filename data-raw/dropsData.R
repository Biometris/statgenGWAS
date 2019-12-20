## Load data from zip files.
dropsMap <- read.csv(unz(description = "./data-raw/infoSNP.zip",
                         filename = "7b-InfoSNP_50K_41722.csv"))

dropsMarkers <- read.csv(unz(description = "./data-raw/genotyping.zip",
                             filename = "7a-Genotyping_50K_41722.csv"),
                         check.names = FALSE)

dropsPheno <- read.csv(unz(description = "./data-raw/grainYield_components_BLUEs.zip", 
                           filename = "2b-GrainYield_components_BLUEs_level.csv"))
## Remove observations from 2011 from dropsPheno.
dropsPheno <- dropsPheno[substring(dropsPheno[["Experiment"]], 
                                   first = 4, last = 5) != "11" &
                           dropsPheno[["Experiment"]] != "Gra13W", ]
dropsPheno <- droplevels(dropsPheno)

usethis::use_data(dropsMap, dropsMarkers, dropsPheno, overwrite = TRUE)
