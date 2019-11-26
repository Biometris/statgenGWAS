## Load data from zip files.
dropsMap <- read.csv(unz(description = system.file("extdata", "infoSNP.zip", 
                                                   package = "statgenGWAS"),
                         filename = "7b-InfoSNP_50K_41722.csv"))

markerFile <- system.file("extdata", "genotyping.zip", package = "statgenGWAS")
dropsMarkers <- read.csv(unz(description = markerFile,
                             filename = "7a-Genotyping_50K_41722.csv"),
                         check.names = FALSE)

phenoFile <- system.file("extdata", "grainYield_components_BLUEs.zip", 
                         package = "statgenGWAS")
dropsPheno <- read.csv(unz(description = phenoFile, 
                           filename = "2b-GrainYield_components_BLUEs_level.csv"))
## Remove observations from 2011 from dropsPheno.
dropsPheno <- dropsPheno[substring(dropsPheno[["Experiment"]], 
                                   first = 4, last = 5) != "11", ]
dropsPheno <- droplevels(dropsPheno)

usethis::use_data(dropsMap, dropsMarkers, dropsPheno, overwrite = TRUE)
