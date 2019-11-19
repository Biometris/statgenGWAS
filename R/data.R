#' Data for vignette
#'
#' Data used in vignette
"gDataRiceCodedBeagle"

#' DROPS data sets
#' 
#' This dataset comes from the European Union project DROPS (DROught-tolerant 
#' yielding PlantS). A panel of 256 maize hybrids was grown with two water 
#' regimes (irrigated or rainfed), in seven fields in 2012 and 2013, 
#' respectively, spread along a climatic transect from western to eastern 
#' Europe, plus one site in Chile in 2013. This resulted in 29 experiments 
#' defined as the combination of one year, one site and one water regime, with 
#' two and three repetitions for rainfed and irrigated treatments, respectively.
#' A detailed environmental characterisation was carried out, with hourly 
#' records of micrometeorological data and soil water status, and associated
#' with precise measurement of phenology. Grain yield and its components were
#' measured at the end of the experiment. The main purpose of this dataset 
#' consists in using the environmental characterisation to quantify the genetic
#' variability of maize grain yield in response to the environmental drivers
#' for genotype-by-environment interaction. For instance, allelic effects at
#' QTLs identified over the field network are consistent within a scenario but
#' largely differ between scenarios. \cr\cr 
#' The data is split in three separate data.frames. 
#' \describe{ 
#' \item{\strong{dropsMarkers}}{This data.frame contains the 50K genotyping 
#' matrix coded in allelic dose (012) filtered and imputed. Genotyping of 
#' 41,722 loci on 246 parental lines were obtained using 50K Illumina Infinium 
#' HD arrays (Ganal et al., 2011). Genotype were coded in allelic dose with 0 
#' for the minor allele, 1 for the heterozygote, and 2 for the major allele. 
#' Genotype were filtered (MAF>1\%) and missing data imputed using Beagle v3.\cr
#' A data.frame with 246 rows and 41723 columns.
#' \describe{
#' \item{Ind}{name of the genotype}
#' \item{SYMN83 to PZE-110111485}{coded QTLs}
#' }
#' }  
#' \item{\strong{dropsMap}}{This data.frame contains the description of the 
#' 41,722 loci genotyped by 50K Illumina Infinium Array on the 246 lines.\cr
#' A data.frame with 41722 rows and 5 columns.
#' \describe{
#'   \item{SNP.names}{name of the SNP}
#'   \item{Chromosome}{number of the B73 reference genome V2}
#'   \item{Position}{position on the B73 reference genome V2 in centimorgan}
#'   \item{allele1}{first original allele (A, T, G or C)}
#'   \item{allele2}{second original allele (A, T, G or C)}
#' }
#' }
#' \item{\strong{dropsPheno}}{This data.frame contains the genotypic means 
#' (Best Linear Unbiased Estimation, BLUEs), with one value per experiment 
#' (Location × year × water regime) per genotype.\cr
#' A data.frane with 7390 rows and 14 columns.\cr
#' \describe{
#' \item{Experiment}{experiments ID described by the three first letters of the
#' city’s name followed by the year of experiment and the water regime with W 
#' for watered and R for rain-fed.}
#' \item{parent1}{identifier of donor dent line}
#' \item{Code_ID}{identifier of the genotype}
#' \item{Variety_ID}{identifier of the genotype}
#' \item{Accession_ID}{identifier of the genotype}
#' \item{geno.panel}{project in which the genetic material was generated} 
#' \item{grain.yield}{genotypic mean for yield adjusted at 15% grain moisture, 
#' in ton per hectare (t ha-1)}
#' \item{grain.number}{genotypic mean for number of grain per square meter}
#' \item{grain.weight}{genotypic mean for individual grain weight in milligram
#'  (mg)}
#' \item{anthesis}{genotypic mean for male flowering (pollen shed), in thermal
#' time cumulated since emergence (d20°C)}
#' \item{silking}{genotypic mean for female flowering (silking emergence), in 
#' thermal time cumulated since emergence (d20°C)}
#' \item{plant.height}{genotypic mean for plant height, from ground level to 
#' the base of the flag leaf (highest) leaf in centimeter (cm)} 
#' \item{tassel.height}{genotypic mean for plant height including tassel, from
#' ground level to the highest point of the tassel in centimeter (cm)}
#' \item{ear.height}{genotypic mean for ear insertion height, from ground level
#' to ligule of the highest ear leaf in centimeter (cm)}
#' }
#' }
#' }
#' 
#' @usage NULL
#' @format NULL
#' 
#' @source \url{https://data.inra.fr/dataset.xhtml?persistentId=doi:10.15454/IASSTN}
#' 
#' @references Millet, E. J., Pommier, C., et al. (2019). A multi-site 
#' experiment in a network of European fields for assessing the maize yield 
#' response to environmental scenarios [Data set]. 
#' \url{https://doi.org/10.15454/IASSTN}
#' @references Ganal MW, et al. (2011) A Large Maize (Zea mays L.) SNP 
#' Genotyping Array: Development and Germplasm Genotyping, and Genetic Mapping
#' to Compare with the B73 Reference Genome. PLoS ONE 6(12): e28334. 
#' \url{https://doi.org/10.1371/journal.pone.0028334}
#' 
#' @name dropsData
NULL

#' @usage NULL
#' @format NULL
#' @rdname dropsData
"dropsMap"

#' @usage NULL
#' @format NULL
#' @rdname dropsData
"dropsMarkers"

#' @usage NULL
#' @format NULL
#' @rdname dropsData 
"dropsPheno"
