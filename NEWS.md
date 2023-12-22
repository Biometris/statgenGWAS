# statgenGWAS 1.0.9.9000

* The QTL plot now has a legend explaining the meaning of the colors used.
* Functions no longer rely on soft-deprecated ggplot2 functions.
* The support for tibbles with phenotypic data is improved.

# statgenGWAS 1.0.9

* The kinship function now has an argument MAF allowing for excluding SNPs with a low minor allele frequency from the calculations.
* Small change in tests to accommodate future ggplot2 release.

# statgenGWAS 1.0.8

* A function readPLINK is added that reads file in PLINK format and converts them directly to a gData object.
* The kinship function now has an option MAF allowing for excluding markers with a MAF below a specified threshold.
* The gData object now has a plot function allowing plotting of the genetic map.
* A bug in gData summaries when only one trait is present in the phenotypic data is fixed.
* A bug in runSingleTraitGwas is fixed now allowing for using traits that are in pheno as covariates in covar.
* Sorting of chromosomes in map is now done correctly for numeric values stored as character.

# statgenGWAS 1.0.7

* No user visual changes.

# statgenGWAS 1.0.6

* The plot function now has an option title.
* The codeMarkers function now has an option MAC equivalent to the same option in runSingleTraitGwas.
* Bug when using tibbles instead of data.frames as input for covariates is fixed.
* The kinship function now has an option for returning a named identity matrix as output.
* An extra check for the presence of . in a character marker matrix is added. 
* The correct SNPs are now excluded when using MAF/MAC.
* The function documentation and vignette have been revised and extra clarification is added where needed.
* The beagle software used for imputation is updated to version 5.2.
* Bug in the computation of the genomic inflation factor is fixed. It was incorrectly computed as 1 / genomic inflation factor.
* The data in the vignette is updated such that the environments included match those in the other statgen packages.
* The manhattan plot function gained arguments for setting the start and end position when plotting results for a single chromosome.

# statgenGWAS 1.0.5

* An option fdr (false discovery rate) is added for selecting significant SNPs when running single trait GWAS. The procedure for this is based on Brzyski et al. (2017). See the function documentation and vignette for an extensive explanation.
* The summary function for objects of class GWAS now has an option traits for restricting the number of traits for which the summary is printed.
* Problem where codeMarkers was getting very slow for large character matrix inputs is fixed by moving part of the code to c++.
* Full ggplot2 is no longer imported.
* Removed use of deprecated `ggplot2::expand_scale`.

# statgenGWAS 1.0.4

* stringsAsFactors = TRUE added where applicable to comply with new defaults in R 4.0.
* Bug in plotting of GWAS objects fixed. The first trait was always plotted in case more than one trait was present.

# statgenGWAS 1.0.3

* Dependency on deprecated `rvg::ph_with_vg_at` removed.
* OMP_THREAD_LIMIT is now respected in parallel code.

# statgenGWAS 1.0.2

* Fixed problem with compilation on Solaris

# statgenGWAS 1.0.1

* Initial CRAN version
