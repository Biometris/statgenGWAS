# statgenGWAS 1.0.5.1

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
