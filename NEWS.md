# statgenGWAS develop

* Removed use of deprecated `ggplot2::expand_scale`.

# statgenGWAS 1.0.4

* StringsAsFactors = TRUE added where applicable to comply with new defaults in R 4.0.
* Bug in plotting of GWAS objects fixed. The first trait was always plotted in case more than one trait was present.

# statgenGWAS 1.0.3

* Dependency on deprecated `rvg::ph_with_vg_at` removed.
* OMP_THREAD_LIMIT is now respected in parallel code.

# statgenGWAS 1.0.2

* Fixed problem with compilation on Solaris

# statgenGWAS 1.0.1

* Initial CRAN version
