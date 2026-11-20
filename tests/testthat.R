# Unit tests of the pipeline's functions. They run in a few seconds and fit no
# models. From the repository root:
#
#     Rscript tests/testthat.R
#
# They complement, and do not replace, a run of the pipeline: they check the
# pieces of arithmetic and bookkeeping that the results rest on (how a
# trapping period's weather is averaged, how observations are matched to
# collection periods, the clustered standard error, how models are paired, the
# logger clock) on small inputs whose right answers can be worked out by hand.
testthat::test_dir(file.path("tests", "testthat"), stop_on_failure = TRUE)
