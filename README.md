binscatter: Binned Scatterplots
===============================

binscatter is a Stata program which generates binned scatterplots.

### Description ###

Binned scatterplots are a non-parametric method of plotting the conditional expectation function (which describes the average y-value for each x-value). The resulting graph is a convenient way of depicting the relationship between two variables, and visualizing OLS regressions.

To generate a binned scatterplot, binscatter groups the x-axis variable into equal-sized bins, computes the mean of the x-axis and y-axis variables within each bin, then creates a scatterplot of these data points.  By default, binscatter also plots a linear fit line using OLS, which represents the best linear approximation to the conditional expectation function.

binscatter provides built-in options to control for covariates before plotting the relationship, and can automatically plot regression discontinuities.  All procedures in binscatter are optimized for speed in large datasets.

### Installation ###

Install binscatter in Stata from the SSC repository: `ssc install binscatter`

If already installed, check for updates to binscatter: `adoupdate`

Verify which version of binscatter is presently being used: `which binscatter`

### Documentation ###

After installing binscatter, run `help binscatter`.
