##BINSCATTER: Stata module to generate binned scatterplots


binscatter is a Stata program which generates binned scatterplots, and is optimized for speed in large datasets.

Binned scatterplots provide a non-parametric visualization of the conditional expectation function.  With a large number of observations, a scatterplot that plots every data point would become too crowded to interpret visually.  binscatter groups the x-axis variable into equal-sized bins, computes the mean of the x-axis and y-axis variables within each bin, then creates a scatterplot of these data points.  The resulting graph is a convenient way of visualizing the relationship between two variables.

binscatter provides built-in options to control for covariates before plotting the relationship. Additionally, binscatter will plot fit lines based on the underlying data, and can automatically handle regression discontinuities.

### Installation

You can install binscatter in Stata from the SSC repository using `ssc install binscatter`.

If you have already installed binscatter but need to update to the latest version, run `adoupdate`.  To verify which version is presently being used, run `which binscatter`.

### Documentation

After installing binscatter, run `help binscatter`.
