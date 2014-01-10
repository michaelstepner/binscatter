binscatter: Binned Scatterplots
===============================

binscatter is a Stata program which generates binned scatterplots.  These are a convenient way of observing the relationship between two variables, or visualizing OLS regressions.  They are especially useful when working with large datasets.

### Installation ###

Install binscatter in Stata from the SSC repository: `ssc install binscatter`

After installing binscatter, you can read the documentation by running `help binscatter`.  The **Examples** section of the help file contains a clickable walk-through of binscatter's various features.

### Technical Description ###

Binned scatterplots are a non-parametric method of plotting the conditional expectation function (which describes the average y-value for each x-value).

To generate a binned scatterplot, binscatter groups the x-axis variable into equal-sized bins, computes the mean of the x-axis and y-axis variables within each bin, then creates a scatterplot of these data points.  By default, binscatter also plots a linear fit line using OLS, which represents the best linear approximation to the conditional expectation function.

binscatter provides built-in options to control for covariates before plotting the relationship, and can automatically plot regression discontinuities.  All procedures in binscatter are optimized for speed in large datasets.

### Example ###

The following graph shows the relationship between quality of teaching in elementary or middle school and a student's earnings at age 28.

![Binscatter example](http://michaelstepner.github.io/img/binscatter_va.jpg "Binscatter example")

This graph is a visual representation of a multivariate regression with 650,965 observations.  The regression finds that after controlling for a number of characteristics that affect student achievement (like class size and parental income), a 1 unit increase in Normalized Teacher Value Added is associated with a $350 increase in Earnings at Age 28.

The graph was created using binscatter:

* binscatter first regressed the y- and x-axis variables on the set of control variables, and generated the residuals from those regressions.  (Note that this is the first step of a partitioned regression.  We could regress the y-residuals on the x-residuals and obtain the coefficient from the full multivariate regression.)

* binscatter then grouped the residualized x-variable into 20 equal-sized bins, computed the mean of the x-variable and y-variable residuals within each bin, and created a scatterplot of these 20 data points. Each dot shows the average "Earnings at Age 28" for a given level of "Teacher Value Added", holding the controls constant.

* Finally, binscatter plotted the best linear fit line, constructed from an OLS regession of the y-residuals on the x-residuals.  The slope of the fit line matches the coefficient of the multivariate regression. 

> **Source:** Raj Chetty, John N. Friedman, and Jonah E. Rockoff. 2013. "Measuring the Impacts of Teacher I: Evaluating Bias in Teacher Value-Added Estimates." NBER Working Paper 19423. http://obs.rc.fas.harvard.edu/chetty/w19423.pdf
>
> The graph above is Figure 2a in the paper, and the associated regression is reported in Table 3, Column 1.
