# MSP-DSC: multispecimen paleointensity plotter and analyser

An interactive **R Shiny** application to analyse multispecimen (MSP) paleointensity data. It implements the MSP-DB, MSP-FC and MSP-DSC protocols (Dekkers & Böhnel, 2006; Fabian & Leonhardt, 2010) and, for MSP-DSC, five ways of handling the multidomain correction parameter alpha, including per-specimen optimization with Ridge regularization.

The application accompanies the paper:

> Camps P., Cych B., Antonio P., Raffray A. (2026) Optimization of the domain-state correction in Multiple-Specimen Paleointensity: Balancing prior and data. *Studia Geophysica et Geodaetica*, XX(X), pp-pp. [DOI to be added]

**Main features**

* MSP-DB, MSP-FC and MSP-DSC protocols, with error-weighted regression anchored at (0,-1) for FC and DSC.
* MSP-DSC with five alpha strategies (constant, single optimized, per-specimen with uniform bounds, per-specimen with individual bounds, per-specimen with a Gaussian prior).
* Detection of influential specimens (Cook's distance) and bootstrapped 95% confidence intervals.
* Paired-bootstrap comparison of the five DSC cases (Holm-corrected p-values).
* Plots and tables that can be downloaded (pdf, png, tiff, csv).

### Requirements and installation

R (a recent version) and the following packages:

```r
install.packages(c("shiny", "shinythemes", "markdown", "DT", "withr", "dplyr",
                   "ggplot2", "pracma", "boot", "htmlwidgets", "yaml",
                   "ggsci", "ggrepel"))
```

Clone or download the repository, then start the application from R or RStudio:

```r
shiny::runApp("path/to/the/repository")
```

or open `ui.R` or `server.R` in RStudio and click on *Run App*.

The repository must contain the following files:

| File | Role |
|---|---|
| `ui.R` | user interface |
| `server.R` | server logic (it loads `MSP_func.R`) |
| `MSP_func.R` | functions and constants: alpha bounds [0, 0.8], alpha of case 1 (0.5), lambda search grid (0.5 to 2, step 0.05), minimum sigma of case 5 (0.05), default k (4), alteration thresholds, confidence levels |
| `ReadMe.md` | this file, also displayed in the *Readme* tab of the application |
| `www/logoGM.png` | logo displayed in the application |
| `data/example.csv` | example of input file, available with the button *Get example of inputfile* |

### Format for the INPUT data files

The program requires a `.csv` file (comma separated values) with a header line exactly as follows:

```
name, Field, m0, D0, I0, m1, D1, I1, m2, D2, I2, m3, D3, I3, m4, D4, I4, Temp, a_max, a_min
```

where:

* `name` is the specimen name;
* `Field` is the laboratory field (in µT);
* `m`, `D` and `I` are the magnetization intensity, the declination and the inclination (in degrees), respectively;
* `0`, `1`, `2`, `3` and `4` are the successive steps of the protocol, as described in Fabian & Leonhardt (2010);
* `Temp` is the dwell temperature;
* `a_max` and `a_min` (optional) are the upper and lower bounds of the parameter alpha of each specimen. They are only used by MSP-DSC cases 4 and 5 (see below).

To use MSP-DSC cases 4 and 5, the two columns `a_max` and `a_min` must be added at the end of the header, as shown above: they must be the **last two columns** of the file. The columns are ignored by MSP-DB, MSP-FC and by MSP-DSC cases 1, 2 and 3. If the columns are not named `a_min` and `a_max`, the last column is read as `a_min` and the second-to-last as `a_max` (order of the header above). If a file has no such columns, the other protocols and cases work as usual.

Use the character `#` to comment a line if you want to discard an entry. Blank lines are allowed in the data file and are automatically skipped.

### Notes on the MSP-DB protocol

The paleointensity is estimated with ordinary least squares regression (OLS), i.e. all residuals have the same variance and are uncorrelated.

```r
# R code for:
# Ordinary least squares regression with possible outliers
    B.QDB <- data.frame(Name = meas()$name)
    B.QDB$x <- B
    B.QDB$y <- QDB
    model_DB <- lm(y ~ x, data = B.QDB)
```

### Notes on the MSP-FC protocol

m2 is reversed if the angle between vector m2 and vector m0 is greater than pi/2.

```r
# Reverse m2 if the angle between m2 and m0 is greater than pi/2
    mm2 <- meas()$m2
    ind_reverse <- which(angle_m0m2 > pi / 2)
    mm2[ind_reverse] <- -mm2[ind_reverse]

# Eq. 16 in Fabian and Leonhardt, 2010.
    QFC <- 2 * (meas()$m1 - meas()$m0) / (2 * meas()$m0 - meas()$m1 - mm2)
```

The paleointensity is estimated with OLS (all residuals have the same variance and are uncorrelated). For physical reasons, i.e. no pTRM is acquired if a null field is applied during the heating/cooling cycle, the linear fit is anchored to the point (0,-1).

```r
# R code for:
# OLS anchored to the point (0,-1), excluding the influential data and the data with an inadequate fraction of NRM
    model_FC <- lm((y + 1) ~ x - 1, data = subset(B.QFC, quality == "Good"))
```

### Notes on the MSP-DSC protocol

The paleointensity is estimated by weighted least squares (WLS), i.e. all residuals have different variances and are uncorrelated. For the same physical reason as above, the linear fit is anchored to the point (0,-1). The error dQ of each measurement is computed from the alteration and from the difference m1 - m3 (Eqs. 19-22 of Fabian & Leonhardt, 2010), and the weights are derived from it.

```r
# R code for:
# WLS on selected points
    model_DSC <- lm((y + 1) ~ x - 1, data = subset(B.QDSC, quality == "Good"), weights = 1/dQ)
```

#### The parameter alpha: five ways of handling it (MSP-DSC)

In the DSC protocol, the corrected quantity of each specimen is the FC quantity plus a correction term built from the step m3 (the term depends on m1 - m3). The parameter alpha is the weight given to this correction term: alpha = 0 gives no DSC correction (the result is the same as MSP-FC), and the larger alpha, the stronger the correction. The paleointensity is the field for which Q = 0, i.e. the inverse of the slope of the regression anchored to (0,-1).

The menu *MSP-DSC: alpha strategy* (left panel, visible when the DSC protocol is selected) offers five cases:

| Case | alpha | Bounds of alpha | Regularization |
|---|---|---|---|
| 1 | constant, 0.5, the same for all specimens | none (fixed) | none |
| 2 | one single value, optimized, common to all specimens | [0, 0.8] | none |
| 3 | one value per specimen, optimized | [0, 0.8] for all specimens | ridge toward 0.5, strength lambda |
| 4 | one value per specimen, optimized | [a_min, a_max] of each specimen, read in the data file | ridge toward the middle of each interval, strength lambda |
| 5 | one value per specimen, optimized | [0, 0.8] for all specimens | Gaussian prior centred on the middle of [a_min, a_max] of each specimen, strength lambda_i = 1/sigma_i^2 with sigma_i = (a_max - a_min)/k |

**Case 1: constant alpha.** alpha is fixed to 0.5 for all specimens. Nothing is optimized. The value is displayed in the left panel. The regression is a weighted least squares anchored to (0,-1), with weights equal to the inverse of dQ.

**Case 2: one optimized alpha for all specimens.** A single alpha, common to all specimens, is searched between 0 and 0.8. It is the value that minimizes the weighted sum of squared residuals of the anchored regression. The optimization is done twice: a first time on all the selected specimens, to detect the influential points with Cook's distance; a second time, after removing these points, to obtain the final alpha and the final slope. The optimized alpha is displayed in the left panel.

**Case 3: one alpha per specimen, same bounds for all.** Each specimen has its own alpha, which can take any value between 0 and 0.8. The slope and the alpha values are estimated together, each alpha being adjusted to bring its specimen close to the regression line. Left alone, the individual alpha values would absorb most of the residuals and the regression would look better than it really is. They are therefore pulled toward 0.5 by a ridge penalty, whose strength is the parameter lambda (see *Choice of lambda* below).

**Case 4: one alpha per specimen, individual bounds.** Same principle as case 3, but each specimen has its own interval [a_min, a_max], read in the data file (the last two columns, see *Format for the INPUT data files*). The intervals typically come from independent constraints on each specimen (e.g. a rock-magnetic proxy of the domain state such as the median destructive field). The ridge penalty pulls each alpha toward the middle of its own interval. The program checks that the values are numerical and that 0 <= a_min <= a_max <= 1; if not, a message is displayed instead of the results.

**Case 5: one alpha per specimen, prior centred on [a_min, a_max].** Same estimation as cases 3 and 4 (slope and alpha values estimated together, ridge penalty), but the intervals of the data file are used as soft information instead of hard bounds. Each alpha is searched in [0, 0.8] and is pulled toward the middle of the interval [a_min, a_max] of its specimen. The penalty of specimen i is lambda_i * (alpha_i - centre_i)^2. In Bayesian terms (lambda_i = 1/sigma_i^2) it is a Gaussian prior of mean centre_i and standard deviation sigma_i = (a_max - a_min)/k, so that the interval is k standard deviations wide (k = 4 by default: the interval is +/- 2 sigma). k can be changed with the slider of the left panel (2 to 8). To avoid an infinite penalty for very narrow intervals, sigma_i is never smaller than 0.05. lambda_i is therefore fixed by the data file and by the choice of k, and is not selected by LOOCV: the tab *Lambda* is not available in case 5. The rest of the procedure (selection, Cook's distance, bootstrap with the optimization repeated for each sample) is the same as in cases 3 and 4. The range of the values lambda_i is displayed in the left panel.

Reading the prior: a large lambda_i (narrow interval) keeps alpha_i close to the middle of its interval; a small lambda_i leaves alpha_i free in [0, 0.8]. Compared with case 4, an alpha can leave its interval if the data require it, at a cost that increases with the distance to the centre. Compared with case 3, the information of the intervals is used. Unlike the LOOCV lambda of cases 3 and 4, which is often below 1 and then gives an almost flat penalty on [0, 0.8], the lambda_i of case 5 are of the order of 40 to 400 for k = 4, and the prior is informative.

**What is common to the five cases.**

* The specimens are selected with the same criteria in all cases (fraction of NRM and angle between the NRM and the pTRM, both set in the left panel), so that the cases can be compared. Two optional alteration thresholds, on (m4 - m3)/m4 and on (m1 - m3)/m1, are defined in `MSP_func.R` (`m4m3threshold` and `m1m3threshold`). They are set to -1, i.e. no alteration threshold is applied: the alteration enters only through the error dQ, which down-weights the altered specimens.
* In cases 2 to 5, the error dQ is computed with alpha = 0.5 and the weights are the inverse of dQ squared (in case 1 they are the inverse of dQ).
* In cases 2 to 5, Cook's distance is calculated again after a first optimization, with the same threshold (3 times the mean), and the final fit is done without the influential points.
* The confidence intervals are obtained by bootstrap. In cases 2 to 5, the optimization of alpha is repeated for each bootstrap sample (lambda is kept fixed), and the interval is read directly in the percentiles of the bootstrapped paleointensities.
* In the tables and in the plots, specimens that are not used in the final fit (*Bad* or *Outlier*) are shown with the global alpha (case 2) or with the ridge target (cases 3, 4 and 5); their alpha has not been optimized.

#### Choice of lambda (cases 3 and 4)

Because cases 3 and 4 have as many alpha values as specimens, they need a limit on how far each alpha may move away from its reference value (0.5 in case 3, the middle of [a_min, a_max] in case 4). This is the role of the ridge penalty, and lambda is its strength:

* a large lambda keeps all the alpha values close to their reference value (the result tends to a constant alpha);
* a small lambda lets the alpha values move freely toward the bounds of their interval, with the risk of fitting the noise.

Lambda is not chosen by hand but by leave-one-out cross-validation (LOOCV), which measures the capacity of the model to predict a specimen that it has not seen:

1. For each candidate value of lambda (from 0.5 to 2, every 0.05), each specimen is left out in turn.
2. The slope and the alpha values are estimated on the other specimens.
3. The alpha of the specimen left out is calculated from this slope (with the same ridge and the same bounds), and the weighted squared difference between its Q and the regression line is computed.
4. The sum of these errors over all the specimens is the LOOCV error of this lambda.
5. The lambda with the smallest LOOCV error is retained.

The search is done separately for case 3 and for case 4, since they do not have the same bounds nor the same reference values, so the two cases generally get two different lambda values. It uses the specimens that pass the selection criteria (before the removal of the influential points). The retained value is displayed in the left panel and used for the final fit and for all the bootstrap samples.

The tab *Lambda* shows the LOOCV error as a function of lambda, with a vertical line at the retained value:

* a clear minimum means that the choice of lambda is well defined;
* a flat curve means that the result depends little on lambda;
* if the minimum is on the first or on the last value of the range, the message "edge of the search range" appears in the left panel: the optimum may lie outside the range, and the result must be taken with caution.

The value of lambda is relative to the weights (inverse of dQ squared), so it depends on the data set: compare lambda values only within the same file.

#### Results specific to cases 3, 4 and 5

When case 3, 4 or 5 is selected, the following tabs are added to the menu *Analysis and Plots*:

* **Alpha values**: table with, for each specimen, the sample name, the laboratory field, a_min, a_max, the alpha value and the quality (Good, Outlier or Bad). In case 3, a_min and a_max are 0 and 0.8. In case 5, they are the intervals of the data file that define the prior (alpha itself is searched in [0, 0.8]). The table can be saved in a csv file.
* **Case 4 vs case 3** (cases 3 and 4): alpha values of case 4 as a function of those of case 3, specimen by specimen, with the 1:1 line and a color for each class of laboratory field. Only the specimens retained in the final fit of both cases are shown. The plot can be downloaded.
* **Case 5 vs case 3** (cases 3 and 5): same format as *Case 4 vs case 3*, with the alpha values of case 5 (prior on the individual intervals) as a function of those of case 3 (full range [0, 0.8]).
* **Lambda** (cases 3 and 4): plot of the search of lambda, described above. The plot can be downloaded.

#### Comparing the five cases (MSP-DSC)

The tab *Compare cases* (visible with the DSC protocol) tests whether the paleointensities of the five cases differ. The five cases are computed on the same specimens, so their paleointensities are strongly correlated: the overlap of the confidence intervals, or a test for independent estimates (t-test, ANOVA), is not appropriate. The comparison is a paired case-resampling bootstrap:

1. the same specimens are drawn with replacement for the five cases at each replicate (the specimens passing the selection criteria are the common set);
2. for each case, the slope is estimated again on the resampled specimens (in cases 2 to 5 the alpha values are optimized again; lambda of cases 3 and 4, the weights, the prior of case 5 and the exclusion of the influential specimens are kept fixed, as in the bootstrap of each case);
3. for the 10 pairs of cases, the difference PI(later case) - PI(earlier case) is recorded. The table gives the difference of the paleointensities obtained with the complete data set, its 95% percentile interval, and a two-sided bootstrap p-value, corrected for the 10 comparisons by the Holm method. A difference is declared significant when the Holm p-value is below 0.05 (the interval of the difference then does not include 0 in most cases).

The test is run only on demand: the five cases must have been computed first (select each case in *MSP-DSC: alpha strategy*; the line *Cases computed* shows which ones are available), then click on *Run the comparison of the five cases*. If the data file, the selection criteria, the number of bootstrap samples (or the width k of case 5, for case 5) are changed, the cases must be computed again. The number of paired replicates is the number of bootstrap samples of the left panel (a few tens of seconds for 10000). The tab also gives a plot of the differences and a csv file. A significant difference says that two cases give different paleointensities, not which one is better: the quality of a case must be judged with other criteria.

#### Reading the plots

* In the plots of Q as a function of the laboratory field, the paleointensity (Q = 0) is the solid vertical line and the colored band is its 95% confidence interval. The symbols show the status of each specimen: *Good* (blue circle), *Bad*, i.e. rejected by the selection criteria (orange square), and *Outlier*, i.e. rejected by Cook's distance (red triangle). The cross is the anchor point (0,-1) (MSP-FC and MSP-DSC).
* In the plot of the bootstrapped paleointensities, the solid line is the retained paleointensity and the colored area under the curve is the 95% confidence interval. The x-axis range can be changed with the slider under the figure, in order to compare the plots of the different protocols and cases.

### Influential points (outliers)

The influential data are detected and eliminated by means of Cook's distance. The threshold value is 3 times the mean value of Cook's distance.

```r
B.QDB$Cook <- cooks.distance(model_DB)
    threshold_Cook <- 3 * mean(B.QDB$Cook)
    DB_influential.data <- which(B.QDB$Cook > threshold_Cook)
    quality_DB <- rep("Good", length(B))
    quality_DB[DB_influential.data] <- "Outlier"
    B.QDB$quality <- quality_DB
    nDB_influential <- length(DB_influential.data)
```

### Confidence intervals

All the regression models are represented with bootstrapped 95% confidence intervals.

```r
# Bootstrapping for uncertainty
    data_sel <- subset(B.QDB, quality == "Good", select = c(x, y))
    reg.model <- boot(data_sel, reg.boot, R = R.rep())
    tt <- boot.ci(reg.model, conf = conf_level, type = "norm")
```

The number of bootstrap samples is set in the left panel (from 5000 to 15000). For MSP-DB, MSP-FC and MSP-DSC case 1, the interval is the normal-approximation bootstrap interval; for MSP-DSC cases 2 to 5, it is the interval between the percentiles of the bootstrapped paleointensities.

In the menu *Summary*, the 68% confidence interval is also provided.

### How to cite

Papers reporting results obtained with this code may cite:

> Camps P., Cych B., Antonio P., Raffray A. (2026) Optimization of the domain-state correction in Multiple-Specimen Paleointensity: Balancing prior and data. *Studia Geophysica et Geodaetica*, XX(X), pp-pp.

The source code and the test data set are archived on github.

### Cautionary notes

This software is provided on an "AS IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or implied.

### License

MIT.

### Contact

Code created by Pierre Camps, CNRS and Univ. Montpellier (Géosciences Montpellier). Please report bugs by opening an issue on this repository or by e-mail to pierre.camps [at] umontpellier.fr.

### References

* Dekkers M.J., Böhnel H.N. (2006) Reliable absolute palaeointensities independent of magnetic domain state. *Earth Planet. Sci. Lett.* 248, 508-517. https://doi.org/10.1016/j.epsl.2006.05.040
* Fabian K., Leonhardt R. (2010) Multiple-specimen absolute paleointensity determination: An optimal protocol including pTRM normalization, domain-state correction, and alteration test. *Earth Planet. Sci. Lett.* 297, 84-94. https://doi.org/10.1016/j.epsl.2010.06.006
* Monster M.W.L., de Groot L.V., Dekkers M.J. (2015) MSP-tool: A VBA-based software tool for the analysis of multispecimen paleointensity data. *Front. Earth Sci.* 3. https://doi.org/10.3389/feart.2015.00086

### Version history

**Shiny version 1.0.0, February 2018.** MSP-DSC Pierre Camps.

**Version 1.0.1.** Themes. Possibility to download the Cook and bootstrap figures.

**Version 2.0.0, September 2026.** MSP-DSC with four ways of handling alpha (cases 1 to 4); lambda selected by leave-one-out cross-validation for cases 3 and 4; new tabs *Alpha values*, *Case 4 vs case 3* and *Lambda*; new graphical style of the plots; number of bootstrap samples between 5000 and 15000; adjustable x-range of the bootstrap plot; Shapiro-Wilk test removed.

**Version 2.1.0, October 2026.** MSP-DSC case 5 (alpha per specimen in [0, 0.8], Gaussian prior centred on the middle of [a_min, a_max], lambda_i = 1/sigma_i^2 fixed by the intervals, adjustable width k); tab *Alpha values* extended to case 5; the R-squared of MSP-DSC case 1 is now the R-squared (it was the adjusted R-squared); a_min and a_max read by position follow the order of the header (a_max, a_min).

**Version 2.2.0, October 2026.** Tab *Compare cases*: paired bootstrap comparison of the five MSP-DSC cases (differences, 95% intervals, Holm-corrected p-values, plot and csv file), run on demand once the five cases have been computed.

**Version 2.2.1, October 2026.** New tab *Case 5 vs case 3* (alpha values of case 5 as a function of case 3); the plot of the tab *Compare cases* now has the graphic style and the font sizes of the bootstrap plot.
