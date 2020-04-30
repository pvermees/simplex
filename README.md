# simplex

**simplex** is an **R** package for SIMS data processing that takes
  into account the *compositional* nature of isotopic data. This means
  that only the relative abundances of the isotopes carry the useful
  information.

## Prerequisites

You must have R installed on your system (see
[http://r-project.org](http://r-project.org)).  Additionally, to
install simplex from Github, you also need the **devtools** package.
This can be installed by typing the following code at the R command
line prompt:

```
install.packages('devtools')
```

## Installation

To install the current development version of simplex from Github, type:

```
library(devtools)
install_github('pvermees/isoplotr')
install_github('pvermees/simplex')
```

## Examples

Enter the following commands at the R prompt to start simplex and load
some example data into memory:

```
library(simplex)
data(Cameca,package="simplex")
```

View the raw time resolved mass spectrometer data of the first SIMS
spot:

```
plot_timeresolved(Cameca[[1]])
```

To view further information about the **plot_timeresolved** function:

```
?plot_timeresolved
```

Plot a calibration curve:

```
Ples <- standards(dat=Cameca,prefix='Plesovice',tst=c(337.13,0.18))
cal <- calibration(stand=Ples,oxide='UO2')
calplot(stand=Ples,fit=cal)
```

Apply the calibration to the Qinghu zircon measurements

```
samp <- unknowns(dat=Cameca,prefix='Qinghu')
UPb <- calibrate(samp,fit=cal)
```

Calculate the Pb-Pb ratios:

```
PbPb <- getPbLogRatios(samp)
```

Merge the results into one data object:

```
results <- mergecal(calUsamp,calPbsamp)
```

Export the results to a data table:

```
tab <- data2table(results)
write.csv(tab,file='Qinghu.csv',row.names=FALSE)
```

## Author

[Pieter Vermeesch](http://ucl.ac.uk/~ucfbpve)

## License

This project is licensed under the GPL-3 License
