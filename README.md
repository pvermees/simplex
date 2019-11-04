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
data(package="simplex")
```

View the raw time resolved mass spectrometer data of the first SIMS
spot:

```
plot_timeresolved(Cameca[[1]],fit=TRUE)
```

To view further information about the **plot_timeresolved** function:

```
?plot_timeresolved
```

Plot a calibration curve:

```
stand <- subset_samples(dat=Cameca,prefix='Plesovice')
lr_stand <- logratios(dat=stand,c64=18.7)
fit <- calibration(lr_stand,dat=stand,plot=TRUE,disp=TRUE,omit=2)
```

where the second aliquot is omitted from the regression. Process
the samples and format as a flat data table:

```
samp <- subset_samples(dat=Cameca,prefix='Qinghu')
lr_samp <- logratios(dat=samp)
PbU <- calibrate(lr_samp,fit=fit,dat=samp,tst=c(337.13,0.18))
tab <- simplex2isoplotr(logPbU)
write.csv(tab,file='~/Desktop/Qinghu.csv')
```

## Author

[Pieter Vermeesch](http://ucl.ac.uk/~ucfbpve)

## License

This project is licensed under the GPL-3 License
