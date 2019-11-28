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
Ples <- subset_samples(dat=Cameca,prefix='Plesovice')
cal <- calibration(stand=Ples,oxide='UO2',tst=c(337.13,0.18))
calplot(dat=Ples,fit=cal)
```

Plot the time resolved data again, but now showing the data fit:

```
plot_timeresolved(Cameca[[1]],fit=cal)
```

Apply the calibration to the Qinghu zircon measurements

```
samp <- subset_samples(Cameca,prefix='Qinghu')
UPb <- calibrate(samp,fit=cal)
```

Export the results to a data table:

```
tab <- data2table(UPb)
write.csv(tab,file='Qinghu.csv',row.names=FALSE)
```

To plot the data on a concordia diagram with [IsoplotR](https://github.com/pvermees/IsoplotR):

```
library(IsoplotR)
Qinghu <- read.data('Qinghu.csv',method='U-Pb',format=5)
concordia(Qinghu)
```

Note that simplex is able to propagate the systematic uncertainties
associated with the calibration curve. However this feature has been
disabled in this release, until the necessary post-processing
functions have been added to IsoplotR.

## Author

[Pieter Vermeesch](http://ucl.ac.uk/~ucfbpve)

## License

This project is licensed under the GPL-3 License
