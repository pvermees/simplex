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

Load the **simplex** package into memory:

```
library(simplex)
```

The `process` function groups all the main data reduction steps,
including the drift correction, logratio calculation, and calibration:

```
m <- method('GA-UPb')
s <- standard(preset="Temora",prefix='TEM')
cd <- process(f='SHRIMP.pd',method=m,stand=s)
plot.calibration(cd)
```

Extracting the results for 91500 zircon, saving the results as a
`.csv` file and plotting in `IsoplotR`:

```
samp <- subset(cd,prefix='915')
tab <- data2table(samp)
write.csv(tab,file='~/Desktop/91500.csv',row.names=FALSE)
UPb <- IsoplotR::read.data('91500.csv',method='U-Pb',format=5)
IsoplotR::concordia(UPb,type=2,show.age=1)
```

Stable isotope analysis of a built-in oxygen dataset:

```
m <- method('IGG-O')
s <- standard(preset="NBS28")
cd <- process(f='*.asc',method=m,stand=s)
del <- delta(cd)
tab <- data2table(del)
```

Built-in help can be obtained at the command prompt:

```
?process
```

## Author

[Pieter Vermeesch](http://ucl.ac.uk/~ucfbpve/)

## License

This project is licensed under the GPL-3 License
