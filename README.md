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
install.packages('remotes')
```

## Installation

To install the current development version of simplex from Github, type:

```
remotes::install_github('tim-band/shinylight')
remotes::install_github('pvermees/isoplotr')
remotes::install_github('pvermees/simplex')
```

## Graphical User Interface (GUI)

The easiest way to run **simplex** is to start its GUI in a browser
tab, using the following **R** command:

```
simplex::simplex()
```

## Command Line Interface (CLI)

Load the **simplex** package into memory:

```
library(simplex)
```

The `process` function groups all the main data reduction steps,
including the drift correction, logratio calculation, and calibration:

```
dat <- read_data(f='SHRIMP.pd',m=method('GA-UPb'))
lr <- logratios(dat)
stand <- standard('Temora')
paired <- pairing(lr,stand=stand)
cal <- calibration(lr,stand=stand,pairing=paired,prefix="TEM")
result <- calibrate(cal,exterr=TRUE)
plot(result)
```

Extracting the results for 91500 zircon, saving the results as a
`.csv` file and plotting in `IsoplotR`:

```
samp <- subset(result,prefix='915')
tab <- data2table(samp)
write.csv(tab,file='~/Desktop/91500.csv',row.names=FALSE)
```

Stable isotope analysis of a built-in oxygen dataset:

```
m <- method('IGG-O')
s <- standard(preset="NBS28")
cd <- process(f='*.asc',m=m,stand=s)
del <- delta(cd)
tab <- data2table(del)
```

Built-in help can be obtained at the command prompt:

```
?process
```

To view a full list of functions:

```
help(package='simplex')
```

## Author

[Pieter Vermeesch](http://ucl.ac.uk/~ucfbpve/)

## License

This project is licensed under the GPL-3 License
