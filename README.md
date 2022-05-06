# simplex

**simplex** is an **R** package for SIMS data processing that takes
  into account the *compositional* nature of isotopic data. This means
  that only the relative abundances of the isotopes carry the useful
  information.

**simplex** can be used online, offline and from the command line. The
easiest way to use the program is by visiting
[http://isoplotr.es.ucl.ac.uk/simplex/](http://isoplotr.es.ucl.ac.uk/simplex/).
If you would like to install your own online mirror of **simplex**,
the please follow the instructions provided [HERE](git.md). The
remainder of this README page provide instructions for offline use,
either using the browser-based Graphical User Interface, or the **R**
command line.

## Prerequisites

You must have **R** installed on your system (see
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
stand <- standard('Temora-t')
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

Stable isotope analysis of oxygen data:

```
dat <- read_data(f='*.asc',m=method('IGG-O'))
lr <- logratios(dat)
stand <- standard('NBS28')
cal <- calibration(lr,stand=stand,prefix="NBS28")
Cameca_oxygen <- calibrate(cal,exterr=TRUE)
tab <- data2table(Cameca_oxygen)
```

Built-in help can be obtained at the command prompt:

```
?logratios
```

To view a full list of functions:

```
help(package='simplex')
```

## Author

[Pieter Vermeesch](http://ucl.ac.uk/~ucfbpve/)

## License

This project is licensed under the GPL-3 License
