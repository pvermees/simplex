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
remainder of this README page provides instructions for offline use,
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

## Vagrant

Launch **simplex** in a virtual machine with `vagrant up`. Stop it
with `vagrant halt`. If you have VirtualBox installed, everything should
work. If you have some other virtualization technology, you might
have to alter the `Vagrantfile` to make sure that your virtual
machine has enough memory to build `Rcpp`. The virtual machine
makes **simplex** available (with eight parallel instances)
on [https://localhost:8080/simplex/].

## .deb installation file

Those running a Debian-based Linux distribution (such as Ubuntu)
can install **simplex** with the `simplex_1.0-1.deb` installation
file.

Firstly you can run:

```sh
sudo apt install ./simplex_1.0-1.deb
```

This is not necessarily the latest simplex, but its first action is to
install the latest version from github.

**simplex** should now be available on [http://localhost/simplex/].
If not, your `nginx` installation is probably not including its location
block, which is installed in `/etc/nginx/app.d/simplex.conf`. Edit
the file `/etc/nginx/sites-available/default` and find the `server`
block containing the line `listen 80 default_server` (or, if you have
set up your own server, you can add this line to whichever one
you like). Add the line `include /etc/nginx/app.d/*.conf;`, so that
the block looks a bit like this:

```
server {
        listen 80 default_server;
        listen [::]:80 default_server;

        # vvv ADDED THIS LINE HERE TO ENABLE SIMPLEX
        include /etc/nginx/app.d/*.conf;

        root /var/www/html;

        # Add index.php to the list if you are using PHP
        index index.html index.htm index.nginx-debian.html;

        server_name _;

        location / {
                # First attempt to serve request as file, then
                # as directory, then fall back to displaying a 404.
                try_files $uri $uri/ =404;
        }
}
```

Now restart nginx with `sudo systemctl restart nginx` and see if
**simplex** appears at the above URL. You can turn **simplex**
off with `sudo simplexctl stop` and back on again with
`sudo simplexctl start`. Stop it from starting on boot with 
`sudo simplexctl disable` and enable starting on boot with
`sudo simplexctl enable`. You can see how all the instances are
doing with `sudo simplexctl status`.

If you want a newer offline version of the `.deb` file, you can
update it with `sudo dpkg -b simplex_1.0-1/`.

## Author

[Pieter Vermeesch](http://ucl.ac.uk/~ucfbpve/)

## License

This project is licensed under the GPL-3 License
