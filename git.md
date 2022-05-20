# Setting up your own *simplex* server

Here are three ways to set up your own **simplex** server:

* On a VirtualBox (or other) virtual machine with Vagrant
* On a Debian-based Linux distribution with an installation file
* By hand

## Vagrant

Vagrant is a technology for setting up a virtual machine simply.

Launch **simplex** in a virtual machine with `vagrant up`. Stop it
with `vagrant halt`. If you have VirtualBox installed, everything should
work. If you have some other virtualization technology, you might
have to alter the `Vagrantfile` to make sure that your virtual
machine has enough memory to build `Rcpp`. The virtual machine
makes **simplex** available (with eight parallel instances)
on [https://localhost:8080/simplex/].

## .deb installation file

Those running a Debian-based Linux distribution (such as Ubuntu)
can install **simplex** with the `simplex.deb` installation
file.

Firstly you can run:

```sh
sudo apt install ./simplex.deb
```

to install **simplex** on your machine with all its dependencies.

### ensuring availability

**simplex** should now be available on [http://localhost/simplex/].
If not, your `nginx` installation is probably not including its location
block, which is installed in `/etc/nginx/app.d/simplex.conf`. Edit
the file `/etc/nginx/sites-available/default` and find the `server`
block containing the line `listen 80 default_server;` (or, if you have
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

### configuring simplex

**simplex** can support many users at once. However, R itself is
single-threaded, so if many users try to run their calculations at the
same time, each calculation will run in turn. This means that some
users might end up waiting a long time, even if more processing power
is available on the server.

This can be mitigated by running many different instances of
**simplex** in parallel.

You can configure **simplex** to run any number of instances
between 1 and 99 with the following commands:

```sh
sudo configureSimplex.sh 8
sudo simplexctl start
```

Here we are requesting 8 instances. The `configureSimplex.sh`
script stops simplex, so in the next line we start it again.

By default, **simplex** takes over port numbers 39101,
39102, 39103... (one port per instance you have). If these
interfere with other processes on your system, you can configure
**simplex** to use a different starting port number like this:

```sh
sudo configureSimplex.sh 8 4000
sudo simplexctl start
```

Here we are specifying ports 4000 to 4007.

### uninstallation

If installed with this method, you can uninstall **simplex** with
`sudo apt remove simplex`.

### updating the package (for developers)

The `.deb` file is defined by the files in the `simplex` directory.
If you change these files, you can update the `.deb` file with
`sudo dpkg -b simplex/`.

## By hand

Here is a way to set up a mirror on a Linux machine using the
following ingredients:

- Ubuntu
- nginx
- R
- git
- crontab

Instructions for offline use are provided in the main
[README](README.md) file.

### Install *nginx*, *R* and *git*

If these packages are not installed on your system already, then you
can add them with the following commands:

```sh
sudo apt-get install nginx git r-base r-base-dev
```

### Create a user to run *simplex*

It can be advantageous to have a non-human user running **simplex**
over the web so as to limit any damage should one behave badly. For
our purposes we will create one called `wwwrunner`:

```sh
sudo useradd -mr wwwrunner
```

### Set up *simplex* for this user

The version of **simplex** that gets run will be the version that our
new user `wwwrunner` has installed.

Install **simplex** for this user:

```sh
sudo -Hu wwwrunner sh -c "mkdir ~/R"
sudo -Hu wwwrunner sh -c "echo R_LIBS_USER=~/R > ~/.Renviron"
sudo -Hu wwwrunner Rscript -e "install.packages(pkgs='remotes',lib='~/R')"
sudo -Hu wwwrunner Rscript -e \
     "remotes::install_github(repo=c('tim-band/shinylight','pvermees/simplex'),lib='~/R')"
```

### Create a systemd service for *simplex*

Copy the file `simplex/DEBIAN/systemd/system/simplex@.service` to
`/systemd/system/`.

#### Copy the scripts onto your path

Copy the files in `simplex/DEBIAN/usr/local/sbin/` to
`/usr/local/sbin/`.

Now we can control **simplex**:

```sh
sudo simplexctl enable
sudo simplexctl start
sudo simplexctl status
```

### Expose *simplex* with *nginx*

Copy the file `simplex/DEBIAN/etc/nginx/app.d/simplex.conf` to
`/etc/nginx/app.d/`, and `simplex/DEBIAN/etc/nginx/conf.d/simplex.conf` to `/etc/nginx/conf.d/`.

As above in the section **ensuring availability**, ensure that your
default `server` block in `/etc/nginx/sites-available/default` contains the
line: `include /etc/nginx/app.d/*.conf;`. Now restart nginx:

```sh
sudo systemctl restart nginx
```

### Auto-updating

Copy the file `simplex/DEBIAN/etc/cron.weekly/simplex` to
`/etc/cron.weekly`.

This will automatically synchronise **shinylight** and **simplex** with **GitHub** on every Sunday.

You can force an update yourself by running the script as the `root` user:

```sh
sudo /usr/local/sbin/updateSimplex.sh
```

### Configuring simplex

See the section **configuring simplex** under **.deb installation file**
above to find out how to horizontally scale **simplex** for more parallelism.
