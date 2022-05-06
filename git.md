## Setting up your own online *simplex* server with *git*

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

Copy following into a new file `/etc/systemd/system/simplex.service`:

```
[Unit]
Description=simplex
After=network.target

[Service]
Type=simple
User=wwwrunner
ExecStart=/usr/bin/Rscript -e simplex::daemon(2829)
Restart=always

[Install]
WantedBy=multi-user.target
```

Note we are setting `User=wwwrunner` to use our new user and we are
running it on port 2829.

Then to make **simplex** start on system boot type:

```sh
sudo systemctl enable simplex
```

Of course you can use other `systemctl` commands such as `start`, `stop`
and `restart` (to control whether it is running), and `disable` (to stop it
from running automatically on boot).

### Expose *simplex* with *nginx*

Ubuntu encourages you to put your configuration files in the
directory `/etc/nginx/sites-enabled`. If this directory is present
(and to be sure, you can check for a line saying `include
/etc/nginx/sites-enabled/*;` in the file `/etc/nginx/nginx.conf`) then
you need to add a file called `/etc/nginx/sites-enabled/default` with
the following contents:

```
server {
    listen 80 default_server;
    listen [::]:80 default_server;

    root /var/www/html;

    index index.html;

    server_name _;

    location /simplex/ {
        proxy_pass http://127.0.0.1:2829/;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
    }
}
```

If you already have a file called `/etc/nginx/sites-enabled/default`,
you will need to copy just the `location {...}` block into the
appropriate `server {...}` block in the existing file.

You can restart nginx to take the changes to its configuration we
made above with:

```sh
sudo systemctl restart nginx
```

If you need to start **simplex** now, call:

```sh
sudo systemctl start simplex
```

### Set up auto-updating

To ensure that **simplex** is up-to-date, it is a good idea to set up
auto-updating.

Put the following in a script `/usr/local/sbin/updateSimplex.sh`:

```sh
sudo -Hu wwwrunner Rscript -e \
     "remotes::install_github(repo=c('tim-band/shinylight','pvermees/simplex'),force=TRUE,lib='~/R')"
systemctl restart simplex
```

Ensure it is executable with:

```sh
sudo chmod a+rx /usr/local/sbin/updateSimplex.sh
```

One way to ensure that this script is regularly run is with **crontab**. First enter `sudo crontab -e` at the command prompt and then enter:

```
# Minute    Hour   Day of Month    Month            Day of Week           Command
# (0-59)   (0-23)    (1-31)    (1-12 or Jan-Dec) (0-6 or Sun-Sat)
    0        0         *             *                  0        /usr/local/sbin/updateSimplex.sh | /usr/bin/logger
```

which will automatically synchronise **shinylight** and **simplex** with **GitHub** on every Sunday.

You can force an update yourself by running the script as the `root` user:

```sh
sudo /usr/local/sbin/updateSimplex.sh
```

### Maintenance

You can view the logs from the various processes mentioned here
as follows:

Process | command for accessing logs
-----|-----
cron (including the update script) | `journalctl -eu cron`
systemD | `journalctl -e _PID=1`
simplex | `journalctl -eu simplex`
nginx | `journalctl -eu nginx`
nginx detail | logs are written into the `/var/log/nginx` directory

`journalctl` has many interesting options; for example `-r` to see
the most recent messages first, `-k` to see messages only from this
boot, or `-f` to show messages as they come in. The `-e` option
we have been using scrolls to the end of the log so that you are
looking at the most recent entries immediately.

If you need to set a custom timeout (say, to 6.5 seconds in this
example), change the `ExecStart` line in `isoplotr.service` like this:

```sh
ExecStart=/usr/bin/Rscript -e simplex::daemon(2829, timeout=6.5)
```