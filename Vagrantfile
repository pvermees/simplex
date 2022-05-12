# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.
Vagrant.configure("2") do |config|
  # The most common configuration options are documented and commented below.
  # For a complete reference, please see the online documentation at
  # https://docs.vagrantup.com.

  # Every Vagrant development environment requires a box. You can search for
  # boxes at https://vagrantcloud.com/search.
  config.vm.box = "ubuntu/focal64"

  # Disable automatic box update checking. If you disable this, then
  # boxes will only be checked for updates when the user runs
  # `vagrant box outdated`. This is not recommended.
  # config.vm.box_check_update = false

  # Create a forwarded port mapping which allows access to a specific port
  # within the machine from a port on the host machine. In the example below,
  # accessing "localhost:8080" will access port 80 on the guest machine.
  # NOTE: This will enable public access to the opened port
  config.vm.network "forwarded_port", guest: 80, host: 8080

  # Create a forwarded port mapping which allows access to a specific port
  # within the machine from a port on the host machine and only allow access
  # via 127.0.0.1 to disable public access
  # config.vm.network "forwarded_port", guest: 80, host: 8080, host_ip: "127.0.0.1"

  # Create a private network, which allows host-only access to the machine
  # using a specific IP.
  # config.vm.network "private_network", ip: "192.168.33.10"

  # Create a public network, which generally matched to bridged network.
  # Bridged networks make the machine appear as another physical device on
  # your network.
  # config.vm.network "public_network"

  # Share an additional folder to the guest VM. The first argument is
  # the path on the host to the actual folder. The second argument is
  # the path on the guest to mount the folder. And the optional third
  # argument is a set of non-required options.
  # config.vm.synced_folder "../data", "/vagrant_data"

  # Provider-specific configuration so you can fine-tune various
  # backing providers for Vagrant. These expose provider-specific options.
  # Example for VirtualBox:
  #
  # config.vm.provider "virtualbox" do |vb|
  #   # Display the VirtualBox GUI when booting the machine
  #   vb.gui = true
  #
  #   # Customize the amount of memory on the VM:
  #   vb.memory = "1024"
  # end
  #
  # View the documentation for the provider you are using for more
  # information on available options.

  # Enable provisioning with a shell script. Additional provisioners such as
  # Puppet, Chef, Ansible, Salt, and Docker are also available. Please see the
  # documentation for more information about their specific syntax and use.
  config.vm.provision "shell", inline: <<-"SHELL"
    useradd -mr wwwrunner
    apt-get update
    apt-get install -y nginx git r-base r-base-dev
    sudo -Hu wwwrunner sh <<- "RPACKAGES"
    mkdir ~/R
    echo R_LIBS_USER=~/R > ~/.Renviron
    Rscript -e "install.packages(pkgs='remotes',lib='~/R')"
    Rscript -e "remotes::install_github(repo=c('tim-band/shinylight','pvermees/simplex'),lib='~/R')"
RPACKAGES
    cat <<- "SYSTEMD" > /etc/systemd/system/simplex@.service
    [Unit]
    Description=simplex
    After=network.target
    [Service]
    Type=simple
    User=wwwrunner
    ExecStart=/usr/bin/Rscript -e "simplex::simplex(%i)"
    Restart=always
    [Install]
    WantedBy=multi-user.target
SYSTEMD
    cat <<- "SIMPLEXCTL" > /usr/local/sbin/simplexctl
    cmd=$1
    shift
    for p in $(seq 3901 3908)
    do systemctl --no-pager ${cmd} simplex@${p} $@
    done
SIMPLEXCTL
    chmod a+x /usr/local/sbin/simplexctl
    cat <<- "UPDATE" > /usr/local/sbin/updateSimplex.sh
    sudo -Hu wwwrunner Rscript -e "remotes::install_github(repo=c('tim-band/shinylight','pvermees/simplex'),force=TRUE,lib='~/R')"
    simplexctl restart
UPDATE
    chmod a+rx /usr/local/sbin/updateSimplex.sh
    simplexctl enable
    simplexctl start
    crontab <<- "CRONTAB"
    0 0 * * 0 /usr/local/sbin/updateSimplex.sh | /usr/bin/logger
CRONTAB
    cat <<- "NGINX" > /etc/nginx/sites-enabled/default
    upstream simplex {
      least_conn;
      server 127.0.0.1:3901;
      server 127.0.0.1:3902;
      server 127.0.0.1:3903;
      server 127.0.0.1:3904;
      server 127.0.0.1:3905;
      server 127.0.0.1:3906;
      server 127.0.0.1:3907;
      server 127.0.0.1:3908;
    }
    server {
      listen 80 default_server;
      listen [::]:80 default_server;
      root /var/www/html;
      index index.html;
      server_name _;
      location / {
        proxy_pass http://simplex/;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
      }
    }
NGINX
    systemctl restart nginx
  SHELL
end
