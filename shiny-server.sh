#!/bin/sh
# Make sure the directory for individual app logs exists
P=$(id -u) ; printf "run_as u$P;" >> /etc/shiny-server/shiny-server.conf
mkdir -p /var/log/shiny-server
chmod 666 /var/log/shiny-server
#exec shiny-server >> /var/log/shiny-server.log 2>&1
shiny-server