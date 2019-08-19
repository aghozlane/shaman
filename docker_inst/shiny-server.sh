#!/bin/sh
mkdir -p /srv/shiny-server/www/masque/todo /srv/shiny-server/www/masque/doing /srv/shiny-server/www/masque/error /srv/shiny-server/www/masque/done
chown shiny.shiny /srv/shiny-server/www/masque/todo /srv/shiny-server/www/masque/doing /srv/shiny-server/www/masque/error /srv/shiny-server/www/masque/done
exec nohup /usr/bin/python3 /usr/bin/shaman_bioblend/shaman_bioblend.py  -w /srv/shiny-server/www/masque/ -s -d &

# Make sure the directory for individual app logs exists
mkdir -p /var/log/shiny-server
chown shiny.shiny /var/log/shiny-server

exec shiny-server >> /var/log/shiny-server.log 2>&1

exec nohup /usr/bin/Rscript /usr/bin/run_kronarshy.R&
