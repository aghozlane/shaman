#!/bin/sh
exec nohup /usr/bin/python3 /usr/bin/shaman_bioblend/shaman_bioblend.py  -w /srv/shiny-server/www/masque/ -s -d &
exec shiny-server >> /var/log/shiny-server.log 2>&1
exec nohup /usr/bin/Rscript /usr/bin/run_kronarshy.R&