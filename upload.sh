#!/bin/bash
rsync -av --progress --delete .  lsst-dev.ncsa.illinois.edu:"/datasets/public_html/sim-data/beta_slair_surveys/"
#