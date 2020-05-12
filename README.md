# Tick abundance and disease prevalence using NEON data

## About
### Authors: Matt Bitters, Melissa Chen, Brendan Hobart, Wynne Moss
A space for data downloading, EDA, & model fitting for Tick project.

Use the `Issues` tab in this repository to post questions about data, make TO DO lists etc.

The goal of this project is to use NEON's tick abundance and tick pathogen data to forecast Lyme disease risk.

## Repo structure

**Important:** To view final project; clone the repo `Final_Proj_Content`. This folder is stand-alone and you will be able to reproduce our final analyses. 

Other Subdirectories:

* `data_raw` - raw data downloaded from NEON. We downloaded small mammal abundances, tick pathogen testing, and tick drag survey data. See scripts for `TickPathogen_01`,  `SMamAbundance_01`, `TickAbundance_01` for how data were downloaded and stacked. 

* `data_derived` - derived datasets (e.g. cleaned, sampled, reworked)

* `source` - any custom functions sourced by scripts

* `output` - temporary figures, tables, etc that you want to save

* `docs` - documentation and resources

* `Model_objects` - fitted model objects saved as .Rdata 

* `scratch` - leftovers and unfinished/broken scripts

## Scripts 
Our main folder contains a variety of scripts. Our data downloading, cleanup and EDA for each separate NEON dataset are labeled by the dataset (e.g. TickPathogen, SMam, or TickAbundance) and numbered in order of how they should be run. Tick abundance and pathogen data were merged (TickMerged).  Other scripts are not necessary but reflect other analyses we ran.