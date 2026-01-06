# README

This repository contains the code used to generate the data products and figures for the manuscript: "Entrained debris records regrowth of the Greenland Ice Sheet after the last interglacial"

It is organized into two main blocks:

## 1. MATLAB example scripts (OPR + MUSIC products)

The `RadarProcessing_Scripts/` directory contains an example script and parameter files that use MATLAB to run the Open Polar Radar (OPR) Toolbox and generate MUSIC radar data volumes.

These scripts are intended to:
- demonstrate the processing workflow used to create the MUSIC 3D volumes included in the associated Dataverse deposit
- provide a starting point for running OPR on similar datasets
- document the key processing steps and parameters used in this project

In general, these scripts assume that:
- you have MATLAB installed
- you have the Open Polar Radar (OPR) Toolbox available locally
- you have access to the radar input files required by OPR (available when running on KU cluster, which all OPR users have access to)

OPR Toolbox wiki (installation, product descriptions, workflows):
https://gitlab.com/openpolarradar/opr/-/wikis/home

## 2. Python notebooks (figure reproduction)

The `Manuscript_FigureGeneration/` directory contains Python notebooks used to reproduce the figures from the manuscript.

These notebooks are designed to:
- load the processed radar products (MUSIC volumes, standard products, and/or derived outputs)
- run analysis steps used in the paper
- regenerate each figure in a reproducible way (typically one notebook per figure or figure group)

The notebooks assume that you have access to the processed data files referenced in the manuscript. These data are available in the accompanying Harvard Dataverse dataset (see below).

## Data availability

The radar products used by these scripts and notebooks are archived in the (Amherst College) Harvard Dataverse dataset associated with the manuscript:

doi.org/10.7910/DVN/9K2J6R

