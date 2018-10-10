PepBedR
=====================
<img src="https://img.shields.io/badge/R%3D-3.4.2-6666ff.svg" alt="minimal R version">

## Introduction

PepBedR: is a R package to process bed and bigBed files from peptide evidences.

## Pre-installation

For some functionalities, PepBedR runs internally some tools from UCSC. To install
the UCSC utilities, follow the instructions [here](https://github.com/ENCODE-DCC/kentUtils).

Alternatively, you could run from terminal the commands below to get installed 
the requiered binaries (tested on Ubuntu 16.04):

      cd $HOME
      
      # make directory
      mkdir kent-utils
      
      cd kent-utils
      
      # download binaries from UCSC
      curl -O http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
      curl -O http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedSummary
      curl -O http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedInfo
      
      # set exec permission
      chmod -R a+x .
      
      # Add the follow line to the .profile file. 
      # You could need to reboot your system to take effect of the changes.
      export PATH=$PATH:$HOME/kent-utils


Note that for Mac OS system, the binaries can be download from [here](http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/).

## Installing pepBed

First we need to install `devtools`:  

      install.packages("devtools")
      library(devtools)
   
Then, we can install the package using: 

      install_github("bigbio/PepBedR")
      library(PepBedR)


## Examples

Some examples to test the package are provided [here](https://github.com/bigbio/PepBedR/blob/master/docs/pepbed_examples_human.pdf). (additional work in process...)




