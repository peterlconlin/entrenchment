# Entrenchment

This repository holds the AvidaMT virtual machine setup, analyses scripts, and assembled supplemental material for the paper "Division of Labor Promotes the Entrenchment of Multicellularity" by Peter L. Conlin, Heather J. Goldsby, Eric Libby, Katherine G. Skocelas, William C. Ratcliff, Charles Ofria, and Benjamin Kerr.
The supplemental material can be viewed in its compiled form here: https://peterlconlin.github.io/entrenchment/index.html


## Navigation

This repository contains three main folders:

- **AvidaMT-VM**: Contains the files necessary to build a Singularity Vagrant Box containing AvidaMT and its dependencies.
- **AnalysisScripts**: Contains all R scripts used to analyze experiment data and generate figures.
- **SupplementalMaterials**: Holds the source files for the supplemental material, but the material itself is best viewed compiled: https://peterlconlin.github.io/entrenchment/index.html


## AvidaMT-VM

The Avida-VM folder contains all necessary files to build a [Singularity](https://sylabs.io/) [Vagrant Box](https://developer.hashicorp.com/vagrant/docs/boxes) capable of running [AvidaMT](https://github.com/kgskocelas/AvidaMT).

AvidaMT is a command line Linux program for conducting computational evolution experiments examining the conditions that promote the entrenchment of multicellularity and its mechanistic basis. First, the program reads in an extensive, user-defined set of environmental and organismal conditions. Next, it evolves a population of digital organisms under these conditions, recording the organisms’ activity, reproduction, and line of descent. This process is highly resource-intensive and designed to be run on a high-performance computing cluster.

For more details see: https://github.com/kgskocelas/AvidaMT


### Why Use a Singularity Vagrant Box?

AvidaMT requires a number of legacy dependencies and configuration with administrator privileges. Full-scale experiments with it are best run on a high-performance computing cluster, where most users are not administrators. As such, the easiest way to install and use it is via a virtual machine built on your local computer and then transferred to the computing cluster of your choice. 


### How to Make and Use Your AvidaMT Singularity Vagrant Box

Step-by-step directions for making and using your AvidaMT Singularity Vagrant Box are available at https://github.com/peterlconlin/entrenchment/AvidaMT-VM/INSTALL.md


## DOI for Supplemental Material

References to the supplemental material should be pointed to our DOI on Zenodo: [INSERT WHEN MADE]
