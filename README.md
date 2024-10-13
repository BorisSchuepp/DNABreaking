# Supplementary Code for "Sequence-specific, mechanophore-free mechanochemistry of DNA"

## Publication by: Johannes Hahmann, Boris N. Schüpp, Aman Ishaqat, Arjuna Selvakumar, Robert Göstl, Frauke Gräter, and Andreas Herrmann

## GitHub maintained by: Boris N. Schüpp, current version 1.1, last updated 11/10/2024 



# How to install?
Download by using\
**git clone https://github.com/BorisSchuepp/DNABreaking** \
in the desired directory.

# Hardware requirements
The software runs on a standard computer with a sufficient amount of RAM. The Next Generation Sequencing code runs for the 
major dataset in less then 10 minutes.  

# Software requirements
**Next Generation Sequencing**: Software will run in a common python3 enviroment (tested with Python 3.12). \
**Required python packages**: scipy, numpy, matplotlib, statistics, time, datetime.\
**Molecular Dynamics**: Software can be run in a Linux terminal (only code with .sh) with GROMACS (tested on version 2021.4) and python3 installed. 
Any other code (.py) runs with in a common python3 enviroment (tested with Python 3.12). \
**Required python packages**: numpy, os, sys, matplotlib, statistics, statsmodels

# What is included?

## Next Generation Sequencing Analysis Pipeline
The subdirectory NextGenerationSequencing contains the code to analize the results of ultrasonic breaking experiments on nicked dsDNA.
The oringal sequencing data on which the manuscript is based, can be found on Zenodo (https://doi.org/10.5281/zenodo.10683294). In order to use the code please copy any .fastq file from Zenodo to the directory RawData/704_AT_GC_400_1500_ATMi/. The script FilteringAlgorithm.py performs the analysis that yields the sequencing result presented in the manuscript. Any relevant processed data is written to ProcessedData/704_AT_GC_400_1500_ATMi. Statistical analysis of these datasets is performed by using the StatisticalAnalysis script and any results are written to ProcessedData/704_AT_GC_400_1500_ATMi_Statistics. Additionally, graphics/plots are created in the procedure and can be created on demand using PlottingNGS.py. Graphics can be found in the gaphics directory. 

To use the FilteringAlgorithm.py pipeline on different data, create a subdirectory in Graphics, RawData and ProcessedData. Afterwards 
save your .fastq files to the RawData directory. Create a run configuration in the RunConfigurations directory and change FilteringAlgorithm.py to 
refer to this run configuration. Refer to the README.txt file in Runconfiguration in order to create a valid run configuration. Plotting ultilies can be used
standalone by modifying the main function in PlottingNGS.py. Statistical analysis is hardcoded for our dataset and needs to be modified further to be used on a different set of data. 

## Molecular Dynamics
The molecular dynamics subdirectory includes scripts for the following four parts of our project:
1. Generation of nicked dsDNA structures. The script NickCreation.py is able to create nicked DNA structures on demand from an input .pdb or .top/.itp file. With this the input strcutures found in /MolecularDynamics/StrucutreGeneration/ have been created. 
2. Step-by-step instructions on how the MD simulations have been performed are found in the /MolecularDynamics/SimulationProcedure/ subdirectory. The resutling .xtc files for our dataset can be found on Zenodo.
3. Data exraction from .xtc trajectories. The script ProcessTrajectories.sh is able to extract relevant bond distances from the trajetories. The script can only be used in a Linux bash terminal mit GROMACS and python3 installed. \
Use the script by excecuting \
bash ProcessTrajectories.sh ../RawData/ ../StructureGeneration/InputStructures/ ../ProcessedData/GCContent/ \
or
bash ProcessTrajectoriesForceVariation.sh ../RawData/ ../StructureGeneration/InputStructures/ForceVariation/ ../ProcessedData/ForceVariation/
after saving the .xtc files from Zenodo in the RawData/ subdirectory.
4. Data analysis and visiualisation: The script DataAnalysis.py will perform in depth analyis on the data obtained in the previous step and can be run in any standard python3 distribution. Note, that this is very hardcoded to our dataset and likely won't be usefull on different MD data. 

# Questions
For any questions regarding the provided software please contact Boris N. Schüpp (boris.schuepp@mtl.maxplanckschools.de).
