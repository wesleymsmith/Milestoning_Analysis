# Milestoning_Analysis
This repository contains notebooks detailing implementation, testing, and application of 1D milestoning methods. See "Margaliano, L. et. al 'Free Energy and Kinetics of Conformational Transitions from Voronoi Tessellated Milestoning with Restraining Potential' JCTC 2009. 5; DOI:10.1021/ct900279z" for an in depth discussion of the background and methodology.

Users who want to run the notebooks through google colab or other online notebook servers will need to clone the repo locally. They will also likely need to locally install the f90nml and bokeh packages for the notebooks that use them.

1. <B>Calculate_Bin_Escape_Matrix_1D.ipynb</B> - Tests binning procedure on a random data set and computes a corresponding 'escape matrix'

2. <B>Well_Test_Sim.ipynb</B> - Runs two testing simulations as milestoning test cases. The first is a non-stochastic simulation of a single particle in a simple harmonic well using a leapfrog algorithm to propogate dynamics. The second simulates the same single particle in a harmonic well but propogates using langevin dynamics to simulate adding thermal coupling and dampening. Running this will create the files Test_Simulation_Milestoning_Data.csv for the non-stochastic simulation and Test_Simulation_Milestoning_Data.lvd.csv for the langevin dynamics simulation.

3. <B>Analyze_Well_Test_Simulation.ipynb</B> - analyzes the results of the two test simulations. The coordinate data is first binned then the corresponding escape matrix is constructed. This escape matrix is then used to reconstruct the equilibrium probability for being in each window. This process is also applied to the MD data set at the end (this has been split into its own set of notebooks now for clarity)

4. <B>Extract_MD_Milestoning_Data.ipynb</B> - Loads the data from the MD milestoning simulation contained in test_md_data. The data is then collected and cleaned to produce a pair of data tables, 'Simulation_Milestone_Coordinate_Data.csv' and 'Simulation_Milestone_Restraint_Data.csv' for the milestoning coordinate data and milestoning softwall restraint parameters respectively. The data within is taken from the first four milestoning windows of a milestoning simulation for permeation of cyclic AMP through a connexin channel. The data tables produced by this notebook can then be analyzed with the 'Analyze_MD_Milestoning_Data.ipynb' notebook

5. <B>Analyze_MD_Milestoning_Data.ipynb</B> - Analyzes the MD simulation milestoning data to produce an escape matrix and corresponding equilibrium probabilities of being in each milestone window.

6. <B>analysis_functions.py</B> - contains a function that takes a data table for milestoning data containing the columns 'Window', 'Time', and 'X' along with the left and right edges of the milestoning restraints for each window. The function returns the corresponding escape matrix and equilibrium probability distribtuion. For a discussion of the theory behind this, see the analysis notebooks. An in depth discussion can also be found in the reference paper listed there in.

7. <B>fort_src</B> - fortran source code for binning 1D coordinate data and constructing a corresponding escape matrix. This code was provided by the author of the milestoning methodology papers referenced in the analysis notebooks

8. <B>test_md_data</B> - raw restraint coordinate and restraint parameter data from a set of molecular dynamics (MD) based milestoning simulations for permeation of a connexin channel by cyclic AMP
