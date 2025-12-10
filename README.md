# Design of Permanent Magnet Synchronous Motor
In this repository, the MATLAB codes related to the design, FEM simulation, and post-processing analysis of a surface-mounted permanent magnet synchronous motor (PMSM) are provided. These scripts were developed as part of an academic study focusing on analytical sizing, electromagnetic field simulation, and torque/EMF quality evaluation of a PMSM.

## Machine Design

The machine is designed analytically based on rated power, speed, magnet properties, and stator/rotor geometry relations. The design outputs are validated through FEMM simulations under rated loading conditions.

The 2D finite element model is generated in FEMM, including stator slots, PM-rotor structure, materials, winding placement, and boundary definition. The FEM simulation provides electromagnetic torque and flux linkage evaluation over electrical angle.

## Matlab Scripts

In this folder, the following MATLAB scripts are provided:

- **s1_geometry_calc.m:** Performs the analytical sizing of the surface-mounted PMSM from user-defined ratings and material data. Computes rotor/stator dimensions, slot geometry, number of turns, yoke and PM dimensions, d–q inductances, leakage inductances, losses (iron, copper, mechanical, PM eddy, additional), rated operating point (current, load angle, efficiency, power factor), and plots the B–H curve, coefficient–flux-density curve, and torque vs. load angle. User interaction is limited to editing the 'Initial Data' block at the top of the script.

- **s2_FEMM_drawing.m:** Uses the analytical design calculations from Script 1 to automatically build the 2D PMSM model in FEMM, including stator slots, rotor with surface-mounted magnets, materials, boundary conditions, and three-phase windings (via dq→abc currents and slot mapping), then runs the magnetostatic solution, plots the flux density and reports the computed electromagnetic torque. User interaction is limited to the parameter block at the top (FEMM path, file name, problem settings, shaft diameter, electrical angle and slot label offsets).

- **s3_torque_vs_thetae.m:** Runs a current-angle sweep (rotor fixed) in FEMM using the analytical operating-point currents from Script 1 to evaluate electromagnetic torque vs. electrical injection angle. Modifies the phase currents in each step, solves the magnetostatic model, extracts torque, identifies the maximum-torque injection angle, saves the results to a .mat file, and plots the torque-angle curve. User interaction is limited to the parameter block at the top (sweep range, step size, output file settings).

- **s4_master_sweep_and_postproc.m:** Performs a rotor-position sweep in FEMM (using the optimal current angle from Script 3) to evaluate torque ripple, flux linkages, and back-EMF over one electrical period. Automatically rotates the rotor, updates phase currents (dq→abc), solves the magnetostatic model at each step, stores torque and flux data, computes back-EMF via numerical differentiation, saves all results to a .mat file, and prints sweep progress. User interaction is limited to the parameter block at the top (step size, sweep range, output file settings).

- **s5_torque_analysis.m:** Loads the FEMM sweep results from Script 4 and computes key torque ripple metrics (average, peak-to-peak, RMS, percentage values). Generates the torque vs. mechanical angle plot and performs FFT-based spectrum analysis to identify dominant ripple orders. Saves the torque KPIs and spectrum data to .mat and .csv files, with an optional automatic labeling of the highest peaks in the spectrum. User interaction is limited to the file names defined at the top.

- **s6_fluxlinkage_emf_analysis.m:** Loads the flux-linkage and EMF sweep results from Script 4 and computes phase and line-to-line RMS back-EMF, fundamental EMF magnitude and THD, and the EMF harmonic spectrum. Saves KPIs and spectrum data to .mat and .csv files and produces plots of flux linkage, phase EMF, line-to-line EMF, and EMF order spectrum. User interaction is limited to file names and peak-labeling settings at the top of the script.

## Installation

To run these scripts, you need:

- MATLAB (R2020a or newer recommended)
- FEMM 4.2 installed
- FEMM MATLAB API added to path

Example: addpath C:\femm42\mfiles

## Usage

run('s1_geometry_calc.m')
run('s2_FEMM_Geometry.m')
run('s3_torque_vs_thetae.m')
run('s4_master_sweep_and_postproc.m')
run('s5_torque_analysis.m')
run('s6_fluxlinkage_emf_analysis.m')

## Contact

The codes were developed by Md Tanbir Siddik Injam (mdtanbir.injam@aalto.fi) as part of an academic machine-design research work.
