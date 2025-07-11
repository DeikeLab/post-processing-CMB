# post-processing-clara

This repository contains code for the post-processing of Wind Wave simulation data. The results presented in the graphs correspond to simulations initialized with a broadbanded wave field. Each script is designed to handle different aspects of the post-processing, ranging from spectral analysis to the comparison of various parameters such as Reynolds numbers and velocities.

## File Descriptions

Some notebooks create DataFrames where I storage the data so it's easier to work with it each time. For most of them, the directory is loaded using:

'work_dir = f'/projects/DEIKE/nscapin/broadband_reorder/re{reA}_bo0{Bo}_P{kpHs}_uoc{uoc}_reW{reW}_L{maxLevel}/' 

So it's necessary to define 
* reA: Air Reynolds number
* Bo: Bond numer
* kpHs: initial condition wave stepness
* uoc: ratio friction velocity and phase velocity
* reW: Water Reynolds numer
* L: Max Level Refinement

- **2D_spectra.ipynb**
  
  
- **2Dspectra.ipynb**  
  Compares branches and includes automatic 2D analysis of spectra. Calculates
  - $E(k_x,k_y)$
  - $E(k,\theta)$
  - Energy spectrum $E(k)$.
  - Wave growth rate $\beta(k)$.
    (I think I will remove it once **2Dspectra_newpath.ipynb** is completely working)
- **2Dspectra_newpath.ipynb**  
  Focuses on correcting the path in 2D spectral data so we can use the path $ work_dir = f'/projects/DEIKE/nscapin/broadband_reorder/re{reA}_bo0{Bo}_P{kpHs}_uoc{uoc}_reW{reW}_L{maxLevel}/'$
  This eta_series is a post-processed eta series where we remove discontinuities in the surface due to drops and bubbles, the creation of this file can be found in **Total3Dspectra.ipynb**.
  
  - $E(k_x,k_y)$
  - $E(k,\theta)$
  - Energy spectrum $E(k)$.
  - Wave growth rate $\beta(k)$.
  - Energy behavior each mode.
  - Caracteristic time where decrease/increase change to increase/decrease
  - $\beta(k)$ for times smaller than this one.
  - Files created here:
    -  'eta_series_re{reA}_bo0{Bo}_P{kpHs}_uoc{uoc}_reW{reW}_L{maxLevel}_energies.csv'**, names=['i', 'k', 'F_integral_interval', 'time']'
    -  df_beta = ([df_growth, df_decay]): '/projects/DEIKE/cmartinb/betas/betas_re{reA}_bo0{Bo}_P{kpHs}_uoc{uoc}_reW{reW}_L{maxLevel}.csv' df_growth = pd.DataFrame({'Tipo': 'Growth', 'Beta': beta_growth, 'k': k_growth}) df_decay = pd.DataFrame({'Tipo': 'Decay','Beta': beta_decay,'k': k_decay })
    - df_datos = pd.DataFrame({'kp_minimos': kp_minimos, 'tiempos_minimos': tiempos_minimos}) saves in df_datos.to_csv(f'/projects/DEIKE/cmartinb/datos/datos_re{reA}_bo0{Bo}_P{kpHs}_uoc{uoc}_reW{reW}_L{maxLevel}.csv.csv', index=False).
    - df_growth_below_time = pd.DataFrame({'Tipo': 'Growth', 'Beta': beta_growth_below_time, 'k': k_growth_below_time}) df_decay_below_time = pd.DataFrame({'Tipo': 'Decay', 'Beta': beta_decay_below_time,'k': k_decay_below_time }), df_beta_below_time = pd.concat([df_growth_below_time, df_decay_below_time], axis=0).reset_index(drop=True) df_beta_below_time.to_csv(f'/projects/DEIKE/cmartinb/betas/betas_re{reA}_bo0{Bo}_P{kpHs}_uoc{uoc}_reW{reW}_L{maxLevel}_below_{tiempo_caracteristico}.csv', index=False)



- **3D_spectra.ipynb**
- 
- **3Dspectra.ipynb**  
  Energy spectra in 3D in time. We divide the data simulation in intervals of longitude 4 periods of the wave. Then we obtain the
  - $E(t,\omega, k_x, k_y)$
  - $E(t,\omega, k, \theta)$
  - $E(t,\omega, k)$ Dispersion relation space.
  - $E(t,\omega, k_fixed)$
  - 3D $E(t,\omega, k, \theta)$ graphs
- **3Dspectra_amplitudes.ipynb**  
  It does the 3D spectra analysis adding the calcul of the amplitudes RMS over the different intervals and wave numbers. 

- **3Dspectra_velocities.ipynb**  

- **Total3Dspectra.ipynb**  
  Similar to `3Dspectra.ipynb`, but not dividing the spectra in different time intervals, so it's just the 3D spectra analysis in $(k,\omega,t)$. The aim is for obtaining the $\omega$ so we can normalize the $\beta(k)$ from the 2D spectra analysis. This frequency omega is saved as a modified last column in the file 'betas_re{reA}_bo0{Bo}_P{kpHs}_uoc{uoc}_reW{reW}_L{maxLevel}.csv'

- **Totalenergy_all.ipynb**  

- **compare_Reynolds.ipynb**  
  Compares Reynolds numbers across all beta values using 3D total spectra.
  
- **compare_beta_diffRE.ipynb**  
  Handles proper path corrections in 2D data and compares different beta values, particularly for different Reynolds numbers.

- **compare_branches.ipynb**  
  Compares branches for different cases.

- **compare_growthcont.ipynb**  
  Post-processing of growth continuation data for different cases.
- **funciones.py**  
  Contain functions used along the notebooks.

- **high_RE.ipynb**  
  Similar to `3Dspectra.ipynb`, but focusing on changing to also high Reynolds, the pathused is the broadbanded_reorder.
- **kx_study.ipynb**  
  Directional spectra analysis in the $k_x$ direction.
- **spectrum.ipynb**  

- **velocity_profiles.ipynb**  
  Automatically generates velocity profiles from simulation data, extraction of the velocity for Air and Water, and its evolution in time. Extraction mean velocity in each time interval at the surface.
## How to Use

Each notebook is self-contained and can be run independently, depending on the specific analysis you wish to conduct. The scripts are designed to handle large datasets typical of Wind Wave simulations. Ensure that all necessary dependencies are installed and that paths to data files are correctly configured before running the notebooks. Some notebooks produce DataFrames that are used in other notebooks for compare different cases.
## Requirements

- Python 3.x
- Jupyter Notebook
