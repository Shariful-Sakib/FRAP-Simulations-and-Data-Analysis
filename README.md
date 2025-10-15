# FRAP-Simulations-and-Data-Analysis
- In this repository are Python programs that simulate Fluorescence Recovery After Photobleaching (FRAP) in spherical, cylindrical, and helical compartments with visuals and data analysis. 
  - Three versions of the simulation are given: standard, batch, and animated
    - The standard version performs one simulation per run. The batch version runs replicates of simulations for data collection. The animated version gives a visual version of the standard version.
- The compartment is defined by discretizing the shape into an array of 3D mid points
  - The spacing, S, between each mid point is calculated by determining the maximum tolerance on the error for any point particle such that the error, E = 4.6 nm (approximate length of a GFP molecule)
    - S = 2*sqrt((E+r)^2 - r^2)
      - This works well for very small S
  - The user can set the values for the different dimensions. If the shape variable is set to 0 then a sphere will be produced, shape= 1 is a cylinder/rod, and shape= 2 is a helical cell. Shape 1 defaults are preset for          _E. coli_ and shape 2 defaults are preset for _(Para)Magnetospirillum magneticum_, also known as AMB-1
  - The parametric equations for a helix is used to generate all the shapes where X(t) = t, Y(t) = A cos(2πt/λ), and Z(t) = A sin(2πt/λ).

- Point particles representing fluorophores are spawned randomly within the compartment within the **internal_radius** around the mid points

- The particles diffuse every timestep with the mean displacement centered at 0 microns and the standard deviation equal to sqrt(2Dt), where D is the **diffusion_constant** (typically 15 um^2/s) and t is the time step (1 ms)
  - There can be a hetergenous mix of immobile, slow, and fast particles if the user wishes

- During each time step, the particles' distance from the compartment is calculated by finding the closest mid point and if the distance from the closest mid point exceeds the internal_radius, then a collision is detected
  - A collision is handled by calculating the following d′ = d − 2(n − r)(n_unit/n) such that the new position is calculated by taking the current uncorrected collision and subtracting 2 times the distance from the membrane      (n-r) multiplied by the unit vector normal to the membrane at the collision point (n_unit/n).

- Microscopy data acquisition is simulated by analyzing the fluorescence (number of point particles) periodically with frame intervals (typically 30 ms). Two methods are used to analyze the fluorescence
  - Bleach region analysis where the count of fluorophores in the bleached and unbleached regions are tracked every frame interval. At the end of the simulation, the characteristic half-time of recovery is calculated by          fitting with an exponential decay function
  - Fourier profile analysis where the count of fluorophores along columns of the compartment (along the x-axis for convenience in this simulation) sectioned by a resolution of 85 nm is observed every frame interval. The         value is then fit with a cosine function i(x, t) = a_1(t) cos(πx/L), then at the end of the simulation, the amplitudes of that first Fourier mode, a_1(t), is then fit to an exponential function, a_{1,0}e^{−t/τ} + a_{1,i}  

# Automated Zeiss microscopy image processing and data analysis using PyImageJ 
Also included in the repository is a script that uses PyImageJ for automating Zeiss microscopy image processing and automated data analysis using PyImageJ.
-The functions and analysis used in the script can be adapted for any micrscope


Author: Shariful Sakib
