# HETEAC-Flex: optimal estimation method for aerosol typing

## 1. Structure of the input file  
The input file (.txt) should contain the following aerosol intensive properties:  
- Particle linear depolarization ratio at 355 nm along with the corresponding uncertainty 
- Lidar ratio at 355 nm along with the corresponding uncertainty 
- Extinction-related Ångström exponent along with the corresponding uncertainty
- Particle linear depolarization ratio at 532 nm along with the corresponding uncertainty
- Lidar ratio at 532 nm along with the corresponding uncertainty
- Backscatter-related color ratio (at the wavelength pair of 532/1064 nm) along with the corresponding uncertainty
 
 It is crucial that the order of the aforementioned properties is respected, in order for the algorithm to run properly.
 
 Example input files and a template can be found in <code><b>git_example_cases</b></code>.  The left column includes the measurements and the right column the corresponding uncertaint.

## 2. OEM routine 

The OEM routine (<code><b>HETEAC_Flex_main.m</b></code>) starts with an initial guess for the state vector. The initial guess is assigned automatically and it is the result of a decision tree (<code><b>des_tree.m</b></code>). Configuration settings can be easily adjusted (<code><b>config.m</b></code>).

## Further details in:
Floutsi, A. A., Baars, H., and Wandinger, U.: HETEAC-Flex: An optimal estimation method for aerosol typing based on lidar-derived intensive optical properties, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-1880, 2023. 

 ## Used in:
 - Heese, B., Floutsi, A. A., Baars, H., Althausen, D., Hofer, J., Herzog, A., Mewes, S., Radenz, M., and Schechner, Y. Y.: The vertical aerosol type distribution above Israel – 2 years of lidar observations at the coastal city of Haifa, Atmos. Chem. Phys., 22, 1633–1648, https://doi.org/10.5194/acp-22-1633-2022, 2022. 


