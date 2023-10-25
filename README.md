# HETEAC-Flex

## 1. Structure of the input file  
The input file (.txt) should contain the following intensive properties:  
- Particle linear depolarization ratio at 355 nm along with the corresponding uncertainty 
- Lidar ratio at 355 nm along with the corresponding uncertainty 
- Extinction-related Angstroem exponent along with the corresponding uncertainty
- Particle linear depolarization ratio at 532 nm along with the corresponding uncertainty
- Lidar ratio at 532 nm along with the corresponding uncertainty
- Backscatter-related color ratio (at the wavelength pair of 532/1064 nm) along with the corresponding uncertainty
 
 It is crucial that the order of the aforementioned properties is respected, in order for the algorithm to run properly.
 
 Example input files and a template can be found in <code><b>git_example_cases</b></code>.  The left column includes the measurements and the left column the corresponding uncertaint.

## 2. OEM routine 

The OEM routine starts with an initial guess for the state vector. The initial guess is assigned automatically and it is the result of a decision tree. 
