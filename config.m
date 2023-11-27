% Configuration file for HETEAC-Flex v4_2 % 

% Maximum number of iterations %
NumbIt = 30;

% A priori %
ap_355  = [10.701 0.09123 0.024; 0.88604 0.05089 0.015; 9.61220 0.15778 0.033; 0.93219 0.01609 0.24]; % Saharan dust - HETEAC-Flex
% ap_355  = [10.701 0.09123 0.024; 0.88604 0.05089 0.015; 9.61220 0.15778 0.033; 0.93219 0.0215 0.25]; % Asian dust - HETEAC-Flex
% ap_355  = [10.701 0.09123 0.0; 0.88604 0.05089 0.0; 9.61220 0.15778 0.0; 0.93219 0.01609 0.24]; % Saharan dust - HETEAC
ap_532  = [6.4455 0.06875 0.024; 0.93604 0.04873 0.015; 5.0313 0.08476 0.033; 0.97321 0.0177 0.33]; % Saharan dust - HETEAC-Flex
% ap_532  = [6.4455 0.06875 0.024; 0.93604 0.04873 0.035; 5.0313 0.08476 0.033; 0.97321 0.0243 0.28]; % Asian dust - HETEAC-Flex
% ap_532  = [6.4455 0.06875 0.0; 0.93604 0.04873 0.0; 5.0313 0.08476 0.0; 0.97321 0.0177 0.33]; % Saharan dust - HETEAC
ap_670  = [4.13820 0.05799 0; 0.97774 0.04398 0; 3.13820 0.06128 0; 1.00740 0.04388 0.29450]; 
ap_1064 = [1.7638 0.03662 0; 1.0618 0.02984 0; 0.99217 0.03106 0; 1.0893 0.04799 0.3197];

for i = 1:length(ap_355)
    ap_355(i,4)  = ap_355(i,2)/(1+ap_355(i,3));
    ap_532(i,4)  = ap_532(i,2)/(1+ap_532(i,3));
    ap_670(i,4)  = ap_670(i,2)/(1+ap_670(i,3));
    ap_1064(i,4) = ap_1064(i,2)/(1+ap_1064(i,3));
end

% Coviariance matrix of a priori %
Sa_val = 0.05; 
cst    = 10^6;
h      = 0.001; % pertubation for numerical differentiation
den    = log(532/355);         

% Additional products 
a_355_meas = [128.8926; 41.0326; 116.9831]; % user defined - change accordingly

sigma_fine   = 0.53; % mode width - fine particles
sigma_coarse = 0.6; % mode width - coarse particles 
r0v_fine     = 0.1626; % mode radius/volume size distribbution - fine particles
r0v_coarse   = 2.32; % mode radius/volume size distribbution - coarse particles 
r0n_fine     = 0.07; % mode radius/number size distribbution - fine particles
r0n_coarse   = 0.788; % mode radius/number size distribbution - coarse particles
rA_fine      = 0.12276786;
rA_coarse    = 1.61889337;
mR(1,1)      = 1.50; % FSA - refractive index/real
mR(2,1)      = 1.45; % CS - refractive index/real
mR(3,1)      = 1.36; % FSNA - refractive index/real
mR(4,1)      = 1.54; % CNS - refractive index/real
mI(1,1)      = 0.043; % FSA - refractive index/imaginary 
mI(2,1)      = 1e-05; % CS - refractive index/imaginary 
mI(3,1)      = 1e-08; % FSNA - refractive index/imaginary 
mI(4,1)      = 0.006; % CNS  - refractive index/imaginary 