%%%%%%%%%%%% HETEAC-Flex %%%%%%%%%%%%
% Version 4.3: November 2023
% Athena Augusta Floutsi (floutsi@tropos.de)
% Update to the decision tree 
% Include new ext/depol a priori for asian dust 
% General optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: Read input file and initial set up 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

config

tic
filePattern = fullfile('./', '*.txt'); % fulfilment criterion for input files 
input_file = dir(filePattern); 
NumbMeas  = length(input_file);  % number of input files

for k = 1:NumbMeas
    baseFileName = input_file(k).name;
    fullFileName = fullfile('./', baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    fid = fopen(fullFileName,'rt');
    fmtstring = repmat('%f', 1, 4);
    data_cells = textscan(fid,fmtstring);
    fclose(fid);
    
    temp = data_cells{1};
    y_m(1:length(temp),k)= data_cells{1}; % read measurements from the input file
    err_m(1:length(temp),k)= data_cells{2}; % read corresponding measurement errors from the input file
         
    y(:,k)   = y_m(~isnan(y_m(:,k)),k);  
    err(:,k) = err_m(~isnan(err_m(:,k)),k);
    NumbParam(:,k) = length(y(:,k)); % number of available intensive optical parameters
    Se(:,:,k) = zeros(NumbParam(:,k),NumbParam(:,k)); % covariance matrix related to the lidar measurement errors
        for i = 1:NumbParam(k) 
            Se(i,i,k) = err(i,k)^2; 
        end
    
    x       = NaN(4,NumbMeas,NumbIt); % state vector
    x_tot   = NaN(NumbIt,NumbMeas); % sum of state vector elements
    F       = NaN(NumbParam(:,k),NumbMeas,NumbIt); % forward model 
   
    Sa(:,:,k)    = [Sa_val 0 0 0; 0 Sa_val 0 0; 0 0 Sa_val 0; 0 0 0 Sa_val]; % coviariance matrix of a priori -diagonal elements are set to zero
    tot_a355  = NaN(NumbIt,NumbMeas); % ext at 355nm (for forward model)
    tot_a532  = NaN(NumbIt,NumbMeas); % " " " 532nm " " "
    tot_a670  = NaN(NumbIt,NumbMeas); % " " " 670nm " " "
    tot_a1064 = NaN(NumbIt,NumbMeas); % " " " 1064nm " " "
    tot_b355  = NaN(NumbIt,NumbMeas); % bsc at 355nm (for forward model)
    tot_b532  = NaN(NumbIt,NumbMeas); % " " " 532nm " " "
    tot_b1064 = NaN(NumbIt,NumbMeas); % " " " 1064nm " " "
    
    J(:,k)    = {NaN(NumbParam(:,k),4,NumbIt)};
    Sdy(:,k)  = {NaN(NumbParam(:,k),NumbParam(:,k),NumbIt)}; % covariance between measurement and F(x)
    S_ret(:,k)= {NaN(4,4,NumbIt)}; % retrieval covariance 
    con2      = NaN(4,NumbMeas,NumbIt); % constraint

    crit(:,k)   = NaN(NumbIt,1); % convergence criterion
    C(:,k)      = NaN(NumbIt,1); % total costs
    C_apr(:,k)  = NaN(NumbIt,1); % a priori costs
    C_obs(:,k)  = NaN(NumbIt,1); % observational costs
    C_con2(:,k) = NaN(NumbIt,1); % constraint costs 

    C2_der1(:,k) = {zeros(4,1,NumbIt)}; % first derivative of 2nd constraint
    C2_der2(:,k) = {zeros(4,4,NumbIt)}; % second derivative of 2nd constraint
    gamma(:,k)   = NaN(NumbIt,1); % Levenberg-Marquardt parameter
    gamma(1,k)   = 2;    % initial value of gamma parameter 
end

for k = 1:NumbMeas
    for i = 1
        x(:,k,1) = des_tree(y, k, i); % state vector initialization based on decision tree
        x_tot(i,k) = sum(abs(x(:,k,i))); % total sum of components     
        x(:,k,i) = x(:,k,i)/x_tot(i,k);  % normalization of the state vector 

        % extinction and backscatter @ initial guess 
        tot_a355(i,k) = (x(i,k,i)*ap_355(i,1))+(x(i+1,k,i)*ap_355(i+1,1))+(x(i+2,k,i)*ap_355(i+2,1))+(x(i+3,k,i)*ap_355(i+3,1));
        tot_a532(i,k) = (x(i,k,i)*ap_532(i,1))+(x(i+1,k,i)*ap_532(i+1,1))+(x(i+2,k,i)*ap_532(i+2,1))+(x(i+3,k,i)*ap_532(i+3,1));
        tot_a670(i,k) = (x(i,k,i)*ap_670(i,1))+(x(i+1,k,i)*ap_670(i+1,1))+(x(i+2,k,i)*ap_670(i+2,1))+(x(i+3,k,i)*ap_670(i+3,1));
        tot_1064(i,k) = (x(i,k,i)*ap_1064(i,1))+(x(i+1,k,i)*ap_1064(i+1,1))+(x(i+2,k,i)*ap_1064(i+2,1))+(x(i+3,k,i)*ap_1064(i+3,1));
        tot_b355(i,k) = (x(i,k,i)*ap_355(i,2)+x(i+1,k,i)*ap_355(i+1,2)+x(i+2,k,i)*ap_355(i+2,2)+x(i+3,k,i)*ap_355(i+3,2));
        tot_b532(i,k) = (x(i,k,i)*ap_532(i,2)+x(i+1,k,i)*ap_532(i+1,2)+x(i+2,k,i)*ap_532(i+2,2)+x(i+3,k,i)*ap_532(i+3,2));
        tot_b1064(i,k) = (x(i,k,i)*ap_1064(i,2)+x(i+1,k,i)*ap_1064(i+1,2)+x(i+2,k,i)*ap_1064(i+2,2)+x(i+3,k,i)*ap_1064(i+3,2));
    
        for j = 1:4 % constrain 1 
            if x(j,k,i)<0
            x(j,k,i)=0;
            end
    
            if (x(j,k,i)>=0) && (x(j,k,i)<=1) % constrain 2 
            con2(j,k,i) = 0;
            else con2(j,k,i) = cst*abs(x(j,k,i))^3; 
            C2_der1{1,k}(j,k,i) = ((cst*abs(x(j,k,i)+h)^3)-(cst*abs(x(j,k,i)-h)^3))/2*h;
            C2_der2{1,k}(j,k,i) = ((cst*abs(x(j,k,i)+h)^3)-2*(cst*abs(x(j,k,i))^3)+(cst*abs(x(j,k,i)-h)^3))/h^2;
            end
        end 

        % Retrieval mode selection & initialization
        if NumbParam(:,k) == 2 && isnan(y_m(3,k)) && isnan(y_m(4,k)) && isnan(y_m(5,k)) && isnan(y_m(6,k))
            disp ('Performing basic 355 nm retrieval mode (1)') 
            retr_code(:,k) = 1;
            F(i,k,i)   = ((x(i,k,i)*ap_355(i,3)*ap_355(i,4))+(x(i+1,k,i)*ap_355(i+1,3)*ap_355(i+1,4))+(x(i+2,k,i)*ap_355(i+2,3)*ap_355(i+2,4))+(x(i+3,k,i)*ap_355(i+3,3)*ap_355(i+3,4)))/((x(i,k,i)*ap_355(i,4))+(x(i+1,k,i)*ap_355(i+1,4))+(x(i+2,k,i)*ap_355(i+2,4))+(x(i+3,k,i)*ap_355(i+3,4))); % depol355
            F(i+1,k,i) = (x(i,k,i)*ap_355(i,1)+x(i+1,k,i)*ap_355(i+1,1)+x(i+2,k,i)*ap_355(i+2,1)+x(i+3,k,i)*ap_355(i+3,1))/(x(i,k,i)*ap_355(i,2)+x(i+1,k,i)*ap_355(i+1,2)+x(i+2,k,i)*ap_355(i+2,2)+x(i+3,k,i)*ap_355(i+3,2)); %  lidar ratio 355
            
            J{1,k}(:,:,1) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
            Sdy{1,k}(:,:,i)  = Se(:,:,k)*(J{1,k}(:,:,i)*Sa(:,:,k)*J{1,k}(:,:,i)'+Se(:,:,k))^-1*Se(:,:,k); % Covariance of measurement and F(xa)

        elseif NumbParam(:,k) == 2 && isnan(y_m(1,k)) && isnan(y_m(2,k)) && isnan(y_m(3,k)) && isnan(y_m(6,k))
            disp ('Performing basic 532 nm retrieval mode (2)')
            retr_code(:,k) = 2;
            F(i,k,i)   = ((x(i,k,i)*ap_532(i,3)*ap_532(i,4))+(x(i+1,k,i)*ap_532(i+1,3)*ap_532(i+1,4))+(x(i+2,k,i)*ap_532(i+2,3)*ap_532(i+2,4))+(x(i+3,k,i)*ap_532(i+3,3)*ap_532(i+3,4)))/((x(i,k,i)*ap_532(i,4))+(x(i+1,k,i)*ap_532(i+1,4))+(x(i+2,k,i)*ap_532(i+2,4))+(x(i+3,k,i)*ap_532(i+3,4))); % depol532
            F(i+1,k,i) = (x(i,k,i)*ap_532(i,1)+x(i+1,k,i)*ap_532(i+1,1)+x(i+2,k,i)*ap_532(i+2,1)+x(i+3,k,i)*ap_532(i+3,1))/(x(i,k,i)*ap_532(i,2)+x(i+1,k,i)*ap_532(i+1,2)+x(i+2,k,i)*ap_532(i+2,2)+x(i+3,k,i)*ap_532(i+3,2)); % lidar ratio 532
    
            J{1,k}(:,:,1) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
            Sdy{1,k}(:,:,i)  = Se(:,:,k)*(J{1,k}(:,:,i)*Sa(:,:,k)*J{1,k}(:,:,i)'+Se(:,:,k))^-1*Se(:,:,k); % Covariance of measurement and F(xa)

        elseif NumbParam(:,k) == 3 && isnan(y_m(4,k)) && isnan(y_m(5,k)) && isnan(y_m(6,k))
            disp ('Performing enhanced 355 nm retrieval mode (3)')
            retr_code(:,k) = 3;
            F(i,k,i)   = ((x(i,k,i)*ap_355(i,3)*ap_355(i,4))+(x(i+1,k,i)*ap_355(i+1,3)*ap_355(i+1,4))+(x(i+2,k,i)*ap_355(i+2,3)*ap_355(i+2,4))+(x(i+3,k,i)*ap_355(i+3,3)*ap_355(i+3,4)))/((x(i,k,i)*ap_355(i,4))+(x(i+1,k,i)*ap_355(i+1,4))+(x(i+2,k,i)*ap_355(i+2,4))+(x(i+3,k,i)*ap_355(i+3,4))); % depol355
            F(i+1,k,i) = (x(i,k,i)*ap_355(i,1)+x(i+1,k,i)*ap_355(i+1,1)+x(i+2,k,i)*ap_355(i+2,1)+x(i+3,k,i)*ap_355(i+3,1))/(x(i,k,i)*ap_355(i,2)+x(i+1,k,i)*ap_355(i+1,2)+x(i+2,k,i)*ap_355(i+2,2)+x(i+3,k,i)*ap_355(i+3,2)); % lidar ratio 355
            F(i+2,k,i) = log(tot_a355(i,k)/tot_a532(i,k))/den; % extinction-related Angstrom exponent 
            
            J{1,k}(:,:,1) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
            Sdy{1,k}(:,:,i)  = Se(:,:,k)*(J{1,k}(:,:,i)*Sa(:,:,k)*J{1,k}(:,:,i)'+Se(:,:,k))^-1*Se(:,:,k); % Covariance of measurement and F(xa)

        elseif NumbParam(:,k) == 3 && isnan(y_m(1,k)) && isnan(y_m(2,k)) && isnan(y_m(3,k))
            disp ('Performing enhanced 532 nm retrieval mode (4)')
            retr_code(:,k) = 4;
            F(i,k,i) = ((x(i,k,i)*ap_532(i,3)*ap_532(i,4))+(x(i+1,k,i)*ap_532(i+1,3)*ap_532(i+1,4))+(x(i+2,k,i)*ap_532(i+2,3)*ap_532(i+2,4))+(x(i+3,k,i)*ap_532(i+3,3)*ap_532(i+3,4)))/((x(i,k,i)*ap_532(i,4))+(x(i+1,k,i)*ap_532(i+1,4))+(x(i+2,k,i)*ap_532(i+2,4))+(x(i+3,k,i)*ap_532(i+3,4))); % depol532
            F(i+1,k,i) = (x(i,k,i)*ap_532(i,1)+x(i+1,k,i)*ap_532(i+1,1)+x(i+2,k,i)*ap_532(i+2,1)+x(i+3,k,i)*ap_532(i+3,1))/(x(i,k,i)*ap_532(i,2)+x(i+1,k,i)*ap_532(i+1,2)+x(i+2,k,i)*ap_532(i+2,2)+x(i+3,k,i)*ap_532(i+3,2)); % lidar ratio 532
            F(i+2,k,i) = tot_b532(i,k)/tot_b1064(i,k); % backscatter-related color ratio 
 
            J{1,k}(:,:,1) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
            Sdy{1,k}(:,:,i)  = Se(:,:,k)*(J{1,k}(:,:,i)*Sa(:,:,k)*J{1,k}(:,:,i)'+Se(:,:,k))^-1*Se(:,:,k); % Covariance of measurement and F(xa)

        elseif NumbParam(:,k) == 4 && isnan(y_m(3,k)) && isnan(y_m(6,k))
            disp ('Performing simultaneous 355/532 nm retrieval mode (5)')
            retr_code(:,k) = 5;
            F(i,k,i)   = ((x(i,k,i)*ap_355(i,3)*ap_355(i,4))+(x(i+1,k,i)*ap_355(i+1,3)*ap_355(i+1,4))+(x(i+2,k,i)*ap_355(i+2,3)*ap_355(i+2,4))+(x(i+3,k,i)*ap_355(i+3,3)*ap_355(i+3,4)))/((x(i,k,i)*ap_355(i,4))+(x(i+1,k,i)*ap_355(i+1,4))+(x(i+2,k,i)*ap_355(i+2,4))+(x(i+3,k,i)*ap_355(i+3,4))); % depol355
            F(i+1,k,i) = (x(i,k,i)*ap_355(i,1)+x(i+1,k,i)*ap_355(i+1,1)+x(i+2,k,i)*ap_355(i+2,1)+x(i+3,k,i)*ap_355(i+3,1))/(x(i,k,i)*ap_355(i,2)+x(i+1,k,i)*ap_355(i+1,2)+x(i+2,k,i)*ap_355(i+2,2)+x(i+3,k,i)*ap_355(i+3,2)); % lidar ratio 355  
            F(i+2,k,i) = ((x(i,k,i)*ap_532(i,3)*ap_532(i,4))+(x(i+1,k,i)*ap_532(i+1,3)*ap_532(i+1,4))+(x(i+2,k,i)*ap_532(i+2,3)*ap_532(i+2,4))+(x(i+3,k,i)*ap_532(i+3,3)*ap_532(i+3,4)))/((x(i,k,i)*ap_532(i,4))+(x(i+1,k,i)*ap_532(i+1,4))+(x(i+2,k,i)*ap_532(i+2,4))+(x(i+3,k,i)*ap_532(i+3,4))); % depol532
            F(i+3,k,i) = (x(i,k,i)*ap_532(i,1)+x(i+1,k,i)*ap_532(i+1,1)+x(i+2,k,i)*ap_532(i+2,1)+x(i+3,k,i)*ap_532(i+3,1))/(x(i,k,i)*ap_532(i,2)+x(i+1,k,i)*ap_532(i+1,2)+x(i+2,k,i)*ap_532(i+2,2)+x(i+3,k,i)*ap_532(i+3,2)); % lidar ratio 532

            J{1,k}(:,:,1) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
            Sdy{1,k}(:,:,i)  = Se(:,:,k)*(J{1,k}(:,:,i)*Sa(:,:,k)*J{1,k}(:,:,i)'+Se(:,:,k))^-1*Se(:,:,k); % Covariance of measurement and F(xa)

        elseif NumbParam(:,k) == 6 
            disp ('Performing enhanced multiwavelength retrieval mode (6)')
            retr_code(:,k) = 6;
            F(i,k,i)   = ((x(i,k,i)*ap_355(i,3)*ap_355(i,4))+(x(i+1,k,i)*ap_355(i+1,3)*ap_355(i+1,4))+(x(i+2,k,i)*ap_355(i+2,3)*ap_355(i+2,4))+(x(i+3,k,i)*ap_355(i+3,3)*ap_355(i+3,4)))/((x(i,k,i)*ap_355(i,4))+(x(i+1,k,i)*ap_355(i+1,4))+(x(i+2,k,i)*ap_355(i+2,4))+(x(i+3,k,i)*ap_355(i+3,4))); % depol355
            F(i+1,k,i) = (x(i,k,i)*ap_355(i,1)+x(i+1,k,i)*ap_355(i+1,1)+x(i+2,k,i)*ap_355(i+2,1)+x(i+3,k,i)*ap_355(i+3,1))/(x(i,k,i)*ap_355(i,2)+x(i+1,k,i)*ap_355(i+1,2)+x(i+2,k,i)*ap_355(i+2,2)+x(i+3,k,i)*ap_355(i+3,2)); % lidar ratio 355    
            F(i+2,k,i) = log(tot_a355(i,k)/tot_a532(i,k))/den; % extinction-related Angstrom exponent
            F(i+3,k,i) = ((x(i,k,i)*ap_532(i,3)*ap_532(i,4))+(x(i+1,k,i)*ap_532(i+1,3)*ap_532(i+1,4))+(x(i+2,k,i)*ap_532(i+2,3)*ap_532(i+2,4))+(x(i+3,k,i)*ap_532(i+3,3)*ap_532(i+3,4)))/((x(i,k,i)*ap_532(i,4))+(x(i+1,k,i)*ap_532(i+1,4))+(x(i+2,k,i)*ap_532(i+2,4))+(x(i+3,k,i)*ap_532(i+3,4))); % depol532
            F(i+4,k,i) = (x(i,k,i)*ap_532(i,1)+x(i+1,k,i)*ap_532(i+1,1)+x(i+2,k,i)*ap_532(i+2,1)+x(i+3,k,i)*ap_532(i+3,1))/(x(i,k,i)*ap_532(i,2)+x(i+1,k,i)*ap_532(i+1,2)+x(i+2,k,i)*ap_532(i+2,2)+x(i+3,k,i)*ap_532(i+3,2)); % lidar ratio 532   
            F(i+5,k,i) = tot_b532(i,k)/tot_b1064(i,k); % backscatter-related color ratio
            
            J{1,k}(:,:,1) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
            Sdy{1,k}(:,:,i)  = Se(:,:,k)*(J{1,k}(:,:,i)*Sa(:,:,k)*J{1,k}(:,:,i)'+Se(:,:,k))^-1*Se(:,:,k); % Covariance of measurement and F(xa)
        end
        
        C(i,k)      = ((x(:,k,i)-x(:,k,i))'*Sa(:,:,k)^-1*(x(:,k,i)-x(:,k,i)))+((y(:,k)-F(:,k,i))'*Se(:,:,k)^-1*(y(:,k)-F(:,k,i)))+con2(1,k,i)+con2(2,k,i)+con2(3,k,i)+con2(4,k,i);
        C_apr(i,k)   = ((x(:,k,i)-x(:,k,1))'*Sa(:,:,k)^-1*(x(:,k,i)-x(:,k,1)));
        C_obs(i,k)   = ((y(:,k)-F(:,k,i))'*Se(:,:,k)^-1*(y(:,k)-F(:,k,i)));
        C_con2(i,k)  = con2(1,k,i)+con2(2,k,i)+con2(3,k,i)+con2(4,k,i);
        df(i,k)  = trace(J{1,k}(:,:,i)*Sa(:,:,k)*J{1,k}(:,:,i)'*(J{1,k}(:,:,i)*Sa(:,:,k)*J{1,k}(:,:,i)'+Se(:,:,k))^-1)+trace(Se(:,:,k)*(J{1,k}(:,:,i)*Sa(:,:,k)*J{1,k}(:,:,i)'+Se(:,:,k))^-1); % Degree of freedom of the measurement (df_signal+df_noise)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: Iterations, Convergence & Optimal solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method = 0; % Gauss-Newton (not yet included)
method = 1; % Levenberg-Marquardt (default) 

if method == true 
    for k = 1:NumbMeas
        for i = 2:NumbIt
            x(:,k,i) = x(:,k,i-1)+(((1+gamma(i-1,k))*Sa(:,:,k)^-1)+J{1,k}(:,:,i-1)'*Se(:,:,k)^-1*J{1,k}(:,:,i-1)+C2_der2{1,k}(:,:,i-1))^-1*(J{1,k}(:,:,i-1)'*Se(:,:,k)^-1*(y(:,k)-F(:,k,i-1))-(Sa(:,:,k)^-1*(x(:,k,i-1)-x(:,k,1)+C2_der1{1,k}(:,:,i-1)))); %optimal solution 
            x_tot(i,k) = sum(abs(x(:,k,i)));
            x(:,k,i) = x(:,k,i)/x_tot(i,k); % normalization 
                for j = 1:4 % constrains check up 
                    if x(j,k,i)<0
                        x(j,k,i)=0;
                    end 
   
                    if (x(j,k,i)>=0) && (x(j,k,i)<=1)
                        con2(j,k,i) = 0;
                    else con2(j,k,i) = cst*abs(x(j,k,i))^3;  
                        C2_der1{1,k}(j,k,i) = ((cst*abs(x(j,k,i)+h)^3)-(cst*abs(x(j,k,i)-h)^3))/2*h;
                        C2_der2{1,k}(j,k,i) = ((cst*abs(x(j,k,i)+h)^3)-2*(cst*abs(x(j,k,i))^3)+(cst*abs(x(j,k,i)-h)^3))/h^2;
                    end
                end 
            
                for j = 1
                    if retr_code(:,k) == 1
                        disp ('Performing basic 355 nm retrieval mode (1)')
                        J{1,k}(:,:,i) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
                
                        tot_a355(i,k) = (x(j,k,i)*ap_355(j,1))+(x(j+1,k,i)*ap_355(j+1,1))+(x(j+2,k,i)*ap_355(j+2,1))+(x(j+3,k,i)*ap_355(j+3,1));
                        tot_a532(i,k) = (x(j,k,i)*ap_532(j,1))+(x(j+1,k,i)*ap_532(j+1,1))+(x(j+2,k,i)*ap_532(j+2,1))+(x(j+3,k,i)*ap_532(j+3,1));
                        tot_a670(i,k) = (x(j,k,i)*ap_670(j,1))+(x(j+1,k,i)*ap_670(j+1,1))+(x(j+2,k,i)*ap_670(j+2,1))+(x(j+3,k,i)*ap_670(j+3,1));
                        tot_1064(i,k) = (x(j,k,i)*ap_1064(j,1))+(x(j+1,k,i)*ap_1064(j+1,1))+(x(j+2,k,i)*ap_1064(j+2,1))+(x(j+3,k,i)*ap_1064(j+3,1));
                        tot_b355(i,k) = (x(j,k,i)*ap_355(j,2)+x(j+1,k,i)*ap_355(j+1,2)+x(j+2,k,i)*ap_355(j+2,2)+x(j+3,k,i)*ap_355(j+3,2));
                        tot_b532(i,k) = (x(j,k,i)*ap_532(j,2)+x(j+1,k,i)*ap_532(j+1,2)+x(j+2,k,i)*ap_532(j+2,2)+x(j+3,k,i)*ap_532(j+3,2));
                        tot_b1064(i,k) = (x(j,k,i)*ap_1064(j,2)+x(j+1,k,i)*ap_1064(j+1,2)+x(j+2,k,i)*ap_1064(j+2,2)+x(j+3,k,i)*ap_1064(j+3,2));
        
                        F(j,k,i)     = ((x(j,k,i)*ap_355(j,3)*ap_355(j,4))+(x(j+1,k,i)*ap_355(j+1,3)*ap_355(j+1,4))+(x(j+2,k,i)*ap_355(j+2,3)*ap_355(j+2,4))+(x(j+3,k,i)*ap_355(j+3,3)*ap_355(j+3,4)))/((x(j,k,i)*ap_355(j,4))+(x(j+1,k,i)*ap_355(j+1,4))+(x(j+2,k,i)*ap_355(j+2,4))+(x(j+3,k,i)*ap_355(j+3,4)));
                        F(j+1,k,i)   = (x(j,k,i)*ap_355(j,1)+x(j+1,k,i)*ap_355(j+1,1)+x(j+2,k,i)*ap_355(j+2,1)+x(j+3,k,i)*ap_355(j+3,1))/(x(j,k,i)*ap_355(j,2)+x(j+1,k,i)*ap_355(j+1,2)+x(j+2,k,i)*ap_355(j+2,2)+x(j+3,k,i)*ap_355(j+3,2));           
                
                    elseif retr_code(:,k) == 2
                        disp ('Performing basic 532 nm retrieval mode (2)')
                        J{1,k}(:,:,i) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
                    
                        tot_a355(i,k) = (x(j,k,i)*ap_355(j,1))+(x(j+1,k,i)*ap_355(j+1,1))+(x(j+2,k,i)*ap_355(j+2,1))+(x(j+3,k,i)*ap_355(j+3,1));
                        tot_a532(i,k) = (x(j,k,i)*ap_532(j,1))+(x(j+1,k,i)*ap_532(j+1,1))+(x(j+2,k,i)*ap_532(j+2,1))+(x(j+3,k,i)*ap_532(j+3,1));
                        tot_a670(i,k) = (x(j,k,i)*ap_670(j,1))+(x(j+1,k,i)*ap_670(j+1,1))+(x(j+2,k,i)*ap_670(j+2,1))+(x(j+3,k,i)*ap_670(j+3,1));
                        tot_1064(i,k) = (x(j,k,i)*ap_1064(j,1))+(x(j+1,k,i)*ap_1064(j+1,1))+(x(j+2,k,i)*ap_1064(j+2,1))+(x(j+3,k,i)*ap_1064(j+3,1));
                        tot_b355(i,k) = (x(j,k,i)*ap_355(j,2)+x(j+1,k,i)*ap_355(j+1,2)+x(j+2,k,i)*ap_355(j+2,2)+x(j+3,k,i)*ap_355(j+3,2));
                        tot_b532(i,k) = (x(j,k,i)*ap_532(j,2)+x(j+1,k,i)*ap_532(j+1,2)+x(j+2,k,i)*ap_532(j+2,2)+x(j+3,k,i)*ap_532(j+3,2));
                        tot_b1064(i,k) = (x(j,k,i)*ap_1064(j,2)+x(j+1,k,i)*ap_1064(j+1,2)+x(j+2,k,i)*ap_1064(j+2,2)+x(j+3,k,i)*ap_1064(j+3,2));
        
                        F(j,k,i) = ((x(j,k,i)*ap_532(j,3)*ap_532(j,4))+(x(j+1,k,i)*ap_532(j+1,3)*ap_532(j+1,4))+(x(j+2,k,i)*ap_532(j+2,3)*ap_532(j+2,4))+(x(j+3,k,i)*ap_532(j+3,3)*ap_532(j+3,4)))/((x(j,k,i)*ap_532(j,4))+(x(j+1,k,i)*ap_532(j+1,4))+(x(j+2,k,i)*ap_532(j+2,4))+(x(j+3,k,i)*ap_532(j+3,4)));
                        F(j+1,k,i) = (x(j,k,i)*ap_532(j,1)+x(j+1,k,i)*ap_532(j+1,1)+x(j+2,k,i)*ap_532(j+2,1)+x(j+3,k,i)*ap_532(j+3,1))/(x(j,k,i)*ap_532(j,2)+x(j+1,k,i)*ap_532(j+1,2)+x(j+2,k,i)*ap_532(j+2,2)+x(j+3,k,i)*ap_532(j+3,2));    
            
                    elseif retr_code(:,k) == 3
                        disp ('Performing enhanced 355 nm retrieval mode (3)')
                        J{1,k}(:,:,i) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
                    
                        tot_a355(i,k) = (x(j,k,i)*ap_355(j,1))+(x(j+1,k,i)*ap_355(j+1,1))+(x(j+2,k,i)*ap_355(j+2,1))+(x(j+3,k,i)*ap_355(j+3,1));
                        tot_a532(i,k) = (x(j,k,i)*ap_532(j,1))+(x(j+1,k,i)*ap_532(j+1,1))+(x(j+2,k,i)*ap_532(j+2,1))+(x(j+3,k,i)*ap_532(j+3,1));
                        tot_a670(i,k) = (x(j,k,i)*ap_670(j,1))+(x(j+1,k,i)*ap_670(j+1,1))+(x(j+2,k,i)*ap_670(j+2,1))+(x(j+3,k,i)*ap_670(j+3,1));
                        tot_1064(i,k) = (x(j,k,i)*ap_1064(j,1))+(x(j+1,k,i)*ap_1064(j+1,1))+(x(j+2,k,i)*ap_1064(j+2,1))+(x(j+3,k,i)*ap_1064(j+3,1));
                        tot_b355(i,k) = (x(j,k,i)*ap_355(j,2)+x(j+1,k,i)*ap_355(j+1,2)+x(j+2,k,i)*ap_355(j+2,2)+x(j+3,k,i)*ap_355(j+3,2));
                        tot_b532(i,k) = (x(j,k,i)*ap_532(j,2)+x(j+1,k,i)*ap_532(j+1,2)+x(j+2,k,i)*ap_532(j+2,2)+x(j+3,k,i)*ap_532(j+3,2));
                        tot_b1064(i,k) = (x(j,k,i)*ap_1064(j,2)+x(j+1,k,i)*ap_1064(j+1,2)+x(j+2,k,i)*ap_1064(j+2,2)+x(j+3,k,i)*ap_1064(j+3,2));
        
                        F(j,k,i)     = ((x(j,k,i)*ap_355(j,3)*ap_355(j,4))+(x(j+1,k,i)*ap_355(j+1,3)*ap_355(j+1,4))+(x(j+2,k,i)*ap_355(j+2,3)*ap_355(j+2,4))+(x(j+3,k,i)*ap_355(j+3,3)*ap_355(j+3,4)))/((x(j,k,i)*ap_355(j,4))+(x(j+1,k,i)*ap_355(j+1,4))+(x(j+2,k,i)*ap_355(j+2,4))+(x(j+3,k,i)*ap_355(j+3,4)));
                        F(j+1,k,i)   = (x(j,k,i)*ap_355(j,1)+x(j+1,k,i)*ap_355(j+1,1)+x(j+2,k,i)*ap_355(j+2,1)+x(j+3,k,i)*ap_355(j+3,1))/(x(j,k,i)*ap_355(j,2)+x(j+1,k,i)*ap_355(j+1,2)+x(j+2,k,i)*ap_355(j+2,2)+x(j+3,k,i)*ap_355(j+3,2));
                        F(j+2,k,i)   = log(tot_a355(i,k)/tot_a532(i,k))/den;
            
                    elseif retr_code(:,k) == 4
                        disp ('Performing enhanced 532 nm retrieval mode (4)')
                        J{1,k}(:,:,i) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
                    
                        tot_a355(i,k) = (x(j,k,i)*ap_355(j,1))+(x(j+1,k,i)*ap_355(j+1,1))+(x(j+2,k,i)*ap_355(j+2,1))+(x(j+3,k,i)*ap_355(j+3,1));
                        tot_a532(i,k) = (x(j,k,i)*ap_532(j,1))+(x(j+1,k,i)*ap_532(j+1,1))+(x(j+2,k,i)*ap_532(j+2,1))+(x(j+3,k,i)*ap_532(j+3,1));
                        tot_a670(i,k) = (x(j,k,i)*ap_670(j,1))+(x(j+1,k,i)*ap_670(j+1,1))+(x(j+2,k,i)*ap_670(j+2,1))+(x(j+3,k,i)*ap_670(j+3,1));
                        tot_1064(i,k) = (x(j,k,i)*ap_1064(j,1))+(x(j+1,k,i)*ap_1064(j+1,1))+(x(j+2,k,i)*ap_1064(j+2,1))+(x(j+3,k,i)*ap_1064(j+3,1));
                        tot_b355(i,k) = (x(j,k,i)*ap_355(j,2)+x(j+1,k,i)*ap_355(j+1,2)+x(j+2,k,i)*ap_355(j+2,2)+x(j+3,k,i)*ap_355(j+3,2));
                        tot_b532(i,k) = (x(j,k,i)*ap_532(j,2)+x(j+1,k,i)*ap_532(j+1,2)+x(j+2,k,i)*ap_532(j+2,2)+x(j+3,k,i)*ap_532(j+3,2));
                        tot_b1064(i,k) = (x(j,k,i)*ap_1064(j,2)+x(j+1,k,i)*ap_1064(j+1,2)+x(j+2,k,i)*ap_1064(j+2,2)+x(j+3,k,i)*ap_1064(j+3,2));
        
                        F(j,k,i) = ((x(j,k,i)*ap_532(j,3)*ap_532(j,4))+(x(j+1,k,i)*ap_532(j+1,3)*ap_532(j+1,4))+(x(j+2,k,i)*ap_532(j+2,3)*ap_532(j+2,4))+(x(j+3,k,i)*ap_532(j+3,3)*ap_532(j+3,4)))/((x(j,k,i)*ap_532(j,4))+(x(j+1,k,i)*ap_532(j+1,4))+(x(j+2,k,i)*ap_532(j+2,4))+(x(j+3,k,i)*ap_532(j+3,4)));
                        F(j+1,k,i) = (x(j,k,i)*ap_532(j,1)+x(j+1,k,i)*ap_532(j+1,1)+x(j+2,k,i)*ap_532(j+2,1)+x(j+3,k,i)*ap_532(j+3,1))/(x(j,k,i)*ap_532(j,2)+x(j+1,k,i)*ap_532(j+1,2)+x(j+2,k,i)*ap_532(j+2,2)+x(j+3,k,i)*ap_532(j+3,2));    
                        F(j+2,k,i) = tot_b532(i,k)/tot_b1064(i,k);
                
                    elseif retr_code(:,k) == 5
                        disp ('Performing simultaneous 355/532 nm retrieval mode (5)')
                        J{1,k}(:,:,i) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
                    
                        tot_a355(i,k) = (x(j,k,i)*ap_355(j,1))+(x(j+1,k,i)*ap_355(j+1,1))+(x(j+2,k,i)*ap_355(j+2,1))+(x(j+3,k,i)*ap_355(j+3,1));
                        tot_a532(i,k) = (x(j,k,i)*ap_532(j,1))+(x(j+1,k,i)*ap_532(j+1,1))+(x(j+2,k,i)*ap_532(j+2,1))+(x(j+3,k,i)*ap_532(j+3,1));
                        tot_a670(i,k) = (x(j,k,i)*ap_670(j,1))+(x(j+1,k,i)*ap_670(j+1,1))+(x(j+2,k,i)*ap_670(j+2,1))+(x(j+3,k,i)*ap_670(j+3,1));
                        tot_1064(i,k) = (x(j,k,i)*ap_1064(j,1))+(x(j+1,k,i)*ap_1064(j+1,1))+(x(j+2,k,i)*ap_1064(j+2,1))+(x(j+3,k,i)*ap_1064(j+3,1));
                        tot_b355(i,k) = (x(j,k,i)*ap_355(j,2)+x(j+1,k,i)*ap_355(j+1,2)+x(j+2,k,i)*ap_355(j+2,2)+x(j+3,k,i)*ap_355(j+3,2));
                        tot_b532(i,k) = (x(j,k,i)*ap_532(j,2)+x(j+1,k,i)*ap_532(j+1,2)+x(j+2,k,i)*ap_532(j+2,2)+x(j+3,k,i)*ap_532(j+3,2));
                        tot_b1064(i,k) = (x(j,k,i)*ap_1064(j,2)+x(j+1,k,i)*ap_1064(j+1,2)+x(j+2,k,i)*ap_1064(j+2,2)+x(j+3,k,i)*ap_1064(j+3,2));
        
                        F(j,k,i)     = ((x(j,k,i)*ap_355(j,3)*ap_355(j,4))+(x(j+1,k,i)*ap_355(j+1,3)*ap_355(j+1,4))+(x(j+2,k,i)*ap_355(j+2,3)*ap_355(j+2,4))+(x(j+3,k,i)*ap_355(j+3,3)*ap_355(j+3,4)))/((x(j,k,i)*ap_355(j,4))+(x(j+1,k,i)*ap_355(j+1,4))+(x(j+2,k,i)*ap_355(j+2,4))+(x(j+3,k,i)*ap_355(j+3,4)));
                        F(j+1,k,i)   = (x(j,k,i)*ap_355(j,1)+x(j+1,k,i)*ap_355(j+1,1)+x(j+2,k,i)*ap_355(j+2,1)+x(j+3,k,i)*ap_355(j+3,1))/(x(j,k,i)*ap_355(j,2)+x(j+1,k,i)*ap_355(j+1,2)+x(j+2,k,i)*ap_355(j+2,2)+x(j+3,k,i)*ap_355(j+3,2));
                        F(j+2,k,i) = ((x(j,k,i)*ap_532(j,3)*ap_532(j,4))+(x(j+1,k,i)*ap_532(j+1,3)*ap_532(j+1,4))+(x(j+2,k,i)*ap_532(j+2,3)*ap_532(j+2,4))+(x(j+3,k,i)*ap_532(j+3,3)*ap_532(j+3,4)))/((x(j,k,i)*ap_532(j,4))+(x(j+1,k,i)*ap_532(j+1,4))+(x(j+2,k,i)*ap_532(j+2,4))+(x(j+3,k,i)*ap_532(j+3,4)));
                        F(j+3,k,i) = (x(j,k,i)*ap_532(j,1)+x(j+1,k,i)*ap_532(j+1,1)+x(j+2,k,i)*ap_532(j+2,1)+x(j+3,k,i)*ap_532(j+3,1))/(x(j,k,i)*ap_532(j,2)+x(j+1,k,i)*ap_532(j+1,2)+x(j+2,k,i)*ap_532(j+2,2)+x(j+3,k,i)*ap_532(j+3,2));            
                
                    elseif retr_code(:,k) == 6
                        disp ('Performing enhanced multiwavelength retrieval mode (6)')
                        J{1,k}(:,:,i) = jacobian(retr_code,k, h, den, ap_355, ap_532, ap_1064, x, i);
                    
                        tot_a355(i,k) = (x(j,k,i)*ap_355(j,1))+(x(j+1,k,i)*ap_355(j+1,1))+(x(j+2,k,i)*ap_355(j+2,1))+(x(j+3,k,i)*ap_355(j+3,1));
                        tot_a532(i,k) = (x(j,k,i)*ap_532(j,1))+(x(j+1,k,i)*ap_532(j+1,1))+(x(j+2,k,i)*ap_532(j+2,1))+(x(j+3,k,i)*ap_532(j+3,1));
                        tot_a670(i,k) = (x(j,k,i)*ap_670(j,1))+(x(j+1,k,i)*ap_670(j+1,1))+(x(j+2,k,i)*ap_670(j+2,1))+(x(j+3,k,i)*ap_670(j+3,1));
                        tot_1064(i,k) = (x(j,k,i)*ap_1064(j,1))+(x(j+1,k,i)*ap_1064(j+1,1))+(x(j+2,k,i)*ap_1064(j+2,1))+(x(j+3,k,i)*ap_1064(j+3,1));
                        tot_b355(i,k) = (x(j,k,i)*ap_355(j,2)+x(j+1,k,i)*ap_355(j+1,2)+x(j+2,k,i)*ap_355(j+2,2)+x(j+3,k,i)*ap_355(j+3,2));
                        tot_b532(i,k) = (x(j,k,i)*ap_532(j,2)+x(j+1,k,i)*ap_532(j+1,2)+x(j+2,k,i)*ap_532(j+2,2)+x(j+3,k,i)*ap_532(j+3,2));
                        tot_b1064(i,k) = (x(j,k,i)*ap_1064(j,2)+x(j+1,k,i)*ap_1064(j+1,2)+x(j+2,k,i)*ap_1064(j+2,2)+x(j+3,k,i)*ap_1064(j+3,2));
        
                        F(j,k,i)     = ((x(j,k,i)*ap_355(j,3)*ap_355(j,4))+(x(j+1,k,i)*ap_355(j+1,3)*ap_355(j+1,4))+(x(j+2,k,i)*ap_355(j+2,3)*ap_355(j+2,4))+(x(j+3,k,i)*ap_355(j+3,3)*ap_355(j+3,4)))/((x(j,k,i)*ap_355(j,4))+(x(j+1,k,i)*ap_355(j+1,4))+(x(j+2,k,i)*ap_355(j+2,4))+(x(j+3,k,i)*ap_355(j+3,4)));
                        F(j+1,k,i)   = (x(j,k,i)*ap_355(j,1)+x(j+1,k,i)*ap_355(j+1,1)+x(j+2,k,i)*ap_355(j+2,1)+x(j+3,k,i)*ap_355(j+3,1))/(x(j,k,i)*ap_355(j,2)+x(j+1,k,i)*ap_355(j+1,2)+x(j+2,k,i)*ap_355(j+2,2)+x(j+3,k,i)*ap_355(j+3,2));
                        F(j+2,k,i) = log(tot_a355(i,k)/tot_a532(i,k))/den;
                        F(j+3,k,i) = ((x(j,k,i)*ap_532(j,3)*ap_532(j,4))+(x(j+1,k,i)*ap_532(j+1,3)*ap_532(j+1,4))+(x(j+2,k,i)*ap_532(j+2,3)*ap_532(j+2,4))+(x(j+3,k,i)*ap_532(j+3,3)*ap_532(j+3,4)))/((x(j,k,i)*ap_532(j,4))+(x(j+1,k,i)*ap_532(j+1,4))+(x(j+2,k,i)*ap_532(j+2,4))+(x(j+3,k,i)*ap_532(j+3,4)));
                        F(j+4,k,i) = (x(j,k,i)*ap_532(j,1)+x(j+1,k,i)*ap_532(j+1,1)+x(j+2,k,i)*ap_532(j+2,1)+x(j+3,k,i)*ap_532(j+3,1))/(x(j,k,i)*ap_532(j,2)+x(j+1,k,i)*ap_532(j+1,2)+x(j+2,k,i)*ap_532(j+2,2)+x(j+3,k,i)*ap_532(j+3,2));    
                        F(j+5,k,i) = tot_b532(i,k)/tot_b1064(i,k);
                    end
                end
            
 % Total costs & cost breakdown
C(i,k)= ((x(:,k,i)-x(:,k,1))'*Sa(:,:,k)^-1*(x(:,k,i)-x(:,k,1)))+((y(:,k)-F(:,k,i))'*Se(:,:,k)^-1*(y(:,k)-F(:,k,i)))+con2(1,k,i)+con2(2,k,i)+con2(3,k,i)+con2(4,k,i); 
C_apr(i,k) = ((x(:,k,i)-x(:,k,1))'*Sa(:,:,k)^-1*(x(:,k,i)-x(:,k,1)));
C_obs(i,k) = ((y(:,k)-F(:,k,i))'*Se(:,:,k)^-1*(y(:,k)-F(:,k,i)));
C_con2(i,k) = con2(1,k,i)+con2(2,k,i)+con2(3,k,i)+con2(4,k,i);


        if C(i,k)>= C(i-1,k) % criterion for gamma change 
            gamma(i,k) = gamma(i-1,k)*10;
        elseif C(i,k)< C(i-1,k)
            gamma(i,k) = gamma(i-1,k)/2;
        end
 
  
Sdy{1,k}(:,:,i)  = Se(:,:,k)*(J{1,k}(:,:,i)*Sa(:,:,k)*J{1,k}(:,:,i)'+Se(:,:,k))^-1*Se(:,:,k); 
S_ret{1,k}(:,:,i) = (J{1,k}(:,:,i)'*Se(:,:,k)^-1*J{1,k}(:,:,i)+Sa(:,:,k)^-1)^-1; 
crit(i,k) =  (F(:,k,i)-F(:,k,i-1))'*Sdy{1,k}(:,:,i)^-1*(F(:,k,i)-F(:,k,i-1));

        if crit(i,k) <= 1/10*df(1,k) % stop iteration once criterion is met   
            it_ind(:,k) = i;
            fprintf('[Levenberg - Marquardt] Just finished iteration #%d\n', i);
            disp ('Convergence met')
            break    
        end           
            
        end
    end
elseif  method == false % GN
   disp ('Update the GN section')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: Statistical significance of retrieval (Chi-squared test) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:NumbMeas %chi_crit values for a 95% confidence interval
    for i = it_ind(k)
        chi_squared(:,k) = (F(:,k,i)-y(:,k)).'*Sdy{1,k}(:,:,i)^-1*(F(:,k,i)-y(:,k));        
        if abs(df(k)-2) <= 0.001
            chi_crit = 5.99; % 95%
        elseif abs(df(k)-3) <= 0.001
            chi_crit = 7.82; % 95%
        elseif abs(df(k)-4) <= 0.001
            chi_crit = 9.49; % 95%
        elseif abs(df(k)-6) <= 0.001
            chi_crit = 11.1; % 95%
        end 
        
        if chi_squared(k) > chi_crit 
            disp('The estimated solution IS NOT statistically significant at the 95% level. Consider a new initial guess')
            stat_sign(k,1) = 0;
        else 
            disp('The estimated solution IS statistically significant at the 95% level.')
            stat_sign(k,1) = 1;
        end 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 4: Important info - quick output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1: NumbMeas
    for w = it_ind(n)
        x_fin(n,:) = x(:,n,w); % optimal solution
        F_fin(n,:) = F(:,n,w); % values of forward model operator @ optimal solution
        S_ret_fin(:,:,n) = sqrt(S_ret{1,n}(:,:,w));
        S_ret_diag_fin(n,:) = diag(S_ret_fin(:,:,n)); % diagonal elements of retrieval error matrix
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 5: Additional products 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT: Define the measured extinction in the configuration file
for k = 1:NumbMeas
        for i = it_ind(k) 
            deno(:,k) = (x(1,k,i)*ap_355(1,1))+(x(2,k,i)*ap_355(2,1))+(x(3,k,i)*ap_355(3,1))+(x(4,k,i)*ap_355(4,1));
            V_abs(1,k) = x(1,k,i)*(a_355_meas(k)/deno(k)); % volume concentration - FSA
            V_abs(2,k) = x(2,k,i)*(a_355_meas(k)/deno(k)); % volume concentration - CS
            V_abs(3,k) = x(3,k,i)*(a_355_meas(k)/deno(k)); % volume concentration - FSNA
            V_abs(4,k) = x(4,k,i)*(a_355_meas(k)/deno(k)); % volume concentration - CNS

            N_abs(1,k) = V_abs(1,k)/(4/3*pi*r0v_fine^3*exp(-9/2*sigma_fine^2)); % number concentration - FSA
            N_abs(2,k) = V_abs(2,k)/(4/3*pi*(r0v_coarse)^3*exp(-9/2*(sigma_coarse)^2)); % number concentration - CS
            N_abs(3,k) = V_abs(3,k)/(4/3*pi*(r0v_fine)^3*exp(-9/2*(sigma_fine)^2)); % number concentration - FSNA
            N_abs(4,k) = V_abs(4,k)/(4/3*pi*(r0v_coarse)^3*exp(-9/2*(sigma_coarse)^2)); % number concentration - CNS
        
            A_fine(:,k)   = 4*pi*rA_fine^2*exp(-2*sigma_fine^2)*(N_abs(1,k)+N_abs(3,k)); % particle surface area - fine 
            A_coarse(:,k) = 4*pi*rA_coarse^2*exp(-2*sigma_coarse^2)*(N_abs(2,k)+N_abs(4,k)); % particle surface area - coarse 
            R_eff(:,k) = 3*(V_abs(1,k)+V_abs(2,k)+V_abs(3,k)+V_abs(4,k))/(A_fine(:,k)+A_coarse(:,k)); % effective radius 
            
            mR_mix(:,k) = (x(1,k,i)*mR(1,1)+x(2,k,i)*mR(2,1)+x(3,k,i)*mR(3,1)+x(4,k,i)*mR(4,1))/(x(1,k,i)+x(2,k,i)+x(3,k,i)+x(4,k,i)); % refractive index/imaginary 
            mI_mix(:,k) =(x(1,k,i)*mI(1,1)+x(2,k,i)*mI(2,1)+x(3,k,i)*mI(3,1)+x(4,k,i)*mI(4,1))/(x(1,k,i)+x(2,k,i)+x(3,k,i)+x(4,k,i)); % refractive index/real 
        end
end
% Error propagation for Monte Carlo
ns = 50000; % number of MC simulations (default)  

for k = 1:NumbMeas
        for i = it_ind(k)
            err_range(1,:,k) = [-S_ret_diag_fin(k,1) S_ret_diag_fin(k,1)]; % error range for FSA component
            err_range(2,:,k) = [-S_ret_diag_fin(k,2) S_ret_diag_fin(k,2)]; % ... CS
            err_range(3,:,k) = [-S_ret_diag_fin(k,3) S_ret_diag_fin(k,3)]; % ... FSNA
            err_range(4,:,k) = [-S_ret_diag_fin(k,4) S_ret_diag_fin(k,4)]; % ... CNS     
        
            rand_err_distr(1,:,k) = rand(ns,1)*range(err_range(1,:,k))+min(err_range(1,:,k));
            rand_err_distr(2,:,k) = rand(ns,1)*range(err_range(2,:,k))+min(err_range(2,:,k));
            rand_err_distr(3,:,k) = rand(ns,1)*range(err_range(3,:,k))+min(err_range(3,:,k));
            rand_err_distr(4,:,k) = rand(ns,1)*range(err_range(4,:,k))+min(err_range(4,:,k));
    
            MC(1,:,k) = rand_err_distr(1,:,k)+x(1,k,i); % ... FSA
            MC(2,:,k) = rand_err_distr(2,:,k)+x(2,k,i); % ... CS
            MC(3,:,k) = rand_err_distr(3,:,k)+x(3,k,i); % ... FSNA
            MC(4,:,k) = rand_err_distr(4,:,k)+x(4,k,i); % ... CNS
 
                for n = 1:ns % Filtering
                    if MC(1,n,k)<0 || MC(2,n,k)<0 || MC(3,n,k)<0 || MC(4,n,k)<0 % ensure that the MC simulated state vector is within the acceptable range 
                    MC(1,n,k) = 0;
                    MC(2,n,k) = 0;
                    MC(3,n,k) = 0;
                    MC(4,n,k) = 0;
                else 
                    MC(1,n,k) = MC(1,n,k);
                    MC(2,n,k) = MC(2,n,k);
                    MC(3,n,k) = MC(3,n,k);
                    MC(4,n,k) = MC(4,n,k);
                    end
                        rel_V_tot_MC(1,n,k) = MC(1,n,k)+MC(2,n,k)+MC(3,n,k)+MC(4,n,k); % ensure that the sum of the MC simulated state vector is not exceed 1
                    if rel_V_tot_MC(1,n,k)>1
                        MC(1,n,k) = 0;
                        MC(2,n,k) = 0;
                        MC(3,n,k) = 0;
                        MC(4,n,k) = 0;
                    end
                end 
        end
end
            
for k = 1:NumbMeas 
    for n = 1:ns 
            deno_MC(:,n,k) = (MC(1,n,k)*ap_355(1,1))+(MC(2,n,k)*ap_355(2,1))+(MC(3,n,k)*ap_355(3,1))+(MC(4,n,k)*ap_355(4,1));
            V_abs_MC(1,n,k) = MC(1,n,k)*(a_355_meas(k)/deno_MC(1,n,k)); % Volume & number concentration errors
            V_abs_MC(2,n,k) = MC(2,n,k)*(a_355_meas(k)/deno_MC(1,n,k));
            V_abs_MC(3,n,k) = MC(3,n,k)*(a_355_meas(k)/deno_MC(1,n,k));
            V_abs_MC(4,n,k) = MC(4,n,k)*(a_355_meas(k)/deno_MC(1,n,k));

            N_abs_MC(1,n,k) = V_abs_MC(1,n,k)/(4/3*pi*r0v_fine^3*exp(-9/2*sigma_fine^2));
            N_abs_MC(2,n,k) = V_abs_MC(2,n,k)/(4/3*pi*(r0v_coarse)^3*exp(-9/2*(sigma_coarse)^2));
            N_abs_MC(3,n,k) = V_abs_MC(3,n,k)/(4/3*pi*(r0v_fine)^3*exp(-9/2*(sigma_fine)^2));
            N_abs_MC(4,n,k) = V_abs_MC(4,n,k)/(4/3*pi*(r0v_coarse)^3*exp(-9/2*(sigma_coarse)^2));
   

            A_fine_MC(1,n,k)   = 4*pi*rA_fine^2*exp(-2*sigma_fine^2)*(N_abs_MC(1,n,k)+N_abs_MC(3,n,k)); % Effective radius errors
            A_coarse_MC(1,n,k) = 4*pi*rA_coarse^2*exp(-2*sigma_coarse^2)*(N_abs_MC(2,n,k)+N_abs_MC(4,n,k));
            R_eff_MC(1,n,k) = 3*(V_abs_MC(1,n,k)+V_abs_MC(2,n,k)+V_abs_MC(3,n,k)+V_abs_MC(4,n,k))/(A_fine_MC(1,n,k)+A_coarse_MC(1,n,k));

            mR_mix_MC(1,n,k) = (MC(1,n,k)*mR(1,1)+MC(2,n,k)*mR(2,1)+MC(3,n,k)*mR(3,1)+MC(4,n,k)*mR(4,1))/(MC(1,n,k)+MC(2,n,k)+MC(3,n,k)+MC(4,n,k)); % refractive index
            mI_mix_MC(1,n,k) = (MC(1,n,k)*mI(1,1)+MC(2,n,k)*mI(2,1)+MC(3,n,k)*mI(3,1)+MC(4,n,k)*mI(4,1))/(MC(1,n,k)+MC(2,n,k)+MC(3,n,k)+MC(4,n,k));
    end
end

for k = 1:NumbMeas 
    V_abs_MC_stat(1,1,k) = nanmean(V_abs_MC(1,:,k)); %FSA mean
    V_abs_MC_stat(2,1,k) = nanmean(V_abs_MC(2,:,k)); %CS mean
    V_abs_MC_stat(3,1,k) = nanmean(V_abs_MC(3,:,k)); %FSNA mean
    V_abs_MC_stat(4,1,k) = nanmean(V_abs_MC(4,:,k)); %CNS mean
    V_abs_MC_stat(1,2,k) = nanstd(V_abs_MC(1,:,k)); %FSA std
    V_abs_MC_stat(2,2,k) = nanstd(V_abs_MC(2,:,k)); %CS std
    V_abs_MC_stat(3,2,k) = nanstd(V_abs_MC(3,:,k)); %FSNA std
    V_abs_MC_stat(4,2,k) = nanstd(V_abs_MC(4,:,k)); %CNS std

    N_abs_MC_stat(1,1,k) = nanmean(N_abs_MC(1,:,k)); %FSA mean
    N_abs_MC_stat(2,1,k) = nanmean(N_abs_MC(2,:,k)); %CS mean
    N_abs_MC_stat(3,1,k) = nanmean(N_abs_MC(3,:,k)); %FSNA mean
    N_abs_MC_stat(4,1,k) = nanmean(N_abs_MC(4,:,k)); %CNS mean
    N_abs_MC_stat(1,2,k) = nanstd(N_abs_MC(1,:,k)); %FSA std
    N_abs_MC_stat(2,2,k) = nanstd(N_abs_MC(2,:,k)); %CS std
    N_abs_MC_stat(3,2,k) = nanstd(N_abs_MC(3,:,k)); %FSNA std
    N_abs_MC_stat(4,2,k) = nanstd(N_abs_MC(4,:,k)); %CNS std
                            
                            
    A_fine_MC_stat(1,1,k) = nanmean(A_fine_MC(1,:,k)); % mean 
    A_fine_MC_stat(1,2,k) = nanstd(A_fine_MC(1,:,k)); % std
    A_coarse_MC_stat(1,1,k) = nanmean(A_coarse_MC(1,:,k)); % mean 
    A_coarse_MC_stat(1,2,k) = nanstd(A_coarse_MC(1,:,k)); % std
    R_eff_MC_stat(1,1,k) = nanmean(R_eff_MC(1,:,k)); % mean 
    R_eff_MC_stat(1,2,k) = nanstd(R_eff_MC(1,:,k)); % std

    mR_mix_MC_stat(1,1,k) = nanmean(mR_mix_MC(1,:,k));
    mR_mix_MC_stat(1,2,k) = nanstd(mR_mix_MC(1,:,k));
    mI_mix_MC_stat(1,1,k) = nanmean(mI_mix_MC(1,:,k));
    mI_mix_MC_stat(1,2,k) = nanstd(mI_mix_MC(1,:,k));
end


toc
