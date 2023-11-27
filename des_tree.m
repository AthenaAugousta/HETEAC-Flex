function [x] = des_tree(y, k, i)
%Decision tree - input depol and lidar ratio 
%              - output initial guess of state vector 
     if (y(i,k)>0) && (y(i,k)<0.11) && (y(i+1,k)>0) && (y(i+1,k)<40.1)
            disp ('(CS*)SEA SALT dominated mix')
            x = [0.05; 0.85; 0.05; 0.05];

     elseif (y(i,k)>0) && (y(i,k)<0.071) && (y(i+1,k)>60) 
            disp ('(FSA*)STR. ABS dominated mix')
            x = [0.85; 0.05; 0.05; 0.05];

     elseif (y(i,k)>0) && (y(i,k)<0.071) && (y(i+1,k)>39.9) && (y(i+1,k)<60.1)
            disp ('(FSNA*)WEAK. ABS dominated mix')
            x = [0.05; 0.05; 0.85; 0.05]; %default

     elseif (y(i,k)>0.18) && (y(i,k)<0.33) && (y(i+1,k)>10) && (y(i+1,k)<90) 
            disp ('(CNS*)DESERT DUST dominated mix')
            x = [0.0; 0.0; 0.0; 1.0 ]; %default

     elseif (y(i,k)>0.05) && (y(i,k)<0.20) && (y(i+1,k)>0) && (y(i+1,k)<40) 
            disp ('(CNS*/CS*)DESERT DUST-SEA SALT mix')
             x = [0.0; 0.7; 0.0; 0.3;]; %default

     elseif (y(i,k)>0.07) && (y(i,k)<0.19) && (y(i+1,k)>60) 
            disp ('(CNS*/FSA*)DESERT DUST-STR. ABS mix')
            x = [0.7; 0.0; 0.0; 0.3]; %default

     elseif (y(i,k)>0.07) && (y(i,k)<0.19) && (y(i+1,k)>39.9) && (y(i+1,k)<60.1) 
            disp ('(CNS*/FSNA*)DESERT DUST-WEAK. ABS mix')
            x = [0.0; 0.0; 0.7; 0.3]; %default

     elseif (y(i,k)>0) && (y(i,k)<0.051) && (y(i+1,k)>55) && (y(i+1,k)<65.1)
            disp ('(FSA*/FSNA*)WEAK. ABS-STR. ABS mix')
            x = [0.5; 0.0; 0.0; 0.5]; %default

     elseif (y(i,k)>0) && (y(i,k)<0.051) && (y(i+1,k)>40.1) && (y(i+1,k)<50)
            disp ('(FSNA*/CS*)WEAK. ABS-SEA SALT mix')
            x = [0.0; 0.5; 0.5; 0.0]; %default
     end     
end 