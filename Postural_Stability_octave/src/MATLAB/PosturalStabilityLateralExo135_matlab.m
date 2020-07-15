%% Postural Stability in Lateral Box Transfer task (Exo Trials): COM deviation from centre of Base of Support

% Based on "Investigation and Analysis of the Effects of Manual Lifting and Carrying Activities on Postural and Gait Stability in Normal Subjects", Mohammed Alamoudi, University of Miami, 2017. 

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.

% For trials: 56, 62, 54, 55. 
% Code for Phases 1 (idxO1:idxF1),3 (idxF2:idxF3) and 5 (idxF4:idxF5)

% Five phases shown in figure 1:
    % 1: Subject picks up the box in the sagittal plane (idxO1:idxF1)
    % 2: Subject takes some steps to rotate to the frontal plane (idxF1:idxF2)
    % 3: Subject deposits the box in the frontal plane (idxF2:idxF3)
    % 4: Subject takes some steps to rotate back to the sagittal plane (idxF3:idxF4)
    % 5: Subject deposits the box in the sagittal plane (idxF4:idxF5)

% In general, AP, which describes fore-aft movement, is shown in x axis and
% ML, which corresponds to sideways movement, is reflected in y axis.
% This occurs in phases 1 and 5. 
% In phase 3 occurs the opposite. 
% In phases 2 and 4, components are calculated as distances between points (other script available to do those computations). 

% What to change to calculate PS in a different phase? 
    % interval of indexes within signal is cropped
    % for: until what number is the i. It should the substraction of the interval extremes + 1
    % phase 3*: extra changes are needed as fore-aft movement is now reflected in y axis are sideways movement is showed in x axis
        % Swap AP_dev and ML_dev: AP_dev = CBoS_ml - CoM_ml ; ML_dev = CBoS_ap - CBoS_ap
                    
clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica56_B.mat') 

% LTOE RTOE LHEE RHEE LFOO2 RFOO2 marker trajectories
    LTOE_x = LTOE(:,1)'; 
    LTOE_y = LTOE(:,2)';
    LTOE_z = LTOE(:,3)';
    
    RTOE_x = RTOE(:,1)';
    RTOE_y = RTOE(:,2)';
    RTOE_z = RTOE(:,3)';
    
    LHEE_x = LHEE(:,1)';
    LHEE_y = LHEE(:,2)';
    LHEE_z = LHEE(:,3)';
    
    RHEE_x = RHEE(:,1)';
    RHEE_y = RHEE(:,2)';
    RHEE_z = RHEE(:,3)';
    
    LFOO2_x = LFOO2(:,1)';
    LFOO2_y = LFOO2(:,2)';
    
    RFOO2_x = RFOO2(:,1)';
    RFOO2_y = RFOO2(:,2)';
    
% Centre of mass 
    % A/P fore-aft direction, thus x axis
    % M/L sideways direction, thus y axis
    
    CoM_ap = ModelData.Raw.(ModelOutput{13})(1,:);  
        
    CoM_ml = ModelData.Raw.(ModelOutput{13})(2,:);
    
%% Segmentation using LTOE, RTOE signals

% Left Toe

[pks, locs] = findpeaks(LTOE_z,'minPeakProminence',10,'MinPeakHeight',80);
TF1 = islocalmin(LTOE_z, 'FlatSelection', 'first');
idx = find(TF1);
flat = idx < locs(1);
idx_flat = find(flat);

flat2 = idx < locs(3)& idx > locs(2);
idx_flat2 = find(flat2);

flat3 = idx > locs(4);
idx_flat3 = find(flat3);

%e = locs < 500 | locs > 800;
%f = find(e);
%locs = locs(f);
%pks = pks(f);

% Right Toe 

[pks2, locs2] = findpeaks(RTOE_z,'minPeakProminence',10);

TF2 = islocalmin(RTOE_z, 'FlatSelection', 'first');
idx2 = find(TF2);
flat4 = idx2 < locs2(1);
idx_flat4 = find(flat4);

flat5 = idx2 > locs2(length(locs2));
idx_flat5 = find(flat5);
flat6 = idx2 > locs2(3);
idx_flat6 = find(flat6);

% Indexes definiton for signal segmentation
idxO1 = 1;
idxF1 = idx2(idx_flat4(length(idx_flat4))); %idxO2
idxF2 = idx2(idx_flat6(1)); %idxO3
idxF3 = idx(idx_flat2(length(idx_flat2))); %idxO4
idxF4 = idx(idx_flat3(3)); %idxO5
idxF5 = length(LTOE_z);

% Visualization of the segmentation
figure(1)
t= 0:(length(LTOE)-1);
plot(t, zscore(LTOE_z));
hold on
plot(t,zscore(RTOE_z));
plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxO1);
plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF2);
plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF3);
plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF5);
plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxO1);
plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF1);
plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF4);
plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF5);
hold off
xlabel('Time in seconds');
ylabel('Angular displacement in degrees');
legend('Left Toe Z', 'Right Toe Z');
title('Segmentation Exoskeleton Trials')


%% Phase...

    LTOE_x = LTOE_x(idxF2:idxF3);
    LTOE_y = LTOE_y(idxF2:idxF3);
    
    RTOE_x = RTOE_x(idxF2:idxF3);
    RTOE_y = RTOE_y(idxF2:idxF3);
    
    LHEE_x = LHEE_x(idxF2:idxF3);
    LHEE_y = LHEE_y(idxF2:idxF3);
    
    RHEE_x = RHEE_x(idxF2:idxF3);
    RHEE_y = RHEE_y(idxF2:idxF3);
    
    LFOO2_x = LFOO2_x(idxF2:idxF3);
    LFOO2_y = LFOO2_y(idxF2:idxF3);
    
    RFOO2_x = RFOO2_x(idxF2:idxF3);
    RFOO2_y = RFOO2_y(idxF2:idxF3);
    
    CoM_ap = CoM_ap(idxF2:idxF3);
        
    CoM_ml = CoM_ml(idxF2:idxF3);

%% Vertex definition using feet markers coordinates in counterclockwise order: 
% RFOO2 RTOE LTOE LFOO2 LHEE RHEE 

X = [];
for i = 1:(idxF3-idxF2+1) % should be changed when calculating another phase
    for j = 1:6 % number of markers that define the BoS
        switch j
            case 1
            X(i,j)= RFOO2_x(i);
            case 2
            X(i,j)= RTOE_x(i);
            case 3
            X(i,j)= LTOE_x(i);
            case 4
            X(i,j)= LFOO2_x(i);
            case 5
            X(i,j)= LHEE_x(i);
            case 6
            X(i,j)= RHEE_x(i);
        end
    end
end

Y = [];

for i = 1:(idxF3-idxF2+1) % should be changed when calculating another phase
    for j = 1:6 % number of markers that define the BoS
        switch j
            case 1
            Y(i,j) = RFOO2_y(i);
            case 2
            Y(i,j)= RTOE_y(i);
            case 3
            Y(i,j)= LTOE_y(i);
            case 4
            Y(i,j)= LFOO2_y(i);
            case 5
            Y(i,j)= LHEE_y(i);
            case 6
            Y(i,j)= RHEE_y(i);
        end
    end
end

%% Centre of Base of Support calculation

for i= 1:(idxF3-idxF2+1)   % should be changed when calculating another phase
    sumaN = 0;
    sumaD= 0;
    sumaN2 = 0;
    
    pgon = polyshape([X(i,3) X(i,5) X(i,6)], [Y(i,3) Y(i,5) Y(i,6)]);
    pgon2 = polyshape([X(i,6) X(i,3) X(i,2)], [Y(i,6) Y(i,3) Y(i,2)]);
    pgon3 = polyshape([X(i,3) X(i,5) X(i,4)], [Y(i,3) Y(i,5) Y(i,4)]);
    pgon4 = polyshape([X(i,6) X(i,2) X(i,1)], [Y(i,6) Y(i,2) Y(i,1)]);
     
    area1 = area(pgon); % triangle 1
    [centroid1_x, centroid1_y] = centroid(pgon);

    area2 = area(pgon2); % triangle 2
    [centroid2_x, centroid2_y] = centroid(pgon2);
 
    area3 = area(pgon3); % triangle 3
    [centroid3_x, centroid3_y] = centroid(pgon3);

    area4 = area(pgon4); % triangle 4
    [centroid4_x, centroid4_y] = centroid(pgon4);

    Area = [area1 area2 area3 area4];
    Centroid = [centroid1_x centroid1_y; centroid2_x centroid2_y; centroid3_x centroid3_y; centroid4_x centroid4_y];
       
    for k = 1:4
        sumaN = sumaN + (Centroid(k,1)*Area(k));  %x
        sumaD = sumaD + Area(k);
    end

    for k = 1:4
        sumaN2 = sumaN2 + (Centroid(k,2)*Area(k));  %x
    end

    CBoS_ap(1,i) = sumaN/sumaD;   %x 
    CBoS_ml(1,i) = sumaN2/sumaD; %y
end

%% "Postural stability in the A/P and M/L directions are measured as the normalized distance of the deviation of the CoM with respect to the vertical and 
% horizontal distance between the CBoS and the edge of the BoS,respectively."

% V: Distance from CBoS to the straight line formed by the points LTOE and RTOE
        % A(x1,y1) is LTOE_x(i) and LTOE_y(i), where i indicates the frame
        % B(x2, y2) is RTOE_x(i) and RTOE_y(i)
             
        % CBoS_ap(1,i) is the x component
        % CBos_ml(1,i) is the y component
        
V =[];
for i = 1:(idxF3-idxF2+1)  % should be changed when calculating another phase
    x1 = LTOE_x(i); 
    y1 = LTOE_y(i);
    x2 = RTOE_x(i);
    y2 = RTOE_y(i);
    
    A = -((y2-y1)/(x2-x1));        % calculating general or implicit equation of the line
    B = 1;
    C = -y1 + ((y2-y1)/(x2-x1)*x1);

    if A<0
    A = A*(-1);
    B = B*(-1);
    C = C*(-1);
    end
    
    V(1,i) = (abs(A*CBoS_ap(1,i) + B*CBoS_ml(1,i) + C) / sqrt(A^2 + B^2));
end
   

% H: Distance from CBoS to the straight line formed by the points RFOO2 and RHEE
H =[];
for i = 1:(idxF3-idxF2+1) % should be changed when calculating another phase
    x3 = RFOO2_x(i); 
    y3 = RFOO2_y(i);
    x4 = RHEE_x(i);
    y4 = RHEE_y(i);
    
    D = -((y4-y3)/(x4-x3));             % calculating general or implicit equation of the line
    E = 1;
    F = -y3 + ((y4-y3)/(x4-x3)*x3);

    if D<0
    D = D*(-1);
    E = E*(-1);
    F = F*(-1);
    end
    
    H(1,i) = (abs(D*CBoS_ap(1,i) + E*CBoS_ml(1,i) + F) / sqrt(D^2 + E^2));
end
     
     
%% Postural Stability Measures

AP_dev = [];
ML_dev = [];
Total_dev =[];

for i =1:(idxF3-idxF2+1)      % should be changed when calculating another phase
    ML_dev(1,i) = abs(CBoS_ap(1,i) - CoM_ap(1,i));  % The deviation of the CoM from CBoS in the fore-aft/ anterior-posterior direction at frame i. 

    AP_dev(1,i) = abs(CBoS_ml(1,i) - CoM_ml(1,i));  % The deviation of the CoM from CBoS in the medio-lateral direction at frame i.

    Total_dev(1,i) = sqrt((AP_dev(1,i))^2 + (ML_dev(1,i))^2);   % The total deviation of the CoM from the CBoS at frame i.
end

% "At the end of a trial of an experiment, all the values for AP_dev calculated at each frame should be averaged by dividing it by the number of 
% frames in order to get a value that represents the stability for that trial."

    % Mean postural stability of the trial 
    AP_dev_mean = mean(AP_dev(1,:)); 
    
    ML_dev_mean = mean(ML_dev(1,:));
   
    Total_dev_mean = mean(Total_dev(1,:));
    
%% Postural Stability values

 PS_AP = [];
 PS_ML = [];
 PS_Total = [];
 for i =1:(idxF3-idxF2+1)
 PS_AP(1,i) = (AP_dev(1,i)/V(1,i));
 PosturalStability_AP = mean(PS_AP(1,:))*100;
 end
 for i =1:(idxF3-idxF2+1)
 PS_ML(1,i) = (ML_dev(1,i)/H(1,i));
 PosturalStability_ML = mean(PS_ML(1,:))*100;
 end
 

%% BoS visualization at last frame of the phase

figure(4)
plot(pgon,'FaceColor','green')
hold on
plot(centroid1_x, centroid1_y,'g*')
hold on
plot(pgon2,'FaceColor','yellow')
hold on
plot(centroid2_x, centroid2_y,'y*')
hold on
plot(pgon3,'FaceColor','red')
hold on
plot(centroid3_x, centroid3_y,'r*','linewidth', 1.25)
hold on
plot(pgon4,'FaceColor','cyan')
hold on
plot(centroid4_x, centroid4_y,'c*','linewidth', 1.25)
hold on
plot(CBoS_ap(1,length(CBoS_ap)), CBoS_ml(1,length(CBoS_ap)),'ko','linewidth', 1.25)
hold on
plot(CoM_ap(1,length(CBoS_ap)), CoM_ml(1,length(CBoS_ap)), 'b+', 'linewidth', 1.25)
hold on
xlabel('Millimetres');
ylabel('Millimetres');
title('Example of Base of Support Representation');

PS = table([3; AP_dev_mean;  ML_dev_mean;  Total_dev_mean; PosturalStability_AP; PosturalStability_ML],'VariableNames',{'Postural_Stability_Measures'},'RowNames',{'Phase','Average A/P Deviation (mm)','Average M/L Deviation (mm)','Average Total Deviation (mm)','Postural Stability (A/P)','Postural Stability(M/L)'})
