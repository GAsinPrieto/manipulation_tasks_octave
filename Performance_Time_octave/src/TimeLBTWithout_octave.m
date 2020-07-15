%% Total time to perform the task in Lateral Box Transfer (without exoskeleton segmentation)

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

pkg load signal

clear all % Clear variables
close all % Close figures
clc

load('../tests/data/input/dinamica42_B.mat')
             
Ts = 1/Fs;

%% Total time to perform the task 

t_total = (double(frames)*Ts);

%% Events of interest in feet markers trajectories

% z LHEE RHEE marker trajectories
LHEE_z = LHEE(:,3)';  
RHEE_z = RHEE(:,3)'; 

% First turn

%[pks, locs] = findpeaks(LHEE_z,'minPeakProminence',10);
[pks, locs] = findpeaks(LHEE_z,'DoubleSided');   
%locs=locs(pks>80);                 
%pks=pks(pks>80);
for count=1:1:length(pks)
    if count==1
        signal_aux_L=LHEE_z(1:locs(count));
        signal_aux_R=LHEE_z(locs(count):locs(count+1));
    else if count==length(pks)
            signal_aux_L=LHEE_z(locs(count-1):locs(count));
            signal_aux_R=LHEE_z(locs(count):end);
        else
            signal_aux_L=LHEE_z(locs(count-1):locs(count));
            signal_aux_R=LHEE_z(locs(count):locs(count+1));
        end
    end
    try
      min_L_aux=findpeaks(1.01*max(signal_aux_L)-signal_aux_L,'DoubleSided');
      min_L=min(1.01*max(signal_aux_L)-min_L_aux);
      if isempty(min_L)
        min_L=min(signal_aux_L);
      end
    catch
      min_L=min(signal_aux_L);
    end
    try
      min_R_aux=findpeaks(1.01*max(signal_aux_R)-signal_aux_R,'DoubleSided');
      min_R=min(1.01*max(signal_aux_R)-min_R_aux);
      if isempty(min_R)
        min_R=min(signal_aux_R);
      end
    catch
      min_R=min(signal_aux_R);
    end
    real_min=min(min_L,min_R);
    prominence(count)=pks(count)-real_min;
end
pks=pks(prominence>10);    
locs=locs(prominence>10);


%TF1 = islocalmin(LHEE_z, 'FlatSelection', 'first');
[TF1,locs_TF1] = findpeaks(1.01*max(LHEE_z)-LHEE_z,'DoubleSided');
TF1=1.01*max(LHEE_z)-TF1;
TF1_aux=zeros(1,length(TF1));
for count=1:1:length(TF1)-1
  if TF1(count+1)==TF1(count)
    TF1_aux(count+1)=1;
  end
  for count_aux=1:1:length(pks)
    if TF1(count)==pks(count_aux)
      TF1_aux(count)=1;
    end
  end
end
TF1(find(TF1_aux))=[];
locs_TF1(find(TF1_aux))=[];


%TF1 = islocalmin(LHEE_z, 'FlatSelection','first');
[TF1,locs_TF1] = findpeaks(1.01*max(LHEE_z(1:locs(1)))-LHEE_z(1:locs(1)),'DoubleSided');
TF1=1.01*max(LHEE_z(1:locs(1)))-TF1;

TF1_aux=zeros(1,length(TF1));
for count=1:1:length(TF1)-1
  if TF1(count+1)==TF1(count)
    TF1_aux(count+1)=1;
  end
end
TF1(find(TF1_aux))=[];
locs_TF1(find(TF1_aux))=[]; 


for count_aux=1:1:length(pks)-1
  %TF1 = islocalmin(LHEE_z, 'FlatSelection','first');
  [TF1_2,locs_TF1_2] = findpeaks(1.01*max(LHEE_z(locs(count_aux):locs(count_aux+1)))-LHEE_z(locs(count_aux):locs(count_aux+1)),'DoubleSided');
  TF1_2=1.01*max(LHEE_z(locs(count_aux):locs(count_aux+1)))-TF1_2;

  TF1_aux_2=zeros(1,length(TF1_2));
  for count=1:1:length(TF1_2)-1
    if TF1_2(count+1)==TF1_2(count)
      TF1_aux_2(count+1)=1;
    end
  end
  TF1_2(find(TF1_aux_2))=[];
  locs_TF1_2(find(TF1_aux_2))=[];
  TF1 = [TF1,TF1_2];
  locs_TF1 = [locs_TF1,locs_TF1_2+locs(count_aux)];
endfor


%TF1 = islocalmin(LHEE_z, 'FlatSelection','first');
[TF1_2,locs_TF1_2] = findpeaks(1.01*max(LHEE_z(locs(end):end))-LHEE_z(locs(end):end),'DoubleSided');
TF1_2=1.01*max(LHEE_z(locs(end):end))-TF1_2;

TF1_aux_2=zeros(1,length(TF1_2));
for count=1:1:length(TF1_2)-1
  if TF1_2(count+1)==TF1_2(count)
    TF1_aux_2(count+1)=1;
  end
end
TF1_2(find(TF1_aux_2))=[];
locs_TF1_2(find(TF1_aux_2))=[]; 

TF1 = [TF1,TF1_2];
locs_TF1 = [locs_TF1,locs_TF1_2+locs(end)];

for count=1:1:length(TF1)-1
  for count_aux=1:1:length(pks)
      if TF1(count)==pks(count_aux)
        TF1_aux(count)=1;
      end
  end
end

TF1(find(TF1_aux))=[];
locs_TF1(find(TF1_aux))=[];









##idx = find(TF1);
idx = locs_TF1;



flat2 = idx > locs(1);
idx_flat2 = find(flat2);

%[pks2, locs2] = findpeaks(RHEE_z,'minPeakProminence',10);
clear prominence; 
[pks2, locs2] = findpeaks(RHEE_z,'DoubleSided');   
%locs2=locs2(pks2>80);                 
%pks2=pks2(pks2>80);
for count=1:1:length(pks2)
    if count==1
        signal_aux_L=RHEE_z(1:locs2(count));
        signal_aux_R=RHEE_z(locs2(count):locs2(count+1));
    else if count==length(pks2)
            signal_aux_L=RHEE_z(locs2(count-1):locs2(count));
            signal_aux_R=RHEE_z(locs2(count):end);
        else
            signal_aux_L=RHEE_z(locs2(count-1):locs2(count));
            signal_aux_R=RHEE_z(locs2(count):locs2(count+1));
        end
    end
    try
      min_L_aux=findpeaks(1.01*max(signal_aux_L)-signal_aux_L,'DoubleSided');
      min_L=min(1.01*max(signal_aux_L)-min_L_aux);
      if isempty(min_L)
        min_L=min(signal_aux_L);
      end
    catch
      min_L=min(signal_aux_L);
    end
    try
      min_R_aux=findpeaks(1.01*max(signal_aux_R)-signal_aux_R,'DoubleSided');
      min_R=min(1.01*max(signal_aux_R)-min_R_aux);
      if isempty(min_R)
        min_R=min(signal_aux_R);
      end
    catch
      min_R=min(signal_aux_R);
    end
    real_min=min(min_L,min_R);
    prominence(count)=pks2(count)-real_min;
end
pks2=pks2(prominence>10);    
locs2=locs2(prominence>10);



##%TF2 = islocalmin(RHEE_z, 'FlatSelection', 'first');
##[TF2,locs_TF2] = findpeaks(1.01*max(RHEE_z)-RHEE_z,'DoubleSided');
##TF2=1.01*max(RHEE_z)-TF2;
##TF2_aux=zeros(1,length(TF2));
##for count=1:1:length(TF2)-1
##  if TF2(count+1)==TF2(count)
##    TF2_aux(count+1)=1;
##  end
##  for count_aux=1:1:length(pks2)
##    if TF2(count)==pks2(count_aux)
##      TF2_aux(count)=1;
##    end
##  end
##end
##TF2(find(TF2_aux))=[];
##locs_TF2(find(TF2_aux))=[];


[TF2,locs_TF2] = findpeaks(1.01*max(RHEE_z(1:locs(1)))-RHEE_z(1:locs(1)),'DoubleSided');
TF2=1.01*max(RHEE_z(1:locs(1)))-TF2;

TF2_aux=zeros(1,length(TF2));
for count=1:1:length(TF2)-1
  if TF2(count+1)==TF2(count)
    TF2_aux(count+1)=1;
  end
end
TF2(find(TF2_aux))=[];
locs_TF2(find(TF2_aux))=[]; 


for count_aux=1:1:length(pks)-1
  %TF2 = islocalmin(RHEE_z, 'FlatSelection','first');
  [TF2_2,locs_TF2_2] = findpeaks(1.01*max(RHEE_z(locs(count_aux):locs(count_aux+1)))-RHEE_z(locs(count_aux):locs(count_aux+1)),'DoubleSided');
  TF2_2=1.01*max(RHEE_z(locs(count_aux):locs(count_aux+1)))-TF2_2;

  TF2_aux_2=zeros(1,length(TF2_2));
  for count=1:1:length(TF2_2)-1
    if TF2_2(count+1)==TF2_2(count)
      TF2_aux_2(count+1)=1;
    end
  end
  TF2_2(find(TF2_aux_2))=[];
  locs_TF2_2(find(TF2_aux_2))=[];
  TF2 = [TF2,TF2_2];
  locs_TF2 = [locs_TF2,locs_TF2_2+locs(count_aux)];
endfor


%TF2 = islocalmin(RHEE_z, 'FlatSelection','first');
[TF2_2,locs_TF2_2] = findpeaks(1.01*max(RHEE_z(locs(end):end))-RHEE_z(locs(end):end),'DoubleSided');
TF2_2=1.01*max(RHEE_z(locs(end):end))-TF2_2;

TF2_aux_2=zeros(1,length(TF2_2));
for count=1:1:length(TF2_2)-1
  if TF2_2(count+1)==TF2_2(count)
    TF2_aux_2(count+1)=1;
  end
end
TF2_2(find(TF2_aux_2))=[];
locs_TF2_2(find(TF2_aux_2))=[]; 

TF2 = [TF2,TF2_2];
locs_TF2 = [locs_TF2,locs_TF2_2+locs(end)];

for count=1:1:length(TF2)-1
  for count_aux=1:1:length(pks2)
      if TF2(count)==pks2(count_aux)
        TF2_aux(count)=1;
      end
  end
end

TF2(find(TF2_aux))=[];
locs_TF2(find(TF2_aux))=[];





##idx2 = find(TF2);
idx2 = locs_TF2;


flat3 = idx2 < locs2(1);  % 42 es <
idx_flat3 = find(flat3);

% Second turn

flat5 = idx < locs(2);      
idx_flat5 = find(flat5);

flat5 = idx > locs(2);
idx_flat5 = find(flat5);

flat8 = idx2 > locs2(2);
idx_flat8 = find(flat8);

% Indexes definition for signal segmentation
idxO1 = 1;
idxF1 = idx2(idx_flat3(length(idx_flat3))); 
idxF2 = idx(idx_flat2(1)); 
idxF3 = idx(idx_flat5(length(idx_flat5))); 
idxF4 = idx2(idx_flat8(1)); 
idxF5 = length(LHEE_z);
   
% Visualization of the segmentation

figure(1)
t= 0:(length(LHEE)-1);

RHEE_z_z = zscore(RHEE_z);
LHEE_z_z = zscore(LHEE_z);

plot(t,LHEE_z_z);
hold on
plot(t,RHEE_z_z);

plot(t(idxO1),LHEE_z_z(idxO1),'r*');
plot(t(idxF2),LHEE_z_z(idxF2),'r*');
plot(t(idxF5),LHEE_z_z(idxF5),'r*');

plot(t(idxO1),RHEE_z_z(idxO1),'r*');
plot(t(idxF4),RHEE_z_z(idxF4),'r*');
plot(t(idxF5),RHEE_z_z(idxF5),'r*');

hold off
legend('Left Heel Z', 'Right Heel Z');
title('Segmentation Without Exoskeleton Trials')

% Three phases:
    % 1: Subject picks up the box in the sagittal plane and takes a step to rotate to the frontal plane (idxO1:idxF2)
    % 2: Subject deposits the box in the frontal plane and takes a step to rotate back to the sagittal plane(idxF2:idxF4)
    % 3: Subject deposits the box in the sagittal plane (idxF4:idxF5)
    

%% Duration of each phase

% phase 1 
phase1 = LHEE_z(idxO1:idxF2);
t_phase1= double(length(phase1))*Ts;

% phase 2
phase2 = LHEE_z(idxF2:idxF4);
t_phase2= double(length(phase2)-1)*Ts;

% phase 3
phase3 = LHEE_z(idxF4:idxF5);
t_phase3= double(length(phase3)-1)*Ts;


%% Verification
Sum = t_phase1 + t_phase2 + t_phase3;

if round(Sum - t_total) == 0
    disp('Está bien hecho!');
else 
    disp('Algo ha fallado');
end

##Time = table([t_phase1; t_phase2; t_phase3; Sum; t_total],'VariableNames',{'Seconds'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Sum', 'Complete Trial'})
Time_seconds = struct;
Time_seconds.Phase_1=t_phase1;
Time_seconds.Phase_2=t_phase2; 
Time_seconds.Phase_3=t_phase3; 
Time_seconds.Sum=Sum;
Time_seconds.Complete_Trial=t_total;
  
Time_seconds



