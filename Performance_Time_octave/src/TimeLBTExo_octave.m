%% Total time to perform the task in Lateral Box Transfer (exoskeleton segmentation)

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

pkg load signal

clear all % Clear variables
close all % Close figures
clc

load('../tests/data/input/dinamica56_B.mat')
             
Ts = 1/Fs;

%% Total time to perform the task 

t_total = (double(frames)*Ts);

%% Events of interest in feet markers trajectories

% Feet markers for signal segmentation
LTOE_z = LTOE(:,3)'; 
RTOE_z = RTOE(:,3)';   

% Left toe:
 
% [pks, locs] = findpeaks(LTOE_z,'minPeakProminence',10,'MinPeakHeight',80);                    
[pks, locs] = findpeaks(LTOE_z,'DoubleSided');   
locs=locs(pks>80);                 
pks=pks(pks>80);
for count=1:1:length(pks)
    if count==1
        signal_aux_L=LTOE_z(1:locs(count));
        signal_aux_R=LTOE_z(locs(count):locs(count+1));
    else if count==length(pks)
            signal_aux_L=LTOE_z(locs(count-1):locs(count));
            signal_aux_R=LTOE_z(locs(count):end);
        else
            signal_aux_L=LTOE_z(locs(count-1):locs(count));
            signal_aux_R=LTOE_z(locs(count):locs(count+1));
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
    
% TF1 = islocalmin(LTOE_z, 'FlatSelection', 'first');
%[TF1,locs_TF1] = findpeaks(1.01*max(LTOE_z)-LTOE_z,'DoubleSided');
%TF1=1.01*max(LTOE_z)-TF1;
%TF1_aux=zeros(1,length(TF1));
%for count=1:1:length(TF1)-1
%  if TF1(count+1)==TF1(count)
%    TF1_aux(count+1)=1;
%  end
%  for count_aux=1:1:length(pks)
%    if TF1(count)==pks(count_aux)
%      TF1_aux(count)=1;
%    end
%  end
%end
%TF1(find(TF1_aux))=[];
%locs_TF1(find(TF1_aux))=[];





%TF1 = islocalmin(LTOE_z, 'FlatSelection','first');
[TF1,locs_TF1] = findpeaks(1.01*max(LTOE_z(1:locs(1)))-LTOE_z(1:locs(1)),'DoubleSided');
TF1=1.01*max(LTOE_z(1:locs(1)))-TF1;

TF1_aux=zeros(1,length(TF1));
for count=1:1:length(TF1)-1
  if TF1(count+1)==TF1(count)
    TF1_aux(count+1)=1;
  end
end
TF1(find(TF1_aux))=[];
locs_TF1(find(TF1_aux))=[]; 


for count_aux=1:1:length(pks)-1
  %TF1 = islocalmin(LTOE_z, 'FlatSelection','first');
  [TF1_2,locs_TF1_2] = findpeaks(1.01*max(LTOE_z(locs(count_aux):locs(count_aux+1)))-LTOE_z(locs(count_aux):locs(count_aux+1)),'DoubleSided');
  TF1_2=1.01*max(LTOE_z(locs(count_aux):locs(count_aux+1)))-TF1_2;

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


%TF1 = islocalmin(LTOE_z, 'FlatSelection','first');
[TF1_2,locs_TF1_2] = findpeaks(1.01*max(LTOE_z(locs(end):end))-LTOE_z(locs(end):end),'DoubleSided');
TF1_2=1.01*max(LTOE_z(locs(end):end))-TF1_2;

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

flat = idx < locs(1);
idx_flat = find(flat);
 
flat2 = idx < locs(3)& idx > locs(2);
idx_flat2 = find(flat2);
 
flat3 = idx > locs(4);
idx_flat3 = find(flat3);
 
 
% Right Toe: 
 
% [pks2, locs2] = findpeaks(RTOE_z,'minPeakProminence',10); 
clear prominence; 
[pks2, locs2] = findpeaks(RTOE_z,'DoubleSided');   
locs2=locs2(pks2>80);                 
pks2=pks2(pks2>80);
for count=1:1:length(pks2)
    if count==1
        signal_aux_L=RTOE_z(1:locs2(count));
        signal_aux_R=RTOE_z(locs2(count):locs2(count+1));
    else if count==length(pks2)
            signal_aux_L=RTOE_z(locs2(count-1):locs2(count));
            signal_aux_R=RTOE_z(locs2(count):end);
        else
            signal_aux_L=RTOE_z(locs2(count-1):locs2(count));
            signal_aux_R=RTOE_z(locs2(count):locs2(count+1));
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



##% TF2 = islocalmin(RTOE_z, 'FlatSelection', 'first');
##[TF2,locs_TF2] = findpeaks(1.01*max(RTOE_z)-RTOE_z,'DoubleSided');
##TF2=1.01*max(RTOE_z)-TF2;
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





[TF2,locs_TF2] = findpeaks(1.01*max(RTOE_z(1:locs(1)))-RTOE_z(1:locs(1)),'DoubleSided');
TF2=1.01*max(RTOE_z(1:locs(1)))-TF2;

TF2_aux=zeros(1,length(TF2));
for count=1:1:length(TF2)-1
  if TF2(count+1)==TF2(count)
    TF2_aux(count+1)=1;
  end
end
TF2(find(TF2_aux))=[];
locs_TF2(find(TF2_aux))=[]; 


for count_aux=1:1:length(pks)-1
  %TF2 = islocalmin(RTOE_z, 'FlatSelection','first');
  [TF2_2,locs_TF2_2] = findpeaks(1.01*max(RTOE_z(locs(count_aux):locs(count_aux+1)))-RTOE_z(locs(count_aux):locs(count_aux+1)),'DoubleSided');
  TF2_2=1.01*max(RTOE_z(locs(count_aux):locs(count_aux+1)))-TF2_2;

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


%TF2 = islocalmin(RTOE_z, 'FlatSelection','first');
[TF2_2,locs_TF2_2] = findpeaks(1.01*max(RTOE_z(locs(end):end))-RTOE_z(locs(end):end),'DoubleSided');
TF2_2=1.01*max(RTOE_z(locs(end):end))-TF2_2;

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
idxF4 = idx(idx_flat3(4)); %idxO5
idxF5 = length(LTOE_z);
 
figure(1)
t= 0:(length(LTOE)-1);

RTOE_z_z = zscore(RTOE_z);
LTOE_z_z = zscore(LTOE_z);

plot(t,LTOE_z_z);
hold on
plot(t,RTOE_z_z);


##plot(t,zscore(LTOE_z),'r-');
plot(t(idxO1),LTOE_z_z(idxO1),'r*');

##plot(t,zscore(RTOE_z),'r-');
plot(t(idxF2),RTOE_z_z(idxF2),'r*');

##plot(t,zscore(LTOE_z),'r-');
plot(t(idxF3),LTOE_z_z(idxF3),'r*');

##plot(t,zscore(LTOE_z),'r-');
plot(t(idxF5),LTOE_z_z(idxF5),'r*');

##plot(t,zscore(RTOE_z),'r-');
plot(t(idxO1),RTOE_z_z(idxO1),'r*');

##plot(t,zscore(RTOE_z),'r-');
plot(t(idxF1),RTOE_z_z(idxF1),'r*');

##plot(t,zscore(LTOE_z),'r-');
plot(t(idxF4),LTOE_z_z(idxF4),'r*');

##plot(t,zscore(RTOE_z),'r-');
plot(t(idxF5),RTOE_z_z(idxF5),'r*');


hold off
legend('Left Toe Z', 'Right Toe Z');
title('Segmentation Exoskeleton Trials in 5 phases')


% Five phases shown in figure 1:
    % 1: Subject picks up the box in the sagittal plane (idxO1:idxF1)
    % 2: Subject takes some steps to rotate to the frontal plane (idxF1:idxF2)
    % 3: Subject deposits the box in the frontal plane (idxF2:idxF3)
    % 4: Subject takes some steps to rotate back to the sagittal plane (idxF3:idxF4)
    % 5: Subject deposits the box in the sagittal plane (idxF4:idxF5)
    

%% Duration of each phase

% phase 1 
phase1 = LTOE_z(idxO1:idxF1);
t_phase1= double(length(phase1))*Ts;

% phase 2
phase2 = LTOE_z(idxF1:idxF2);
t_phase2= double(length(phase2)-1)*Ts;

% phase 3
phase3 = LTOE_z(idxF2:idxF3);
t_phase3= double(length(phase3)-1)*Ts;

% phase 4
phase4 = LTOE_z(idxF3:idxF4);
t_phase4= double(length(phase4)-1)*Ts;

% phase 5
phase5 = LTOE_z(idxF4:idxF5);
t_phase5= double(length(phase5)-1)*Ts;


%% Verification
Sum = t_phase1 + t_phase2 + t_phase3 + t_phase4 + t_phase5;

if round(Sum - t_total) == 0
    disp('Está bien hecho!');
else 
    disp('Algo ha fallado');
end

##Time = table([t_phase1; t_phase2; t_phase3; t_phase4; t_phase5; Sum; t_total],'VariableNames',{'Seconds'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4', 'Phase 5', 'Sum', 'Complete Trial'})
Time_seconds = struct;
Time_seconds.Phase_1=t_phase1;
Time_seconds.Phase_2=t_phase2; 
Time_seconds.Phase_3=t_phase3; 
Time_seconds.Phase_4=t_phase4; 
Time_seconds.Phase_5=t_phase5;
Time_seconds.Sum=Sum;
Time_seconds.Complete_Trial=t_total;
  
Time_seconds