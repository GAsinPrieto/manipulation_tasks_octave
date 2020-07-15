%% Total time to perform the task in Sagittal Lifting

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

pkg load signal

clear all % Clear variables
close all % Close figures
clc

load('../tests/data/input/dinamica44_B.mat')
             
Ts = 1/Fs;

%% Total time to perform the task 

t_total = (double(frames)*Ts);

%% Events of interest in vertical box trajectory identification
 
position_Z = SIDEBOX3(:,3)';
 
%[pks, locs] = findpeaks(position_Z,'minPeakProminence',20);
[pks, locs] = findpeaks(position_Z,'DoubleSided');   
for count=1:1:length(pks)
    if count==1
        signal_aux_L=position_Z(1:locs(count));
        signal_aux_R=position_Z(locs(count):locs(count+1));
    else if count==length(pks)
            signal_aux_L=position_Z(locs(count-1):locs(count));
            signal_aux_R=position_Z(locs(count):end);
        else
            signal_aux_L=position_Z(locs(count-1):locs(count));
            signal_aux_R=position_Z(locs(count):locs(count+1));
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
pks=pks(prominence>20);    
locs=locs(prominence>20);


t = (0:1:(length(position_Z)-1))*Ts;
 
 
%TF1 = islocalmin(position_Z, 'FlatSelection','first');
[TF1,locs_TF1] = findpeaks(1.01*max(position_Z(1:locs(1)))-position_Z(1:locs(1)),'DoubleSided');
TF1=1.01*max(position_Z(1:locs(1)))-TF1;

TF1_aux=zeros(1,length(TF1));
for count=1:1:length(TF1)-1
  if TF1(count+1)==TF1(count)
    TF1_aux(count+1)=1;
  end
end
TF1(find(TF1_aux))=[];
locs_TF1(find(TF1_aux))=[]; 


for count_aux=1:1:length(pks)-1
  %TF1 = islocalmin(position_Z, 'FlatSelection','first');
  [TF1_2,locs_TF1_2] = findpeaks(1.01*max(position_Z(locs(count_aux):locs(count_aux+1)))-position_Z(locs(count_aux):locs(count_aux+1)),'DoubleSided');
  TF1_2=1.01*max(position_Z(locs(count_aux):locs(count_aux+1)))-TF1_2;

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


%TF1 = islocalmin(position_Z, 'FlatSelection','first');
[TF1_2,locs_TF1_2] = findpeaks(1.01*max(position_Z(locs(end):end))-position_Z(locs(end):end),'DoubleSided');
TF1_2=1.01*max(position_Z(locs(end):end))-TF1_2;

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



%for i = 1:frames                               % solo para dinamica 03
 %   TF2(1,i) = position_Z(1,i) > 0;
%end
%idx = find(TF1);
idx = locs_TF1;

%idx_03 = find(TF2);
flat = idx < locs(1);
idx_flat = find(flat);
flat2 = idx > locs(1);
idx_flat2 = find(flat2);
 
flat3 = idx < locs(2);
idx_flat3 = find(flat3);
 
flat4 = idx > locs(2);
idx_flat4 = find(flat4);
 
% Indexes definiton for signal segmentation
idxO1 = 1;
idxF1 = idx(idx_flat(length(idx_flat)));
%idxF1 = idx_03(1)-1;                               % solo para dinamica 03
idxF2 = idx(idx_flat2(1)); 
idxO3 = idx(idx_flat3(length(idx_flat3))); 
idxF3 = idx(idx_flat4(1));                          % para dinamica 59 es un 4, para din 03 din 61 es un 2 y para el resto es un 1          
idxF4 = length(position_Z)


figure(1)
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
plot(t,position_Z,t(locs),position_Z(locs),'o');
hold on

%plot(t,position_Z,'r*','MarkerIndices',idxO1);
plot(t(idxO1),position_Z(idxO1),'r*');

%plot(t,position_Z,'r*','MarkerIndices',idxF1);
plot(t(idxF1),position_Z(idxF1),'r*');

%plot(t,position_Z,'r*','MarkerIndices',idxF2);
plot(t(idxF2),position_Z(idxF2),'r*');

%plot(t,position_Z,'r*','MarkerIndices',idxO3);
plot(t(idxO3),position_Z(idxO3),'r*');

%plot(t,position_Z,'r*','MarkerIndices',idxF3);
plot(t(idxF3),position_Z(idxF3),'r*');

%plot(t,position_Z,'r*','MarkerIndices',idxF4);
plot(t(idxF4),position_Z(idxF4),'r*');

hold off
xlabel('Time in seconds'); 
ylabel('Millimetres');
title('Box signal segmentation in 5 phases');

% phase 1: Lowering without box (idxO1:idxF1)
% phase 2: Holding and lifting the box (idxF1:idxF2)
% phase 3: Flat section where subject goes back to anatomic position and afterwards, subject prepares to hold the box again (idxF2:idxO3)
% phase 4: Holding the box again and lowering it (idxO3:idxF3)
% phase 5: Raising up without box (idxF3:idxF4)        

%% Duration of each phase

% phase 1 
phase1 = position_Z(idxO1:idxF1);
t_phase1= double(length(phase1))*Ts;

% phase 2
phase2 = position_Z(idxF1:idxF2);
t_phase2= double(length(phase2)-1)*Ts;

% phase 3
phase3 = position_Z(idxF2:idxO3);
t_phase3= double(length(phase3)-1)*Ts;

% phase 4
phase4 = position_Z(idxO3:idxF3);
t_phase4= double(length(phase4)-1)*Ts;

% phase 5
phase5 = position_Z(idxF3:idxF4);
t_phase5= double(length(phase5)-1)*Ts;

%% Verification
Sum = t_phase1 + t_phase2 + t_phase3 + t_phase4 + t_phase5;

if round(Sum - t_total) == 0
    disp('Está bien hecho!');
else 
    disp('Algo ha fallado');
end

%Time = table([t_phase1; t_phase2; t_phase3; t_phase4; t_phase5; Sum; t_total],'VariableNames',{'Seconds'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4', 'Phase 5', 'Sum', 'Complete Trial'})
Time_seconds = struct;
Time_seconds.Phase_1=t_phase1;
Time_seconds.Phase_2=t_phase2; 
Time_seconds.Phase_3=t_phase3; 
Time_seconds.Phase_4=t_phase4; 
Time_seconds.Phase_5=t_phase5;
Time_seconds.Sum=Sum;
Time_seconds.Complete_Trial=t_total;
  
Time_seconds
