%% Range of movement of Ankle in Lateral Box Transfer

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

% Exoskeleton trials segmentation: 54, 55, 56, 62.
% Five phases shown in figure 1:
    % 1: Subject picks up the box in the sagittal plane (idxO1:idxF1)
    % 2: Subject takes some steps to rotate to the frontal plane (idxF1:idxF2)
    % 3: Subject deposits the box in the frontal plane (idxF2:idxF3)
    % 4: Subject takes some steps to rotate back to the sagittal plane (idxF3:idxF4)
    % 5: Subject deposits the box in the sagittal plane (idxF4:idxF5)


pkg load signal

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica56_B.mat')
           
Ts = 1/Fs;
t_total = (double(frames)*Ts);
t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));

% Normal Ranges of Ankle Motion
for i =1:(length(t_1000))
    flex_ankle(1,i) = 45; % flexion
    ext_ankle(1,i) = -20;   % extension
    supination_ankle(1,i) = 20;  % supination
    pronation_ankle(1,i) = -20;  % pronation
end

% Feet markers for signal segmentation
LTOE_z = LTOE(:,3)'; 
RTOE_z = RTOE(:,3)';   
%LHEE_z = LHEE(:,3)';  
%RHEE_z = RHEE(:,3)'; 


t= 0:(length(LTOE)-1);

%% SEGMENTATION LEFT TOE: 

%[pks, locs] = findpeaks(LTOE_z,'minPeakProminence',10,'MinPeakHeight',80);                     
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



%TF1 = islocalmin(LTOE_z, 'FlatSelection', 'first');
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




%idx = find(TF1);
idx = locs_TF1;


flat = idx < locs(1);
idx_flat = find(flat);

flat2 = idx < locs(3)& idx > locs(2);
idx_flat2 = find(flat2);

flat3 = idx > locs(4);
idx_flat3 = find(flat3);


figure(1)
plot(t, LTOE_z, t(locs), LTOE_z(locs),'o')
hold on
%plot(t, LTOE_z,'r*','MarkerIndices',idx(idx_flat(length(idx_flat))));
plot(t(idx(idx_flat(length(idx_flat)))),LTOE_z(idx(idx_flat(length(idx_flat)))),'r*');

%plot(t, LTOE_z,'r*','MarkerIndices',idx(idx_flat2(3)));                    
plot(t(idx(idx_flat2(3))),LTOE_z(idx(idx_flat2(3))),'r*');

%plot(t, LTOE_z,'r*','MarkerIndices',idx(idx_flat2(length(idx_flat2))));
plot(t(idx(idx_flat2(length(idx_flat2)))),LTOE_z(idx(idx_flat2(length(idx_flat2)))),'r*');

%plot(t, LTOE_z,'r*','MarkerIndices',idx(idx_flat3(4)));                     
plot(t(idx(idx_flat3(4))),LTOE_z(idx(idx_flat3(4))),'r*');

hold off
title('Left Toe');

%% SEGMENTATION RIGHT TOE: 

%[pks2, locs2] = findpeaks(RTOE_z,'minPeakProminence',10);
[pks2, locs2] = findpeaks(RTOE_z,'DoubleSided');   
%locs2=locs2(pks2>80);                 
%pks2=pks2(pks2>80);
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





%TF2 = islocalmin(RTOE_z, 'FlatSelection', 'first');
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
    for count_aux=1:1:length(pks)
        if TF2(count)==pks(count_aux)
          TF2_aux(count)=1;
        end
    end
end

TF2(find(TF2_aux))=[];
locs_TF2(find(TF2_aux))=[];









%idx2 = find(TF2);
idx2 = locs_TF2;


flat4 = idx2 < locs2(1);
idx_flat4 = find(flat4);

flat5 = idx2 > locs2(length(locs2));
idx_flat5 = find(flat5);
flat6 = idx2 > locs2(3);
idx_flat6 = find(flat6);

figure(2)
plot(t, RTOE_z, t(locs2), RTOE_z(locs2),'o')
hold on

%plot(t, RTOE_z,'r*','MarkerIndices',idx2(idx_flat4(length(idx_flat4))));
plot(t(idx2(idx_flat4(length(idx_flat4)))),RTOE_z(idx2(idx_flat4(length(idx_flat4)))),'r*');

%plot(t, RTOE_z,'r*','MarkerIndices',idx2(idx_flat6(1)));
plot(t(idx2(idx_flat6(1))),RTOE_z(idx2(idx_flat6(1))),'r*');

%plot(t, RTOE_z,'r*','MarkerIndices',idx2(idx_flat5(2))); 
plot(t(idx2(idx_flat5(2))),RTOE_z(idx2(idx_flat5(2))),'r*');

title('Right Toe');

% Indexes definiton for signal segmentation
idxO1 = 1;
idxF1 = idx2(idx_flat4(length(idx_flat4))); %idxO2
idxF2 = idx2(idx_flat6(1)); %idxO3
idxF3 = idx(idx_flat2(length(idx_flat2))); %idxO4
idxF4 = idx(idx_flat3(3)); %idxO5
idxF5 = length(LTOE_z);

figure(3)

RTOE_z_z = zscore(RTOE_z);
LTOE_z_z = zscore(LTOE_z);

plot(t,LTOE_z_z);
hold on
plot(t,RTOE_z_z);

%plot(t, zscore(LTOE_z));
%hold on
%plot(t,zscore(RTOE_z));

%plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxO1);
plot(t(idxO1),LTOE_z_z(idxO1),'r*');

%plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF2);
plot(t(idxF2),RTOE_z_z(idxF2),'r*');

%plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF3);
plot(t(idxF3),LTOE_z_z(idxF3),'r*');

%plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF5);
plot(t(idxF5),LTOE_z_z(idxF5),'r*');

%plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxO1);
plot(t(idxO1),RTOE_z_z(idxO1),'r*');

%plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF1);
plot(t(idxF1),RTOE_z_z(idxF1),'r*');

%plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF4);
plot(t(idxF4),LTOE_z_z(idxF4),'r*');

%plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF5);
plot(t(idxF5),RTOE_z_z(idxF5),'r*');

hold off
xlabel('Time in frames');
ylabel('Milimeters');
legend('Left Toe Z', 'Right Toe Z');
legend('Left Toe Z', 'Right Toe Z');
title('Lateral Box Transfer Segmentation for Exoskeleton Trials')


%% Import model data for Ankle Relative Angles

% Ankle  
        RAnkleAngles_x = ModelData.Raw.(ModelOutput{5})(1,:); 
        RAnkleAngles_y = ModelData.Raw.(ModelOutput{5})(2,:);
        RAnkleAngles_z = ModelData.Raw.(ModelOutput{5})(3,:);
        

        LAnkleAngles_x = ModelData.Raw.(ModelOutput{6})(1,:);
        LAnkleAngles_y = ModelData.Raw.(ModelOutput{6})(2,:);
        LAnkleAngles_z = ModelData.Raw.(ModelOutput{6})(3,:);
        

%% Ankle signals visualization before going through segmentation

figure(4)
frames = length(LTOE_z);
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,RAnkleAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_ankle,'--',t_1000, ext_ankle,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-30 50]);
title('Right Ankle angles: Flexo-Extension');

figure(5)
plot(t_1000, RAnkleAngles_y(1:(frames/1000):frames));
hold on 
plot(t_1000, supination_ankle,'--',t_1000, pronation_ankle,'--'); % ver cual es la positiva
legend({'y', 'Normal range of supination-pronation'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-30 30]);
title('Right Ankle angles: Supination-Pronation');

figure(6)
plot(t_1000, RAnkleAngles_z(1:(frames/1000):frames));
legend({'z'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 50]);
title('Right Ankle angles: Internal-External Rotation');

figure(7)
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,LAnkleAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_ankle,'--',t_1000, ext_ankle,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-30 50]);
title('Left Ankle angles: Flexo-Extension');

figure(8)
plot(t_1000, LAnkleAngles_y(1:(frames/1000):frames));
hold on 
plot(t_1000, supination_ankle,'--',t_1000, pronation_ankle,'--'); % ver cual es la positiva
legend({'y', 'Normal range of supination-pronation'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-30 30]);
title('Left Ankle angles: Supination-Pronation');

figure(9)
plot(t_1000, LAnkleAngles_z(1:(frames/1000):frames));
legend({'z'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 80]);
title('Left Ankle angles: Internal-External Rotation');

%figure(3)
%plot(t ,RAnkleAngles_x, t, RAnkleAngles_y, t, RAnkleAngles_z);
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Right Ankle angles');

%figure(4)
%plot(t,LAnkleAngles_x, t, LAnkleAngles_y, t, LAnkleAngles_z);   
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Left Ankle angles');

%% Signal Segmentation

% Right Ankle

phase1_RAnkle_x = RAnkleAngles_x(idxO1:idxF1);  
phase1_RAnkle_y = RAnkleAngles_y(idxO1:idxF1);
phase1_RAnkle_z = RAnkleAngles_z(idxO1:idxF1);

phase2_RAnkle_x = RAnkleAngles_x(idxF1:idxF2); 
phase2_RAnkle_y = RAnkleAngles_y(idxF1:idxF2);
phase2_RAnkle_z = RAnkleAngles_z(idxF1:idxF2);

phase3_RAnkle_x = RAnkleAngles_x(idxF2:idxF3);
phase3_RAnkle_y = RAnkleAngles_y(idxF2:idxF3);
phase3_RAnkle_z = RAnkleAngles_z(idxF2:idxF3);

phase4_RAnkle_x = RAnkleAngles_x(idxF3:idxF4); 
phase4_RAnkle_y = RAnkleAngles_y(idxF3:idxF4);
phase4_RAnkle_z = RAnkleAngles_z(idxF3:idxF4);

phase5_RAnkle_x = RAnkleAngles_x(idxF4:idxF5); 
phase5_RAnkle_y = RAnkleAngles_y(idxF4:idxF5);
phase5_RAnkle_z = RAnkleAngles_z(idxF4:idxF5);

% Left Ankle

phase1_LAnkle_x = LAnkleAngles_x(idxO1:idxF1);
phase1_LAnkle_y = LAnkleAngles_y(idxO1:idxF1);
phase1_LAnkle_z = LAnkleAngles_z(idxO1:idxF1);

phase2_LAnkle_x = LAnkleAngles_x(idxF1:idxF2);
phase2_LAnkle_y = LAnkleAngles_y(idxF1:idxF2);
phase2_LAnkle_z = LAnkleAngles_z(idxF1:idxF2);

phase3_LAnkle_x = LAnkleAngles_x(idxF2:idxF3);
phase3_LAnkle_y = LAnkleAngles_y(idxF2:idxF3);
phase3_LAnkle_z = LAnkleAngles_z(idxF2:idxF3);

phase4_LAnkle_x = LAnkleAngles_x(idxF3:idxF4);
phase4_LAnkle_y = LAnkleAngles_y(idxF3:idxF4);
phase4_LAnkle_z = LAnkleAngles_z(idxF3:idxF4);

phase5_LAnkle_x = LAnkleAngles_x(idxF4:idxF5);
phase5_LAnkle_y = LAnkleAngles_y(idxF4:idxF5);
phase5_LAnkle_z = LAnkleAngles_z(idxF4:idxF5);

%% time vectors definition

t2 = 0:1:(length(phase3_RAnkle_x)-1);
t3 = 0:1:(length(phase1_RAnkle_x)-1);
t4 = 0:1:(length(phase2_RAnkle_x)-1);
t5 = 0:1:(length(phase4_RAnkle_x)-1);
t6 = 0:1:(length(phase5_RAnkle_x)-1);

% Right Ankle figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(10)
subplot(5,1,1)
plot (t3, phase1_RAnkle_x, t3, phase1_RAnkle_y, '--',t3,phase1_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('RIGHT ANGLES: Phase 1'); % Picking up the box sagittal plane
subplot(5,1,2)
plot(t4, phase2_RAnkle_x,t4,phase2_RAnkle_y, '--',t4,phase2_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('Phase 2'); % Taking steps from sagittal to frontal
subplot(5,1,3)
plot(t2, phase3_RAnkle_x, t2, phase3_RAnkle_y,'--', t2, phase3_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3'); % Placing the box frontal plane
subplot(5,1,4)
plot(t5, phase4_RAnkle_x, t5, phase4_RAnkle_y,'--', t5, phase4_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('Phase 4'); % Taking steps from frontal to sagittal
subplot(5,1,5)
plot(t6, phase5_RAnkle_x, t6, phase5_RAnkle_y,'--', t6, phase5_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('Phase 5'); % Placing the box sagittal plane again

% Left Ankle figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(11)
subplot(5,1,1)
plot (t3, phase1_LAnkle_x, t3, phase1_LAnkle_y, '--',t3,phase1_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('LEFT ANGLES: Phase 1'); % Picking up the box sagittal plane
subplot(5,1,2)
plot(t4, phase2_LAnkle_x,t4, phase2_LAnkle_y, '--',t4, phase2_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 2'); % Taking steps from sagittal to frontal
subplot(5,1,3)
plot(t2, phase3_LAnkle_x, t2, phase3_LAnkle_y,'--', t2, phase3_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3'); % Placing the box frontal plane
subplot(5,1,4)
plot(t5, phase4_LAnkle_x, t5, phase4_LAnkle_y,'--', t5, phase4_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('Phase 4'); % Taking steps from frontal to sagittal
subplot(5,1,5)
plot(t6, phase5_LAnkle_x, t6, phase5_LAnkle_y,'--', t6, phase5_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 5'); % Placing the box sagittal plane again


%% ROM calculation

% Right Ankle

RAnkle_ph1(:,1) = phase1_RAnkle_x;
RAnkle_ph1(:,2) = phase1_RAnkle_y;
RAnkle_ph1(:,3) = phase1_RAnkle_z;

RAnkle_ph2(:,1) = phase2_RAnkle_x;
RAnkle_ph2(:,2) = phase2_RAnkle_y;
RAnkle_ph2(:,3) = phase2_RAnkle_z;

RAnkle_ph3(:,1) = phase3_RAnkle_x;
RAnkle_ph3(:,2) = phase3_RAnkle_y;
RAnkle_ph3(:,3) = phase3_RAnkle_z;

RAnkle_ph4(:,1) = phase4_RAnkle_x;
RAnkle_ph4(:,2) = phase4_RAnkle_y;
RAnkle_ph4(:,3) = phase4_RAnkle_z;

RAnkle_ph5(:,1) = phase5_RAnkle_x;
RAnkle_ph5(:,2) = phase5_RAnkle_y;
RAnkle_ph5(:,3) = phase5_RAnkle_z;


RAnkle_ph1_max = max(RAnkle_ph1);
RAnkle_ph1_min = min(RAnkle_ph1);

RAnkle_ph2_max = max(RAnkle_ph2);
RAnkle_ph2_min = min(RAnkle_ph2);

RAnkle_ph3_max = max(RAnkle_ph3);
RAnkle_ph3_min = min(RAnkle_ph3);

RAnkle_ph4_max = max(RAnkle_ph4);
RAnkle_ph4_min = min(RAnkle_ph4);

RAnkle_ph5_max = max(RAnkle_ph5);
RAnkle_ph5_min = min(RAnkle_ph5);

for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph1_ROM(1,i) = RAnkle_ph1_max(i)- RAnkle_ph1_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph2_ROM(1,i) = RAnkle_ph2_max(i)- RAnkle_ph2_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph3_ROM(1,i) = RAnkle_ph3_max(i)- RAnkle_ph3_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph4_ROM(1,i) = RAnkle_ph4_max(i)- RAnkle_ph4_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph5_ROM(1,i) = RAnkle_ph5_max(i)- RAnkle_ph5_min(i);
end

% Left Ankle

LAnkle_ph1(:,1) = phase1_LAnkle_x;
LAnkle_ph1(:,2) = phase1_LAnkle_y;
LAnkle_ph1(:,3) = phase1_LAnkle_z;

LAnkle_ph2(:,1) = phase2_LAnkle_x;
LAnkle_ph2(:,2) = phase2_LAnkle_y;
LAnkle_ph2(:,3) = phase2_LAnkle_z;

LAnkle_ph3(:,1) = phase3_LAnkle_x;
LAnkle_ph3(:,2) = phase3_LAnkle_y;
LAnkle_ph3(:,3) = phase3_LAnkle_z;

LAnkle_ph4(:,1) = phase4_LAnkle_x;
LAnkle_ph4(:,2) = phase4_LAnkle_y;
LAnkle_ph4(:,3) = phase4_LAnkle_z;

LAnkle_ph5(:,1) = phase5_LAnkle_x;
LAnkle_ph5(:,2) = phase5_LAnkle_y;
LAnkle_ph5(:,3) = phase5_LAnkle_z;


LAnkle_ph1_max = max(LAnkle_ph1);
LAnkle_ph1_min = min(LAnkle_ph1);

LAnkle_ph2_max = max(LAnkle_ph2);
LAnkle_ph2_min = min(LAnkle_ph2);

LAnkle_ph3_max = max(LAnkle_ph3);
LAnkle_ph3_min = min(LAnkle_ph3);

LAnkle_ph4_max = max(LAnkle_ph4);
LAnkle_ph4_min = min(LAnkle_ph4);

LAnkle_ph5_max = max(LAnkle_ph5);
LAnkle_ph5_min = min(LAnkle_ph5);

for i = 1:length(RAnkle_ph1_max)
    LAnkle_ph1_ROM(1,i) = LAnkle_ph1_max(i)- LAnkle_ph1_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    LAnkle_ph2_ROM(1,i) = LAnkle_ph2_max(i)- LAnkle_ph2_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    LAnkle_ph3_ROM(1,i) = LAnkle_ph3_max(i)- LAnkle_ph3_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    LAnkle_ph4_ROM(1,i) = LAnkle_ph4_max(i)- LAnkle_ph4_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    LAnkle_ph5_ROM(1,i) = LAnkle_ph5_max(i)- LAnkle_ph5_min(i);
end

% Total ROM Ankle Angles

RAnkle_total(:,1) = RAnkleAngles_x;
RAnkle_total(:,2) = RAnkleAngles_y;
RAnkle_total(:,3) = RAnkleAngles_z;

RAnkle_total_max = max(RAnkle_total);
RAnkle_total_min = min(RAnkle_total);

for i = 1:length(RAnkle_total_max)
    RAnkle_total_ROM(1,i) = RAnkle_total_max(i)- RAnkle_total_min(i);
end

LAnkle_total(:,1) = LAnkleAngles_x;
LAnkle_total(:,2) = LAnkleAngles_y;
LAnkle_total(:,3) = LAnkleAngles_z;

LAnkle_total_max = max(LAnkle_total);
LAnkle_total_min = min(LAnkle_total);

for i = 1:length(LAnkle_total_max)
    LAnkle_total_ROM(1,i) = LAnkle_total_max(i)- LAnkle_total_min(i);
end

%ROM_Right_Ankle = table([RAnkle_ph1_ROM(1);RAnkle_ph2_ROM(1);RAnkle_ph3_ROM(1);RAnkle_ph4_ROM(1); RAnkle_ph5_ROM(1)],[RAnkle_ph1_ROM(2);RAnkle_ph2_ROM(2);RAnkle_ph3_ROM(2);RAnkle_ph4_ROM(2); RAnkle_ph5_ROM(2)],[RAnkle_ph1_ROM(3);RAnkle_ph2_ROM(3);RAnkle_ph3_ROM(3);RAnkle_ph4_ROM(3); RAnkle_ph5_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4', 'Phase 5'})
ROM_Right_Ankle = struct;
ROM_Right_Ankle.Phase1_X=RAnkle_ph1_ROM(1);
ROM_Right_Ankle.Phase1_Y=RAnkle_ph1_ROM(2);
ROM_Right_Ankle.Phase1_Z=RAnkle_ph1_ROM(3);
ROM_Right_Ankle.Phase2_X=RAnkle_ph2_ROM(1);
ROM_Right_Ankle.Phase2_Y=RAnkle_ph2_ROM(2);
ROM_Right_Ankle.Phase2_Z=RAnkle_ph2_ROM(3);
ROM_Right_Ankle.Phase3_X=RAnkle_ph3_ROM(1);
ROM_Right_Ankle.Phase3_Y=RAnkle_ph3_ROM(2);
ROM_Right_Ankle.Phase3_Z=RAnkle_ph3_ROM(3);
ROM_Right_Ankle.Phase4_X=RAnkle_ph4_ROM(1);
ROM_Right_Ankle.Phase4_Y=RAnkle_ph4_ROM(2);
ROM_Right_Ankle.Phase4_Z=RAnkle_ph4_ROM(3);   
ROM_Right_Ankle.Phase5_X=RAnkle_ph5_ROM(1);
ROM_Right_Ankle.Phase5_Y=RAnkle_ph5_ROM(2);
ROM_Right_Ankle.Phase5_Z=RAnkle_ph5_ROM(3);

ROM_Right_Ankle

%ROM_Left_Ankle = table([LAnkle_ph1_ROM(1);LAnkle_ph2_ROM(1);LAnkle_ph3_ROM(1);LAnkle_ph4_ROM(1); LAnkle_ph5_ROM(1)],[LAnkle_ph1_ROM(2);LAnkle_ph2_ROM(2);LAnkle_ph3_ROM(2);LAnkle_ph4_ROM(2); LAnkle_ph5_ROM(2)],[LAnkle_ph1_ROM(3);LAnkle_ph2_ROM(3);LAnkle_ph3_ROM(3);LAnkle_ph4_ROM(3); LAnkle_ph5_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4', 'Phase 5'})
ROM_Left_Ankle = struct;
ROM_Left_Ankle.Phase1_X=LAnkle_ph1_ROM(1);
ROM_Left_Ankle.Phase1_Y=LAnkle_ph1_ROM(2);
ROM_Left_Ankle.Phase1_Z=LAnkle_ph1_ROM(3);
ROM_Left_Ankle.Phase2_X=LAnkle_ph2_ROM(1);
ROM_Left_Ankle.Phase2_Y=LAnkle_ph2_ROM(2);
ROM_Left_Ankle.Phase2_Z=LAnkle_ph2_ROM(3);
ROM_Left_Ankle.Phase3_X=LAnkle_ph3_ROM(1);
ROM_Left_Ankle.Phase3_Y=LAnkle_ph3_ROM(2);
ROM_Left_Ankle.Phase3_Z=LAnkle_ph3_ROM(3);
ROM_Left_Ankle.Phase4_X=LAnkle_ph4_ROM(1);
ROM_Left_Ankle.Phase4_Y=LAnkle_ph4_ROM(2);
ROM_Left_Ankle.Phase4_Z=LAnkle_ph4_ROM(3);
ROM_Left_Ankle.Phase5_X=LAnkle_ph5_ROM(1);
ROM_Left_Ankle.Phase5_Y=LAnkle_ph5_ROM(2);
ROM_Left_Ankle.Phase5_Z=LAnkle_ph5_ROM(3);

ROM_Left_Ankle

%TotalROM_Ankle = table([RAnkle_total_ROM(1); LAnkle_total_ROM(1)],[RAnkle_total_ROM(2); LAnkle_total_ROM(2)],[RAnkle_total_ROM(3); LAnkle_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})
TotalROM_Ankle = struct;
TotalROM_Ankle.Right_X=RAnkle_total_ROM(1);
TotalROM_Ankle.Right_Y=RAnkle_total_ROM(2);
TotalROM_Ankle.Right_Z=RAnkle_total_ROM(3);
TotalROM_Ankle.Left_X=LAnkle_total_ROM(1);
TotalROM_Ankle.Left_Y=LAnkle_total_ROM(2);
TotalROM_Ankle.Left_Z=LAnkle_total_ROM(3);

TotalROM_Ankle

