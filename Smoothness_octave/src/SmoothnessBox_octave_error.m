%% Smoothness of box trajectory marker

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

% Proposed only for Sagittal Lifting tasks
% It is observed in Vicon captures that when the subject deposits a loaded box 
% it bounces. Thus, measuring its smoothness can show 
% differences regarding precision or delicacy between placing a loaded box
% and a unloaded one. 

pkg load signal

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica59_B.mat')

Ts = 1/Fs;
parameters = [0.015; Fs/2; 4];   % threshold between [0.015,0.020]


%% Vertical SIDEBOX1 box marker trajectory

    position_Z = SIDEBOX1(:,3)';
    t = (0:1:(length(position_Z)-1))*Ts;
    
%% Angular Velocity Calculation(rad/s)

    Box_angVelocity_z = dv(t,position_Z);
    Box_angVelocity_z = Box_angVelocity_z'/Ts;
    
    
%% Events of interest in vertical box trajectory identification
 
position_Z = SIDEBOX3(:,3)';
 
%[pks, locs] = findpeaks(position_Z,'minPeakProminence',4);
[pks, locs] = findpeaks(position_Z,'DoubleSided');   
%locs=locs(pks>80);                 
%pks=pks(pks>80);
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
pks=pks(prominence>4);    
locs=locs(prominence>4);





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
end


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
 
 
% Indexes definiton for signal segmentation
idxOLF = idx(idx_flat(length(idx_flat)));
idxFLF = idx(idx_flat2(1));

 
figure(1)
plot(t,position_Z,t(locs(1)),position_Z(locs(1)),'o');
hold on
%plot(t,position_Z,'r*','MarkerIndices',idxOLF);
plot(t(idxOLF),position_Z(idxOLF),'r*');

hold on
%plot(t,position_Z,'r*','MarkerIndices',idxFLF);
plot(t(idxFLF),position_Z(idxFLF),'r*');

hold off
xlabel('Time in seconds');
ylabel('Millimetres');
title('Events of interest in box trajectory during lifting for SPARC');

% Lowering phase

flat3 = idx < locs(2);
idx_flat3 = find(flat3);

flat4 = idx > locs(2);
idx_flat4 = find(flat4);

idxOLW = idx(idx_flat3(length(idx_flat3)));
idxFLW = idx(idx_flat4(3));

figure(2)
plot(t,position_Z,t(locs(2)),position_Z(locs(2)),'o');
hold on
%plot(t,position_Z,'r*','MarkerIndices',idxOLW );
plot(t(idxOLW),position_Z(idxOLW),'r*');

hold on
%plot(t,position_Z,'r*','MarkerIndices',idxFLW);
plot(t(idxFLW),position_Z(idxFLW),'r*');

hold off
xlabel('Time in seconds');
ylabel('Millimetres');
title('Events of interest in box trajectory during lowering for SPARC');


%% Crop relative angular velocities to lifting phase

    Box_angVelocity_z_lifting_sparc =  Box_angVelocity_z(idxOLF:idxFLF);
    
    Box_angVelocity_z_lowering_sparc =  Box_angVelocity_z(idxOLW:idxFLW);

    %t2 = (0:1:(length(Box_angVelocity_z_lifting_sparc)-1))*Ts;
    
    %t3 = (0:1:(length(Box_angVelocity_z_lowering_sparc)-1))*Ts;
    
%% SPARC algorithm: smoothness calculation

    speed_lifting = Box_angVelocity_z_lifting_sparc(1,:);
    S_lifting_SPARC = SpectralArcLength(speed_lifting, Ts, parameters);
    
    speed_lowering = Box_angVelocity_z_lowering_sparc(1,:);
    S_lowering_SPARC = SpectralArcLength(speed_lowering, Ts, parameters);
    
%% NP
    
    speed_lf_np = Box_angVelocity_z_lifting_sparc(1,:);
    [S_lifting_NP, locs_np] = np(t2,speed_lf_np);
   
    speed_lw_np = Box_angVelocity_z_lowering_sparc(1,:);
    [S_lowering_NP, locs2_np] = np(t3,speed_lw_np);

%% LDLJ
    
    speed_lf_ldlj = Box_angVelocity_z_lifting_sparc(1,:);
    [Dlj S_lifting_LDLJ] = ldlj(t2,speed_lf_ldlj);
       
    speed_lw_ldlj = Box_angVelocity_z_lowering_sparc(1,:);
    [Dlj S_lowering_LDLJ] = ldlj(t3,speed_lw_ldlj);


%% Table

%Smoothness_Lifting_Box = table([S_lifting_SPARC],[S_lifting_NP],[S_lifting_LDLJ],'VariableNames',{'SPARC','NP','LDLJ'},'RowNames',{'Box Z'})
Smoothness_Lifting_Box = struct;
Smoothness_Lifting_Box.Right_X_SPARC=S_lifting_SPARC{1};
Smoothness_Lifting_Box.Right_X_NP=S_lifting_NP{1};
Smoothness_Lifting_Box.Right_X_LDLJ=S_lifting_LDLJ{1};
Smoothness_Lifting_Box.Left_X_SPARC=S_lifting_SPARC{2};
Smoothness_Lifting_Box.Left_X_NP=S_lifting_NP{2};
Smoothness_Lifting_Box.Left_X_LDLJ=S_lifting_LDLJ{2};
Smoothness_Lifting_Box.Right_Y_SPARC=S_lifting_SPARC{3};
Smoothness_Lifting_Box.Right_Y_NP=S_lifting_NP{3};
Smoothness_Lifting_Box.Right_Y_LDLJ=S_lifting_LDLJ{3};
Smoothness_Lifting_Box.Left_Y_SPARC=S_lifting_SPARC{4};
Smoothness_Lifting_Box.Left_Y_NP=S_lifting_NP{4};
Smoothness_Lifting_Box.Left_Y_LDLJ=S_lifting_LDLJ{4};
Smoothness_Lifting_Box.Right_Z_SPARC=S_lifting_SPARC{5};
Smoothness_Lifting_Box.Right_Z_NP=S_lifting_NP{5};
Smoothness_Lifting_Box.Right_Z_LDLJ=S_lifting_LDLJ{5};
Smoothness_Lifting_Box.Left_Z_SPARC=S_lifting_SPARC{6};
Smoothness_Lifting_Box.Left_Z_NP=S_lifting_NP{6};
Smoothness_Lifting_Box.Left_Z_LDLJ=S_lifting_LDLJ{6};

Smoothness_Lifting_Box

%Smoothness_Lowering_Box = table([S_lowering_SPARC],[S_lowering_NP],[S_lowering_LDLJ],'VariableNames',{'SPARC','NP','LDLJ'},'RowNames',{'Box Z'})
Smoothness_Lowering_Box = struct;
Smoothness_Lowering_Box.Right_X_SPARC=S_lowering_SPARC{1};
Smoothness_Lowering_Box.Right_X_NP=S_lowering_NP{1};
Smoothness_Lowering_Box.Right_X_LDLJ=S_lowering_LDLJ{1};
Smoothness_Lowering_Box.Left_X_SPARC=S_lowering_SPARC{2};
Smoothness_Lowering_Box.Left_X_NP=S_lowering_NP{2};
Smoothness_Lowering_Box.Left_X_LDLJ=S_lowering_LDLJ{2};
Smoothness_Lowering_Box.Right_Y_SPARC=S_lowering_SPARC{3};
Smoothness_Lowering_Box.Right_Y_NP=S_lowering_NP{3};
Smoothness_Lowering_Box.Right_Y_LDLJ=S_lowering_LDLJ{3};
Smoothness_Lowering_Box.Left_Y_SPARC=S_lowering_SPARC{4};
Smoothness_Lowering_Box.Left_Y_NP=S_lowering_NP{4};
Smoothness_Lowering_Box.Left_Y_LDLJ=S_lowering_LDLJ{4};
Smoothness_Lowering_Box.Right_Z_SPARC=S_lowering_SPARC{5};
Smoothness_Lowering_Box.Right_Z_NP=S_lowering_NP{5};
Smoothness_Lowering_Box.Right_Z_LDLJ=S_lowering_LDLJ{5};
Smoothness_Lowering_Box.Left_Z_SPARC=S_lowering_SPARC{6};
Smoothness_Lowering_Box.Left_Z_NP=S_lowering_NP{6};
Smoothness_Lowering_Box.Left_Z_LDLJ=S_lowering_LDLJ{6};

Smoothness_Lowering_Box



