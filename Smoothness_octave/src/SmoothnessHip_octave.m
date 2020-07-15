%% Smoothness of hip angular velocity signal: SPARC, NP an LDLJ

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

% Proposed only for Sagittal Lifting tasks
% Calculated for two phases named:
    % Lifting phase: Subject holds and deposits the box (idxOLF:idxFLF)
    % Lowering phase: Subject holds and deposits the box in its original location (idxOLW:idxFLW)

pkg load signal

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica59_B.mat')

Ts = 1/Fs;
parameters = [0.025; Fs/2; 4];  % Threshold 0.025 or less, between [0.015, 0.025]

%% Import model data for Hip Relative Angles

% Hip: 3 dofs

        RHipAngles_x = ModelData.Raw.(ModelOutput{1})(1,:); 
        RHipAngles_y = ModelData.Raw.(ModelOutput{1})(2,:);
        RHipAngles_z = ModelData.Raw.(ModelOutput{1})(3,:);
       

        LHipAngles_x = ModelData.Raw.(ModelOutput{2})(1,:);
        LHipAngles_y = ModelData.Raw.(ModelOutput{2})(2,:);
        LHipAngles_z = ModelData.Raw.(ModelOutput{2})(3,:);
        
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

figure(2)
plot(t, position_Z, t, RHipAngles_x, t, RHipAngles_y, t, RHipAngles_z, t, LHipAngles_x, t, LHipAngles_y, t, LHipAngles_z);
legend('Z', 'RHipAngles_x', 'RHipAngles_y', 'RHipAngles_z', 'LHipAngles_x', 'LHipAngles_y', 'LHipAngles_z');
xlabel('Time in seconds');
title('Vertical Box Trajectory and Hip Relative Angles');


%% Angular Velocities Calculation(rad/s)

    RHip_angVelocity_x = dv(t,RHipAngles_x);
    RHip_angVelocity_x = RHip_angVelocity_x'/Ts;
    RHip_angVelocity_y = dv(t,RHipAngles_y);
    RHip_angVelocity_y = RHip_angVelocity_y'/Ts;
    RHip_angVelocity_z = dv(t,RHipAngles_z);
    RHip_angVelocity_z = RHip_angVelocity_z'/Ts;

    LHip_angVelocity_x = dv(t,LHipAngles_x);
    LHip_angVelocity_x = LHip_angVelocity_x'/Ts;
    LHip_angVelocity_y = dv(t,LHipAngles_y);
    LHip_angVelocity_y = LHip_angVelocity_y'/Ts;
    LHip_angVelocity_z = dv(t,LHipAngles_z);
    LHip_angVelocity_z = LHip_angVelocity_z'/Ts;  
      
    
%% Crop Relative angular velocities to lifting phase
   
    RHip_angVelocity_x_lifting_sparc =  RHip_angVelocity_x(idxOLF:idxFLF);
    RHip_angVelocity_y_lifting_sparc =  RHip_angVelocity_y(idxOLF:idxFLF);
    RHip_angVelocity_z_lifting_sparc =  RHip_angVelocity_z(idxOLF:idxFLF);
    
    LHip_angVelocity_x_lifting_sparc =  LHip_angVelocity_x(idxOLF:idxFLF);
    LHip_angVelocity_y_lifting_sparc =  LHip_angVelocity_y(idxOLF:idxFLF);
    LHip_angVelocity_z_lifting_sparc =  LHip_angVelocity_z(idxOLF:idxFLF);
    
figure (3) 
t2 = (0:1:(length(RHip_angVelocity_x_lifting_sparc)-1))*Ts;
subplot (2,1,1)
plot(t2, RHip_angVelocity_x_lifting_sparc, t2, RHip_angVelocity_y_lifting_sparc, t2, RHip_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Right Hip Angular Velocities Lifting Phase') 
subplot (2,1,2)
plot(t2, LHip_angVelocity_x_lifting_sparc, t2, LHip_angVelocity_y_lifting_sparc, t2, LHip_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Left Hip Angular Velocities Lifting Phase')


%% Smoothness Lifting Hip: SPARC, NP and LDLJ

S_lifting_Hip = {'S_R_lifting_X', 'S_R_lifting_Y', 'S_R_lifting_Z', 'S_L_lifting_X', 'S_L_lifting_Y', 'S_L_lifting_Z'};
    Hip_angVel(1,:) =  RHip_angVelocity_x_lifting_sparc;
    Hip_angVel(2,:) =  RHip_angVelocity_y_lifting_sparc;
    Hip_angVel(3,:) =  RHip_angVelocity_z_lifting_sparc;
    Hip_angVel(4,:) =  LHip_angVelocity_x_lifting_sparc;
    Hip_angVel(5,:) =  LHip_angVelocity_y_lifting_sparc;
    Hip_angVel(6,:) =  LHip_angVelocity_z_lifting_sparc;

for i = 1:length(S_lifting_Hip)
    
    speed = Hip_angVel(i,:);
    %SPARC
    S_sparc = SpectralArcLength(speed, Ts, parameters);
    [S_lifting_Hip_SPARC{i}] = S_sparc;
    
    %NP
    [S_np, locs_np] = np_octave(t2,speed);
    [S_lifting_Hip_NP{i}] = S_np;
    [peaks_lifting{i}] = locs_np;
    
    %LDLJ
    [Dlj S_ldlj] = ldlj(t2,speed);
    [S_lifting_Hip_LDLJ{i}] = S_ldlj;
end 


%% Events of interest in vertical box trajectory identification

% Lowering phase

flat3 = idx < locs(2);
idx_flat3 = find(flat3);

flat4 = idx > locs(2);
idx_flat4 = find(flat4);

idxOLW = idx(idx_flat3(length(idx_flat3)));
idxFLW = idx(idx_flat4(3));                 % para din 59 es 4, para din 03 y 61 es 2 y para el resto es 1

figure(4)
plot(t,position_Z,t(locs(2)),position_Z(locs(2)),'o');
hold on
%plot(t,position_Z,'r*','MarkerIndices',idxOLW);
plot(t(idxOLW),position_Z(idxOLW),'r*');

hold on
%plot(t,position_Z,'r*','MarkerIndices',idxFLW);
plot(t(idxFLW),position_Z(idxFLW),'r*');

hold off
xlabel('Time in seconds');
ylabel('Millimetres');
title('Events of interest in box trajectory during lowering for SPARC');


%% Crop Relative angular velocities to lowering phase


    RHip_angVelocity_x_lowering_sparc = RHip_angVelocity_x(idxOLW:idxFLW);
    RHip_angVelocity_y_lowering_sparc = RHip_angVelocity_y(idxOLW:idxFLW);
    RHip_angVelocity_z_lowering_sparc = RHip_angVelocity_z(idxOLW:idxFLW);

    LHip_angVelocity_x_lowering_sparc =  LHip_angVelocity_x(idxOLW:idxFLW);
    LHip_angVelocity_y_lowering_sparc =  LHip_angVelocity_y(idxOLW:idxFLW);
    LHip_angVelocity_z_lowering_sparc =  LHip_angVelocity_z(idxOLW:idxFLW);
    
figure (5) 
t3 = 0:1:(length(RHip_angVelocity_x_lowering_sparc)-1);
subplot (2,1,1)
plot(t3, RHip_angVelocity_x_lowering_sparc, t3, RHip_angVelocity_y_lowering_sparc, t3, RHip_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Right Hip Angular Velocities Lowering Phase') 
subplot (2,1,2)
plot(t3, LHip_angVelocity_x_lowering_sparc, t3, LHip_angVelocity_y_lowering_sparc, t3, LHip_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Left Hip Angular Velocities Lowering Phase')

%% Smoothness Lowering Hip: SPARC, NP and LDLJ
   
S_lowering_Hip = {'S_R_lowering_X', 'S_R_lowering_Y', 'S_R_lowering_Z', 'S_L_lowering_X', 'S_L_lowering_Y', 'S_L_lowering_Z'};
    Hip_angVel_low(1,:) =  RHip_angVelocity_x_lowering_sparc;
    Hip_angVel_low(2,:) =  RHip_angVelocity_y_lowering_sparc;
    Hip_angVel_low(3,:) =  RHip_angVelocity_z_lowering_sparc;
    Hip_angVel_low(4,:) =  LHip_angVelocity_x_lowering_sparc;
    Hip_angVel_low(5,:) =  LHip_angVelocity_y_lowering_sparc;
    Hip_angVel_low(6,:) =  LHip_angVelocity_z_lowering_sparc;

for i = 1:length(S_lowering_Hip)
    
    speed = Hip_angVel_low(i,:);
    %SPARC
     S_sparc = SpectralArcLength(speed, Ts, parameters);
    [S_lowering_Hip_SPARC{i}] = S_sparc;
    
    %NP
    [S_np, locs_np] = np_octave(t3,speed);
    [S_lowering_Hip_NP{i}] = S_np;
    [peaks_lowering{i}] = locs_np;
    
    %LDLJ
    [Dlj S_ldlj] = ldlj(t3,speed);
    [S_lowering_Hip_LDLJ{i}] = S_ldlj;
end 
  

%Smoothness_Lifting_Hip = table([S_lifting_Hip_SPARC(1); S_lifting_Hip_SPARC(4); S_lifting_Hip_SPARC(2); S_lifting_Hip_SPARC(5);S_lifting_Hip_SPARC(3);S_lifting_Hip_SPARC(6)],[S_lifting_Hip_NP(1); S_lifting_Hip_NP(4); S_lifting_Hip_NP(2); S_lifting_Hip_NP(5);S_lifting_Hip_NP(3);S_lifting_Hip_NP(6)],[S_lifting_Hip_LDLJ(1); S_lifting_Hip_LDLJ(4); S_lifting_Hip_LDLJ(2); S_lifting_Hip_LDLJ(5);S_lifting_Hip_LDLJ(3);S_lifting_Hip_LDLJ(6)],'VariableNames',{'SPARC','NP','LDLJ'},'RowNames',{'Right X','Left  X','Right Y','Left  Y','Right Z','Left  Z'})
Smoothness_Lifting_Hip = struct;
Smoothness_Lifting_Hip.Right_X_SPARC=S_lifting_Hip_SPARC{1};
Smoothness_Lifting_Hip.Right_X_NP=S_lifting_Hip_NP{1};
Smoothness_Lifting_Hip.Right_X_LDLJ=S_lifting_Hip_LDLJ{1};
Smoothness_Lifting_Hip.Left_X_SPARC=S_lifting_Hip_SPARC{2};
Smoothness_Lifting_Hip.Left_X_NP=S_lifting_Hip_NP{2};
Smoothness_Lifting_Hip.Left_X_LDLJ=S_lifting_Hip_LDLJ{2};
Smoothness_Lifting_Hip.Right_Y_SPARC=S_lifting_Hip_SPARC{3};
Smoothness_Lifting_Hip.Right_Y_NP=S_lifting_Hip_NP{3};
Smoothness_Lifting_Hip.Right_Y_LDLJ=S_lifting_Hip_LDLJ{3};
Smoothness_Lifting_Hip.Left_Y_SPARC=S_lifting_Hip_SPARC{4};
Smoothness_Lifting_Hip.Left_Y_NP=S_lifting_Hip_NP{4};
Smoothness_Lifting_Hip.Left_Y_LDLJ=S_lifting_Hip_LDLJ{4};
Smoothness_Lifting_Hip.Right_Z_SPARC=S_lifting_Hip_SPARC{5};
Smoothness_Lifting_Hip.Right_Z_NP=S_lifting_Hip_NP{5};
Smoothness_Lifting_Hip.Right_Z_LDLJ=S_lifting_Hip_LDLJ{5};
Smoothness_Lifting_Hip.Left_Z_SPARC=S_lifting_Hip_SPARC{6};
Smoothness_Lifting_Hip.Left_Z_NP=S_lifting_Hip_NP{6};
Smoothness_Lifting_Hip.Left_Z_LDLJ=S_lifting_Hip_LDLJ{6};

Smoothness_Lifting_Hip

%Smoothness_Lowering_Hip = table([S_lowering_Hip_SPARC(1); S_lowering_Hip_SPARC(4); S_lowering_Hip_SPARC(2); S_lowering_Hip_SPARC(5);S_lowering_Hip_SPARC(3);S_lowering_Hip_SPARC(6)],[S_lowering_Hip_NP(1); S_lowering_Hip_NP(4); S_lowering_Hip_NP(2); S_lowering_Hip_NP(5);S_lowering_Hip_NP(3);S_lowering_Hip_NP(6)],[S_lowering_Hip_LDLJ(1); S_lowering_Hip_LDLJ(4); S_lowering_Hip_LDLJ(2); S_lowering_Hip_LDLJ(5);S_lowering_Hip_LDLJ(3);S_lowering_Hip_LDLJ(6)],'VariableNames',{'SPARC','NP','LDLJ'},'RowNames',{'Right X','Left  X','Right Y','Left  Y','Right Z','Left  Z'})
Smoothness_Lowering_Hip = struct;
Smoothness_Lowering_Hip.Right_X_SPARC=S_lowering_Hip_SPARC{1};
Smoothness_Lowering_Hip.Right_X_NP=S_lowering_Hip_NP{1};
Smoothness_Lowering_Hip.Right_X_LDLJ=S_lowering_Hip_LDLJ{1};
Smoothness_Lowering_Hip.Left_X_SPARC=S_lowering_Hip_SPARC{2};
Smoothness_Lowering_Hip.Left_X_NP=S_lowering_Hip_NP{2};
Smoothness_Lowering_Hip.Left_X_LDLJ=S_lowering_Hip_LDLJ{2};
Smoothness_Lowering_Hip.Right_Y_SPARC=S_lowering_Hip_SPARC{3};
Smoothness_Lowering_Hip.Right_Y_NP=S_lowering_Hip_NP{3};
Smoothness_Lowering_Hip.Right_Y_LDLJ=S_lowering_Hip_LDLJ{3};
Smoothness_Lowering_Hip.Left_Y_SPARC=S_lowering_Hip_SPARC{4};
Smoothness_Lowering_Hip.Left_Y_NP=S_lowering_Hip_NP{4};
Smoothness_Lowering_Hip.Left_Y_LDLJ=S_lowering_Hip_LDLJ{4};
Smoothness_Lowering_Hip.Right_Z_SPARC=S_lowering_Hip_SPARC{5};
Smoothness_Lowering_Hip.Right_Z_NP=S_lowering_Hip_NP{5};
Smoothness_Lowering_Hip.Right_Z_LDLJ=S_lowering_Hip_LDLJ{5};
Smoothness_Lowering_Hip.Left_Z_SPARC=S_lowering_Hip_SPARC{6};
Smoothness_Lowering_Hip.Left_Z_NP=S_lowering_Hip_NP{6};
Smoothness_Lowering_Hip.Left_Z_LDLJ=S_lowering_Hip_LDLJ{6};

Smoothness_Lowering_Hip
