vicon = ViconNexus;
SubjectName = vicon.GetSubjectNames;
Fs = vicon.GetFrameRate; % Sampling frequency
Ts = 1/Fs;
frames = vicon.GetFrameCount;



% Import SIDEBOX3 trajectory
    % Input = Subject name, marker name
    % Output = Trajectory data X, Y, Z & E
[SIDEBOX3(:,1),SIDEBOX3(:,2),SIDEBOX3(:,3),SIDEBOX3(:,4)] = vicon.GetTrajectory(SubjectName{1},'SIDEBOX3');

% Import model data for right wrist, interest in X and Y due to wrist's dofs
ModelOutput = {'RHipAngles','LHipAngles','RKneeAngles','LKneeAngles','RAnkleAngles','LAnkleAngles','RElbowAngles', 'LElbowAngles'};

for i = 1:length(ModelOutput)
    [ModelData.Raw.(ModelOutput{i}), ModelData.Exists.(ModelOutput{i})] = vicon.GetModelOutput(SubjectName{1},ModelOutput{i});
end    
% Hip 1,2. 3 dofs PENSAR HACERLO CON UN SWITCH CASE
        RHipAngles_x = ModelData.Raw.(ModelOutput{1})(1,:); %CUIDADO PORQUE PARA SPARC QUIERO VELOCIDAD
        RHipAngles_x = RHipAngles_x(8:end);
        RHipAngles_y = ModelData.Raw.(ModelOutput{1})(2,:);
        RHipAngles_y = RHipAngles_y(8:end);
        RHipAngles_z = ModelData.Raw.(ModelOutput{1})(3,:);
        RHipAngles_z = RHipAngles_z(8:end);

        LHipAngles_x = ModelData.Raw.(ModelOutput{2})(1,:);
        LHipAngles_x = LHipAngles_x(8:end);
        LHipAngles_y = ModelData.Raw.(ModelOutput{2})(2,:);
        LHipAngles_y = LHipAngles_y(8:end);
        LHipAngles_z = ModelData.Raw.(ModelOutput{2})(3,:);
        LHipAngles_z = LHipAngles_z(8:end);
        
% Knee 3,4. 2 dofs: flexoext y rotacion, transversal(y) y longitudinal (Z)
        RKneeAngles_x = ModelData.Raw.(ModelOutput{3})(1,:); 
        RKneeAngles_x = RKneeAngles_x(8:end);
        RKneeAngles_y = ModelData.Raw.(ModelOutput{3})(2,:);
        RKneeAngles_y = RKneeAngles_y(8:end);
        RKneeAngles_z = ModelData.Raw.(ModelOutput{3})(3,:);
        RKneeAngles_z = RKneeAngles_z(8:end);

        LKneeAngles_x = ModelData.Raw.(ModelOutput{4})(1,:); 
        LKneeAngles_x = LKneeAngles_x(8:end);
        LKneeAngles_y = ModelData.Raw.(ModelOutput{4})(2,:);
        LKneeAngles_y = LKneeAngles_y(8:end);
        LKneeAngles_z = ModelData.Raw.(ModelOutput{4})(3,:);
        LKneeAngles_z = LKneeAngles_z(8:end);
        
% Ankle 5,6. 1 dof: flexo ext, transversal(y)
        RAnkleAngles_x = ModelData.Raw.(ModelOutput{5})(1,:); 
        RAnkleAngles_x = RAnkleAngles_x(8:end);
        RAnkleAngles_y = ModelData.Raw.(ModelOutput{5})(2,:);
        RAnkleAngles_y = RAnkleAngles_y(8:end);
        RAnkleAngles_z = ModelData.Raw.(ModelOutput{5})(3,:);
        RAnkleAngles_z = RAnkleAngles_z(8:end);

        LAnkleAngles_x = ModelData.Raw.(ModelOutput{6})(1,:); 
        LAnkleAngles_x = LAnkleAngles_x(8:end);
        LAnkleAngles_y = ModelData.Raw.(ModelOutput{6})(2,:);
        LAnkleAngles_y = LAnkleAngles_y(8:end);
        LAnkleAngles_z = ModelData.Raw.(ModelOutput{6})(3,:);
        LAnkleAngles_z = LAnkleAngles_z(8:end);
        
% Elbow 7,8. 2 dofs: flexo ext,pronosupinacion, transversal(y)y longitudinal (z)
        RElbowAngles_x = ModelData.Raw.(ModelOutput{7})(1,:); 
        RElbowAngles_x = RElbowAngles_x(8:end);
        RElbowAngles_y = ModelData.Raw.(ModelOutput{7})(2,:);
        RElbowAngles_y = RElbowAngles_y(8:end);
        RElbowAngles_z = ModelData.Raw.(ModelOutput{7})(3,:);
        RElbowAngles_z = RElbowAngles_z(8:end);

        LElbowAngles_x = ModelData.Raw.(ModelOutput{8})(1,:); 
        LElbowAngles_x = LElbowAngles_x(8:end);
        LElbowAngles_y = ModelData.Raw.(ModelOutput{8})(2,:);
        LElbowAngles_y = LElbowAngles_y(8:end);
        LElbowAngles_z = ModelData.Raw.(ModelOutput{8})(3,:);
        LElbowAngles_z = LElbowAngles_z(8:end);

% ANGULAR VELOCITIES

% Hip Angular Velocities x, y, z
    RHip_angVelocity_x = diff((RHipAngles_x), 1);
    RHip_angVelocity_x = RHip_angVelocity_x/Ts;
    RHip_angVelocity_y = diff((RHipAngles_y), 1);
    RHip_angVelocity_y = RHip_angVelocity_y/Ts;
    RHip_angVelocity_z = diff((RHipAngles_z), 1);
    RHip_angVelocity_z = RHip_angVelocity_z/Ts;

    LHip_angVelocity_x = diff((LHipAngles_x), 1);
    LHip_angVelocity_x = LHip_angVelocity_x/Ts;
    LHip_angVelocity_y = diff((LHipAngles_y), 1);
    LHip_angVelocity_y = LHip_angVelocity_y/Ts;
    LHip_angVelocity_z = diff((LHipAngles_z), 1);
    LHip_angVelocity_z = LHip_angVelocity_z/Ts;

% Knee Angular Velocities x, y, z
    RKnee_angVelocity_x = diff((RKneeAngles_x), 1);
    RKnee_angVelocity_x = RKnee_angVelocity_x/Ts;
    RKnee_angVelocity_y = diff((RKneeAngles_y), 1);
    RKnee_angVelocity_y = RKnee_angVelocity_y/Ts;
    RKnee_angVelocity_z = diff((RKneeAngles_z), 1);
    RKnee_angVelocity_z = RKnee_angVelocity_z/Ts;

    LKnee_angVelocity_x = diff((LKneeAngles_x), 1);
    LKnee_angVelocity_x = LKnee_angVelocity_x/Ts;
    LKnee_angVelocity_y = diff((LKneeAngles_y), 1);
    LKnee_angVelocity_y = LKnee_angVelocity_y/Ts;
    LKnee_angVelocity_z = diff((LKneeAngles_z), 1);
    LKnee_angVelocity_z = LKnee_angVelocity_z/Ts;

% Ankle Angular Velocities x, y, z
    RAnkle_angVelocity_x = diff((RAnkleAngles_x), 1);
    RAnkle_angVelocity_x = RAnkle_angVelocity_x/Ts;
    RAnkle_angVelocity_y = diff((RAnkleAngles_y), 1);
    RAnkle_angVelocity_y = RAnkle_angVelocity_y/Ts;
    RAnkle_angVelocity_z = diff((RAnkleAngles_z), 1);
    RAnkle_angVelocity_z = RAnkle_angVelocity_z/Ts;

    LAnkle_angVelocity_x = diff((LAnkleAngles_x), 1);
    LAnkle_angVelocity_x = LAnkle_angVelocity_x/Ts;
    LAnkle_angVelocity_y = diff((LAnkleAngles_y), 1);
    LAnkle_angVelocity_y = LAnkle_angVelocity_y/Ts;
    LAnkle_angVelocity_z = diff((LAnkleAngles_z), 1);
    LAnkle_angVelocity_z = LAnkle_angVelocity_z/Ts;

% Elbow Angular Velocities x, y, z
    RElbow_angVelocity_x = diff((RElbowAngles_x), 1);
    RElbow_angVelocity_x = RElbow_angVelocity_x/Ts;
    RElbow_angVelocity_y = diff((RElbowAngles_y), 1);
    RElbow_angVelocity_y = RElbow_angVelocity_y/Ts;
    RElbow_angVelocity_z = diff((RElbowAngles_z), 1);
    RElbow_angVelocity_z = RElbow_angVelocity_z/Ts;

    LElbow_angVelocity_x = diff((LElbowAngles_x), 1);
    LElbow_angVelocity_x = LElbow_angVelocity_x/Ts;
    LElbow_angVelocity_y = diff((LElbowAngles_y), 1);
    LElbow_angVelocity_y = LElbow_angVelocity_y/Ts;
    LElbow_angVelocity_z = diff((LElbowAngles_z), 1);
    LElbow_angVelocity_z = LElbow_angVelocity_z/Ts;

% Vertical box marker trajectories
position_Z = SIDEBOX3(:,3)';
position_Z = position_Z(8:end);

% SEGMENTATION LIFTING: Encontrar minimos y maximos

% PREVIO 
[pks, locs] = findpeaks(position_Z,'minPeakProminence',20);
t = 0:1:(length(position_Z)-1);

TF1 = islocalmin(position_Z, 'FlatSelection', 'first');
idx = find(TF1);
flat = idx < locs(1);
idx_flat = find(flat);
flat2 = idx > locs(1);
idx_flat2 = find(flat2);

figure(1)
plot(t,position_Z,t(locs(1)),position_Z(locs(1)),'o');
hold on
plot(t,position_Z,'r*','MarkerIndices',idx(idx_flat(length(idx_flat))));
hold on
plot(t,position_Z,'r*','MarkerIndices',idx(idx_flat2(1)));
hold off
title('Events of interest  for SPARC pi while lifting in box trajectories');

idxOLF = idx(idx_flat(length(idx_flat)));
idxFLF = idx(idx_flat2(1));

figure(2)
subplot (2,2,1)
plot(t, position_Z, t, RHipAngles_x, t, RHipAngles_y, t, RHipAngles_z, t, LHipAngles_x, t, LHipAngles_y, t, LHipAngles_z);
legend('Z', 'RHipAngles_x', 'RHipAngles_y', 'RHipAngles_z', 'LHipAngles_x', 'LHipAngles_y', 'LHipAngles_z');
title(' Vertical Box Trajectory and Hip Relative Angles');

subplot (2,2,2)
plot(t, position_Z, t, RKneeAngles_x, t, RKneeAngles_y, t, RKneeAngles_z, t, LKneeAngles_x, t, LKneeAngles_y, t, LKneeAngles_z);
legend('Z', 'RKneeAngles_x', 'RKneeAngles_y', 'RKneeAngles_z', 'LKneeAngles_x', 'LKneeAngles_y', 'LKneeAngles_z');
title('Vertical Box Trajectory and Knee Relative Angles');

subplot (2,2,3)
plot(t, position_Z, t, RAnkleAngles_x, t, RAnkleAngles_y, t, RAnkleAngles_z, t, LAnkleAngles_x, t, LAnkleAngles_y, t, LAnkleAngles_z);
legend('Z', 'RAnkleAngles_x', 'RAnkleAngles_y', 'RAnkleAngles_z', 'LAnkleAngles_x', 'LAnkleAngles_y', 'LAnkleAngles_z');
title('Vertical Box Trajectory and Ankle Relative Angles');

subplot (2,2,4)
plot(t, position_Z, t, RElbowAngles_x, t, RElbowAngles_y, t, RElbowAngles_z, t, LElbowAngles_x, t, LElbowAngles_y, t, LElbowAngles_z);
legend('Z', 'RElbowAngles_x', 'RElbowAngles_y', 'RElbowAngles_z', 'LElbowAngles_x', 'LElbowAngles_y', 'LElbowAngles_z');
title('Vertical Box Trajectory and Elbow Relative Angles');

% Crop Relative angular velocities to lift and deposit box phase during lifting
% from index first flat to index flat signal again
    % HIP 
    RHip_angVelocity_x_lifting_sparc =  RHip_angVelocity_x(idxOLF:idxFLF);
    RHip_angVelocity_y_lifting_sparc =  RHip_angVelocity_y(idxOLF:idxFLF);
    RHip_angVelocity_z_lifting_sparc =  RHip_angVelocity_z(idxOLF:idxFLF);
    
    LHip_angVelocity_x_lifting_sparc =  LHip_angVelocity_x(idxOLF:idxFLF);
    LHip_angVelocity_y_lifting_sparc =  LHip_angVelocity_y(idxOLF:idxFLF);
    LHip_angVelocity_z_lifting_sparc =  LHip_angVelocity_z(idxOLF:idxFLF);
    
    % KNEE 
    RKnee_angVelocity_x_lifting_sparc =  RKnee_angVelocity_x(idxOLF:idxFLF);
    RKnee_angVelocity_y_lifting_sparc =  RKnee_angVelocity_y(idxOLF:idxFLF);
    RKnee_angVelocity_z_lifting_sparc =  RKnee_angVelocity_z(idxOLF:idxFLF);

    LKnee_angVelocity_x_lifting_sparc =  LKnee_angVelocity_x(idxOLF:idxFLF);
    LKnee_angVelocity_y_lifting_sparc =  LKnee_angVelocity_y(idxOLF:idxFLF);
    LKnee_angVelocity_z_lifting_sparc =  LKnee_angVelocity_z(idxOLF:idxFLF);
    
    % ANKLE 
    RAnkle_angVelocity_x_lifting_sparc =  RAnkle_angVelocity_x(idxOLF:idxFLF);
    RAnkle_angVelocity_y_lifting_sparc =  RAnkle_angVelocity_y(idxOLF:idxFLF);
    RAnkle_angVelocity_z_lifting_sparc =  RAnkle_angVelocity_z(idxOLF:idxFLF);

    LAnkle_angVelocity_x_lifting_sparc =  LAnkle_angVelocity_x(idxOLF:idxFLF);
    LAnkle_angVelocity_y_lifting_sparc =  LAnkle_angVelocity_y(idxOLF:idxFLF);
    LAnkle_angVelocity_z_lifting_sparc =  LAnkle_angVelocity_z(idxOLF:idxFLF);

    % ELBOW
    RElbow_angVelocity_x_lifting_sparc =  RElbow_angVelocity_x(idxOLF:idxFLF);
    RElbow_angVelocity_y_lifting_sparc =  RElbow_angVelocity_y(idxOLF:idxFLF);
    RElbow_angVelocity_z_lifting_sparc =  RElbow_angVelocity_z(idxOLF:idxFLF);

    LElbow_angVelocity_x_lifting_sparc =  LElbow_angVelocity_x(idxOLF:idxFLF);
    LElbow_angVelocity_y_lifting_sparc =  LElbow_angVelocity_y(idxOLF:idxFLF);
    LElbow_angVelocity_z_lifting_sparc =  LElbow_angVelocity_z(idxOLF:idxFLF);
    
figure (3) 
t3 = 0:1:(length(RHip_angVelocity_x_lifting_sparc)-1);
subplot (4,2,1)
plot(t3, RHip_angVelocity_x_lifting_sparc, t3, RHip_angVelocity_y_lifting_sparc, t3, RHip_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
title ('Cropped Right Hip Angular Velocities Lifting Phase') 
subplot (4,2,2)
plot(t3, LHip_angVelocity_x_lifting_sparc, t3, LHip_angVelocity_y_lifting_sparc, t3, LHip_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
title ('Cropped Left Hip Angular Velocities Lifting Phase')
subplot (4,2,3)
plot(t3, RKnee_angVelocity_x_lifting_sparc, t3, RKnee_angVelocity_y_lifting_sparc, t3, RKnee_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
title ('Cropped Right Knee Angular Velocities Lifting Phase')
subplot (4,2,4)
plot(t3, LKnee_angVelocity_x_lifting_sparc, t3, LKnee_angVelocity_y_lifting_sparc, t3, LKnee_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
title ('Cropped Left Knee Angular Velocities Lifting Phase')
subplot (4,2,5)
plot(t3, RAnkle_angVelocity_x_lifting_sparc, t3, RAnkle_angVelocity_y_lifting_sparc, t3, RAnkle_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
title ('Cropped Right Ankle Angular Velocities Lifting Phase')
subplot (4,2,6)
plot(t3, LAnkle_angVelocity_x_lifting_sparc, t3, LAnkle_angVelocity_y_lifting_sparc, t3, LAnkle_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
title ('Cropped Left Ankle Angular Velocities Lifting Phase')
subplot (4,2,7)
plot(t3, RElbow_angVelocity_x_lifting_sparc, t3, RElbow_angVelocity_y_lifting_sparc, t3, RElbow_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
title ('Cropped Right Elbow Angular Velocities Lifting Phase')
subplot (4,2,8)
plot(t3, LElbow_angVelocity_x_lifting_sparc, t3, LElbow_angVelocity_y_lifting_sparc, t3, LElbow_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
title ('Cropped Left Elbow Angular Velocities Lifting Phase')

% SPARC pi
    % HIP 
S_lifting_Hip = {'S_R_lifting_X', 'S_R_lifting_Y', 'S_R_lifting_Z', 'S_L_lifting_X', 'S_L_lifting_Y', 'S_L_lifting_Z'};
    Hip_angVel(1,:) =  RHip_angVelocity_x_lifting_sparc;
    Hip_angVel(2,:) =  RHip_angVelocity_y_lifting_sparc;
    Hip_angVel(3,:) =  RHip_angVelocity_z_lifting_sparc;
    Hip_angVel(4,:) =  LHip_angVelocity_x_lifting_sparc;
    Hip_angVel(5,:) =  LHip_angVelocity_y_lifting_sparc;
    Hip_angVel(6,:) =  LHip_angVelocity_z_lifting_sparc;

for i = 1:length(S_lifting_Hip)
    speed = Hip_angVel(i,:);
    S = SpectralArcLength(speed);
    [S_lifting_Hip{i}] = S;
end    

    % KNEE
S_lifting_Knee = {'S_R_lifting_X', 'S_R_lifting_Y', 'S_R_lifting_Z', 'S_L_lifting_X', 'S_L_lifting_Y', 'S_L_lifting_Z'};
    Knee_angVel(1,:) = RKnee_angVelocity_x_lifting_sparc;
    Knee_angVel(2,:) = RKnee_angVelocity_y_lifting_sparc;
    Knee_angVel(3,:) =  RKnee_angVelocity_z_lifting_sparc;
    Knee_angVel(4,:) =  LKnee_angVelocity_x_lifting_sparc;
    Knee_angVel(5,:) =  LKnee_angVelocity_y_lifting_sparc;
    Knee_angVel(6,:) =  LKnee_angVelocity_z_lifting_sparc;

for i = 1:length(S_lifting_Knee)
    speed = Knee_angVel(i,:);
    S = SpectralArcLength(speed);
    [S_lifting_Knee{i}] = S;
end 

% ANKLE
S_lifting_Ankle = {'S_R_lifting_X', 'S_R_lifting_Y', 'S_R_lifting_Z', 'S_L_lifting_X', 'S_L_lifting_Y', 'S_L_lifting_Z'};
    Ankle_angVel(1,:) = RAnkle_angVelocity_x_lifting_sparc;
    Ankle_angVel(2,:) = RAnkle_angVelocity_y_lifting_sparc;
    Ankle_angVel(3,:) = RAnkle_angVelocity_z_lifting_sparc;
    Ankle_angVel(4,:) = LAnkle_angVelocity_x_lifting_sparc;
    Ankle_angVel(5,:) = LAnkle_angVelocity_y_lifting_sparc;
    Ankle_angVel(6,:) = LAnkle_angVelocity_z_lifting_sparc;

for i = 1:length(S_lifting_Ankle)
    speed = Ankle_angVel(i,:);
    S = SpectralArcLength(speed);
    [S_lifting_Ankle{i}] = S;
end

% ELBOW
S_lifting_Elbow = {'S_R_lifting_X', 'S_R_lifting_Y', 'S_R_lifting_Z', 'S_L_lifting_X', 'S_L_lifting_Y', 'S_L_lifting_Z'};
    Elbow_angVel(1,:) = RElbow_angVelocity_x_lifting_sparc;
    Elbow_angVel(2,:) = RElbow_angVelocity_y_lifting_sparc;
    Elbow_angVel(3,:) = RElbow_angVelocity_z_lifting_sparc;
    Elbow_angVel(4,:) = LElbow_angVelocity_x_lifting_sparc;
    Elbow_angVel(5,:) = LElbow_angVelocity_y_lifting_sparc;
    Elbow_angVel(6,:) = LElbow_angVelocity_z_lifting_sparc;

for i = 1:length(S_lifting_Elbow)
    speed = Elbow_angVel(i,:);
    S = SpectralArcLength(speed);
    [S_lifting_Elbow{i}] = S;
end 

% SEGMENTATION LOWERING: Encontrar minimos y maximos

% Ultimo indice flat antes segundo maximo Z
flat3 = idx < locs(2);
idx_flat3 = find(flat3);
% Primer indice flat despues segundo maximo Z
flat4 = idx > locs(2);
idx_flat4 = find(flat4);

figure(4)
plot(t,position_Z,t(locs(2)),position_Z(locs(2)),'o');
hold on
plot(t,position_Z,'r*','MarkerIndices',idx(idx_flat3(length(idx_flat3))));
hold on
plot(t,position_Z,'r*','MarkerIndices',idx(idx_flat4(1)));
hold off
title('Events of interest  for SPARC pi while lowering in box trajectories');

idxOLW = idx(idx_flat3(length(idx_flat3)));
idxFLW = idx(idx_flat4(1));

% Crop Relative angular velocities to lower and deposit box phase during lifting
% from index first flat to index flat signal again
    % HIP 
    RHip_angVelocity_x_lowering_sparc = RHip_angVelocity_x(idxOLW:idxFLW);
    RHip_angVelocity_y_lowering_sparc = RHip_angVelocity_y(idxOLW:idxFLW);
    RHip_angVelocity_z_lowering_sparc = RHip_angVelocity_z(idxOLW:idxFLW);

    LHip_angVelocity_x_lowering_sparc =  LHip_angVelocity_x(idxOLW:idxFLW);
    LHip_angVelocity_y_lowering_sparc =  LHip_angVelocity_y(idxOLW:idxFLW);
    LHip_angVelocity_z_lowering_sparc =  LHip_angVelocity_z(idxOLW:idxFLW);
    
    % KNEE 
    RKnee_angVelocity_x_lowering_sparc =  RKnee_angVelocity_x(idxOLW:idxFLW);
    RKnee_angVelocity_y_lowering_sparc =  RKnee_angVelocity_y(idxOLW:idxFLW);
    RKnee_angVelocity_z_lowering_sparc =  RKnee_angVelocity_z(idxOLW:idxFLW);

    LKnee_angVelocity_x_lowering_sparc =  LKnee_angVelocity_x(idxOLW:idxFLW);
    LKnee_angVelocity_y_lowering_sparc =  LKnee_angVelocity_y(idxOLW:idxFLW);
    LKnee_angVelocity_z_lowering_sparc =  LKnee_angVelocity_z(idxOLW:idxFLW);
    
    % ANKLE 
    RAnkle_angVelocity_x_lowering_sparc =  RAnkle_angVelocity_x(idxOLW:idxFLW);
    RAnkle_angVelocity_y_lowering_sparc =  RAnkle_angVelocity_y(idxOLW:idxFLW);
    RAnkle_angVelocity_z_lowering_sparc =  RAnkle_angVelocity_z(idxOLW:idxFLW);

    LAnkle_angVelocity_x_lowering_sparc =  LAnkle_angVelocity_x(idxOLW:idxFLW);
    LAnkle_angVelocity_y_lowering_sparc =  LAnkle_angVelocity_y(idxOLW:idxFLW);
    LAnkle_angVelocity_z_lowering_sparc =  LAnkle_angVelocity_z(idxOLW:idxFLW);

    % ELBOW
    RElbow_angVelocity_x_lowering_sparc =  RElbow_angVelocity_x(idxOLW:idxFLW);
    RElbow_angVelocity_y_lowering_sparc =  RElbow_angVelocity_y(idxOLW:idxFLW);
    RElbow_angVelocity_z_lowering_sparc =  RElbow_angVelocity_z(idxOLW:idxFLW);

    LElbow_angVelocity_x_lowering_sparc =  LElbow_angVelocity_x(idxOLW:idxFLW);
    LElbow_angVelocity_y_lowering_sparc =  LElbow_angVelocity_y(idxOLW:idxFLW);
    LElbow_angVelocity_z_lowering_sparc =  LElbow_angVelocity_z(idxOLW:idxFLW);

figure (5) 
t4 = 0:1:(length(RHip_angVelocity_x_lowering_sparc)-1);
subplot (4,2,1)
plot(t4, RHip_angVelocity_x_lowering_sparc, t4, RHip_angVelocity_y_lowering_sparc, t4, RHip_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
title ('Cropped Right Hip Angular Velocities Lowering Phase') 
subplot (4,2,2)
plot(t4, LHip_angVelocity_x_lowering_sparc, t4, LHip_angVelocity_y_lowering_sparc, t4, LHip_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
title ('Cropped Left Hip Angular Velocities Lowering Phase')
subplot (4,2,3)
plot(t4, RKnee_angVelocity_x_lowering_sparc, t4, RKnee_angVelocity_y_lowering_sparc, t4, RKnee_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
title ('Cropped Right Knee Angular Velocities Lowering Phase')
subplot (4,2,4)
plot(t4, LKnee_angVelocity_x_lowering_sparc, t4, LKnee_angVelocity_y_lowering_sparc, t4, LKnee_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
title ('Cropped Left Knee Angular Velocities Lowering Phase')
subplot (4,2,5)
plot(t4, RAnkle_angVelocity_x_lowering_sparc, t4, RAnkle_angVelocity_y_lowering_sparc, t4, RAnkle_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
title ('Cropped Right Ankle Angular Velocities Lifting Phase')
subplot (4,2,6)
plot(t4, LAnkle_angVelocity_x_lowering_sparc, t4, LAnkle_angVelocity_y_lowering_sparc, t4, LAnkle_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
title ('Cropped Left Ankle Angular Velocities Lowering Phase')
subplot (4,2,7)
plot(t4, RElbow_angVelocity_x_lowering_sparc, t4, RElbow_angVelocity_y_lowering_sparc, t4, RElbow_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
title ('Cropped Right Elbow Angular Velocities Lowering Phase')
subplot (4,2,8)
plot(t4, LElbow_angVelocity_x_lowering_sparc, t4, LElbow_angVelocity_y_lowering_sparc, t4, LElbow_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
title ('Cropped Left Elbow Angular Velocities Lowering Phase')

% SPARC pi
    % HIP 
S_lowering_Hip = {'S_R_lowering_X', 'S_R_lowering_Y', 'S_R_lowering_Z', 'S_L_lowering_X', 'S_L_lowering_Y', 'S_L_lowering_Z'};
    Hip_angVel_low(1,:) =  RHip_angVelocity_x_lowering_sparc;
    Hip_angVel_low(2,:) =  RHip_angVelocity_y_lowering_sparc;
    Hip_angVel_low(3,:) =  RHip_angVelocity_z_lowering_sparc;
    Hip_angVel_low(4,:) =  LHip_angVelocity_x_lowering_sparc;
    Hip_angVel_low(5,:) =  LHip_angVelocity_y_lowering_sparc;
    Hip_angVel_low(6,:) =  LHip_angVelocity_z_lowering_sparc;

for i = 1:length(S_lowering_Hip)
    speed = Hip_angVel_low(i,:);
    S = SpectralArcLength(speed);
    [S_lowering_Hip{i}] = S;
end    

    % KNEE
S_lowering_Knee = {'S_R_lowering_X', 'S_R_lowering_Y', 'S_R_lowering_Z', 'S_L_lowering_X', 'S_L_lowering_Y', 'S_L_lowering_Z'};
    Knee_angVel_low(1,:) =  RKnee_angVelocity_x_lowering_sparc;
    Knee_angVel_low(2,:) =  RKnee_angVelocity_y_lowering_sparc;
    Knee_angVel_low(3,:) =  RKnee_angVelocity_z_lowering_sparc;
    Knee_angVel_low(4,:) =  LKnee_angVelocity_x_lowering_sparc;
    Knee_angVel_low(5,:) =  LKnee_angVelocity_y_lowering_sparc;
    Knee_angVel_low(6,:) =  LKnee_angVelocity_z_lowering_sparc;

for i = 1:length(S_lowering_Knee)
    speed = Knee_angVel_low(i,:);
    S = SpectralArcLength(speed);
    [S_lowering_Knee{i}] = S;
end  

% ANKLE
S_lowering_Ankle = {'S_R_lowering_X', 'S_R_lowering_Y', 'S_R_lowering_Z', 'S_L_lowering_X', 'S_L_lowering_Y', 'S_L_lowering_Z'};
    Ankle_angVel_low(1,:) =  RAnkle_angVelocity_x_lowering_sparc;
    Ankle_angVel_low(2,:) =  RAnkle_angVelocity_y_lowering_sparc;
    Ankle_angVel_low(3,:) =  RAnkle_angVelocity_z_lowering_sparc;
    Ankle_angVel_low(4,:) =  LAnkle_angVelocity_x_lowering_sparc;
    Ankle_angVel_low(5,:) =  LAnkle_angVelocity_y_lowering_sparc;
    Ankle_angVel_low(6,:) =  LAnkle_angVelocity_z_lowering_sparc;

for i = 1:length(S_lowering_Ankle)
    speed = Ankle_angVel_low(i,:);
    S = SpectralArcLength(speed);
    [S_lowering_Ankle{i}] = S;
end  

% ELBOW
S_lowering_Elbow = {'S_R_lowering_X', 'S_R_lowering_Y', 'S_R_lowering_Z', 'S_L_lowering_X', 'S_L_lowering_Y', 'S_L_lowering_Z'};
    Elbow_angVel_low(1,:) =  RElbow_angVelocity_x_lowering_sparc;
    Elbow_angVel_low(2,:) =  RElbow_angVelocity_y_lowering_sparc;
    Elbow_angVel_low(3,:) =  RElbow_angVelocity_z_lowering_sparc;
    Elbow_angVel_low(4,:) =  LElbow_angVelocity_x_lowering_sparc;
    Elbow_angVel_low(5,:) =  LElbow_angVelocity_y_lowering_sparc;
    Elbow_angVel_low(6,:) =  LElbow_angVelocity_z_lowering_sparc;

for i = 1:length(S_lowering_Elbow)
    speed = Elbow_angVel_low(i,:);
    S = SpectralArcLength(speed);
    [S_lowering_Elbow{i}] = S;
end  