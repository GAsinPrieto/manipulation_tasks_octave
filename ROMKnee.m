vicon = ViconNexus;
SubjectName = vicon.GetSubjectNames;
Fs = vicon.GetFrameRate;               % Sampling frequency

% Mechanical limitations H2 exoskeleton in degrees
lim_knee = 110;
lim_ankle = 60;

% Import SIDEBOX3 trajectory
    % Input = Subject name, marker name
    % Output = Trajectory data X, Y, Z & E
[SIDEBOX3(:,1),SIDEBOX3(:,2),SIDEBOX3(:,3),SIDEBOX3(:,4)] = vicon.GetTrajectory(SubjectName{1},'SIDEBOX3');

% Events of interest in vertical box trajectory identification
position_Z = SIDEBOX3(:,3)';
position_Z = position_Z(8:end);

[pks, locs] = findpeaks(position_Z,'minPeakProminence',20);
t = 0:1:(length(position_Z)-1);

TF1 = islocalmin(position_Z, 'FlatSelection', 'first');
idx = find(TF1);
flat = idx < locs(1);
idx_flat = find(flat);
flat2 = idx > locs(1);
idx_flat2 = find(flat2);

% Ultimo indice flat antes segundo maximo Z
flat3 = idx < locs(2);
idx_flat3 = find(flat3);
% Primer indice flat despues segundo maximo Z
flat4 = idx > locs(2);
idx_flat4 = find(flat4);

figure(1)
plot(t,position_Z,t(locs),position_Z(locs),'o');
hold on
plot(t,position_Z,'r*','MarkerIndices',idx(idx_flat(length(idx_flat))));
hold on
plot(t,position_Z,'r*','MarkerIndices',idx(idx_flat2(1)));
hold on
plot(t,position_Z,'r*','MarkerIndices',idx(idx_flat3(length(idx_flat3))));
hold on
plot(t,position_Z,'r*','MarkerIndices',idx(idx_flat4(1)));
hold off

% Signal segementation
idxO1 = 1;
idxF1 = idx(idx_flat(length(idx_flat))); % = idxOLF 
idxF2 = idx(idx_flat2(1)); % = idxFLF
idxO3 = idx(idx_flat3(length(idx_flat3))); % = idxOLW
idxF3 = idx(idx_flat4(1)); % = idxFLW
idxF4 = length(position_Z);

idx = [idxO1 idxF1 idxF2 idxO3 idxF3 idxF4];
        
figure(2) % Comprobacion
plot (t, position_Z)
hold on
plot(t,position_Z,'r*', 'MarkerIndices', idx(1));
plot(t,position_Z,'r*', 'MarkerIndices', idx(2));
plot(t,position_Z,'r*', 'MarkerIndices', idx(3));
plot(t,position_Z,'r*', 'MarkerIndices', idx(4));
plot(t,position_Z,'r*', 'MarkerIndices', idx(5));
plot(t,position_Z,'r*', 'MarkerIndices', idx(6));
hold off

% Import model data for Knee Relative Angles
ModelOutput = {'RKneeAngles','LKneeAngles'};

for i = 1:length(ModelOutput)
    [ModelData.Raw.(ModelOutput{i}), ModelData.Exists.(ModelOutput{i})] = vicon.GetModelOutput(SubjectName{1},ModelOutput{i});
end

% Knee 1,2. 2 dofs 
        RKneeAngles_x = ModelData.Raw.(ModelOutput{1})(1,:); 
        RKneeAngles_x = RKneeAngles_x(8:end);
        RKneeAngles_y = ModelData.Raw.(ModelOutput{1})(2,:);
        RKneeAngles_y = RKneeAngles_y(8:end);
        RKneeAngles_z = ModelData.Raw.(ModelOutput{1})(3,:);
        RKneeAngles_z = RKneeAngles_z(8:end);

        LKneeAngles_x = ModelData.Raw.(ModelOutput{2})(1,:);
        LKneeAngles_x = LKneeAngles_x(8:end);
        LKneeAngles_y = ModelData.Raw.(ModelOutput{2})(2,:);
        LKneeAngles_y = LKneeAngles_y(8:end);
        LKneeAngles_z = ModelData.Raw.(ModelOutput{2})(3,:);
        LKneeAngles_z = LKneeAngles_z(8:end);
        
% Visualizacion señal antes de segmentarlas
figure(3)
plot(t,RKneeAngles_x, t, RKneeAngles_y, '--', t, RKneeAngles_z,'.');
legend('x','y','z');
title('Right knee angles');

figure(4)
plot(t,LKneeAngles_x, t, LKneeAngles_y, '--', t, LKneeAngles_z,'.');   
legend('x','y','z');
title('Left knee angles');

% Right Knee      
phase1_RKnee_x = RKneeAngles_x(idxO1:idxF1);  % phase 1
phase1_RKnee_y = RKneeAngles_y(idxO1:idxF1);
phase1_RKnee_z = RKneeAngles_z(idxO1:idxF1);

phase2_RKnee_x = RKneeAngles_x(idxF1:idxF2); % phase 2
phase2_RKnee_y = RKneeAngles_y(idxF1:idxF2);
phase2_RKnee_z = RKneeAngles_z(idxF1:idxF2);

phase3_RKnee_x = RKneeAngles_x(idxO3:idxF3); % phase 3
phase3_RKnee_y = RKneeAngles_y(idxO3:idxF3);
phase3_RKnee_z = RKneeAngles_z(idxO3:idxF3);

phase4_RKnee_x = RKneeAngles_x(idxF3:idxF4); % phase 4
phase4_RKnee_y = RKneeAngles_y(idxF3:idxF4);
phase4_RKnee_z = RKneeAngles_z(idxF3:idxF4);

% Left Knee 
phase1_LKnee_x = LKneeAngles_x(idxO1:idxF1);
phase1_LKnee_y = LKneeAngles_y(idxO1:idxF1);
phase1_LKnee_z = LKneeAngles_z(idxO1:idxF1);

phase2_LKnee_x = LKneeAngles_x(idxF1:idxF2);
phase2_LKnee_y = LKneeAngles_y(idxF1:idxF2);
phase2_LKnee_z = LKneeAngles_z(idxF1:idxF2);

phase3_LKnee_x = LKneeAngles_x(idxO3:idxF3);
phase3_LKnee_y = LKneeAngles_y(idxO3:idxF3);
phase3_LKnee_z = LKneeAngles_z(idxO3:idxF3);

phase4_LKnee_x = LKneeAngles_x(idxF3:idxF4);
phase4_LKnee_y = LKneeAngles_y(idxF3:idxF4);
phase4_LKnee_z = LKneeAngles_z(idxF3:idxF4);

% time vectors definition
t2 = 0:1:length(phase3_RKnee_x)-1;
t3 = 0:1:(length(phase1_RKnee_x)-1);
t4 = 0:1:(length(phase2_RKnee_x)-1);
t5 = 0:1:(length(phase4_RKnee_x)-1);

% Right Knee Representacion Separando X,Y,Z y misma fase dist. peso juntos

% PARA COMPARAR SIN CAJA Y CON CAJA Y COMPARAR SEGUN LA DIRECCION X,Y,Z

% Right Knee Bajada
figure(5) 
subplot(3,1,1)
phase1_RKnee_x2 = phase1_RKnee_x((length(phase1_RKnee_x) - (length(phase3_RKnee_x)-1)): length(phase1_RKnee_x));
plot(t2, phase1_RKnee_x2, t2, phase3_RKnee_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Rodilla Derecha');
subplot(3,1,2)
phase1_RKnee_y2 = phase1_RKnee_y((length(phase1_RKnee_y) - (length(phase3_RKnee_y)-1)): length(phase1_RKnee_y));
plot(t2, phase1_RKnee_y2, t2, phase3_RKnee_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Rodilla Derecha');
subplot(3,1,3)
phase1_RKnee_z2 = phase1_RKnee_z((length(phase1_RKnee_z) - (length(phase3_RKnee_z)-1)): length(phase1_RKnee_z));
plot(t2, phase1_RKnee_z2, t2, phase3_RKnee_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Rodilla Derecha');

% Right Knee Subida

figure(6)
subplot(3,1,1)
if length(phase4_RKnee_x) > length(phase2_RKnee_x)
        phase4_RKnee_x2 = phase4_RKnee_x(1:length(phase2_RKnee_x));
        plot(t4, phase4_RKnee_x2, t4, phase2_RKnee_x, '--');
else if length(phase4_RKnee_x) < length(phase2_RKnee_x)
        phase2_RKnee_x2 = phase2_RKnee_x(1:length(phase4_RKnee_x));
        plot(t5, phase4_RKnee_x, t5, phase2_RKnee_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Rodilla Derecha');
subplot(3,1,2)
if length(phase4_RKnee_y) > length(phase2_RKnee_y)
        phase4_RKnee_y2 = phase4_RKnee_y(1:length(phase2_RKnee_y));
        plot(t4, phase4_RKnee_y2, t4, phase2_RKnee_y, '--');
else if length(phase4_RKnee_y) < length(phase2_RKnee_y)
        phase2_RKnee_y2 = phase2_RKnee_y(1:length(phase4_RKnee_y));
        plot(t5, phase4_RKnee_y, t5, phase2_RKnee_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Rodilla Derecha');
subplot(3,1,3)
if length(phase4_RKnee_z) > length(phase2_RKnee_z)
        phase4_RKnee_z2 = phase4_RKnee_z(1:length(phase2_RKnee_z));
        plot(t4, phase4_RKnee_z2, t4, phase2_RKnee_z, '--');
else if length(phase4_RKnee_z) < length(phase2_RKnee_z)
        phase2_RKnee_z2 = phase2_RKnee_z(1:length(phase4_RKnee_z));
        plot(t5, phase4_RKnee_z, t5, phase2_RKnee_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Rodilla Derecha');

% Left Knee Representacion Separando X,Y,Z y misma fase dist. peso juntos
% Left Knee Bajada
figure(7)
subplot(3,1,1)
phase1_LKnee_x2 = phase1_LKnee_x((length(phase1_LKnee_x) - (length(phase3_LKnee_x)-1)): length(phase1_LKnee_x));
plot(t2, phase1_LKnee_x2, t2, phase3_LKnee_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Rodilla Izquierda');
subplot(3,1,2)
phase1_LKnee_y2 = phase1_LKnee_y((length(phase1_LKnee_y) - (length(phase3_LKnee_y)-1)): length(phase1_LKnee_y));
plot(t2, phase1_LKnee_y2, t2, phase3_LKnee_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Rodilla Izquierda');
subplot(3,1,3)
phase1_LKnee_z2 = phase1_LKnee_z((length(phase1_LKnee_z) - (length(phase3_LKnee_z)-1)): length(phase1_LKnee_z));
plot(t2, phase1_LKnee_z2, t2, phase3_LKnee_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Rodilla Izquierda');

% Left Hip Subida
figure(8)
subplot(3,1,1)
if length(phase4_LKnee_x) > length(phase2_LKnee_x)
        phase4_LKnee_x2 = phase4_LKnee_x(1:length(phase2_LKnee_x));
        plot(t4, phase4_LKnee_x2, t4, phase2_LKnee_x, '--');
else if length(phase4_LKnee_x) < length(phase2_LKnee_x)
        phase2_LKnee_x2 = phase2_LKnee_x(1:length(phase4_LKnee_x));
        plot(t5, phase4_LKnee_x, t5, phase2_LKnee_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Rodilla Izquierda');
subplot(3,1,2)
if length(phase4_LKnee_y) > length(phase2_LKnee_y)
        phase4_LKnee_y2 = phase4_LKnee_y(1:length(phase2_LKnee_y));
        plot(t4, phase4_LKnee_y2, t4, phase2_LKnee_y, '--');
else if length(phase4_LKnee_y) < length(phase2_LKnee_y)
        phase2_LKnee_y2 = phase2_LKnee_y(1:length(phase4_LKnee_y));
        plot(t5, phase4_LKnee_y, t5, phase2_LKnee_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Rodilla Izquierda');
subplot(3,1,3)
if length(phase4_LKnee_z) > length(phase2_LKnee_z)
        phase4_LKnee_z2 = phase4_LKnee_z(1:length(phase2_LKnee_z));
        plot(t4, phase4_LKnee_z2, t4, phase2_LKnee_z, '--');
else if length(phase4_LKnee_z) < length(phase2_LKnee_z)
        phase2_LKnee_z2 = phase2_LKnee_z(1:length(phase4_LKnee_z));
        plot(t5, phase4_LKnee_z, t5, phase2_LKnee_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Rodilla Izquierda');

% Right Knee Representacion Separando 4 fases e X,Y,Z de cada fase juntos

% PARA COMPARAR DERECHA- IZQUIERDA EN LAS 4 FASES

figure(9)
subplot(4,1,1)
plot (t3, phase1_RKnee_x, t3, phase1_RKnee_y, '--',t3,phase1_RKnee_z, '.');
legend('x','y','z')
title('Bajada sin caja Rodilla Derecha');
subplot(4,1,2)
plot(t4, phase2_RKnee_x,t4,phase2_RKnee_y, '--',t4,phase2_RKnee_z, '.');
legend('x','y','z')
title('Subida con caja Rodilla Derecha');
subplot(4,1,3)
plot(t2, phase3_RKnee_x, t2, phase3_RKnee_y,'--', t2, phase3_RKnee_z, '.');
legend('x','y','z')
title('Bajada con caja Rodilla Derecha');
subplot(4,1,4)
plot(t5, phase4_RKnee_x, t5, phase4_RKnee_y,'--', t5, phase4_RKnee_z, '.');
legend('x','y','z')
title('Subida sin caja Rodilla Derecha');

% Left Knee Representacion Separando 4 fases e X,Y,Z de cada fase juntos

figure(10)
subplot(4,1,1)
plot (t3, phase1_LKnee_x, t3, phase1_LKnee_y, '--',t3,phase1_LKnee_z, '.');
legend('x','y','z')
title('Bajada sin caja Rodilla Izquierda');
subplot(4,1,2)
plot(t4, phase2_LKnee_x,t4, phase2_LKnee_y, '--',t4, phase2_LKnee_z, '.');
legend('x','y','z')
title('Subida con caja Rodilla Izquierda');
subplot(4,1,3)
plot(t2, phase3_LKnee_x, t2, phase3_LKnee_y,'--', t2, phase3_LKnee_z, '.');
legend('x','y','z')
title('Bajada con caja Rodilla Izquierda');
subplot(4,1,4)
plot(t5, phase4_LKnee_x, t5, phase4_LKnee_y,'--', t5, phase4_LKnee_z, '.');
legend('x','y','z')
title('Subida sin caja Rodilla Izquierda');

% ROM calculation

% Right Knee 
RKnee_ph1(:,1) = phase1_RKnee_x;
RKnee_ph1(:,2) = phase1_RKnee_y;
RKnee_ph1(:,3) = phase1_RKnee_z;

RKnee_ph2(:,1) = phase2_RKnee_x;
RKnee_ph2(:,2) = phase2_RKnee_y;
RKnee_ph2(:,3) = phase2_RKnee_z;

RKnee_ph3(:,1) = phase3_RKnee_x;
RKnee_ph3(:,2) = phase3_RKnee_y;
RKnee_ph3(:,3) = phase3_RKnee_z;

RKnee_ph4(:,1) = phase4_RKnee_x;
RKnee_ph4(:,2) = phase4_RKnee_y;
RKnee_ph4(:,3) = phase4_RKnee_z;

RKnee_ph1_max = max(RKnee_ph1);
RKnee_ph1_min = min(RKnee_ph1);

RKnee_ph2_max = max(RKnee_ph2);
RKnee_ph2_min = min(RKnee_ph2);

RKnee_ph3_max = max(RKnee_ph3);
RKnee_ph3_min = min(RKnee_ph3);

RKnee_ph4_max = max(RKnee_ph4);
RKnee_ph4_min = min(RKnee_ph4);

for i = 1:length(RKnee_ph1_max)
    RKnee_ph1_ROM(1,i) = RKnee_ph1_max(i)- RKnee_ph1_min(i);
end

for i = 1:length(RKnee_ph1_max)
    RKnee_ph2_ROM(1,i) = RKnee_ph2_max(i)- RKnee_ph2_min(i);
end

for i = 1:length(RKnee_ph1_max)
    RKnee_ph3_ROM(1,i) = RKnee_ph3_max(i)- RKnee_ph3_min(i);
end

for i = 1:length(RKnee_ph1_max)
    RKnee_ph4_ROM(1,i) = RKnee_ph4_max(i)- RKnee_ph4_min(i);
end

% Left Knee 
LKnee_ph1(:,1) = phase1_LKnee_x;
LKnee_ph1(:,2) = phase1_LKnee_y;
LKnee_ph1(:,3) = phase1_LKnee_z;

LKnee_ph2(:,1) = phase2_LKnee_x;
LKnee_ph2(:,2) = phase2_LKnee_y;
LKnee_ph2(:,3) = phase2_LKnee_z;

LKnee_ph3(:,1) = phase3_LKnee_x;
LKnee_ph3(:,2) = phase3_LKnee_y;
LKnee_ph3(:,3) = phase3_LKnee_z;

LKnee_ph4(:,1) = phase4_LKnee_x;
LKnee_ph4(:,2) = phase4_LKnee_y;
LKnee_ph4(:,3) = phase4_LKnee_z;

LKnee_ph1_max = max(LKnee_ph1);
LKnee_ph1_min = min(LKnee_ph1);

LKnee_ph2_max = max(LKnee_ph2);
LKnee_ph2_min = min(LKnee_ph2);

LKnee_ph3_max = max(LKnee_ph3);
LKnee_ph3_min = min(LKnee_ph3);

LKnee_ph4_max = max(LKnee_ph4);
LKnee_ph4_min = min(LKnee_ph4);

for i = 1:length(LKnee_ph1_max)
    LKnee_ph1_ROM(1,i) = LKnee_ph1_max(i)- LKnee_ph1_min(i);
end

for i = 1:length(LKnee_ph1_max)
    LKnee_ph2_ROM(1,i) = LKnee_ph2_max(i)- LKnee_ph2_min(i);
end

for i = 1:length(LKnee_ph1_max)
    LKnee_ph3_ROM(1,i) = LKnee_ph3_max(i)- LKnee_ph3_min(i);
end

for i = 1:length(LKnee_ph1_max)
    LKnee_ph4_ROM(1,i) = LKnee_ph4_max(i)- LKnee_ph4_min(i);
end

% Total ROM Knee Angles 
RKnee_total(:,1) = RKneeAngles_x;
RKnee_total(:,2) = RKneeAngles_y;
RKnee_total(:,3) = RKneeAngles_z;

RKnee_total_max = max(RKnee_total);
RKnee_total_min = min(RKnee_total);

for i = 1:length(RKnee_total_max)
    RKnee_total_max_ROM(1,i) = RKnee_total_max(i)- RKnee_total_min(i);
end

LKnee_total(:,1) = LKneeAngles_x;
LKnee_total(:,2) = LKneeAngles_y;
LKnee_total(:,3) = LKneeAngles_z;

LKnee_total_max = max(LKnee_total);
LKnee_total_min = min(LKnee_total);

for i = 1:length(LKnee_total_max)
    LKnee_total_max_ROM(1,i) = LKnee_total_max(i)- LKnee_total_min(i);
end

ROM_Right_Knee = table([RKnee_ph1_ROM(1);RKnee_ph2_ROM(1);RKnee_ph3_ROM(1);RKnee_ph4_ROM(1)],[RKnee_ph1_ROM(2);RKnee_ph2_ROM(2);RKnee_ph3_ROM(2);RKnee_ph4_ROM(2)],[RKnee_ph1_ROM(3);RKnee_ph2_ROM(3);RKnee_ph3_ROM(3);RKnee_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

ROM_Left_Knee = table([LKnee_ph1_ROM(1);LKnee_ph2_ROM(1);LKnee_ph3_ROM(1);LKnee_ph4_ROM(1)],[LKnee_ph1_ROM(2);LKnee_ph2_ROM(2);LKnee_ph3_ROM(2);LKnee_ph4_ROM(2)],[LKnee_ph1_ROM(3);LKnee_ph2_ROM(3);LKnee_ph3_ROM(3);LKnee_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

TotalROM_Knee = table([RKnee_total_max_ROM(1); LKnee_total_max_ROM(1)],[RKnee_total_max_ROM(2); LKnee_total_max_ROM(2)],[RKnee_total_max_ROM(3); LKnee_total_max_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})