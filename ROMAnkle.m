vicon = ViconNexus;
SubjectName = vicon.GetSubjectNames;
Fs = vicon.GetFrameRate;               % Sampling frequency

% Mechanical limitations H2 exoskeleton in degrees
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
ModelOutput = {'RAnkleAngles','LAnkleAngles'};

for i = 1:length(ModelOutput)
    [ModelData.Raw.(ModelOutput{i}), ModelData.Exists.(ModelOutput{i})] = vicon.GetModelOutput(SubjectName{1},ModelOutput{i});
end

% Ankle 1 dof 
        RAnkleAngles_x = ModelData.Raw.(ModelOutput{1})(1,:); 
        RAnkleAngles_x = RAnkleAngles_x(8:end);
        RAnkleAngles_y = ModelData.Raw.(ModelOutput{1})(2,:);
        RAnkleAngles_y = RAnkleAngles_y(8:end);
        RAnkleAngles_z = ModelData.Raw.(ModelOutput{1})(3,:);
        RAnkleAngles_z = RAnkleAngles_z(8:end);

        LAnkleAngles_x = ModelData.Raw.(ModelOutput{2})(1,:);
        LAnkleAngles_x = LAnkleAngles_x(8:end);
        LAnkleAngles_y = ModelData.Raw.(ModelOutput{2})(2,:);
        LAnkleAngles_y = LAnkleAngles_y(8:end);
        LAnkleAngles_z = ModelData.Raw.(ModelOutput{2})(3,:);
        LAnkleAngles_z = LAnkleAngles_z(8:end);
        
% Visualizacion señal antes de segmentarlas
figure(3)
plot(t,RAnkleAngles_x, t, RAnkleAngles_y, '--', t, RAnkleAngles_z,'.');
legend('x','y','z');
title('Right Ankle angles');

figure(4)
plot(t,LAnkleAngles_x, t, LAnkleAngles_y, '--', t, LAnkleAngles_z,'.');   
legend('x','y','z');
title('Left Ankle angles');

% Right Ankle      
phase1_RAnkle_x = RAnkleAngles_x(idxO1:idxF1);  % phase 1
phase1_RAnkle_y = RAnkleAngles_y(idxO1:idxF1);
phase1_RAnkle_z = RAnkleAngles_z(idxO1:idxF1);

phase2_RAnkle_x = RAnkleAngles_x(idxF1:idxF2); % phase 2
phase2_RAnkle_y = RAnkleAngles_y(idxF1:idxF2);
phase2_RAnkle_z = RAnkleAngles_z(idxF1:idxF2);

phase3_RAnkle_x = RAnkleAngles_x(idxO3:idxF3); % phase 3
phase3_RAnkle_y = RAnkleAngles_y(idxO3:idxF3);
phase3_RAnkle_z = RAnkleAngles_z(idxO3:idxF3);

phase4_RAnkle_x = RAnkleAngles_x(idxF3:idxF4); % phase 4
phase4_RAnkle_y = RAnkleAngles_y(idxF3:idxF4);
phase4_RAnkle_z = RAnkleAngles_z(idxF3:idxF4);

% Left Ankle 
phase1_LAnkle_x = LAnkleAngles_x(idxO1:idxF1);
phase1_LAnkle_y = LAnkleAngles_y(idxO1:idxF1);
phase1_LAnkle_z = LAnkleAngles_z(idxO1:idxF1);

phase2_LAnkle_x = LAnkleAngles_x(idxF1:idxF2);
phase2_LAnkle_y = LAnkleAngles_y(idxF1:idxF2);
phase2_LAnkle_z = LAnkleAngles_z(idxF1:idxF2);

phase3_LAnkle_x = LAnkleAngles_x(idxO3:idxF3);
phase3_LAnkle_y = LAnkleAngles_y(idxO3:idxF3);
phase3_LAnkle_z = LAnkleAngles_z(idxO3:idxF3);

phase4_LAnkle_x = LAnkleAngles_x(idxF3:idxF4);
phase4_LAnkle_y = LAnkleAngles_y(idxF3:idxF4);
phase4_LAnkle_z = LAnkleAngles_z(idxF3:idxF4);

% time vectors definition
t2 = 0:1:length(phase3_RAnkle_x)-1;
t3 = 0:1:(length(phase1_RAnkle_x)-1);
t4 = 0:1:(length(phase2_RAnkle_x)-1);
t5 = 0:1:(length(phase4_RAnkle_x)-1);

% Right Ankle Representacion Separando X,Y,Z y misma fase dist. peso juntos

% PARA COMPARAR SIN CAJA Y CON CAJA Y COMPARAR SEGUN LA DIRECCION X,Y,Z

% Right Ankle Bajada
figure(5) 
subplot(3,1,1)
phase1_RAnkle_x2 = phase1_RAnkle_x((length(phase1_RAnkle_x) - (length(phase3_RAnkle_x)-1)): length(phase1_RAnkle_x));
plot(t2, phase1_RAnkle_x2, t2, phase3_RAnkle_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Tobillo Derecho');
subplot(3,1,2)
phase1_RAnkle_y2 = phase1_RAnkle_y((length(phase1_RAnkle_y) - (length(phase3_RAnkle_y)-1)): length(phase1_RAnkle_y));
plot(t2, phase1_RAnkle_y2, t2, phase3_RAnkle_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Tobillo Derecho');
subplot(3,1,3)
phase1_RAnkle_z2 = phase1_RAnkle_z((length(phase1_RAnkle_z) - (length(phase3_RAnkle_z)-1)): length(phase1_RAnkle_z));
plot(t2, phase1_RAnkle_z2, t2, phase3_RAnkle_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Tobillo Derecho');

% Right Ankle Subida
figure(6)
subplot(3,1,1)
if length(phase4_RAnkle_x) > length(phase2_RAnkle_x)
        phase4_RAnkle_x2 = phase4_RAnkle_x(1:length(phase2_RAnkle_x));
        plot(t4, phase4_RAnkle_x2, t4, phase2_RAnkle_x, '--');
else if length(phase4_RAnkle_x) < length(phase2_RAnkle_x)
        phase2_RAnkle_x2 = phase2_RAnkle_x(1:length(phase4_RAnkle_x));
        plot(t5, phase4_RAnkle_x, t5, phase2_RAnkle_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Tobillo Derecho');
subplot(3,1,2)
if length(phase4_RAnkle_y) > length(phase2_RAnkle_y)
        phase4_RAnkle_y2 = phase4_RAnkle_y(1:length(phase2_RAnkle_y));
        plot(t4, phase4_RAnkle_y2, t4, phase2_RAnkle_y, '--');
else if length(phase4_RAnkle_y) < length(phase2_RAnkle_y)
        phase2_RAnkle_y2 = phase2_RAnkle_y(1:length(phase4_RAnkle_y));
        plot(t5, phase4_RAnkle_y, t5, phase2_RAnkle_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Tobillo Derecho');
subplot(3,1,3)
if length(phase4_RAnkle_z) > length(phase2_RAnkle_z)
        phase4_RAnkle_z2 = phase4_RAnkle_z(1:length(phase2_RAnkle_z));
        plot(t4, phase4_RAnkle_z2, t4, phase2_RAnkle_z, '--');
else if length(phase4_RAnkle_z) < length(phase2_RAnkle_z)
        phase2_RAnkle_z2 = phase2_RAnkle_z(1:length(phase4_RAnkle_z));
        plot(t5, phase4_RAnkle_z, t5, phase2_RAnkle_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Tobillo Derecho');

% Left Ankle Representacion Separando X,Y,Z y misma fase dist. peso juntos
% Left Ankle Bajada
figure(7)
subplot(3,1,1)
phase1_LAnkle_x2 = phase1_LAnkle_x((length(phase1_LAnkle_x) - (length(phase3_LAnkle_x)-1)): length(phase1_LAnkle_x));
plot(t2, phase1_LAnkle_x2, t2, phase3_LAnkle_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Tobillo Izquierdo');
subplot(3,1,2)
phase1_LAnkle_y2 = phase1_LAnkle_y((length(phase1_LAnkle_y) - (length(phase3_LAnkle_y)-1)): length(phase1_LAnkle_y));
plot(t2, phase1_LAnkle_y2, t2, phase3_LAnkle_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Tobillo Izquierdo');
subplot(3,1,3)
phase1_LAnkle_z2 = phase1_LAnkle_z((length(phase1_LAnkle_z) - (length(phase3_LAnkle_z)-1)): length(phase1_LAnkle_z));
plot(t2, phase1_LAnkle_z2, t2, phase3_LAnkle_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Tobillo Izquierdo');

% Left Ankle Subida
figure(8)
subplot(3,1,1)
if length(phase4_LAnkle_x) > length(phase2_LAnkle_x)
        phase4_LAnkle_x2 = phase4_LAnkle_x(1:length(phase2_LAnkle_x));
        plot(t4, phase4_LAnkle_x2, t4, phase2_LAnkle_x, '--');
else if length(phase4_LAnkle_x) < length(phase2_LAnkle_x)
        phase2_LAnkle_x2 = phase2_LAnkle_x(1:length(phase4_LAnkle_x));
        plot(t5, phase4_LAnkle_x, t5, phase2_LAnkle_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Tobillo Izquierdo');
subplot(3,1,2)
if length(phase4_LAnkle_y) > length(phase2_LAnkle_y)
        phase4_LAnkle_y2 = phase4_LAnkle_y(1:length(phase2_LAnkle_y));
        plot(t4, phase4_LAnkle_y2, t4, phase2_LAnkle_y, '--');
else if length(phase4_LAnkle_y) < length(phase2_LAnkle_y)
        phase2_LAnkle_y2 = phase2_LAnkle_y(1:length(phase4_LAnkle_y));
        plot(t5, phase4_LAnkle_y, t5, phase2_LAnkle_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Tobillo Izquierdo');
subplot(3,1,3)
if length(phase4_LAnkle_z) > length(phase2_LAnkle_z)
        phase4_LAnkle_z2 = phase4_LAnkle_z(1:length(phase2_LAnkle_z));
        plot(t4, phase4_LAnkle_z2, t4, phase2_LAnkle_z, '--');
else if length(phase4_LAnkle_z) < length(phase2_LAnkle_z)
        phase2_LAnkle_z2 = phase2_LAnkle_z(1:length(phase4_LAnkle_z));
        plot(t5, phase4_LAnkle_z, t5, phase2_LAnkle_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Tobillo Izquierdo');

% Right Ankle Representacion Separando 4 fases e X,Y,Z de cada fase juntos

% PARA COMPARAR DERECHA- IZQUIERDA EN LAS 4 FASES

figure(9)
subplot(4,1,1)
plot (t3, phase1_RAnkle_x, t3, phase1_RAnkle_y, '--',t3,phase1_RAnkle_z, '.');
legend('x','y','z')
title('Bajada sin caja Tobillo Derecho');
subplot(4,1,2)
plot(t4, phase2_RAnkle_x,t4,phase2_RAnkle_y, '--',t4,phase2_RAnkle_z, '.');
legend('x','y','z')
title('Subida con caja Tobillo Derecho');
subplot(4,1,3)
plot(t2, phase3_RAnkle_x, t2, phase3_RAnkle_y,'--', t2, phase3_RAnkle_z, '.');
legend('x','y','z')
title('Bajada con caja Tobillo Derecho');
subplot(4,1,4)
plot(t5, phase4_RAnkle_x, t5, phase4_RAnkle_y,'--', t5, phase4_RAnkle_z, '.');
legend('x','y','z')
title('Subida sin caja Tobillo Derecho');

% Left Ankle Representacion Separando 4 fases e X,Y,Z de cada fase juntos

figure(10)
subplot(4,1,1)
plot (t3, phase1_LAnkle_x, t3, phase1_LAnkle_y, '--',t3,phase1_LAnkle_z, '.');
legend('x','y','z')
title('Bajada sin caja Tobillo Izquierdo');
subplot(4,1,2)
plot(t4, phase2_LAnkle_x,t4, phase2_LAnkle_y, '--',t4, phase2_LAnkle_z, '.');
legend('x','y','z')
title('Subida con caja Tobillo Izquierdo');
subplot(4,1,3)
plot(t2, phase3_LAnkle_x, t2, phase3_LAnkle_y,'--', t2, phase3_LAnkle_z, '.');
legend('x','y','z')
title('Bajada con caja Tobillo Izquierdo');
subplot(4,1,4)
plot(t5, phase4_LAnkle_x, t5, phase4_LAnkle_y,'--', t5, phase4_LAnkle_z, '.');
legend('x','y','z')
title('Subida sin caja Tobillo Izquierdo');

% ROM calculation

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

RAnkle_ph1_max = max(RAnkle_ph1);
RAnkle_ph1_min = min(RAnkle_ph1);

RAnkle_ph2_max = max(RAnkle_ph2);
RAnkle_ph2_min = min(RAnkle_ph2);

RAnkle_ph3_max = max(RAnkle_ph3);
RAnkle_ph3_min = min(RAnkle_ph3);

RAnkle_ph4_max = max(RAnkle_ph4);
RAnkle_ph4_min = min(RAnkle_ph4);

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

LAnkle_ph1_max = max(LAnkle_ph1);
LAnkle_ph1_min = min(LAnkle_ph1);

LAnkle_ph2_max = max(LAnkle_ph2);
LAnkle_ph2_min = min(LAnkle_ph2);

LAnkle_ph3_max = max(LAnkle_ph3);
LAnkle_ph3_min = min(LAnkle_ph3);

LAnkle_ph4_max = max(LAnkle_ph4);
LAnkle_ph4_min = min(LAnkle_ph4);

for i = 1:length(LAnkle_ph1_max)
    LAnkle_ph1_ROM(1,i) = LAnkle_ph1_max(i)- LAnkle_ph1_min(i);
end

for i = 1:length(LAnkle_ph1_max)
    LAnkle_ph2_ROM(1,i) = LAnkle_ph2_max(i)- LAnkle_ph2_min(i);
end

for i = 1:length(LAnkle_ph1_max)
    LAnkle_ph3_ROM(1,i) = LAnkle_ph3_max(i)- LAnkle_ph3_min(i);
end

for i = 1:length(LAnkle_ph1_max)
    LAnkle_ph4_ROM(1,i) = LAnkle_ph4_max(i)- LAnkle_ph4_min(i);
end

% Total ROM Ankle Angles 
RAnkle_total(:,1) = RAnkleAngles_x;
RAnkle_total(:,2) = RAnkleAngles_y;
RAnkle_total(:,3) = RAnkleAngles_z;

RAnkle_total_max = max(RAnkle_total);
RAnkle_total_min = min(RAnkle_total);

for i = 1:length(RAnkle_total_max)
    RAnkle_total_max_ROM(1,i) = RAnkle_total_max(i)- RAnkle_total_min(i);
end

LAnkle_total(:,1) = LAnkleAngles_x;
LAnkle_total(:,2) = LAnkleAngles_y;
LAnkle_total(:,3) = LAnkleAngles_z;

LAnkle_total_max = max(LAnkle_total);
LAnkle_total_min = min(LAnkle_total);

for i = 1:length(LAnkle_total_max)
    LAnkle_total_max_ROM(1,i) = LAnkle_total_max(i)- LAnkle_total_min(i);
end

ROM_Right_Ankle = table([RAnkle_ph1_ROM(1);RAnkle_ph2_ROM(1);RAnkle_ph3_ROM(1);RAnkle_ph4_ROM(1)],[RAnkle_ph1_ROM(2);RAnkle_ph2_ROM(2);RAnkle_ph3_ROM(2);RAnkle_ph4_ROM(2)],[RAnkle_ph1_ROM(3);RAnkle_ph2_ROM(3);RAnkle_ph3_ROM(3);RAnkle_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

ROM_Left_Ankle = table([LAnkle_ph1_ROM(1);LAnkle_ph2_ROM(1);LAnkle_ph3_ROM(1);LAnkle_ph4_ROM(1)],[LAnkle_ph1_ROM(2);LAnkle_ph2_ROM(2);LAnkle_ph3_ROM(2);LAnkle_ph4_ROM(2)],[LAnkle_ph1_ROM(3);LAnkle_ph2_ROM(3);LAnkle_ph3_ROM(3);LAnkle_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

TotalROM_Ankle = table([RAnkle_total_max_ROM(1); LAnkle_total_max_ROM(1)],[RAnkle_total_max_ROM(2); LAnkle_total_max_ROM(2)],[RAnkle_total_max_ROM(3); LAnkle_total_max_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})
