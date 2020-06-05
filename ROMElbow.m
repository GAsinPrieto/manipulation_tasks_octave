vicon = ViconNexus;
SubjectName = vicon.GetSubjectNames;
Fs = vicon.GetFrameRate;               % Sampling frequency

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
ModelOutput = {'RElbowAngles','LElbowAngles'};

for i = 1:length(ModelOutput)
    [ModelData.Raw.(ModelOutput{i}), ModelData.Exists.(ModelOutput{i})] = vicon.GetModelOutput(SubjectName{1},ModelOutput{i});
end

% Elbow  
        RElbowAngles_x = ModelData.Raw.(ModelOutput{1})(1,:); 
        RElbowAngles_x = RElbowAngles_x(8:end);
        RElbowAngles_y = ModelData.Raw.(ModelOutput{1})(2,:);
        RElbowAngles_y = RElbowAngles_y(8:end);
        RElbowAngles_z = ModelData.Raw.(ModelOutput{1})(3,:);
        RElbowAngles_z = RElbowAngles_z(8:end);

        LElbowAngles_x = ModelData.Raw.(ModelOutput{2})(1,:);
        LElbowAngles_x = LElbowAngles_x(8:end);
        LElbowAngles_y = ModelData.Raw.(ModelOutput{2})(2,:);
        LElbowAngles_y = LElbowAngles_y(8:end);
        LElbowAngles_z = ModelData.Raw.(ModelOutput{2})(3,:);
        LElbowAngles_z = LElbowAngles_z(8:end);
        
% Visualizacion señal antes de segmentarlas
figure(3)
plot(t,RElbowAngles_x, t, RElbowAngles_y, '--', t, RElbowAngles_z,'.');
legend('x','y','z');
title('Right Elbow angles');

figure(4)
plot(t,LElbowAngles_x, t, LElbowAngles_y, '--', t, LElbowAngles_z,'.');   
legend('x','y','z');
title('Left Elbow angles');

% Right Elbow      
phase1_RElbow_x = RElbowAngles_x(idxO1:idxF1);  % phase 1
phase1_RElbow_y = RElbowAngles_y(idxO1:idxF1);
phase1_RElbow_z = RElbowAngles_z(idxO1:idxF1);

phase2_RElbow_x = RElbowAngles_x(idxF1:idxF2); % phase 2
phase2_RElbow_y = RElbowAngles_y(idxF1:idxF2);
phase2_RElbow_z = RElbowAngles_z(idxF1:idxF2);

phase3_RElbow_x = RElbowAngles_x(idxO3:idxF3); % phase 3
phase3_RElbow_y = RElbowAngles_y(idxO3:idxF3);
phase3_RElbow_z = RElbowAngles_z(idxO3:idxF3);

phase4_RElbow_x = RElbowAngles_x(idxF3:idxF4); % phase 4
phase4_RElbow_y = RElbowAngles_y(idxF3:idxF4);
phase4_RElbow_z = RElbowAngles_z(idxF3:idxF4);

% Left Elbow 
phase1_LElbow_x = LElbowAngles_x(idxO1:idxF1);
phase1_LElbow_y = LElbowAngles_y(idxO1:idxF1);
phase1_LElbow_z = LElbowAngles_z(idxO1:idxF1);

phase2_LElbow_x = LElbowAngles_x(idxF1:idxF2);
phase2_LElbow_y = LElbowAngles_y(idxF1:idxF2);
phase2_LElbow_z = LElbowAngles_z(idxF1:idxF2);

phase3_LElbow_x = LElbowAngles_x(idxO3:idxF3);
phase3_LElbow_y = LElbowAngles_y(idxO3:idxF3);
phase3_LElbow_z = LElbowAngles_z(idxO3:idxF3);

phase4_LElbow_x = LElbowAngles_x(idxF3:idxF4);
phase4_LElbow_y = LElbowAngles_y(idxF3:idxF4);
phase4_LElbow_z = LElbowAngles_z(idxF3:idxF4);

% time vectors definition
t2 = 0:1:(length(phase3_RElbow_x)-1);
t3 = 0:1:(length(phase1_RElbow_x)-1);
t4 = 0:1:(length(phase2_RElbow_x)-1);
t5 = 0:1:(length(phase4_RElbow_x)-1);

% Right Elbow Representacion Separando X,Y,Z y misma fase dist. peso juntos

% PARA COMPARAR SIN CAJA Y CON CAJA Y COMPARAR SEGUN LA DIRECCION X,Y,Z

% Right Elbow Bajada
figure(5) 
subplot(3,1,1)
phase1_RElbow_x2 = phase1_RElbow_x((length(phase1_RElbow_x) - (length(phase3_RElbow_x)-1)): length(phase1_RElbow_x));
plot(t2, phase1_RElbow_x2, t2, phase3_RElbow_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Right Elbow');
subplot(3,1,2)
phase1_RElbow_y2 = phase1_RElbow_y((length(phase1_RElbow_y) - (length(phase3_RElbow_y)-1)): length(phase1_RElbow_y));
plot(t2, phase1_RElbow_y2, t2, phase3_RElbow_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Right Elbow');
subplot(3,1,3)
phase1_RElbow_z2 = phase1_RElbow_z((length(phase1_RElbow_z) - (length(phase3_RElbow_z)-1)): length(phase1_RElbow_z));
plot(t2, phase1_RElbow_z2, t2, phase3_RElbow_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Right Elbow');

% Right Elbow Subida
figure(6)
subplot(3,1,1)
if length(phase4_RElbow_x) > length(phase2_RElbow_x)
        phase4_RElbow_x2 = phase4_RElbow_x(1:length(phase2_RElbow_x));
        plot(t4, phase4_RElbow_x2, t4, phase2_RElbow_x, '--');
else if length(phase4_RElbow_x) < length(phase2_RElbow_x)
        phase2_RElbow_x2 = phase2_RElbow_x(1:length(phase4_RElbow_x));
        plot(t5, phase4_RElbow_x, t5, phase2_RElbow_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Right Elbow');
subplot(3,1,2)
if length(phase4_RElbow_y) > length(phase2_RElbow_y)
        phase4_RElbow_y2 = phase4_RElbow_y(1:length(phase2_RElbow_y));
        plot(t4, phase4_RElbow_y2, t4, phase2_RElbow_y, '--');
else if length(phase4_RElbow_y) < length(phase2_RElbow_y)
        phase2_RElbow_y2 = phase2_RElbow_y(1:length(phase4_RElbow_y));
        plot(t5, phase4_RElbow_y, t5, phase2_RElbow_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Right Elbow');
subplot(3,1,3)
if length(phase4_RElbow_z) > length(phase2_RElbow_z)
        phase4_RElbow_z2 = phase4_RElbow_z(1:length(phase2_RElbow_z));
        plot(t4, phase4_RElbow_z2, t4, phase2_RElbow_z, '--');
else if length(phase4_RElbow_z) < length(phase2_RElbow_z)
        phase2_RElbow_z2 = phase2_RElbow_z(1:length(phase4_RElbow_z));
        plot(t5, phase4_RElbow_z, t5, phase2_RElbow_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Right Elbow');

% Left Elbow Representacion Separando X,Y,Z y misma fase dist. peso juntos
% Left Elbow Bajada
figure(7)
subplot(3,1,1)
phase1_LElbow_x2 = phase1_LElbow_x((length(phase1_LElbow_x) - (length(phase3_LElbow_x)-1)): length(phase1_LElbow_x));
plot(t2, phase1_LElbow_x2, t2, phase3_LElbow_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Left Elbow');
subplot(3,1,2)
phase1_LElbow_y2 = phase1_LElbow_y((length(phase1_LElbow_y) - (length(phase3_LElbow_y)-1)): length(phase1_LElbow_y));
plot(t2, phase1_LElbow_y2, t2, phase3_LElbow_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Left Elbow');
subplot(3,1,3)
phase1_LElbow_z2 = phase1_LElbow_z((length(phase1_LElbow_z) - (length(phase3_LElbow_z)-1)): length(phase1_LElbow_z));
plot(t2, phase1_LElbow_z2, t2, phase3_LElbow_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Left Elbow');

% Left Elbow Subida
figure(8)
subplot(3,1,1)
if length(phase4_LElbow_x) > length(phase2_LElbow_x)
        phase4_LElbow_x2 = phase4_LElbow_x(1:length(phase2_LElbow_x));
        plot(t4, phase4_LElbow_x2, t4, phase2_LElbow_x, '--');
else if length(phase4_LElbow_x) < length(phase2_LElbow_x)
        phase2_LElbow_x2 = phase2_LElbow_x(1:length(phase4_LElbow_x));
        plot(t5, phase4_LElbow_x, t5, phase2_LElbow_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Left Elbow');
subplot(3,1,2)
if length(phase4_LElbow_y) > length(phase2_LElbow_y)
        phase4_LElbow_y2 = phase4_LElbow_y(1:length(phase2_LElbow_y));
        plot(t4, phase4_LElbow_y2, t4, phase2_LElbow_y, '--');
else if length(phase4_LElbow_y) < length(phase2_LElbow_y)
        phase2_LElbow_y2 = phase2_LElbow_y(1:length(phase4_LElbow_y));
        plot(t5, phase4_LElbow_y, t5, phase2_LElbow_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Left Elbow');
subplot(3,1,3)
if length(phase4_LElbow_z) > length(phase2_LElbow_z)
        phase4_LElbow_z2 = phase4_LElbow_z(1:length(phase2_LElbow_z));
        plot(t4, phase4_LElbow_z2, t4, phase2_LElbow_z, '--');
else if length(phase4_LElbow_z) < length(phase2_LElbow_z)
        phase2_LElbow_z2 = phase2_LElbow_z(1:length(phase4_LElbow_z));
        plot(t5, phase4_LElbow_z, t5, phase2_LElbow_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Left Elbow');

% Right Elbow Representacion Separando 4 fases e X,Y,Z de cada fase juntos

% PARA COMPARAR DERECHA- IZQUIERDA EN LAS 4 FASES

figure(9)
subplot(4,1,1)
plot (t3, phase1_RElbow_x, t3, phase1_RElbow_y, '--',t3,phase1_RElbow_z, '.');
legend('x','y','z')
title('Bajada sin caja Right Elbow');
subplot(4,1,2)
plot(t4, phase2_RElbow_x,t4,phase2_RElbow_y, '--',t4,phase2_RElbow_z, '.');
legend('x','y','z')
title('Subida con caja Right Elbow');
subplot(4,1,3)
plot(t2, phase3_RElbow_x, t2, phase3_RElbow_y,'--', t2, phase3_RElbow_z, '.');
legend('x','y','z')
title('Bajada con caja Right Elbow');
subplot(4,1,4)
plot(t5, phase4_RElbow_x, t5, phase4_RElbow_y,'--', t5, phase4_RElbow_z, '.');
legend('x','y','z')
title('Subida sin caja Right Elbow');

% Left Elbow Representacion Separando 4 fases e X,Y,Z de cada fase juntos

figure(10)
subplot(4,1,1)
plot (t3, phase1_LElbow_x, t3, phase1_LElbow_y, '--',t3,phase1_LElbow_z, '.');
legend('x','y','z')
title('Bajada sin caja Left Elbow');
subplot(4,1,2)
plot(t4, phase2_LElbow_x,t4, phase2_LElbow_y, '--',t4, phase2_LElbow_z, '.');
legend('x','y','z')
title('Subida con caja Left Elbow');
subplot(4,1,3)
plot(t2, phase3_LElbow_x, t2, phase3_LElbow_y,'--', t2, phase3_LElbow_z, '.');
legend('x','y','z')
title('Bajada con caja Left Elbow');
subplot(4,1,4)
plot(t5, phase4_LElbow_x, t5, phase4_LElbow_y,'--', t5, phase4_LElbow_z, '.');
legend('x','y','z')
title('Subida sin caja Left Elbow');

% ROM calculation

% Right Elbow 
RElbow_ph1(:,1) = phase1_RElbow_x;
RElbow_ph1(:,2) = phase1_RElbow_y;
RElbow_ph1(:,3) = phase1_RElbow_z;

RElbow_ph2(:,1) = phase2_RElbow_x;
RElbow_ph2(:,2) = phase2_RElbow_y;
RElbow_ph2(:,3) = phase2_RElbow_z;

RElbow_ph3(:,1) = phase3_RElbow_x;
RElbow_ph3(:,2) = phase3_RElbow_y;
RElbow_ph3(:,3) = phase3_RElbow_z;

RElbow_ph4(:,1) = phase4_RElbow_x;
RElbow_ph4(:,2) = phase4_RElbow_y;
RElbow_ph4(:,3) = phase4_RElbow_z;

RElbow_ph1_max = max(RElbow_ph1);
RElbow_ph1_min = min(RElbow_ph1);

RElbow_ph2_max = max(RElbow_ph2);
RElbow_ph2_min = min(RElbow_ph2);

RElbow_ph3_max = max(RElbow_ph3);
RElbow_ph3_min = min(RElbow_ph3);

RElbow_ph4_max = max(RElbow_ph4);
RElbow_ph4_min = min(RElbow_ph4);

for i = 1:length(RElbow_ph1_max)
    RElbow_ph1_ROM(1,i) = RElbow_ph1_max(i)- RElbow_ph1_min(i);
end

for i = 1:length(RElbow_ph1_max)
    RElbow_ph2_ROM(1,i) = RElbow_ph2_max(i)- RElbow_ph2_min(i);
end

for i = 1:length(RElbow_ph1_max)
    RElbow_ph3_ROM(1,i) = RElbow_ph3_max(i)- RElbow_ph3_min(i);
end

for i = 1:length(RElbow_ph1_max)
    RElbow_ph4_ROM(1,i) = RElbow_ph4_max(i)- RElbow_ph4_min(i);
end

% Left Elbow
LElbow_ph1(:,1) = phase1_LElbow_x;
LElbow_ph1(:,2) = phase1_LElbow_y;
LElbow_ph1(:,3) = phase1_LElbow_z;

LElbow_ph2(:,1) = phase2_LElbow_x;
LElbow_ph2(:,2) = phase2_LElbow_y;
LElbow_ph2(:,3) = phase2_LElbow_z;

LElbow_ph3(:,1) = phase3_LElbow_x;
LElbow_ph3(:,2) = phase3_LElbow_y;
LElbow_ph3(:,3) = phase3_LElbow_z;

LElbow_ph4(:,1) = phase4_LElbow_x;
LElbow_ph4(:,2) = phase4_LElbow_y;
LElbow_ph4(:,3) = phase4_LElbow_z;

LElbow_ph1_max = max(LElbow_ph1);
LElbow_ph1_min = min(LElbow_ph1);

LElbow_ph2_max = max(LElbow_ph2);
LElbow_ph2_min = min(LElbow_ph2);

LElbow_ph3_max = max(LElbow_ph3);
LElbow_ph3_min = min(LElbow_ph3);

LElbow_ph4_max = max(LElbow_ph4);
LElbow_ph4_min = min(LElbow_ph4);

for i = 1:length(RElbow_ph1_max)
    LElbow_ph1_ROM(1,i) = LElbow_ph1_max(i)- LElbow_ph1_min(i);
end

for i = 1:length(RElbow_ph1_max)
    LElbow_ph2_ROM(1,i) = LElbow_ph2_max(i)- LElbow_ph2_min(i);
end

for i = 1:length(RElbow_ph1_max)
    LElbow_ph3_ROM(1,i) = LElbow_ph3_max(i)- LElbow_ph3_min(i);
end

for i = 1:length(RElbow_ph1_max)
    LElbow_ph4_ROM(1,i) = LElbow_ph4_max(i)- LElbow_ph4_min(i);
end

% Total ROM Elbow Angles 
RElbow_total(:,1) = RElbowAngles_x;
RElbow_total(:,2) = RElbowAngles_y;
RElbow_total(:,3) = RElbowAngles_z;

RElbow_total_max = max(RElbow_total);
RElbow_total_min = min(RElbow_total);

for i = 1:length(RElbow_total_max)
    RElbow_total_max_ROM(1,i) = RElbow_total_max(i)- RElbow_total_min(i);
end

LElbow_total(:,1) = LElbowAngles_x;
LElbow_total(:,2) = LElbowAngles_y;
LElbow_total(:,3) = LElbowAngles_z;

LElbow_total_max = max(LElbow_total);
LElbow_total_min = min(LElbow_total);

for i = 1:length(LElbow_total_max)
    LElbow_total_max_ROM(1,i) = LElbow_total_max(i)- LElbow_total_min(i);
end

ROM_Right_Elbow = table([RElbow_ph1_ROM(1);RElbow_ph2_ROM(1);RElbow_ph3_ROM(1);RElbow_ph4_ROM(1)],[RElbow_ph1_ROM(2);RElbow_ph2_ROM(2);RElbow_ph3_ROM(2);RElbow_ph4_ROM(2)],[RElbow_ph1_ROM(3);RElbow_ph2_ROM(3);RElbow_ph3_ROM(3);RElbow_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

ROM_Left_Elbow = table([LElbow_ph1_ROM(1);LElbow_ph2_ROM(1);LElbow_ph3_ROM(1);LElbow_ph4_ROM(1)],[LElbow_ph1_ROM(2);LElbow_ph2_ROM(2);LElbow_ph3_ROM(2);LElbow_ph4_ROM(2)],[LElbow_ph1_ROM(3);LElbow_ph2_ROM(3);LElbow_ph3_ROM(3);LElbow_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

TotalROM_Elbow = table([RElbow_total_max_ROM(1); LElbow_total_max_ROM(1)],[RElbow_total_max_ROM(2); LElbow_total_max_ROM(2)],[RElbow_total_max_ROM(3); LElbow_total_max_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})

