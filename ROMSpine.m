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
ModelOutput = {'RSpineAngles','LSpineAngles'};

for i = 1:length(ModelOutput)
    [ModelData.Raw.(ModelOutput{i}), ModelData.Exists.(ModelOutput{i})] = vicon.GetModelOutput(SubjectName{1},ModelOutput{i});
end

% Spine  
        RSpineAngles_x = ModelData.Raw.(ModelOutput{1})(1,:); 
        RSpineAngles_x = RSpineAngles_x(8:end);
        RSpineAngles_y = ModelData.Raw.(ModelOutput{1})(2,:);
        RSpineAngles_y = RSpineAngles_y(8:end);
        RSpineAngles_z = ModelData.Raw.(ModelOutput{1})(3,:);
        RSpineAngles_z = RSpineAngles_z(8:end);

        LSpineAngles_x = ModelData.Raw.(ModelOutput{2})(1,:);
        LSpineAngles_x = LSpineAngles_x(8:end);
        LSpineAngles_y = ModelData.Raw.(ModelOutput{2})(2,:);
        LSpineAngles_y = LSpineAngles_y(8:end);
        LSpineAngles_z = ModelData.Raw.(ModelOutput{2})(3,:);
        LSpineAngles_z = LSpineAngles_z(8:end);
        
% Visualizacion señal antes de segmentarlas
figure(3)
plot(t,RSpineAngles_x, t, RSpineAngles_y, '--', t, RSpineAngles_z,'.');
legend('x','y','z');
title('Right Spine angles');

figure(4)
plot(t,LSpineAngles_x, t, LSpineAngles_y, '--', t, LSpineAngles_z,'.');   
legend('x','y','z');
title('Left Spine angles');

% Right Spine      
phase1_RSpine_x = RSpineAngles_x(idxO1:idxF1);  % phase 1
phase1_RSpine_y = RSpineAngles_y(idxO1:idxF1);
phase1_RSpine_z = RSpineAngles_z(idxO1:idxF1);

phase2_RSpine_x = RSpineAngles_x(idxF1:idxF2); % phase 2
phase2_RSpine_y = RSpineAngles_y(idxF1:idxF2);
phase2_RSpine_z = RSpineAngles_z(idxF1:idxF2);

phase3_RSpine_x = RSpineAngles_x(idxO3:idxF3); % phase 3
phase3_RSpine_y = RSpineAngles_y(idxO3:idxF3);
phase3_RSpine_z = RSpineAngles_z(idxO3:idxF3);

phase4_RSpine_x = RSpineAngles_x(idxF3:idxF4); % phase 4
phase4_RSpine_y = RSpineAngles_y(idxF3:idxF4);
phase4_RSpine_z = RSpineAngles_z(idxF3:idxF4);

% Left Spine 
phase1_LSpine_x = LSpineAngles_x(idxO1:idxF1);
phase1_LSpine_y = LSpineAngles_y(idxO1:idxF1);
phase1_LSpine_z = LSpineAngles_z(idxO1:idxF1);

phase2_LSpine_x = LSpineAngles_x(idxF1:idxF2);
phase2_LSpine_y = LSpineAngles_y(idxF1:idxF2);
phase2_LSpine_z = LSpineAngles_z(idxF1:idxF2);

phase3_LSpine_x = LSpineAngles_x(idxO3:idxF3);
phase3_LSpine_y = LSpineAngles_y(idxO3:idxF3);
phase3_LSpine_z = LSpineAngles_z(idxO3:idxF3);

phase4_LSpine_x = LSpineAngles_x(idxF3:idxF4);
phase4_LSpine_y = LSpineAngles_y(idxF3:idxF4);
phase4_LSpine_z = LSpineAngles_z(idxF3:idxF4);

% time vectors definition
t2 = 0:1:length(phase3_RSpine_x)-1;
t3 = 0:1:(length(phase1_RSpine_x)-1);
t4 = 0:1:(length(phase2_RSpine_x)-1);
t5 = 0:1:(length(phase4_RSpine_x)-1);

% Right Spine Representacion Separando X,Y,Z y misma fase dist. peso juntos

% PARA COMPARAR SIN CAJA Y CON CAJA Y COMPARAR SEGUN LA DIRECCION X,Y,Z

% Right Spine Bajada
figure(5) 
subplot(3,1,1)
phase1_RSpine_x2 = phase1_RSpine_x((length(phase1_RSpine_x) - (length(phase3_RSpine_x)-1)): length(phase1_RSpine_x));
plot(t2, phase1_RSpine_x2, t2, phase3_RSpine_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Right Spine');
subplot(3,1,2)
phase1_RSpine_y2 = phase1_RSpine_y((length(phase1_RSpine_y) - (length(phase3_RSpine_y)-1)): length(phase1_RSpine_y));
plot(t2, phase1_RSpine_y2, t2, phase3_RSpine_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Right Spine');
subplot(3,1,3)
phase1_RSpine_z2 = phase1_RSpine_z((length(phase1_RSpine_z) - (length(phase3_RSpine_z)-1)): length(phase1_RSpine_z));
plot(t2, phase1_RSpine_z2, t2, phase3_RSpine_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Right Spine');

% Right Spine Subida
figure(6)
subplot(3,1,1)
if length(phase4_RSpine_x) > length(phase2_RSpine_x)
        phase4_RSpine_x2 = phase4_RSpine_x(1:length(phase2_RSpine_x));
        plot(t4, phase4_RSpine_x2, t4, phase2_RSpine_x, '--');
else if length(phase4_RSpine_x) < length(phase2_RSpine_x)
        phase2_RSpine_x2 = phase2_RSpine_x(1:length(phase4_RSpine_x));
        plot(t5, phase4_RSpine_x, t5, phase2_RSpine_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Right Spine');
subplot(3,1,2)
if length(phase4_RSpine_y) > length(phase2_RSpine_y)
        phase4_RSpine_y2 = phase4_RSpine_y(1:length(phase2_RSpine_y));
        plot(t4, phase4_RSpine_y2, t4, phase2_RSpine_y, '--');
else if length(phase4_RSpine_y) < length(phase2_RSpine_y)
        phase2_RSpine_y2 = phase2_RSpine_y(1:length(phase4_RSpine_y));
        plot(t5, phase4_RSpine_y, t5, phase2_RSpine_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Right Spine');
subplot(3,1,3)
if length(phase4_RSpine_z) > length(phase2_RSpine_z)
        phase4_RSpine_z2 = phase4_RSpine_z(1:length(phase2_RSpine_z));
        plot(t4, phase4_RSpine_z2, t4, phase2_RSpine_z, '--');
else if length(phase4_RSpine_z) < length(phase2_RSpine_z)
        phase2_RSpine_z2 = phase2_RSpine_z(1:length(phase4_RSpine_z));
        plot(t5, phase4_RSpine_z, t5, phase2_RSpine_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Right Spine');

% Left Spine Representacion Separando X,Y,Z y misma fase dist. peso juntos
% Left Spine Bajada
figure(7)
subplot(3,1,1)
phase1_LSpine_x2 = phase1_LSpine_x((length(phase1_LSpine_x) - (length(phase3_LSpine_x)-1)): length(phase1_LSpine_x));
plot(t2, phase1_LSpine_x2, t2, phase3_LSpine_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Left Spine');
subplot(3,1,2)
phase1_LSpine_y2 = phase1_LSpine_y((length(phase1_LSpine_y) - (length(phase3_LSpine_y)-1)): length(phase1_LSpine_y));
plot(t2, phase1_LSpine_y2, t2, phase3_LSpine_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Left Spine');
subplot(3,1,3)
phase1_LSpine_z2 = phase1_LSpine_z((length(phase1_LSpine_z) - (length(phase3_LSpine_z)-1)): length(phase1_LSpine_z));
plot(t2, phase1_LSpine_z2, t2, phase3_LSpine_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Left Spine');

% Left Spine Subida
figure(8)
subplot(3,1,1)
if length(phase4_LSpine_x) > length(phase2_LSpine_x)
        phase4_LSpine_x2 = phase4_LSpine_x(1:length(phase2_LSpine_x));
        plot(t4, phase4_LSpine_x2, t4, phase2_LSpine_x, '--');
else if length(phase4_LSpine_x) < length(phase2_LSpine_x)
        phase2_LSpine_x2 = phase2_LSpine_x(1:length(phase4_LSpine_x));
        plot(t5, phase4_LSpine_x, t5, phase2_LSpine_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Left Spine');
subplot(3,1,2)
if length(phase4_LSpine_y) > length(phase2_LSpine_y)
        phase4_LSpine_y2 = phase4_LSpine_y(1:length(phase2_LSpine_y));
        plot(t4, phase4_LSpine_y2, t4, phase2_LSpine_y, '--');
else if length(phase4_LSpine_y) < length(phase2_LSpine_y)
        phase2_LSpine_y2 = phase2_LSpine_y(1:length(phase4_LSpine_y));
        plot(t5, phase4_LSpine_y, t5, phase2_LSpine_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Left Spine');
subplot(3,1,3)
if length(phase4_LSpine_z) > length(phase2_LSpine_z)
        phase4_LSpine_z2 = phase4_LSpine_z(1:length(phase2_LSpine_z));
        plot(t4, phase4_LSpine_z2, t4, phase2_LSpine_z, '--');
else if length(phase4_LSpine_z) < length(phase2_LSpine_z)
        phase2_LSpine_z2 = phase2_LSpine_z(1:length(phase4_LSpine_z));
        plot(t5, phase4_LSpine_z, t5, phase2_LSpine_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Left Spine');

% Right Spine Representacion Separando 4 fases e X,Y,Z de cada fase juntos

% PARA COMPARAR DERECHA- IZQUIERDA EN LAS 4 FASES

figure(9)
subplot(4,1,1)
plot (t3, phase1_RSpine_x, t3, phase1_RSpine_y, '--',t3,phase1_RSpine_z, '.');
legend('x','y','z')
title('Bajada sin caja Right Spine');
subplot(4,1,2)
plot(t4, phase2_RSpine_x,t4,phase2_RSpine_y, '--',t4,phase2_RSpine_z, '.');
legend('x','y','z')
title('Subida con caja Right Spine');
subplot(4,1,3)
plot(t2, phase3_RSpine_x, t2, phase3_RSpine_y,'--', t2, phase3_RSpine_z, '.');
legend('x','y','z')
title('Bajada con caja Right Spine');
subplot(4,1,4)
plot(t5, phase4_RSpine_x, t5, phase4_RSpine_y,'--', t5, phase4_RSpine_z, '.');
legend('x','y','z')
title('Subida sin caja Right Spine');

% Left Spine Representacion Separando 4 fases e X,Y,Z de cada fase juntos

figure(10)
subplot(4,1,1)
plot (t3, phase1_LSpine_x, t3, phase1_LSpine_y, '--',t3,phase1_LSpine_z, '.');
legend('x','y','z')
title('Bajada sin caja Left Spine');
subplot(4,1,2)
plot(t4, phase2_LSpine_x,t4, phase2_LSpine_y, '--',t4, phase2_LSpine_z, '.');
legend('x','y','z')
title('Subida con caja Left Spine');
subplot(4,1,3)
plot(t2, phase3_LSpine_x, t2, phase3_LSpine_y,'--', t2, phase3_LSpine_z, '.');
legend('x','y','z')
title('Bajada con caja Left Spine');
subplot(4,1,4)
plot(t5, phase4_LSpine_x, t5, phase4_LSpine_y,'--', t5, phase4_LSpine_z, '.');
legend('x','y','z')
title('Subida sin caja Left Spine');

% ROM calculation

% Right Spine 
RSpine_ph1(:,1) = phase1_RSpine_x;
RSpine_ph1(:,2) = phase1_RSpine_y;
RSpine_ph1(:,3) = phase1_RSpine_z;

RSpine_ph2(:,1) = phase2_RSpine_x;
RSpine_ph2(:,2) = phase2_RSpine_y;
RSpine_ph2(:,3) = phase2_RSpine_z;

RSpine_ph3(:,1) = phase3_RSpine_x;
RSpine_ph3(:,2) = phase3_RSpine_y;
RSpine_ph3(:,3) = phase3_RSpine_z;

RSpine_ph4(:,1) = phase4_RSpine_x;
RSpine_ph4(:,2) = phase4_RSpine_y;
RSpine_ph4(:,3) = phase4_RSpine_z;

RSpine_ph1_max = max(RSpine_ph1);
RSpine_ph1_min = min(RSpine_ph1);

RSpine_ph2_max = max(RSpine_ph2);
RSpine_ph2_min = min(RSpine_ph2);

RSpine_ph3_max = max(RSpine_ph3);
RSpine_ph3_min = min(RSpine_ph3);

RSpine_ph4_max = max(RSpine_ph4);
RSpine_ph4_min = min(RSpine_ph4);

for i = 1:length(RSpine_ph1_max)
    RSpine_ph1_ROM(1,i) = RSpine_ph1_max(i)- RSpine_ph1_min(i);
end

for i = 1:length(RSpine_ph1_max)
    RSpine_ph2_ROM(1,i) = RSpine_ph2_max(i)- RSpine_ph2_min(i);
end

for i = 1:length(RSpine_ph1_max)
    RSpine_ph3_ROM(1,i) = RSpine_ph3_max(i)- RSpine_ph3_min(i);
end

for i = 1:length(RSpine_ph1_max)
    RSpine_ph4_ROM(1,i) = RSpine_ph4_max(i)- RSpine_ph4_min(i);
end

% Left Spine
LSpine_ph1(:,1) = phase1_LSpine_x;
LSpine_ph1(:,2) = phase1_LSpine_y;
LSpine_ph1(:,3) = phase1_LSpine_z;

LSpine_ph2(:,1) = phase2_LSpine_x;
LSpine_ph2(:,2) = phase2_LSpine_y;
LSpine_ph2(:,3) = phase2_LSpine_z;

LSpine_ph3(:,1) = phase3_LSpine_x;
LSpine_ph3(:,2) = phase3_LSpine_y;
LSpine_ph3(:,3) = phase3_LSpine_z;

LSpine_ph4(:,1) = phase4_LSpine_x;
LSpine_ph4(:,2) = phase4_LSpine_y;
LSpine_ph4(:,3) = phase4_LSpine_z;

LSpine_ph1_max = max(LSpine_ph1);
LSpine_ph1_min = min(LSpine_ph1);

LSpine_ph2_max = max(LSpine_ph2);
LSpine_ph2_min = min(LSpine_ph2);

LSpine_ph3_max = max(LSpine_ph3);
LSpine_ph3_min = min(LSpine_ph3);

LSpine_ph4_max = max(LSpine_ph4);
LSpine_ph4_min = min(LSpine_ph4);

for i = 1:length(RSpine_ph1_max)
    LSpine_ph1_ROM(1,i) = LSpine_ph1_max(i)- LSpine_ph1_min(i);
end

for i = 1:length(RSpine_ph1_max)
    LSpine_ph2_ROM(1,i) = LSpine_ph2_max(i)- LSpine_ph2_min(i);
end

for i = 1:length(RSpine_ph1_max)
    LSpine_ph3_ROM(1,i) = LSpine_ph3_max(i)- LSpine_ph3_min(i);
end

for i = 1:length(RSpine_ph1_max)
    LSpine_ph4_ROM(1,i) = LSpine_ph4_max(i)- LSpine_ph4_min(i);
end

% Total ROM Spine Angles 
RSpine_total(:,1) = RSpineAngles_x;
RSpine_total(:,2) = RSpineAngles_y;
RSpine_total(:,3) = RSpineAngles_z;

RSpine_total_max = max(RSpine_total);
RSpine_total_min = min(RSpine_total);

for i = 1:length(RSpine_total_max)
    RSpine_total_max_ROM(1,i) = RSpine_total_max(i)- RSpine_total_min(i);
end

LSpine_total(:,1) = LSpineAngles_x;
LSpine_total(:,2) = LSpineAngles_y;
LSpine_total(:,3) = LSpineAngles_z;

LSpine_total_max = max(LSpine_total);
LSpine_total_min = min(LSpine_total);

for i = 1:length(LSpine_total_max)
    LSpine_total_max_ROM(1,i) = LSpine_total_max(i)- LSpine_total_min(i);
end

ROM_Right_Spine = table([RSpine_ph1_ROM(1);RSpine_ph2_ROM(1);RSpine_ph3_ROM(1);RSpine_ph4_ROM(1)],[RSpine_ph1_ROM(2);RSpine_ph2_ROM(2);RSpine_ph3_ROM(2);RSpine_ph4_ROM(2)],[RSpine_ph1_ROM(3);RSpine_ph2_ROM(3);RSpine_ph3_ROM(3);RSpine_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

ROM_Left_Spine = table([LSpine_ph1_ROM(1);LSpine_ph2_ROM(1);LSpine_ph3_ROM(1);LSpine_ph4_ROM(1)],[LSpine_ph1_ROM(2);LSpine_ph2_ROM(2);LSpine_ph3_ROM(2);LSpine_ph4_ROM(2)],[LSpine_ph1_ROM(3);LSpine_ph2_ROM(3);LSpine_ph3_ROM(3);LSpine_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

TotalROM_Spine = table([RSpine_total_max_ROM(1); LSpine_total_max_ROM(1)],[RSpine_total_max_ROM(2); LSpine_total_max_ROM(2)],[RSpine_total_max_ROM(3); LSpine_total_max_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})
