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
ModelOutput = {'RShoulderAngles','LShoulderAngles'};

for i = 1:length(ModelOutput)
    [ModelData.Raw.(ModelOutput{i}), ModelData.Exists.(ModelOutput{i})] = vicon.GetModelOutput(SubjectName{1},ModelOutput{i});
end

% Shoulder  
        RShoulderAngles_x = ModelData.Raw.(ModelOutput{1})(1,:); 
        RShoulderAngles_x = RShoulderAngles_x(8:end);
        RShoulderAngles_y = ModelData.Raw.(ModelOutput{1})(2,:);
        RShoulderAngles_y = RShoulderAngles_y(8:end);
        RShoulderAngles_z = ModelData.Raw.(ModelOutput{1})(3,:);
        RShoulderAngles_z = RShoulderAngles_z(8:end);

        LShoulderAngles_x = ModelData.Raw.(ModelOutput{2})(1,:);
        LShoulderAngles_x = LShoulderAngles_x(8:end);
        LShoulderAngles_y = ModelData.Raw.(ModelOutput{2})(2,:);
        LShoulderAngles_y = LShoulderAngles_y(8:end);
        LShoulderAngles_z = ModelData.Raw.(ModelOutput{2})(3,:);
        LShoulderAngles_z = LShoulderAngles_z(8:end);
        
% Visualizacion señal antes de segmentarlas
figure(3)
plot(t,RShoulderAngles_x, t, RShoulderAngles_y, '--', t, RShoulderAngles_z,'.');
legend('x','y','z');
title('Right Shoulder angles');

figure(4)
plot(t,LShoulderAngles_x, t, LShoulderAngles_y, '--', t, LShoulderAngles_z,'.');   
legend('x','y','z');
title('Left Shoulder angles');

% Right Shoulder     
phase1_RShoulder_x = RShoulderAngles_x(idxO1:idxF1);  % phase 1
phase1_RShoulder_y = RShoulderAngles_y(idxO1:idxF1);
phase1_RShoulder_z = RShoulderAngles_z(idxO1:idxF1);

phase2_RShoulder_x = RShoulderAngles_x(idxF1:idxF2); % phase 2
phase2_RShoulder_y = RShoulderAngles_y(idxF1:idxF2);
phase2_RShoulder_z = RShoulderAngles_z(idxF1:idxF2);

phase3_RShoulder_x = RShoulderAngles_x(idxO3:idxF3); % phase 3
phase3_RShoulder_y = RShoulderAngles_y(idxO3:idxF3);
phase3_RShoulder_z = RShoulderAngles_z(idxO3:idxF3);

phase4_RShoulder_x = RShoulderAngles_x(idxF3:idxF4); % phase 4
phase4_RShoulder_y = RShoulderAngles_y(idxF3:idxF4);
phase4_RShoulder_z = RShoulderAngles_z(idxF3:idxF4);

% Left Shoulder 
phase1_LShoulder_x = LShoulderAngles_x(idxO1:idxF1);
phase1_LShoulder_y = LShoulderAngles_y(idxO1:idxF1);
phase1_LShoulder_z = LShoulderAngles_z(idxO1:idxF1);

phase2_LShoulder_x = LShoulderAngles_x(idxF1:idxF2);
phase2_LShoulder_y = LShoulderAngles_y(idxF1:idxF2);
phase2_LShoulder_z = LShoulderAngles_z(idxF1:idxF2);

phase3_LShoulder_x = LShoulderAngles_x(idxO3:idxF3);
phase3_LShoulder_y = LShoulderAngles_y(idxO3:idxF3);
phase3_LShoulder_z = LShoulderAngles_z(idxO3:idxF3);

phase4_LShoulder_x = LShoulderAngles_x(idxF3:idxF4);
phase4_LShoulder_y = LShoulderAngles_y(idxF3:idxF4);
phase4_LShoulder_z = LShoulderAngles_z(idxF3:idxF4);

% time vectors definition
t2 = 0:1:(length(phase3_RShoulder_x)-1);
t3 = 0:1:(length(phase1_RShoulder_x)-1);
t4 = 0:1:(length(phase2_RShoulder_x)-1);
t5 = 0:1:(length(phase4_RShoulder_x)-1);

% Right Elbow Representacion Separando X,Y,Z y misma fase dist. peso juntos

% PARA COMPARAR SIN CAJA Y CON CAJA Y COMPARAR SEGUN LA DIRECCION X,Y,Z

% Right Shoulder Bajada
figure(5) 
subplot(3,1,1)
phase1_RShoulder_x2 = phase1_RShoulder_x((length(phase1_RShoulder_x) - (length(phase3_RShoulder_x)-1)): length(phase1_RShoulder_x));
plot(t2, phase1_RShoulder_x2, t2, phase3_RShoulder_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Right Shoulder');
subplot(3,1,2)
phase1_RShoulder_y2 = phase1_RShoulder_y((length(phase1_RShoulder_y) - (length(phase3_RShoulder_y)-1)): length(phase1_RShoulder_y));
plot(t2, phase1_RShoulder_y2, t2, phase3_RShoulder_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Right Shoulder');
subplot(3,1,3)
phase1_RShoulder_z2 = phase1_RShoulder_z((length(phase1_RShoulder_z) - (length(phase3_RShoulder_z)-1)): length(phase1_RShoulder_z));
plot(t2, phase1_RShoulder_z2, t2, phase3_RShoulder_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Right Shoulder');

% Right Shoulder Subida
figure(6)
subplot(3,1,1)
if length(phase4_RShoulder_x) > length(phase2_RShoulder_x)
        phase4_RShoulder_x2 = phase4_RShoulder_x(1:length(phase2_RShoulder_x));
        plot(t4, phase4_RShoulder_x2, t4, phase2_RShoulder_x, '--');
else if length(phase4_RShoulder_x) < length(phase2_RShoulder_x)
        phase2_RShoulder_x2 = phase2_RShoulder_x(1:length(phase4_RShoulder_x));
        plot(t5, phase4_RShoulder_x, t5, phase2_RShoulder_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Right Shoulder');
subplot(3,1,2)
if length(phase4_RShoulder_y) > length(phase2_RShoulder_y)
        phase4_RShoulder_y2 = phase4_RShoulder_y(1:length(phase2_RShoulder_y));
        plot(t4, phase4_RShoulder_y2, t4, phase2_RShoulder_y, '--');
else if length(phase4_RShoulder_y) < length(phase2_RShoulder_y)
        phase2_RShoulder_y2 = phase2_RShoulder_y(1:length(phase4_RShoulder_y));
        plot(t5, phase4_RShoulder_y, t5, phase2_RShoulder_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Right Shoulder');
subplot(3,1,3)
if length(phase4_RShoulder_z) > length(phase2_RShoulder_z)
        phase4_RShoulder_z2 = phase4_RShoulder_z(1:length(phase2_RShoulder_z));
        plot(t4, phase4_RShoulder_z2, t4, phase2_RShoulder_z, '--');
else if length(phase4_RShoulder_z) < length(phase2_RShoulder_z)
        phase2_RShoulder_z2 = phase2_RShoulder_z(1:length(phase4_RShoulder_z));
        plot(t5, phase4_RShoulder_z, t5, phase2_RShoulder_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Right Shoulder');

% Left Shoulder Representacion Separando X,Y,Z y misma fase dist. peso juntos
% Left Shoulder Bajada
figure(7)
subplot(3,1,1)
phase1_LShoulder_x2 = phase1_LShoulder_x((length(phase1_LShoulder_x) - (length(phase3_LShoulder_x)-1)): length(phase1_LShoulder_x));
plot(t2, phase1_LShoulder_x2, t2, phase3_LShoulder_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Left Shoulder');
subplot(3,1,2)
phase1_LShoulder_y2 = phase1_LShoulder_y((length(phase1_LShoulder_y) - (length(phase3_LShoulder_y)-1)): length(phase1_LShoulder_y));
plot(t2, phase1_LShoulder_y2, t2, phase3_LShoulder_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Left Shoulder');
subplot(3,1,3)
phase1_LShoulder_z2 = phase1_LShoulder_z((length(phase1_LShoulder_z) - (length(phase3_LShoulder_z)-1)): length(phase1_LShoulder_z));
plot(t2, phase1_LShoulder_z2, t2, phase3_LShoulder_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Left Shoulder');

% Left Shoulder Subida
figure(8)
subplot(3,1,1)
if length(phase4_LShoulder_x) > length(phase2_LShoulder_x)
        phase4_LShoulder_x2 = phase4_LShoulder_x(1:length(phase2_LShoulder_x));
        plot(t4, phase4_LShoulder_x2, t4, phase2_LShoulder_x, '--');
else if length(phase4_LShoulder_x) < length(phase2_LShoulder_x)
        phase2_LShoulder_x2 = phase2_LShoulder_x(1:length(phase4_LShoulder_x));
        plot(t5, phase4_LShoulder_x, t5, phase2_LShoulder_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Left Shoulder');
subplot(3,1,2)
if length(phase4_LShoulder_y) > length(phase2_LShoulder_y)
        phase4_LShoulder_y2 = phase4_LShoulder_y(1:length(phase2_LShoulder_y));
        plot(t4, phase4_LShoulder_y2, t4, phase2_LShoulder_y, '--');
else if length(phase4_LShoulder_y) < length(phase2_LShoulder_y)
        phase2_LShoulder_y2 = phase2_LShoulder_y(1:length(phase4_LShoulder_y));
        plot(t5, phase4_LShoulder_y, t5, phase2_LShoulder_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Left Shoulder');
subplot(3,1,3)
if length(phase4_LShoulder_z) > length(phase2_LShoulder_z)
        phase4_LShoulder_z2 = phase4_LShoulder_z(1:length(phase2_LShoulder_z));
        plot(t4, phase4_LShoulder_z2, t4, phase2_LShoulder_z, '--');
else if length(phase4_LShoulder_z) < length(phase2_LShoulder_z)
        phase2_LShoulder_z2 = phase2_LShoulder_z(1:length(phase4_LShoulder_z));
        plot(t5, phase4_LShoulder_z, t5, phase2_LShoulder_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Left Shoulder');

% Right Shoulder Representacion Separando 4 fases e X,Y,Z de cada fase juntos

% PARA COMPARAR DERECHA- IZQUIERDA EN LAS 4 FASES

figure(9)
subplot(4,1,1)
plot (t3, phase1_RShoulder_x, t3, phase1_RShoulder_y, '--',t3,phase1_RShoulder_z, '.');
legend('x','y','z')
title('Bajada sin caja Right Shoulder');
subplot(4,1,2)
plot(t4, phase2_RShoulder_x,t4,phase2_RShoulder_y, '--',t4,phase2_RShoulder_z, '.');
legend('x','y','z')
title('Subida con caja Right Shoulder');
subplot(4,1,3)
plot(t2, phase3_RShoulder_x, t2, phase3_RShoulder_y,'--', t2, phase3_RShoulder_z, '.');
legend('x','y','z')
title('Bajada con caja Right Shoulder');
subplot(4,1,4)
plot(t5, phase4_RShoulder_x, t5, phase4_RShoulder_y,'--', t5, phase4_RShoulder_z, '.');
legend('x','y','z')
title('Subida sin caja Right Shoulder');

% Left Shoulder Representacion Separando 4 fases e X,Y,Z de cada fase juntos

figure(10)
subplot(4,1,1)
plot (t3, phase1_LShoulder_x, t3, phase1_LShoulder_y, '--',t3,phase1_LShoulder_z, '.');
legend('x','y','z')
title('Bajada sin caja Left Shoulder');
subplot(4,1,2)
plot(t4, phase2_LShoulder_x,t4, phase2_LShoulder_y, '--',t4, phase2_LShoulder_z, '.');
legend('x','y','z')
title('Subida con caja Left Shoulder');
subplot(4,1,3)
plot(t2, phase3_LShoulder_x, t2, phase3_LShoulder_y,'--', t2, phase3_LShoulder_z, '.');
legend('x','y','z')
title('Bajada con caja Left Shoulder');
subplot(4,1,4)
plot(t5, phase4_LShoulder_x, t5, phase4_LShoulder_y,'--', t5, phase4_LShoulder_z, '.');
legend('x','y','z')
title('Subida sin caja Left Shoulder');

% ROM calculation

% Right Shoulder 
RShoulder_ph1(:,1) = phase1_RShoulder_x;
RShoulder_ph1(:,2) = phase1_RShoulder_y;
RShoulder_ph1(:,3) = phase1_RShoulder_z;

RShoulder_ph2(:,1) = phase2_RShoulder_x;
RShoulder_ph2(:,2) = phase2_RShoulder_y;
RShoulder_ph2(:,3) = phase2_RShoulder_z;

RShoulder_ph3(:,1) = phase3_RShoulder_x;
RShoulder_ph3(:,2) = phase3_RShoulder_y;
RShoulder_ph3(:,3) = phase3_RShoulder_z;

RShoulder_ph4(:,1) = phase4_RShoulder_x;
RShoulder_ph4(:,2) = phase4_RShoulder_y;
RShoulder_ph4(:,3) = phase4_RShoulder_z;

RShoulder_ph1_max = max(RShoulder_ph1);
RShoulder_ph1_min = min(RShoulder_ph1);

RShoulder_ph2_max = max(RShoulder_ph2);
RShoulder_ph2_min = min(RShoulder_ph2);

RShoulder_ph3_max = max(RShoulder_ph3);
RShoulder_ph3_min = min(RShoulder_ph3);

RShoulder_ph4_max = max(RShoulder_ph4);
RShoulder_ph4_min = min(RShoulder_ph4);

for i = 1:length(RShoulder_ph1_max)
    RShoulder_ph1_ROM(1,i) = RShoulder_ph1_max(i)- RShoulder_ph1_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    RShoulder_ph2_ROM(1,i) = RShoulder_ph2_max(i)- RShoulder_ph2_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    RShoulder_ph3_ROM(1,i) = RShoulder_ph3_max(i)- RShoulder_ph3_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    RShoulder_ph4_ROM(1,i) = RShoulder_ph4_max(i)- RShoulder_ph4_min(i);
end

% Left Shoulder
LShoulder_ph1(:,1) = phase1_LShoulder_x;
LShoulder_ph1(:,2) = phase1_LShoulder_y;
LShoulder_ph1(:,3) = phase1_LShoulder_z;

LShoulder_ph2(:,1) = phase2_LShoulder_x;
LShoulder_ph2(:,2) = phase2_LShoulder_y;
LShoulder_ph2(:,3) = phase2_LShoulder_z;

LShoulder_ph3(:,1) = phase3_LShoulder_x;
LShoulder_ph3(:,2) = phase3_LShoulder_y;
LShoulder_ph3(:,3) = phase3_LShoulder_z;

LShoulder_ph4(:,1) = phase4_LShoulder_x;
LShoulder_ph4(:,2) = phase4_LShoulder_y;
LShoulder_ph4(:,3) = phase4_LShoulder_z;

LShoulder_ph1_max = max(LShoulder_ph1);
LShoulder_ph1_min = min(LShoulder_ph1);

LShoulder_ph2_max = max(LShoulder_ph2);
LShoulder_ph2_min = min(LShoulder_ph2);

LShoulder_ph3_max = max(LShoulder_ph3);
LShoulder_ph3_min = min(LShoulder_ph3);

LShoulder_ph4_max = max(LShoulder_ph4);
LShoulder_ph4_min = min(LShoulder_ph4);

for i = 1:length(RShoulder_ph1_max)
    LShoulder_ph1_ROM(1,i) = LShoulder_ph1_max(i)- LShoulder_ph1_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    LShoulder_ph2_ROM(1,i) = LShoulder_ph2_max(i)- LShoulder_ph2_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    LShoulder_ph3_ROM(1,i) = LShoulder_ph3_max(i)- LShoulder_ph3_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    LShoulder_ph4_ROM(1,i) = LShoulder_ph4_max(i)- LShoulder_ph4_min(i);
end

% Total ROM Shoulder Angles 
RShoulder_total(:,1) = RShoulderAngles_x;
RShoulder_total(:,2) = RShoulderAngles_y;
RShoulder_total(:,3) = RShoulderAngles_z;

RShoulder_total_max = max(RShoulder_total);
RShoulder_total_min = min(RShoulder_total);

for i = 1:length(RShoulder_total_max)
    RShoulder_total_max_ROM(1,i) = RShoulder_total_max(i)- RShoulder_total_min(i);
end

LShoulder_total(:,1) = LShoulderAngles_x;
LShoulder_total(:,2) = LShoulderAngles_y;
LShoulder_total(:,3) = LShoulderAngles_z;

LShoulder_total_max = max(LShoulder_total);
LShoulder_total_min = min(LShoulder_total);

for i = 1:length(LShoulder_total_max)
    LShoulder_total_max_ROM(1,i) = LShoulder_total_max(i)- LShoulder_total_min(i);
end

ROM_Right_Shoulder = table([RShoulder_ph1_ROM(1);RShoulder_ph2_ROM(1);RShoulder_ph3_ROM(1);RShoulder_ph4_ROM(1)],[RShoulder_ph1_ROM(2);RShoulder_ph2_ROM(2);RShoulder_ph3_ROM(2);RShoulder_ph4_ROM(2)],[RShoulder_ph1_ROM(3);RShoulder_ph2_ROM(3);RShoulder_ph3_ROM(3);RShoulder_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

ROM_Left_Shoulder = table([LShoulder_ph1_ROM(1);LShoulder_ph2_ROM(1);LShoulder_ph3_ROM(1);LShoulder_ph4_ROM(1)],[LShoulder_ph1_ROM(2);LShoulder_ph2_ROM(2);LShoulder_ph3_ROM(2);LShoulder_ph4_ROM(2)],[LShoulder_ph1_ROM(3);LShoulder_ph2_ROM(3);LShoulder_ph3_ROM(3);LShoulder_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

TotalROM_Shoulder = table([RShoulder_total_max_ROM(1); LShoulder_total_max_ROM(1)],[RShoulder_total_max_ROM(2); LShoulder_total_max_ROM(2)],[RShoulder_total_max_ROM(3); LShoulder_total_max_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})




