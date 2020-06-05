vicon = ViconNexus;
SubjectName = vicon.GetSubjectNames;
Fs = vicon.GetFrameRate;               % Sampling frequency

% Mechanical limitations H2 exoskeleton in degrees
lim_hip = 140;
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

% Import model data for right wrist, interest in X and Y due to wrist's dofs
ModelOutput = {'RHipAngles','LHipAngles'};

for i = 1:length(ModelOutput)
    [ModelData.Raw.(ModelOutput{i}), ModelData.Exists.(ModelOutput{i})] = vicon.GetModelOutput(SubjectName{1},ModelOutput{i});
end

% Hip 1,2. 3 dofs PENSAR HACERLO CON UN SWITCH CASE
        RHipAngles_x = ModelData.Raw.(ModelOutput{1})(1,:); 
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
        
% Visualizacion señal antes de segmentarlas
figure(3)
plot(t,RHipAngles_x, t, RHipAngles_y, '--', t, RHipAngles_z,'.');
legend('x','y','z');
title('Right hip angles');

figure(4)
plot(t,LHipAngles_x, t, LHipAngles_y, '--', t, LHipAngles_z,'.');   
legend('x','y','z');
title('Left hip angles');

% Right Hip      
bajada_sincaja_Rhip_x = RHipAngles_x(idxO1:idxF1);  % phase 1
bajada_sincaja_Rhip_y = RHipAngles_y(idxO1:idxF1);
bajada_sincaja_Rhip_z = RHipAngles_z(idxO1:idxF1);

subida_concaja_Rhip_x = RHipAngles_x(idxF1:idxF2); % phase 2
subida_concaja_Rhip_y = RHipAngles_y(idxF1:idxF2);
subida_concaja_Rhip_z = RHipAngles_z(idxF1:idxF2);

bajada_concaja_Rhip_x = RHipAngles_x(idxO3:idxF3); % phase 3
bajada_concaja_Rhip_y = RHipAngles_y(idxO3:idxF3);
bajada_concaja_Rhip_z = RHipAngles_z(idxO3:idxF3);

subida_sincaja_Rhip_x = RHipAngles_x(idxF3:idxF4); % phase 4
subida_sincaja_Rhip_y = RHipAngles_y(idxF3:idxF4);
subida_sincaja_Rhip_z = RHipAngles_z(idxF3:idxF4);

% Left Hip 
bajada_sincaja_Lhip_x = LHipAngles_x(idxO1:idxF1);
bajada_sincaja_Lhip_y = LHipAngles_y(idxO1:idxF1);
bajada_sincaja_Lhip_z = LHipAngles_z(idxO1:idxF1);

subida_concaja_Lhip_x = LHipAngles_x(idxF1:idxF2);
subida_concaja_Lhip_y = LHipAngles_y(idxF1:idxF2);
subida_concaja_Lhip_z = LHipAngles_z(idxF1:idxF2);

bajada_concaja_Lhip_x = LHipAngles_x(idxO3:idxF3);
bajada_concaja_Lhip_y = LHipAngles_y(idxO3:idxF3);
bajada_concaja_Lhip_z = LHipAngles_z(idxO3:idxF3);

subida_sincaja_Lhip_x = LHipAngles_x(idxF3:idxF4);
subida_sincaja_Lhip_y = LHipAngles_y(idxF3:idxF4);
subida_sincaja_Lhip_z = LHipAngles_z(idxF3:idxF4);

% time vectors definition
t2 = 0:1:length(bajada_concaja_Rhip_x)-1;
t3 = 0:1:(length(bajada_sincaja_Rhip_x)-1);
t4 = 0:1:(length(subida_concaja_Rhip_x)-1);
t5 = 0:1:(length(subida_sincaja_Rhip_x)-1);

% Right Hip Representacion Separando X,Y,Z y misma fase dist. peso juntos

% PARA COMPARAR SIN CAJA Y CON CAJA Y COMPARAR SEGUN LA DIRECCION X,Y,Z

% Right Hip Bajada
figure(5) 
subplot(3,1,1)
bajada_sincaja_Rhip_x2 = bajada_sincaja_Rhip_x((length(bajada_sincaja_Rhip_x) - (length(bajada_concaja_Rhip_x)-1)): length(bajada_sincaja_Rhip_x));
plot(t2, bajada_sincaja_Rhip_x2, t2, bajada_concaja_Rhip_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Cadera Derecha');
subplot(3,1,2)
bajada_sincaja_Rhip_y2 = bajada_sincaja_Rhip_y((length(bajada_sincaja_Rhip_y) - (length(bajada_concaja_Rhip_y)-1)): length(bajada_sincaja_Rhip_y));
plot(t2, bajada_sincaja_Rhip_y2, t2, bajada_concaja_Rhip_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Cadera Derecha');
subplot(3,1,3)
bajada_sincaja_Rhip_z2 = bajada_sincaja_Rhip_x((length(bajada_sincaja_Rhip_z) - (length(bajada_concaja_Rhip_z)-1)): length(bajada_sincaja_Rhip_z));
plot(t2, bajada_sincaja_Rhip_z2, t2, bajada_concaja_Rhip_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Cadera Derecha');

% Right Hip Subida

figure(6)
subplot(3,1,1)
if length(subida_sincaja_Rhip_x) > length(subida_concaja_Rhip_x)
        subida_sincaja_Rhip_x2 = subida_sincaja_Rhip_x(1:length(subida_concaja_Rhip_x));
        plot(t4, subida_sincaja_Rhip_x2, t4, subida_concaja_Rhip_x, '--');
else if length(subida_sincaja_Rhip_x) < length(subida_concaja_Rhip_x)
        subida_concaja_Rhip_x2 = subida_concaja_Rhip_x(1:length(subida_sincaja_Rhip_x));
        plot(t5, subida_sincaja_Rhip_x, t5, subida_concaja_Rhip_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Cadera Derecha');
subplot(3,1,2)
if length(subida_sincaja_Rhip_y) > length(subida_concaja_Rhip_y)
        subida_sincaja_Rhip_y2 = subida_sincaja_Rhip_y(1:length(subida_concaja_Rhip_y));
        plot(t4, subida_sincaja_Rhip_y2, t4, subida_concaja_Rhip_y, '--');
else if length(subida_sincaja_Rhip_y) < length(subida_concaja_Rhip_y)
        subida_concaja_Rhip_y2 = subida_concaja_Rhip_y(1:length(subida_sincaja_Rhip_y));
        plot(t5, subida_sincaja_Rhip_y, t5, subida_concaja_Rhip_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Cadera Derecha');
subplot(3,1,3)
if length(subida_sincaja_Rhip_z) > length(subida_concaja_Rhip_z)
        subida_sincaja_Rhip_z2 = subida_sincaja_Rhip_z(1:length(subida_concaja_Rhip_z));
        plot(t4, subida_sincaja_Rhip_z2, t4, subida_concaja_Rhip_z, '--');
else if length(subida_sincaja_Rhip_z) < length(subida_concaja_Rhip_z)
        subida_concaja_Rhip_z2 = subida_concaja_Rhip_z(1:length(subida_sincaja_Rhip_z));
        plot(t5, subida_sincaja_Rhip_z, t5, subida_concaja_Rhip_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Cadera Derecha');

% Left Hip Representacion Separando X,Y,Z y misma fase dist. peso juntos
% Left Hip Bajada
figure(7)
t2 = 0:1:length(bajada_concaja_Lhip_x)-1;
subplot(3,1,1)
bajada_sincaja_Lhip_x2 = bajada_sincaja_Lhip_x((length(bajada_sincaja_Lhip_x) - (length(bajada_concaja_Lhip_x)-1)): length(bajada_sincaja_Lhip_x));
plot(t2, bajada_sincaja_Lhip_x2, t2, bajada_concaja_Lhip_x, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje X Cadera Izquierda');
subplot(3,1,2)
bajada_sincaja_Lhip_y2 = bajada_sincaja_Lhip_y((length(bajada_sincaja_Lhip_y) - (length(bajada_concaja_Lhip_y)-1)): length(bajada_sincaja_Lhip_y));
plot(t2, bajada_sincaja_Lhip_y2, t2, bajada_concaja_Lhip_y, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Y Cadera Izquierda');
subplot(3,1,3)
bajada_sincaja_Lhip_z2 = bajada_sincaja_Lhip_x((length(bajada_sincaja_Lhip_z) - (length(bajada_concaja_Lhip_z)-1)): length(bajada_sincaja_Lhip_z));
plot(t2, bajada_sincaja_Lhip_z2, t2, bajada_concaja_Lhip_z, '--');
legend('sin caja', 'con caja');
title('Bajada sin caja vs Bajada con caja Eje Z Cadera Izquierda');

% Left Hip Subida

figure(8)
subplot(3,1,1)
if length(subida_sincaja_Lhip_x) > length(subida_concaja_Lhip_x)
        subida_sincaja_Lhip_x2 = subida_sincaja_Lhip_x(1:length(subida_concaja_Lhip_x));
        plot(t4, subida_sincaja_Lhip_x2, t4, subida_concaja_Lhip_x, '--');
else if length(subida_sincaja_Lhip_x) < length(subida_concaja_Lhip_x)
        subida_concaja_Lhip_x2 = subida_concaja_Lhip_x(1:length(subida_sincaja_Lhip_x));
        plot(t5, subida_sincaja_Lhip_x, t5, subida_concaja_Lhip_x2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje X Cadera Izquierda');
subplot(3,1,2)
if length(subida_sincaja_Lhip_y) > length(subida_concaja_Lhip_y)
        subida_sincaja_Lhip_y2 = subida_sincaja_Lhip_y(1:length(subida_concaja_Lhip_y));
        plot(t4, subida_sincaja_Lhip_y2, t4, subida_concaja_Lhip_y, '--');
else if length(subida_sincaja_Lhip_y) < length(subida_concaja_Lhip_y)
        subida_concaja_Lhip_y2 = subida_concaja_Lhip_y(1:length(subida_sincaja_Lhip_y));
        plot(t5, subida_sincaja_Lhip_y, t5, subida_concaja_Lhip_y2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Y Cadera Izquierda');
subplot(3,1,3)
if length(subida_sincaja_Lhip_z) > length(subida_concaja_Lhip_z)
        subida_sincaja_Lhip_z2 = subida_sincaja_Lhip_z(1:length(subida_concaja_Lhip_z));
        plot(t4, subida_sincaja_Lhip_z2, t4, subida_concaja_Lhip_z, '--');
else if length(subida_sincaja_Lhip_z) < length(subida_concaja_Lhip_z)
        subida_concaja_Lhip_z2 = subida_concaja_Lhip_z(1:length(subida_sincaja_Lhip_z));
        plot(t5, subida_sincaja_Lhip_z, t5, subida_concaja_Lhip_z2, '--');
    end
end
legend('sin caja', 'con caja');
title('Subida sin caja vs Subida con caja Eje Z Cadera Izquierda');

% Right Hip Representacion Separando 4 fases e X,Y,Z de cada fase juntos

% PARA COMPARAR DERECHA- IZQUIERDA EN LAS 4 FASES

figure(9)
subplot(4,1,1)
plot (t3, bajada_sincaja_Rhip_x, t3, bajada_sincaja_Rhip_y, '--',t3,bajada_sincaja_Rhip_z, '.');
legend('x','y','z')
title('Bajada sin caja Cadera Derecha');
subplot(4,1,2)
plot(t4, subida_concaja_Rhip_x,t4,subida_concaja_Rhip_y, '--',t4,subida_concaja_Rhip_z, '.');
legend('x','y','z')
title('Subida con caja Cadera Derecha');
subplot(4,1,3)
plot(t2, bajada_concaja_Rhip_x, t2, bajada_concaja_Rhip_y,'--', t2, bajada_concaja_Rhip_z, '.');
legend('x','y','z')
title('Bajada con caja Cadera Derecha');
subplot(4,1,4)
plot(t5, subida_sincaja_Rhip_x, t5, subida_sincaja_Rhip_y,'--', t5, subida_sincaja_Rhip_z, '.');
legend('x','y','z')
title('Subida sin caja Cadera Derecha');

% Left Hip Representacion Separando 4 fases e X,Y,Z de cada fase juntos

figure(10)
subplot(4,1,1)
plot (t3, bajada_sincaja_Lhip_x, t3, bajada_sincaja_Lhip_y, '--',t3,bajada_sincaja_Lhip_z, '.');
legend('x','y','z')
title('Bajada sin caja Cadera Izquierda');
subplot(4,1,2)
plot(t4, subida_concaja_Lhip_x,t4,subida_concaja_Lhip_y, '--',t4,subida_concaja_Lhip_z, '.');
legend('x','y','z')
title('Subida con caja Cadera Izquierda');
subplot(4,1,3)
plot(t2, bajada_concaja_Lhip_x, t2, bajada_concaja_Lhip_y,'--', t2, bajada_concaja_Lhip_z, '.');
legend('x','y','z')
title('Bajada con caja Cadera Izquierda');
subplot(4,1,4)
plot(t5, subida_sincaja_Lhip_x, t5, subida_sincaja_Lhip_y,'--', t5, subida_sincaja_Lhip_z, '.');
legend('x','y','z')
title('Subida sin caja Cadera Izquierda');

% ROM calculation

% Right Hip 
Rhip_ph1(:,1) = bajada_sincaja_Rhip_x;
Rhip_ph1(:,2) = bajada_sincaja_Rhip_y;
Rhip_ph1(:,3) = bajada_sincaja_Rhip_z;

Rhip_ph2(:,1) = subida_concaja_Rhip_x;
Rhip_ph2(:,2) = subida_concaja_Rhip_y;
Rhip_ph2(:,3) = subida_concaja_Rhip_z;

Rhip_ph3(:,1) = bajada_concaja_Rhip_x;
Rhip_ph3(:,2) = bajada_concaja_Rhip_y;
Rhip_ph3(:,3) = bajada_concaja_Rhip_z;

Rhip_ph4(:,1) = subida_sincaja_Rhip_x;
Rhip_ph4(:,2) = subida_sincaja_Rhip_y;
Rhip_ph4(:,3) = subida_sincaja_Rhip_z;

Rhip_ph1_max = max(Rhip_ph1);
Rhip_ph1_min = min(Rhip_ph1);

Rhip_ph2_max = max(Rhip_ph2);
Rhip_ph2_min = min(Rhip_ph2);

Rhip_ph3_max = max(Rhip_ph3);
Rhip_ph3_min = min(Rhip_ph3);

Rhip_ph4_max = max(Rhip_ph4);
Rhip_ph4_min = min(Rhip_ph4);

for i = 1:length(Rhip_ph1_max)
    Rhip_ph1_ROM(1,i) = Rhip_ph1_max(i)- Rhip_ph1_min(i);
end

for i = 1:length(Rhip_ph1_max)
    Rhip_ph2_ROM(1,i) = Rhip_ph2_max(i)- Rhip_ph2_min(i);
end

for i = 1:length(Rhip_ph1_max)
    Rhip_ph3_ROM(1,i) = Rhip_ph3_max(i)- Rhip_ph3_min(i);
end

for i = 1:length(Rhip_ph1_max)
    Rhip_ph4_ROM(1,i) = Rhip_ph4_max(i)- Rhip_ph4_min(i);
end

% Left Hip 
Lhip_ph1(:,1) = bajada_sincaja_Lhip_x;
Lhip_ph1(:,2) = bajada_sincaja_Lhip_y;
Lhip_ph1(:,3) = bajada_sincaja_Lhip_z;

Lhip_ph2(:,1) = subida_concaja_Lhip_x;
Lhip_ph2(:,2) = subida_concaja_Lhip_y;
Lhip_ph2(:,3) = subida_concaja_Lhip_z;

Lhip_ph3(:,1) = bajada_concaja_Lhip_x;
Lhip_ph3(:,2) = bajada_concaja_Lhip_y;
Lhip_ph3(:,3) = bajada_concaja_Lhip_z;

Lhip_ph4(:,1) = subida_sincaja_Lhip_x;
Lhip_ph4(:,2) = subida_sincaja_Lhip_y;
Lhip_ph4(:,3) = subida_sincaja_Lhip_z;

Lhip_ph1_max = max(Lhip_ph1);
Lhip_ph1_min = min(Lhip_ph1);

Lhip_ph2_max = max(Lhip_ph2);
Lhip_ph2_min = min(Lhip_ph2);

Lhip_ph3_max = max(Lhip_ph3);
Lhip_ph3_min = min(Lhip_ph3);

Lhip_ph4_max = max(Lhip_ph4);
Lhip_ph4_min = min(Lhip_ph4);

for i = 1:length(Lhip_ph1_max)
    Lhip_ph1_ROM(1,i) = Lhip_ph1_max(i)- Lhip_ph1_min(i);
end

for i = 1:length(Lhip_ph1_max)
    Lhip_ph2_ROM(1,i) = Lhip_ph2_max(i)- Lhip_ph2_min(i);
end

for i = 1:length(Lhip_ph1_max)
    Lhip_ph3_ROM(1,i) = Lhip_ph3_max(i)- Lhip_ph3_min(i);
end

for i = 1:length(Lhip_ph1_max)
    Lhip_ph4_ROM(1,i) = Lhip_ph4_max(i)- Lhip_ph4_min(i);
end

% Total ROM Hip Angles 
Rhip_total(:,1) = RHipAngles_x;
Rhip_total(:,2) = RHipAngles_y;
Rhip_total(:,3) = RHipAngles_z;

Rhip_total_max = max(Rhip_total);
Rhip_total_min = min(Rhip_total);

for i = 1:length(Rhip_total_max)
    Rhip_total_max_ROM(1,i) = Rhip_total_max(i)- Rhip_total_min(i);
end

Lhip_total(:,1) = LHipAngles_x;
Lhip_total(:,2) = LHipAngles_y;
Lhip_total(:,3) = LHipAngles_z;

Lhip_total_max = max(Lhip_total);
Lhip_total_min = min(Lhip_total);

for i = 1:length(Lhip_total_max)
    Lhip_total_max_ROM(1,i) = Lhip_total_max(i)- Lhip_total_min(i);
end

ROM_Right_Hip = table([Rhip_ph1_ROM(1);Rhip_ph2_ROM(1);Rhip_ph3_ROM(1);Rhip_ph4_ROM(1)],[Rhip_ph1_ROM(2);Rhip_ph2_ROM(2);Rhip_ph3_ROM(2);Rhip_ph4_ROM(2)],[Rhip_ph1_ROM(3);Rhip_ph2_ROM(3);Rhip_ph3_ROM(3);Rhip_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

ROM_Left_Hip = table([Lhip_ph1_ROM(1);Lhip_ph2_ROM(1);Lhip_ph3_ROM(1);Lhip_ph4_ROM(1)],[Lhip_ph1_ROM(2);Lhip_ph2_ROM(2);Lhip_ph3_ROM(2);Lhip_ph4_ROM(2)],[Lhip_ph1_ROM(3);Lhip_ph2_ROM(3);Lhip_ph3_ROM(3);Lhip_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

TotalROM_Hip = table([Rhip_total_max_ROM(1); Lhip_total_max_ROM(1)],[Rhip_total_max_ROM(2); Lhip_total_max_ROM(2)],[Rhip_total_max_ROM(3); Lhip_total_max_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})