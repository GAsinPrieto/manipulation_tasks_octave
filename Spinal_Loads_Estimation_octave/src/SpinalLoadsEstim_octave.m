%% Spinal Loads Estimation using simplified regression equations at flexed postures

% Based on "Subject-specific regression equations to estimate lower spinal loads during symmetric and asymmetric static lifting", F. Ghezelbash et al., Journal of Biomechanics, 2020.

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

% Comparison between "No Exo" and "Exo" trials

                            % NO EXO                 % EXO

% FLOOR-KNEE 0              44                      17_02       
% FLOOR-SHOULDER 0          50                      53             
% KNEE-SHOULDER  0          16                      57          
% FLOOR-KNEE 10%            01                      59            
% FLOOR-SHOULDER 10%        03                      58         
% KNEE-SHOULDER  10%        05                      61          
%
% KNEE-KNEE 0               42                      56          
% SHOULDER-SHOULDER 0       21                      62          
% KNEE-KNEE 10%             08                      54          
% SHOULDER-SHOULDER 10%     47                      55          

pkg load signal

clear all % Clear variables
close all % Close figures
clc

load('../tests/data/input/dinamica21_B.mat');
%load('../tests/data/input/dinamica62_B.mat');

% sin numero es no exo
% 2 es con exo

t_total = (double(frames)*Ts);
%t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));
t_100 = 1:1:100;

%% Input Parameters Calculation

% External Load magnitude(kg)
%M = round(10/100*BW); %10 %, pero no siempre es asi, hay trials que son 0
M = 0;

% Shoulder and Box trajectories (mm)
Shoulder = LSHO(:,2)'; % (-1)*
Load = SIDEBOX1(:,2)'; % (-1)*
D_mm = abs(Shoulder - Load);
D_cm = D_mm/10;
%D_cm2 = D_mm/10;

% Trunk flexion or Load Elevation
LThoraxAngles_x = ModelData.Raw.(ModelOutput{14})(1,:); 
F = LThoraxAngles_x;
%F2 = LThoraxAngles_x;


%% Spinal loads estimation

% L4L5 Compression

% no exo
for i = 1:frames
   L4L5_compression(1,i) = 92.21616 + 8.722677*BW -3.24825*F(1,i) -4.71189*D_cm(1,i) + 0.231067*BW*F(1,i) + 0.048462*BW*D_cm(1,i) +0.905165*M*F(1,i) +2.441194*M*D_cm(1,i) + 0.052136*F(1,i)*D_cm(1,i);
end
% exo
%for i = 1:frames
 %   L4L5_compression2(1,i) = 92.21616 + 8.722677*(BW) -3.24825*F(1,i) -4.71189*D_cm(1,i) + 0.231067*(BW)*F(1,i) + 0.048462*(BW)*D_cm(1,i) +0.905165*M*F(1,i) +2.441194*M*D_cm(1,i) + 0.052136*F(1,i)*D_cm(1,i);
%end

% L4L5 Shear

% no exo
for i = 1:frames
    L4L5_shear(1,i) = -70.1973 + 1.314008*BW + 3.130555*M + 0.49279*F(1,i) + 1.049813*D_cm(1,i) -0.02008*BW*M + 0.034287*BW*F(1,i) + 0.108536*M*F(1,i) + 0.302128*M*D_cm(1,i) -0.00393*F(1,i)*D_cm(1,i);
end
% exo
%for i = 1:frames
 %   L4L5_shear2(1,i) = -70.1973 + 1.314008*(BW) + 3.130555*M + 0.49279*F(1,i) + 1.049813*D_cm(1,i) -0.02008*(BW)*M + 0.034287*(BW)*F(1,i) + 0.108536*M*F(1,i) + 0.302128*M*D_cm(1,i) -0.00393*F(1,i)*D_cm(1,i);
%end

% L5S1 Compression
% no exo
for i = 1:frames
    L5S1_compression(1,i) = 71.3142 + 9.654267*BW -5.42554*F(1,i) + 0.069926*BW*M +0.240893*BW*F(1,i) +0.029001*BW*D_cm(1,i) +0.913375*M*F(1,i) +2.538994*M*D_cm(1,i) +0.022737*F(1,i)*D_cm(1,i);
end
% exo
%for i = 1:frames
 %   L5S1_compression2(1,i) = 71.3142 + 9.654267*(BW) -5.42554*F(1,i) + 0.069926*(BW)*M +0.240893*(BW)*F(1,i) +0.029001*(BW)*D_cm(1,i) +0.913375*M*F(1,i) +2.538994*M*D_cm(1,i) +0.022737*F(1,i)*D_cm(1,i);
%end

% L5S1 Shear
% no exo
for i = 1:frames
    L5S1_shear(1,i) = -49.8027 + 4.434174*BW + 3.66045*M +0.068312*BW*F(1,i) + 0.22104*M*F(1,i) +1.068737*M*D_cm(1,i) -0.00484*F(1,i)*D_cm(1,i);
end
% exo
%for i = 1:frames
 %   L5S1_shear2(1,i) = -49.8027 + 4.434174*(BW) + 3.66045*M +0.068312*(BW)*F(1,i) + 0.22104*M*F(1,i) +1.068737*M*D_cm(1,i) -0.00484*F(1,i)*D_cm(1,i);
%end

%% Visualization

% (1:(frames1/100):frames1)  si frames>1000
% (1:(frames1/100):frames1+(frames1/100)) si frames<1000 (dinamica 08)

figure(1)
frames1 =length(L4L5_compression);
frames2=length(L4L5_compression2);

plot(t_100, L4L5_compression(1:(frames1/100):frames1),'r');
hold on
plot(t_100, L4L5_compression2(1:(frames2/100):frames2),'b');
xlabel('% Lifting Cycle');
ylabel('Force (N)');
legend({'No exoskeleton','Exoskeleton'},'Location','best');
%title('Floor-to-Knee 0 kg. L4-L5 compression forces: Without Exo vs Exoskeleton.')
%title('Floor-to-Shoulder 0 kg. L4-L5 compression forces: Without Exo vs Exoskeleton.')     %2
%title('Knee-to-Shoulder 0 kg. L4-L5 compression forces: Without Exo vs Exoskeleton.')      %3
%title('Floor-to-Knee 10% BW kg. L4-L5 compression forces: Without Exo vs Exoskeleton.')    %4
%title('Floor-to-Shoulder 10% BW kg. L4-L5 compression forces: Without Exo vs Exoskeleton.')%5
%title('Knee-to-Shoulder 10% BW kg. L4-L5 compression forces: Without Exo vs Exoskeleton.') %6
%title('Knee-to-Knee 0 kg. L4-L5 compression forces: Without Exo vs Exoskeleton.')          %7
title('Shoulder-to-Shoulder 0 kg. L4-L5 compression forces: Without Exo vs Exoskeleton.')  %8
%title('Knee-to-Knee 10% BW kg. L4-L5 compression forces: Without Exo vs Exoskeleton.')     %9
%title('Shoulder-to-Shoulder 10% BW kg. L4-L5 compression forces: Without Exo vs Exoskeleton.')  %10

figure(2)
frames1 =length(L4L5_compression);
frames2=length(L4L5_compression2);
plot(t_100, L4L5_shear(1:(frames1/100):frames1),'r');
hold on
plot(t_100, L4L5_shear2(1:(frames2/100):frames2),'b');
xlabel('% Lifting Cycle');
ylabel('Force (N)');
legend({'No exoskeleton','Exoskeleton'},'Location','best');
%title('Floor-to-Knee 0 kg. L4-L5 shear forces: Without Exo vs Exoskeleton.')
%title('Floor-to-Shoulder 0 kg. L4-L5 shear forces: Without Exo vs Exoskeleton.')     %2
%title('Knee-to-Shoulder 0 kg. L4-L5 shear forces: Without Exo vs Exoskeleton.')      %3
%title('Floor-to-Knee 10% BW kg. L4-L5 shear forces: Without Exo vs Exoskeleton.')    %4
%title('Floor-to-Shoulder 10% BW kg. L4-L5 shear forces: Without Exo vs Exoskeleton.')%5
%title('Knee-to-Shoulder 10% BW kg. L4-L5 shear forces: Without Exo vs Exoskeleton.') %6
%title('Knee-to-Knee 0 kg. L4-L5 shear forces: Without Exo vs Exoskeleton.')          %7
title('Shoulder-to-Shoulder 0 kg. L4-L5 shear forces: Without Exo vs Exoskeleton.')  %8
%title('Knee-to-Knee 10% BW kg. L4-L5 shear forces: Without Exo vs Exoskeleton.')     %9
%title('Shoulder-to-Shoulder 10% BW kg. L4-L5 shear forces: Without Exo vs Exoskeleton.')  %10

figure(3)
frames1 =length(L4L5_compression);
frames2=length(L4L5_compression2);
plot(t_100, L5S1_compression(1:(frames1/100):frames1),'r');
hold on
plot(t_100, L5S1_compression2(1:(frames2/100):frames2),'b');
xlabel('% Lifting Cycle');
ylabel('Force (N)');
legend({'No exoskeleton','Exoskeleton'},'Location','best');
%title('Floor-to-Knee 0 kg. L5-S1 compression forces: Without Exo vs Exoskeleton.')
%title('Floor-to-Shoulder 0 kg. L5-S1 compression forces: Without Exo vs Exoskeleton.')     %2
%title('Knee-to-Shoulder 0 kg. L5-S1 compression forces: Without Exo vs Exoskeleton.')      %3
%title('Floor-to-Knee 10% BW kg. L5-S1 compression forces: Without Exo vs Exoskeleton.')    %4
%title('Floor-to-Shoulder 10% BW kg. L5-S1 compression forces: Without Exo vs Exoskeleton.')%5
%title('Knee-to-Shoulder 10% BW kg. L5-S1 compression forces: Without Exo vs Exoskeleton.') %6
%title('Knee-to-Knee 0 kg. L5-S1 compression forces: Without Exo vs Exoskeleton.')          %7
title('Shoulder-to-Shoulder 0 kg. L5-S1 compression forces: Without Exo vs Exoskeleton.')  %8
%title('Knee-to-Knee 10% BW kg. L5-S1 compression forces: Without Exo vs Exoskeleton.')     %9
%title('Shoulder-to-Shoulder 10% BW kg. L5-S1 compression forces: Without Exo vs Exoskeleton.')  %10

figure(4)
frames1 =length(L4L5_compression);
frames2=length(L4L5_compression2);
plot(t_100, L5S1_shear(1:(frames1/100):frames1),'r');
hold on
plot(t_100, L5S1_shear2(1:(frames2/100):frames2),'b');
xlabel('% Lifting Cycle');
ylabel('Force (N)');
legend({'No exoskeleton','Exoskeleton'},'Location','best');
%title('Floor-to-Knee 0 kg. L5-S1 shear forces: Without Exo vs Exoskeleton.')
%title('Floor-to-Shoulder 0 kg. L5-S1 shear forces: Without Exo vs Exoskeleton.')     %2
%title('Knee-to-Shoulder 0 kg. L5-S1 shear forces: Without Exo vs Exoskeleton.')      %3
%title('Floor-to-Knee 10% BW kg. L5-S1 shear forces: Without Exo vs Exoskeleton.')    %4
%title('Floor-to-Shoulder 10% BW kg. L5-S1 shear forces: Without Exo vs Exoskeleton.')%5
%title('Knee-to-Shoulder 10% BW kg. L5-S1 shear forces: Without Exo vs Exoskeleton.') %6
%title('Knee-to-Knee 0 kg. L5-S1 shear forces: Without Exo vs Exoskeleton.')          %7
title('Shoulder-to-Shoulder 0 kg. L5-S1 shear forces: Without Exo vs Exoskeleton.')  %8
%title('Knee-to-Knee 10% BW kg. L5-S1 shear forces: Without Exo vs Exoskeleton.')     %9
%title('Shoulder-to-Shoulder 10% BW kg. L5-S1 shear forces: Without Exo vs Exoskeleton.')  %10

figure(5)
frames1 =length(F);
frames2=length(F2);
plot(t_100, F(1:(frames1/100):frames1),'r');
hold on
plot(t_100, F2(1:(frames2/100):frames2),'b');
xlabel('% Lifting Cycle');
ylabel('Thorax Angles (Degrees)');
legend({'No exoskeleton','Exoskeleton'},'Location','best');
title('Knee-to-Shoulder 0 kg. L5-S1 compression forces: Without Exo vs Exoskeleton.')      %3

figure(6)
frames1 =length(F);
frames2=length(F2);
plot(t_100, D_cm(1:(frames1/100):frames1),'r');
hold on
plot(t_100, D_cm2(1:(frames2/100):frames2),'b');
xlabel('% Lifting Cycle');
ylabel('D (cm)');
legend({'No exoskeleton','Exoskeleton'},'Location','best');
title('Knee-to-Shoulder 0 kg. L5-S1 compression forces: Without Exo vs Exoskeleton.')      %3