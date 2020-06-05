function [struc_free, struc_exo, no_struc_free, no_struc_exo]= randomise(seed)
%UNTITLED2 Summary of this function goes here
    % momentaneamente quitada output system
% Pensar si voy a randomizar system o empiezo siempre sin exo
%rng(seed)
%system = randperm(2);
%if system(1) == 1
    %disp('Sin exoesqueleto');
%else 
    %disp('Con exoesqueleto');
%end
% Orden tareas con estructura sin exo: Lower-lifting sagittal plane and Lateral Load Transfer
rng(seed);
struc_free = randperm(10)

% Orden tarea sin estructura sin exo: Isolated basic movements 
rng(seed);
no_struc_free = randperm(3)

% Orden tarea sin estructura con exo: Isolated basic movements 
rng(2*seed);
no_struc_exo = randperm(3)

% Orden tareas con estructura con exo: Lower-lifting sagittal plane and Lateral Load Transfer
rng(2*seed);
struc_exo = randperm(10)




