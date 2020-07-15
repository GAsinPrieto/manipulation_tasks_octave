% Peaks metric
% Adapted to Octave by Guillermo Asín Prieto

function [y,locs_np] = np_octave(time,speed)
	%pkg load signal
    [pks,locs_np] = findpeaks(speed,'DoubleSided');
	y = -length(findpeaks(speed,'DoubleSided'));	
end