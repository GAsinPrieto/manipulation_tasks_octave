% Peaks metric

function [y,locs_np] = np(time,speed)
	%pkg load signal
    [pks,locs_np] = findpeaks(speed);
	y = -length(findpeaks(speed));	
end