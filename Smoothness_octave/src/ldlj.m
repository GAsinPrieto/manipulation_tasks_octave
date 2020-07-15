% Dimensionless jerk (DLJ) and Log Dimensionless Jerk (LDLJ)
%
function [dlj ldlj] = ldlj(time,speed)
	vPeak = max(speed);
	%
	DDspeed2 = ddv(time,speed).^2;
	%	
	dlj = - integ(time,DDspeed2)*( (time(end)-time(1))^3 ) / (vPeak^2);
	ldlj = -log(abs(dlj));
end