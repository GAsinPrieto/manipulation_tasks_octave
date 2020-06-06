% Numerical integration, trapezoidal
% Constant time-step size
%
function y = integ(t,v)
	deltat = t(2)-t(1);
	y = 0;
	for i=1:length(v)-1
		y = y + (v(i)+v(i+1))*deltat/2;
	end
end
