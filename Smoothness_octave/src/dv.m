% First numerical derivative, central differences
% Constant time-step size
%
function y = dv(t,v)
	deltat = t(2)-t(1);
	y = zeros(length(v),1);
	y(1,1) = v(2)-v(1);
	for i=2:length(v)-1
		y(i,1) = ( v(i+1)-v(i-1) ) / 2;
	end
	y(end,1) = v(end)-v(end-1);
	y = y./deltat;
end