function f = fseries(t,x,n)
% Function calculates Fourier series of a periodic function x(t) from one 
% given period of (t,x) data to order n
% t = independent variable data
% x = dependent variable data
% n = order of series
    % Period
    T = t(end)-t(1);
    % Zeroth term
    a0 = 2/T*trapz(t-t(1),x);
    f = @(t) a0/2;
    % Higher order terms
    for i=1:n
        ai = 2/T*trapz(t,x.*cos((t-t(1))*2*pi*i/T));
        bi = 2/T*trapz(t,x.*sin((t-t(1))*2*pi*i/T));
        f = @(t) f(t)+ai*cos((t-t(1))*2*pi*i/T)+bi*sin((t-t(1))*2*pi*i/T);
    end
end