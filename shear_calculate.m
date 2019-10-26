function [shear] = shear_calculate(type,rpm,posx)
%%%% type is 'tapered' or 'straight',
%%%% shera stress formula from DOI: 10.1039/B712486D

if nargin < 3 || isempty(posx)
    posx=0;
end

Q_1rpm=0.024;  %%% ul/min
Q_ml_min=Q_1rpm*rpm;  %%% ul/min
Q= Q_ml_min*(10^-6)/60;%%% m3/s
mu=0.8*10^-3; %%% Pa.s
h=10^-3; %%% height channel in meters
n=2;
if strcmp(type,'tapered'); 
    
    w= posx*(10^-6)/2 +(10^-3);  %%% posx is in um
elseif strcmp(type,'straight');
    w=h;
elseif strcmp(type,'step')
    if posx<1500 ; w=h;
    elseif posx>1500 & posx<4500; w=3*h;
    elseif posx>4500 & posx<8000; w=5*h;
    elseif  posx>8000; w=nan;
    end
end
 %m=1.7+ 0.5*(h/w)^(1.4);

% shear= (2*mu*Q_L_s/(w*(h^2)))*(1+1/m)*(n+1);  %%%% Pa
 shear= (6*mu*Q./(w*(h^2)));  %%%% Pa two plate approssimation
 shear=shear*10; %%%% dyne/cm2   
end

