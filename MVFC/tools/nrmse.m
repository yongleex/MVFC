function y=NRMSE(Vx,Vy,u,v)
% 2d2c data NRMSE
I =  isfinite(Vx+Vy+u+v);
Vx(~I)=0;Vy(~I)=0;u(~I)=0;v(~I)=0;

y=(norm(Vx-u)^2+norm(Vy-v)^2)^0.5/(norm(u)^2+norm(v)^2+eps)^0.5;
end