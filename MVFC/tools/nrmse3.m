function y=nrmse3(Vx,Vy,Vz,u,v,w)
%3D version of nrmse
y=(norm(Vx(:)-u(:))^2+norm(Vy(:)-v(:))^2+norm(Vz(:)-w(:))^2)^0.5/(norm(u(:))^2+norm(v(:))^2+norm(w(:))^2+eps)^0.5;
end