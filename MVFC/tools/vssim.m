function y=vssim(Vx,Vy,u,v)

U  = [Vx,Vy];
U0 = [u,v];
I =  isfinite(U+U0);
U(~I)=0;
L  = max(abs(U(:)))+eps;
y  =ssim(U,U0,'DynamicRange',L);
end