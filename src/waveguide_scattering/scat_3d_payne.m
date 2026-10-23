function integrant=scat_3d_payne(theta,beta,lambda,nclad,ncore,Lc,sigma)
% integrant=scat_3d_payne(theta,beta,lambda,nclad,ncore,Lc,sigma)
% Integrant of Payne's 2D Theory
% Testing 3D scattering losses theory
% Tymon Barwicz NSL, MIT 2002

Phi=1;
epsilon=8.854e-12;
mu=4*pi*1e-7;
c=(epsilon*mu)^(-0.5);
omega=2*pi*c./lambda;
Cni=(ncore^2-nclad^2);
k1=nclad.*omega./c;

integrant=1./4.*Phi.^2.*Cni.^2.*omega.^3.*(epsilon.*mu).^(3./2).*sigma.^2.*Lc./(1+Lc.^2.*(beta-k1.*cos(theta)).^2)./pi.*(epsilon./mu).^(1./2);