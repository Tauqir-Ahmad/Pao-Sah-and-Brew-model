%SEMICONDUCTOR PARAMETER AND CONSTANTS
q=1.6e-19;
epsi_0=8.85e-12;
vkbt=26e-3;
Nsub=-4e17*1e6;             %negative for NMOS 
epsi_si=11.7*epsi_0;
epsi_ox=3.9*epsi_0;
eg=1.1*q;
ni=1.5e10*1e6;
na=abs(Nsub);
w=1e-6;
l=1e-6;
mu_eff=200*1e-4;
chi_si=4.05*q;
tox=10e-9;
cox=epsi_ox/tox;
phi_m=chi_si/q;
phi_b=-sign(Nsub)*vkbt*log(abs(Nsub)/ni);
phi_s=chi_si/q+eg/(2*q)+phi_b;
vfb=phi_m-phi_s;


vg=1.5;
vd=0:.05:2;
n=length(vd);
for i=1:n
psi_s_min=-vg-vd(i)-abs(vfb);
psi_s_max=vg+vd(i)+abs(vfb);
dpsi_s=1e-3;

psi_svec=psi_s_min:dpsi_s:psi_s_max;

idvs=[];
dv=1e-3;
del_psi=10e-3;

for v=0:dv:vd(i)
    f1= @(psi) ni^2/na*exp((psi-v)/vkbt);
    f2= @(psi) (2*vkbt*q*na/epsi_si).^0.5*(psi/vkbt+f1(psi)/na).^0.5;
    f3= @(psi_s) vfb+psi_s+epsi_si/cox*f2(psi_s);
    f1byf2=@(psi) f1(psi)./f2(psi);
    
    vgs=f3(psi_svec);
    psi_s=interp1(real(vgs),real(psi_svec),vg);
    inn_int=integral(f1byf2,del_psi,psi_s);
    idvs=[idvs q*mu_eff*w/l*inn_int];
end

id(i)=sum(idvs)*dv;
end   
%semilogy(vd,id*1e6,'linewidth',2)
%xlabel('Vds (V)','FontSize',15)
%ylabel('log(Id) (uA)','FontSize',15)
%title('log(Id) VS Vds with Vgs=1.5V','FontSize',15)
%figure(2)
plot(vd,abs(id*1e6),'linewidth',2)
xlabel('Vds (V)','FontSize',15)
ylabel('(Id) (uA)','FontSize',15)
title('(Id) VS Vds with Vgs=1.5V','FontSize',15)