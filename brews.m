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

vg=0:.1:5;
vd=2.5;
n=length(vg);
for i=1:n

QT=@(psi_s) cox*(vg(i)-vfb-psi_s);
QD=@(psi_s) (2*epsi_si*q*na*psi_s).^0.5;
QI=@(psi_s) QT(psi_s)-QD(psi_s);
dvpsi_s=@(psi_s) 1+2*vkbt*(cox*QT(psi_s)+epsi_si*q*na)./((QT(psi_s)).^2-(QD(psi_s)).^2);
vgf=@(psi_s,v) vfb+psi_s+1/cox*(2*epsi_si*vkbt*q*na)^0.5*(psi_s/vkbt+ni^2/na^2*exp((psi_s-v)/vkbt)).^0.5;

psi_s_min=-vg(i)-vd-abs(vfb);
psi_s_max=vg(i)+vd+abs(vfb);
dpsi_s=1e-3;
psi_svec=psi_s_min:dpsi_s:psi_s_max;

psi_ss=interp1(real(vgf(psi_svec,0)),real(psi_svec),vg(i));
psi_sd=interp1(real(vgf(psi_svec,vd)),real(psi_svec),vg(i));

intf=@(psi_s) QI(psi_s).*dvpsi_s(psi_s);

id(i)=mu_eff*w/l*integral(intf,psi_ss,psi_sd);
end
figure(1)
semilogy(vg,abs(id*1e6),'linewidth',2)
xlabel('Vgs (V)','FontSize',15)
ylabel('log(Id) (uA)','FontSize',15)
title('log(Id) VS Vgs ','FontSize',15)
%legend('vd=2v','vd=1.5v','vd=1v')
legend('vd=1.8v','vd=2.2v','vd=2.5v')
hold on
figure(2)
plot(vg,abs(id*1e6),'linewidth',2)
xlabel('Vgs (V)','FontSize',15)
ylabel('Id (uA)','FontSize',15)
title('Id VS Vgs with Vds=1V','FontSize',15)
legend('vd=1.8v','vd=2.2v','vd=2.5v')
hold on