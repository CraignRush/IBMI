%% Problem 1
close all;
clear all;
df = 1; % Hz off-resonance.
T1 = 1000; % ms.
T2 = 80; % ms.
TE = 40; % ms.
TR = 1500; % ms.
TI = 300; % ms.
flip1 = pi; % pi y pulse
flip2 = -pi/2;% -pi/2 y pulse
flip3 = pi; % pi x pulse
dT = 1;
Ntr = round(TR/dT);
Nti = round(TI/dT);
Nte = round(TE/dT);
Ntehalf = round(TE/(2*dT));
Nex = 5; % number of RF excitations.
% magnetization vector
M=zeros(3,Nex*Ntr);
% initial magnetization
M(:,1) = [0;0;1];

% Bloch equation matrices
Rflip1 = yrot(flip1);
Rflip2 = yrot(flip2);
Rflip3 = xrot(flip3);
[A1,B1] = freeprecess(dT,T1,T2,df);

%% Simulate Bloch equations
Mcount=1;
% simulate Nex TRs
for n=1:Nex
% pi RF excitation
M(:,Mcount) = Rflip1*M(:,Mcount);
% free-precession and relaxation over period TI
for k=1:Nti
Mcount=Mcount+1;
M(:,Mcount)=A1*M(:,Mcount-1)+B1;
end
% pi/2 RF excitation
M(:,Mcount) = Rflip2*M(:,Mcount);
% free-precession and relaxation over period TR-TI
for k=1:Ntehalf
Mcount=Mcount+1;
M(:,Mcount)=A1*M(:,Mcount-1)+B1;
end
% pi RF excitation
M(:,Mcount) = Rflip3*M(:,Mcount);
% free-precession and relaxation over period TR-TI-TE/2
for k=1:(Ntr-Nti-Ntehalf)
Mcount=Mcount+1;
M(:,Mcount)=A1*M(:,Mcount-1)+B1;
end

end
time = (0:Mcount-1)*dT;
%% a)
% verify formation of echo by plotting angle of transverse magnetization
% for second TR
Mtrans=M(1,:)+ 1i*M(2,:);
figure;
plot(time(Ntr:2*Ntr),angle(Mtrans(Ntr:2*Ntr)),'b-');
hold on;
scatter(Ntr+TI+TE,angle(Mtrans(Ntr+TI+TE+1)), 300, 'red');
xlabel('Time (ms)');
ylabel('Angle of transverse magnetization (rad)');
title('spin echo sequence');
%% b)
% plot magnetization evolution
figure;
time = [0:Mcount-1]*dT;
plot(time,M(1,:),'b-',time,M(2,:),'g-',time,M(3,:),'r-');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;
title('inversion recovery spin echo sequence');
%% c)
% signal at TE of the 5th TR based on numerical Bloch equation simulation
M_echo_5thTR = abs(M(:,4*Ntr+ Nti + Nte));
M_echo_5thTR_x = M_echo_5thTR(1,:)
% signal based on analytical formula
M_anyalytical_x = abs((exp(-TE/T2))*(1-2*exp(-TI/T1)+exp(-TR/T1)))
%% d)
T1fat=360; %ms
TR=1500; %ms
TI=-T1fat*log(.5+.5*exp(-TR/T1fat))
S=(1-2*exp(-TI/T1fat) + exp(-TR/T1fat))