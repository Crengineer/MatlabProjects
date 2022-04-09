%% Project work of Marina Colombo, Elena Bongiorni, Marco Crespi
% Three-rooms building example
% Marcello Farina, 19/12/2018

u=2; %Transmittance between the rooms
U=0.5; %
VA=240; %Volume room A - constant
VBC=24; %Volume room B, volume room C - constant
sr=12; %Wall surface between A and B, wall surface between A and C - constant
sA=84; %Wall surface between A and the environment - constant
sBC=24; %Wall surface between B and the environment, same for room C - constant
c=1.225*1005; %specific heat of the air - constant

%TE=0 / TA=TB=TC=T_const= 20°
%qa=qa_const=sA*U*T_const / qB=qB_const=qC=qC_const=sBC*u*T_const 
%deltaTA= TA - T_const / deltaTB= TB - T_const / deltaTC= TC - T_const
%deltaqA= qA - qA_const / deltaqB= qB - qB_const / deltaqC= qC - qC_const

gamma=sr*u/(c*VA); %higher u, higher gamma
GAMMA=sr*u/(c*VBC); %higher u, higher GAMMA
gammaA=sA*U/(c*VA); %higher U, higher gammaA
gammar=sBC*U/(c*VBC); %higher U, higher gammar

Atot=[-(GAMMA+gammar) GAMMA             0
      gamma           -(2*gamma+gammaA) gamma
      0               GAMMA             -(GAMMA+gammar)];

Bdec{1}=[1/c/VBC;
         0;
         0];
Bdec{2}=[0;
         1/c/VA;
         0];
Bdec{3}=[0;
         0;
         1/c/VBC];
Cdec{1}=[1 0 0];
Cdec{2}=[0 1 0];
Cdec{3}=[0 0 1];

%Definition of starting condition of the system, taking into account: deltaT= T° - 20°
x_0=[-20 -20 -20;]'; %change of sign in initial condition means that the graph is mirrored 
                     %starting temperature T = 0° -> T=TE

%% 1)

%1 System matrices in continuous and discrete time
A=Atot;
B=[Bdec{1}, Bdec{2}, Bdec{3}];
C=[Cdec{1}; Cdec{2}; Cdec{3}];
D=0;
systemCT=ss(A, B, C, D);
i=0;
h=300; %Sample Time
systemDT=c2d(systemCT, h);
[F,G,H,L,Ts]=ssdata(systemDT);
for i=1:3
    Gdec{i}=G(:,i);
    Hdec=Cdec;
end

%1a eigenvalues and spectral abscissa of CT system
eigA=eig(A)
rhoCT=max(real(eig(A)))
%1b eigenvalues and spectral abscissa of DT system
eigF=eig(F)
rhoDT=max(abs(eig(F)))
%the system is asymptoically stable since in CT all eigenvalues
%have strictly negative real part and in DT all eigenvalues are
%inside the unit circle.
%% 2) CENTRALIZED CONTROL STRUCTURE

N=3;
ContStrucC=ones(N,N);
rounding=4;

%2a continuous time fixed modes
CFMct=di_fixed_modes(A,Bdec,Cdec,N,ContStrucC,rounding)

%2b discrete time fixed modes
CFMdt=di_fixed_modes(F,Gdec,Hdec,N,ContStrucC,rounding)

%2c stabilizing CT control gain
[K_c,rho_c,feas_c]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC)

%2d stabilizing DT control gain
[KD_c,rhoD_c,feasD_c]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStrucC)

%2e eigenvalues, stability, CL trajectories in CT e DT
An_c=A+B*K_c

h=300;
k=0;
tt=h/3;
for t=0:0.1:tt; 
    k=k+1;
    yfree_c_ct(:,k)=C*expm(An_c*t)*x_0;
    u_c_ct(:,k)=K_c*yfree_c_ct(:,k);
end

figure(1)
subplot 211
hold on
grid on
plot(0:1:tt,yfree_c_ct(1,1:tt+1),0:1:tt,yfree_c_ct(2,1:tt+1),0:1:tt,yfree_c_ct(3,1:tt+1));
title ('Centralized Control in Continuous Time')
legend('\deltaTB','\deltaTA','\deltaTC')
xlabel(['Interval time from 0 to ', num2str(tt)])
ylabel(['Output free motions for Centralized Control'])

subplot 212
hold on
grid on
plot(0:1:tt,u_c_ct(1,1:tt+1),0:1:tt,u_c_ct(2,1:tt+1),0:1:tt,u_c_ct(3,1:tt+1))
title ('Control action for Centralized Control, Continuous Time')
legend('\deltaqB','\deltaqA','\deltaqC')
xlabel(['Interval time from 0 to ',  num2str(tt)])
ylabel(['Static Output-feedback Control Law u(t)'])


eigen_CT_centralized=eig(An_c)
rho_CT_centralized=max(abs(real(eig(An_c))))

h=300;
i=1;
k=0;
Tfinal=10*h;
steps=[0:Tfinal/h];
yfree_c_dt=zeros(3,length(steps));
u_c_dt=zeros(3,length(steps));
for k=steps  
    yfree_c_dt(:,i)=(F+G*KD_c)^(k)*x_0;
    u_c_dt(:,i)=KD_c*yfree_c_dt(:,k+1);
    i=i+1;
end

figure(2)
subplot 211
hold on
grid on
plot(steps*h,yfree_c_dt(1,:),'ko',steps*h,yfree_c_dt(2,:),'k*',steps*h,yfree_c_dt(3,:),'k.');
title ('Centralized Control in Discrete Time')
legend('\deltaTB','\deltaTA','\deltaTC')
xlabel(['Time instants from 0 to 3000 [s] with sampling time h= ',  num2str(h)])
ylabel(['Output free motions for Centralized Control'])

subplot 212
hold on
grid on
plot(steps*h,u_c_dt(1,:),'ko',steps*h,u_c_dt(2,:),'k*',steps*h,u_c_dt(3,:),'k.');
title ('Control action for Centralized Control, Discrete Time')
legend('\deltaqB','\deltaqA','\deltaqC')
xlabel(['Time instants from 0 to 3000 [s] with sampling time h= ', num2str(h)])
ylabel(['Static Output-feedback Control Law u(kh)'])

eigen_DT_centralized=eig(F+G*KD_c)
rho_DT_centralized=max(abs(real(eigen_DT_centralized)))

%% DECENTRALIZED CONTROL STRUCTURE

N=3;
ContStrucD=eye(3,3);
rounding=4;
%2a continuous time fixed modes
DFMct=di_fixed_modes(A,Bdec,Cdec,N,ContStrucD,rounding)

%2b discrete time fixed modes
DFMdt=di_fixed_modes(F,Gdec,Hdec,N,ContStrucD,rounding)

%2c stabilizing CT control gain
[K_d,rho_d,feas_d]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucD)

%2d stabilizing DT control gain
[KD_d,rhoD_d,feasD_d]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStrucD)

%2e eigenvalues, stability, CL trajectories in CT e DT
An_d=A+B*K_d

h=300;
k=0;
tt=h/3;
for t=0:0.1:tt; 
    k=k+1;
    yfree_de_ct(:,k)=C*expm(An_d*t)*x_0;
    u_de_ct(:,k)=K_d*yfree_de_ct(:,k);
end

figure(3)
subplot 211
hold on
grid on
plot(0:1:tt,yfree_de_ct(1,1:tt+1),0:1:tt,yfree_de_ct(2,1:tt+1),0:1:tt,yfree_de_ct(3,1:tt+1));
title ('Decentralized Control in Continuous Time')
legend('\deltaTB','\deltaTA','\deltaTC')
xlabel(['Interval time from 0 to ', num2str(tt)])
ylabel('Output free motions for Decentralized Control')

subplot 212
hold on
grid on
plot(0:1:tt,u_de_ct(1,1:tt+1),0:1:tt,u_de_ct(2,1:tt+1),0:1:tt,u_de_ct(3,1:tt+1))
title ('Control action for Decentralized Control, Continuous Time')
legend('\deltaqB','\deltaqA','\deltaqC')
xlabel(['Interval time from 0 to ', num2str(tt)])
ylabel('Static Output-feedback Control Law u(t)')

eigen_CT_decentralized=eig(An_d)
rhoCT_decentralized=max(abs(real(eig(An_d))))


h=300;
i=1;
k=0;
Tfinal=10*h;
steps=[0:Tfinal/h];
for k=steps  
    yfree_de_dt(:,i)=(F+G*KD_d)^(k)*x_0;
    u_de_dt(:,i)=KD_d*yfree_de_dt(:,k+1);
    i=i+1;
end

figure(4)
subplot 211
hold on
grid on
plot(steps*h,yfree_de_dt(1,:),'ko',steps*h,yfree_de_dt(2,:),'k*',steps*h,yfree_de_dt(3,:),'k.');
title ('Decentralized Control in Discrete Time')
legend('\deltaTB','\deltaTA','\deltaTC')
xlabel(['Time instants from 0 to 3000 [s] with sampling time h= ',  num2str(h)])
ylabel('Output free motions for Decentralized Control')

subplot 212
hold on
grid on
plot(steps*h,u_de_dt(1,:),'ko',steps*h,u_de_dt(2,:),'k*',steps*h,u_de_dt(3,:),'k.');
title ('Control action in Decentralized Control, Discrete Time')
legend('\deltaqB','\deltaqA','\deltaqC')
xlabel(['Time instants from 0 to 3000 [s] with sampling time h= ',  num2str(h)])
ylabel('Static Output-feedback Control Law u(kh)')

eigen_DT_decentralized=eig(F+G*KD_d)
rho_DT_decentralized=max(abs(real(eigen_DT_decentralized)))

%% DISTRIBUTED CONTROL STRUCTURES
N=3;
rounding=4;

%Cycle Configuration

ContStrucC=[1 0 1;
            1 1 0;
            0 1 1;];
        
%2a continuous time fixed modes
DiFMct_Cycle=di_fixed_modes(A,Bdec,Cdec,N,ContStrucC,rounding)

%2b discrete time fixed modes
DiFMdt_Cycle=di_fixed_modes(F,Gdec,Hdec,N,ContStrucC,rounding)

%2c stabilizing CT control gain
[K_di_Cycle,rho_di_Cycle,feas_di_Cycle]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC)

%2d stabilizing DT control gain
[KD_di_Cycle,rhoD_di_Cycle,feasD_di_Cycle]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStrucC)

%2e eigenvalues, stability, CL trajectories in CT e DT
An_di_Cycle=A+B*K_di_Cycle

h=300;
k=0;
tt=h/3;
for t=0:0.1:tt; 
    k=k+1;
    yfree_di_ct_Cycle(:,k)=C*expm(An_di_Cycle*t)*x_0;
    u_di_ct_Cycle(:,k)=K_di_Cycle*yfree_di_ct_Cycle(:,k);
end

figure(5)
subplot 211
hold on
grid on
plot(0:1:tt,yfree_di_ct_Cycle(1,1:tt+1),0:1:tt,yfree_di_ct_Cycle(2,1:tt+1),0:1:tt,yfree_di_ct_Cycle(3,1:tt+1));
title ('Distributed Control in Continuous Time-Cycle Configuration')
legend('\deltaTB','\deltaTA','\deltaTC')
xlabel(['Interval time from 0 to ', num2str(tt)])
ylabel('Output free motions for Distributed Control')

subplot 212
hold on
grid on
plot(0:1:tt,u_di_ct_Cycle(1,1:tt+1),0:1:tt,u_di_ct_Cycle(2,1:tt+1),0:1:tt,u_di_ct_Cycle(3,1:tt+1))
title ('Control action for Distributed Control-Cycle Configuration, Continuous Time')
legend('\deltaqB','\deltaqA','\deltaqC')
xlabel(['Interval time from 0 to ', num2str(tt)])
ylabel('Static Output-feedback Control Law u(t)')

eigen_CT_distributed_Cycle=eig(An_di_Cycle)
rhoCT_distributed_Cycle=max(abs(real(eig(An_di_Cycle))))


h=300;
i=1;
k=0;
Tfinal=10*h;
steps=[0:Tfinal/h];
for k=steps      
    i=i+1;
    yfree_di_dt_Cycle(:,k+1)=(F+G*KD_di_Cycle)^(k)*x_0;
    u_di_dt_Cycle(:,k+1)=KD_di_Cycle*yfree_di_dt_Cycle(:,k+1);
end

figure(6)
subplot 211
hold on
grid on
plot(steps*h,yfree_di_dt_Cycle(1,:),'ko',steps*h,yfree_di_dt_Cycle(2,:),'k*',steps*h,yfree_di_dt_Cycle(3,:),'k.');
title ('Distributed Control in Discrete Time-Cycle Configuration')
legend('\deltaTB','\deltaTA','\deltaTC')
xlabel(['Time instants from 0 to 3000 [s] with sampling time h= ',  num2str(h)])
ylabel('Output free motions for Distributed Control')

subplot 212
hold on
grid on
plot(steps*h,u_di_dt_Cycle(1,:),'ko',steps*h,u_di_dt_Cycle(2,:),'k*',steps*h,u_di_dt_Cycle(3,:),'k.');
title ('Control action in Distributed Control-Cycle Configuration, Discrete Time')
legend('\deltaqB','\deltaqA','\deltaqC')
xlabel(['Time instants from 0 to 3000 [s] with sampling time h= ',  num2str(h)])
ylabel('Static Output-feedback Control Law u(kh)')

eigen_DT_distributed_Cycle=eig(F+G*KD_di_Cycle)
rho_DT_distributed_Cycle=max(abs(real(eigen_DT_distributed_Cycle)))

%% Star with bidirectional communication - Centre in A

ContStrucC=[1 1 0;
            1 1 1;
            0 1 1;];
        
%2a continuous time fixed modes
DiFMct_sbi=di_fixed_modes(A,Bdec,Cdec,N,ContStrucC,rounding)

%2b discrete time fixed modes
DiFMdt_sbi=di_fixed_modes(F,Gdec,Hdec,N,ContStrucC,rounding)

%2c stabilizing CT control gain
[K_di_sbi,rho_di_sbi,feas_di_sbi]=LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStrucC)

%2d stabilizing DT control gain
[KD_di_sbi,rhoD_di_sbi,feasD_di_sbi]=LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStrucC)

%2e eigenvalues, stability, CL trajectories in CT e DT
An_di_sbi=A+B*K_di_sbi

h=300;
k=0;
tt=h/3;
for t=0:0.1:tt; 
    k=k+1;
    yfree_di_ct_sbi(:,k)=C*expm(An_di_sbi*t)*x_0;
    u_di_ct_sbi(:,k)=K_di_sbi*yfree_di_ct_sbi(:,k);
end

figure(7)
subplot 211
hold on
grid on
plot(0:1:tt,yfree_di_ct_sbi(1,1:tt+1),0:1:tt,yfree_di_ct_sbi(2,1:tt+1),0:1:tt,yfree_di_ct_sbi(3,1:tt+1));
title ('Distributed Control in Continuous Time-Star with bidirectional communication ')
legend('\deltaTB','\deltaTA','\deltaTC')
xlabel(['Interval time from 0 to ', num2str(tt)])
ylabel('Output free motions for Distributed Control')

subplot 212
hold on
grid on
plot(0:1:tt,u_di_ct_sbi(1,1:tt+1),0:1:tt,u_di_ct_sbi(2,1:tt+1),0:1:tt,u_di_ct_sbi(3,1:tt+1))
title ('Control action for Distributed Control-Star with bidirectional communication , Continuous Time')
legend('\deltaqB','\deltaqA','\deltaqC')
xlabel(['Interval time from 0 to ', num2str(tt)])
ylabel('Static Output-feedback Control Law u(t)')

eigen_CT_distributed_sbi=eig(An_di_sbi)
rhoCT_distributed_sbi=max(abs(real(eig(An_di_sbi))))


h=300;
i=1;
k=0;
Tfinal=10*h;
steps=[0:Tfinal/h];
for k=steps      
    i=i+1;
    yfree_di_dt_sbi(:,k+1)=(F+G*KD_di_sbi)^(k)*x_0;
    u_di_dt_sbi(:,k+1)=KD_di_sbi*yfree_di_dt_sbi(:,k+1);
end

figure(8)
subplot 211
hold on
grid on
plot(steps*h,yfree_di_dt_sbi(1,:),'ko',steps*h,yfree_di_dt_sbi(2,:),'k*',steps*h,yfree_di_dt_sbi(3,:),'k.');
title ('Distributed Control in Discrete Time-Star with bidirectional communication ')
legend('\deltaTB','\deltaTA','\deltaTC')
xlabel(['Time instants from 0 to 3000 [s] with sampling time h= ',  num2str(h)])
ylabel('Output free motions for Distributed Control')

subplot 212
hold on
grid on
plot(steps*h,u_di_dt_sbi(1,:),'ko',steps*h,u_di_dt_sbi(2,:),'k*',steps*h,u_di_dt_sbi(3,:),'k.');
title ('Control action in Distributed Control-Star with bidirectional communication , Discrete Time')
legend('\deltaqB','\deltaqA','\deltaqC')
xlabel(['Time instants from 0 to 3000 [s] with sampling time h= ',  num2str(h)])
ylabel('Static Output-feedback Control Law u(kh)')

eigen_DT_distributed_sbi=eig(F+G*KD_di_sbi)
rho_DT_distributed_sbi=max(abs(real(eigen_DT_distributed_sbi)))
