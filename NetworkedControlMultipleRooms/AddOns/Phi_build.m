function Phi=Phi_build(A,B,C,K,h,tau,theta)

% h: sampling time
% tau: delay
% theta=1 if transmission occurs
% theta=0 if packet is lost

n=size(A,2);
p=size(C,1);
Gamma_tau=A\(expm(A*tau)-eye(n));
Gamma_h_m_tau=A\(expm(A*(h-tau))-eye(n));
Phi=[expm(A*h)+theta*Gamma_h_m_tau*B*K*C, expm(A*(h-tau))*Gamma_tau*B*K+(1-theta)*Gamma_h_m_tau*B*K 
    theta*C (1-theta)*eye(p)];