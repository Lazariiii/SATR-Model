clc
clear

% Model parameters
theta = 0.8;
alpha = 0.072;
mu = 0.9;
b1 = 0.002;
b2 = 0.4;% rate of infection
gamma = 0.25; % rate of recovery (try also 0.07)
delta = 0.02; % rate of immunity loss
sigma = 0.009;
wS = 0.21;
wA = 0.2;
wT = 0.2;
wR = 0.02;
%N = 6*10^7; % Total population N = S + A + T + R
%I = 10^4; % initial number of infected
p = 360; % period of 365 days
dt = 6; % time interval of 6 hours (1/4 of a day)
fprintf('Value of parameter R0 is %.2f',(alpha*theta)/(wA*(mu+b1+wA)));
% Calculate the model
alpha = 1/5;
[S,A,T,R] = satr_model(theta,alpha,mu,b1,b2,delta,sigma,wS,wA,wT,wR,gamma,p,dt);
% Plots that display the epidemic outbreak
tt = 0:dt:(p-dt);
% Curve
plot(tt,S,'r',tt,A,'b',tt,T,'g',tt,R,'y','LineWidth',2); grid on;
xlabel('Days'); ylabel('Population');
legend('S','A','T','R');
% Map
%plot(A(1:(p/dt)-1),A(2:p/dt),"LineWidth",1,"Color",'r');
%hold on; grid on;
%plot(A(2),A(1),'ob','MarkerSize',4);
%xlabel('Infected at time t'); ylabel('Infected at time t+1');
%hold on;
%delta = 1/60;
%[S,A,T,R] = satr_model(theta,alpha,mu,b1,b2,delta,sigma,wS,wA,wT,wR,gamma,p,dt);
%plot(tt,S,'r',tt,A,'b',tt,T,'g',tt,R,'y','LineWidth',2); grid on;
%xlabel('Days'); ylabel('Population');
%legend('S','A','T','R');

function [S,A,T,R] = satr_model(theta,alpha,mu,b1,b2,delta,sigma,wS,wA,wT,wR,gamma,p,dt)
    S = zeros(1,p/dt);
    S(1) = 5;
    A = zeros(1,p/dt); 
    A(1) = 2;
    T = zeros(1,p/dt);
    T(1) = 0.5;
    R = zeros(1,p/dt);
    R(1) = 0.2;

    for tt = 1:(p/dt)-1
        % Equations of the model
        dS = theta - (alpha*(S(tt))*(A(tt))) - (wS*S(tt)) +(delta*(T(tt))) + (sigma*(R(tt)));
        dA = (alpha*(S(tt))*(A(tt))) - (mu*(A(tt))) - (b1*(A(tt))) - (wA*(A(tt)));
        dT = (mu*(A(tt))) - (b2*(T(tt))) - (delta*(T(tt))) - (wT*(T(tt)));
        dR = (b1*(A(tt))) + (b2*(T(tt))) - (sigma*(R(tt))) - (wR*(T(tt))) - (gamma*(R(tt)));
        S(tt+1) = S(tt) + dS;
        A(tt+1) = A(tt) + dA;
        T(tt+1) = T(tt) + dT;
        R(tt+1) = R(tt) + dR;
    end
end