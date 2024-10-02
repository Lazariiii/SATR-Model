clc
clear

%Parameter
theta = 0.8;
alpha = 0.072;
mu = 0.9;
b1 = 0.002;
b2 = 0.4;
%gamma = 0.25; 
delta = 0.02; 
sigma = 0.009;

%Variabel
%S = 5;
%A = 2;
%T = 0.5;
%R = 0.2;

%Fungsi
syms wS wA wT wR gamma;
syms S A T R;
dS = theta - (alpha*S*A) - (wS*S) +(delta*T) + (sigma*R);
dA = (alpha*S*A) - (mu*A) - (b1*A) - (wA*A);
dT = (mu*A) - (b2*T) - (delta*T) - (wT*T);
dR = (b1*A) + (b2*T) - (sigma*R) - (wR*T) - (gamma*R);

%Matrik A
A = [diff(dS,S) diff(dS,A) diff(dS,T) diff(dS,R);diff(dA,S) diff(dA,A) diff(dA,T) diff(dA,R);diff(dT,S) diff(dT,A) diff(dT,T) diff(dT,R);diff(dR,S) diff(dR,A) diff(dR,T) diff(dR,R)];
disp('Nilai Matrik A adalah: ');
gamma_val = 0.25;
wS_val = 0.21;
wA_val = 0.2;
wT_val = 0.2;
wR_val = 0.02;
a = subs(A,[wS wA wT wR gamma],[wS_val wA_val wT_val wR_val gamma_val]);
format short;
disp(a);

%Matrik B
B = [diff(dS,wS) diff(dS,gamma);diff(dA,wA) diff(dA,gamma);diff(dT,wT) diff(dT,gamma);diff(dR,wR) diff(dR,gamma)];
disp('Nilai Matrik B adalah:');
disp(B);

%Nilai Eigen dan Vektor Eigen
S_val = 5;
A_val = 2;
a1 = subs(a,[S A],[S_val A_val] );
[V, N] = eig(a1);
disp('Nilai Eigen:');
disp(diag(N));
disp('Vektor Eigen:');
disp(V);
