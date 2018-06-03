clear all
close all

L = 1e-3; % Vs/A
C = 1e-6; % As/V
R = 10;    % V/A

s = tf('s');

figure(1)
F = 1/(1+s*R*C+s^2*L*C);
step(F), grid, hold on;

figure(2);
bode(F), grid;

% Zustandsraum
A = [[-R/L,-1/L];[1/C,0]];
%b = [1/L,0]';
b = [1/L;0];
c = [0,1];
d = 0;

figure(1);
step(ss(A,b,c,d),'.r');

% Regelungsnormalform
An = [[0,1];[-1/L/C,-R/L]];
bn = [0,1/L/C]';
cn = [1,0];
dn = 0;

step(ss(An,bn,cn,dn),'og');

% Jordan Normalform
[M,eigVal] = eig(A);
lamda = diag(eigVal);

T = M^-1;
A_ = T*A*M;
b_ = T*b;
c_ = c*M;

tend = 1.4e-3;
t = 0:1.4e-3/100:1.4e-3;
u0 = 1;
x1 = b_(1)*u0/-lamda(1)*(1-exp(lamda(1)*t));
x2 = b_(2)*u0/-lamda(2)*(1-exp(lamda(2)*t));
y = c_(1)*x1+c_(2)*x2;
%fprintf('y = %1.4fV\n',y);
figure(1);
plot(t,real(y),'sk');

% Diskretisierung
f0 = sqrt(real(lamda(1))^2+imag(lamda(1))^2)/2/pi;
T0 = 1/f0;
Ts = T0/10;

Phi =   eye(size(A));
gamma = eye(size(A));
for i = 1:20
    Phi = Phi + A^i*Ts^i/factorial(i);
    gamma = gamma + A^i*Ts^i/factorial(i+1);
end;
gamma = Ts*gamma * b;

RLC_d = c2d(ss(A,b,c,d),Ts);
% Phi = RLC_d.a;
% gamma = RLC_d.b;

% Darstellung der Sprungantwort
x = [0,0]';
t = 0:Ts:tend;
y_vek = zeros(size(t));

for j=1:length(t)
    y_vek(j) = c*x;
    x = Phi*x+gamma*u0;
end;
figure(1);
stairs(t,y_vek,'m');





