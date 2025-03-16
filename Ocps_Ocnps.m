%% Optimal control - no preview signals and preview signals
clc; clear; %close all;
%% Simple system
%% Simple system
N = 4;
% (1)
A{1} = [-0.6 0.4; -1.2 1.2]; B{1} = [1.1; 0.6]; D{1} = [1; 1];
Q{1} = eye(2); R{1} = eye(1);
% (2)
A{2} = [-0.8 0.4; -1.2 -0.2]; B{2} = [1.3; 0.7]; D{2} = [-1; 1];
Q{2} = eye(2); R{2} = eye(1);
% (3)
A{3} = [0.5 -0.7; 0.3 1.2]; B{3} = [0.5; -0.5]; D{3} = [1; -1];
Q{3} = eye(2); R{3} = eye(1);
% (4)
A{4} = [-0.6 1.4; -1.2 -0.4]; B{4} = [1.1; 0.7]; D{4} = [-1; 1];
Q{4} = eye(2); R{4} = eye(1);
% x_0
x_0 = [0.5 0.5]';
%% Parameter
a = 5; t_end = 100;
n = size(A{1},2); m = size(B{1},2); l = size(D{1},2);
%% Signals
% Disturbance
for i = 1:t_end
        if i >= 30 && i<= 70    d{i} = 0.07*sin((1/15)*pi*i);
        else                    d{i} = 0;
        end
end
%% Optimal control with no preview signals
% Transform system
xn_0 = [zeros(n,1); zeros(n,1); x_0];
An = []; Bn = []; Qn = []; Rn = [];
M0 = eye(n); %(M0,M1,M2) are medium matrix
% Process follow column
for i = 1:N
    M0 = A{i}*M0;
    M1 = eye(n); M2 = [];
    for j = 1:N
        if j < i
            M2 = [M2; zeros(n,m)];
        elseif j == i
            M2 = [M2; B{i}];
        else
            M1 = A{j}*M1;
            M2 = [M2; M1*B{i}];
        end
    end
    An = [An; M0]; Bn = [Bn, M2];
    Qn = blkdiag(Qn,Q{i}); Rn = blkdiag(Rn,R{i});
end
An = [zeros(N*n,(N-1)*n),An];
% Optimal controller
[K,P,~] = dlqr(An,Bn,Qn,Rn);
% Simulation
x_a{1} = x_0;
for i = 1:t_end-1
    k = mod(i,N)+1; % State matrix
    % Control signals
    u_a{i} = -K(k,size(K,2)*(N-1)/N+1:end)*x_a{i};
    % State value
    x_a{i+1} = A{k}*x_a{i} + B{k}*u_a{i} + D{k}*d{i};
end
%% Optimal control with preview signals
Psin = []; Lambdan = []; Qsn = []; Rsn = [];
Ad = zeros((a+N+1)*l); Ad(1:end-l,l+1:end) = eye((a+N)*l);
% Trans1
for i = 1:N
    G{i} = [D{i} zeros(n,(a+N)*l)];
    Psi{i} = [A{i} G{i}; zeros((a+N+1)*l,n) Ad];
    Lambda{i} = [B{i}; zeros((a+N+1)*l,m)];
    Qs{i} = blkdiag(Q{i}, zeros((a+N+1)*l));
    Rs{i} = R{i};
end
% Trans2
M0 = eye(n+(a+N+1)*l);
for i = 1:N
    M0 = Psi{i}*M0;
    M1 = eye(n+(a+N+1)*l); M2 = [];
    for j = 1:N
        if j < i
            M2 = [M2; zeros(n+(a+N+1)*l,m)];
        elseif j == i
            M2 = [M2; Lambda{i}];
        else
            M1 = Psi{j}*M1;
            M2 = [M2; M1*Lambda{i}];
        end
    end
    Psin = [Psin; M0]; Lambdan = [Lambdan, M2];
    Qsn = blkdiag(Qsn,Qs{i}); Rsn = blkdiag(Rsn,Rs{i});
end
Psin = [zeros(N*(n+(a+N+1)*l),(N-1)*(n+(a+N+1)*l)),Psin];
% Optimal controller
[K,P,~] = dlqr(Psin,Lambdan,Qsn,Rsn);
K = K(:,size(K,2)*(N-1)/N+1:end);
% Simulation
x_b{1} = x_0;
for i = 1:t_end-1
    k = mod(i,N)+1; % State matrix
    % Find dot control signal
    ud = zeros(l,1);
    for j = 0:a
        if i+j <= t_end
                ud = ud + K(k,n+j*l+1:n+(j+1)*l)*d{i+j};
        else    ud = ud;
        end
    end
    u_b{i} = - (K(k,1:n)*x_b{i} + ud);
    % State value
    x_b{i+1} = A{k}*x_b{i} + B{k}*u_b{i} + D{k}*d{i};
end

figure
subplot(1, 2, 1);
xxa = cell2mat(x_a);
hold on
plot(0:size(xxa,2)-1,xxa(1,:))
plot(0:size(xxa,2)-1,xxa(2,:))
%xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
legend('$x_1$','$x_2$', 'Interpreter', 'latex', 'FontSize', 10);

subplot(1, 2, 2);
xxb = cell2mat(x_b);
hold on
plot(0:size(xxb,2)-1,xxb(1,:))
plot(0:size(xxb,2)-1,xxb(2,:))
%xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
legend('$x_1$','$x_2$', 'Interpreter', 'latex', 'FontSize', 10);

%{
figure(2)
subplot(1, 2, 1);
uua = cell2mat(u_a);
plot(0:size(uua,2)-1,uua)
subplot(1, 2, 2);
uub = cell2mat(u_b);
plot(0:size(uub,2)-1,uub)
%}
