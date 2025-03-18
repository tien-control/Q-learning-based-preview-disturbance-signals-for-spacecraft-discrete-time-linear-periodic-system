%% LQR co Review control 00 va LQR khong co Review control va LMI 00
clc, clear, close all;

%% Input
% Nhap AA = [A(1); A(2);...; A(N)]
N = 2;
AA = [-0.6 0.4; -1.2 1.2; -0.8 0.4; -1.3 1.2];
BB = [1.1 ; 0.6; 1.3; 0.7];
CC = [-0.7 0.9; -0.6 0.9];
DD = [0.3; 0.2; 0.3; 0.2];
EE = [0; 0];

% Tinh kich co ma tran
n = size(AA,2);
q = size(BB,2);
l = size(DD,2);

% Nhap tham so dieu chinh
t_end = 80;
t = 1:t_end;
Mr = 2;
Md = 2;
e_ts = 1;
e_ts_kRC = 1;

%% Khoi tao ma tran
% Bien trang thai thuc
x_k = zeros(n,1);
u_k = zeros(q,1);
y_k = zeros(q,1);
v_k = zeros(q,1);
dx_k = zeros(n,1);
du_k = zeros(q,1);
dy_k = zeros(q,1);
dv_k = zeros(q,1);

% Ma tran he thong thuc
A = cell(1,N);
B = cell(1,N);
C = cell(1,N);
D = cell(1,N);
E = cell(1,N);

% Ma tran he thong trung gian
Ar = zeros((Mr+N+1)*q);
Aw = zeros((Md+N+1)*l);

An_k = zeros(n+(Mr+N+1)*q+(Md+N+1)*l, n+(Mr+N+1)*q+(l+Md+N));
Bn_k = zeros(n+(Mr+N+1)*q+(Md+N+1)*l,q);
Gxw_k = zeros(1, l+Md+N);
Cn_k = zeros(size(CC,1)/N, size(CC,2)+size(CC,1)/N+Mr+N+size(EE,1)+Md+N);
Cn = cell(1,N);
G1 = zeros(size(CC,1)/N,size(CC,1)/N+Mr+N);
G2_k = zeros(size(EE,1)/N,size(EE,1)/N+Md+N);

% Ma tran he thong moi
Am_k = zeros(n+(Mr+N+1)*q+(Md+N+1)*l+size(CC,1)/N,n+(Mr+N+1)*q+(Md+N+1)*l+l);
Bm_k = zeros(n+(Mr+N+1)*q+(Md+N+1)*l+size(CC,1)/N,q);
Cz_k = zeros(size(CC,1)+(Mr+N+1)*q+(Md+N+1)*l+q,size(CC,2)+(Mr+N+1)*q+(Md+N+1)*l+q);

Am = cell(1,N);
Bm = cell(1,N);
Cz = cell(1,N);

%% Tinh toan ma tran
for o = 1:(Mr+N+1)*q-1
    Ar(o,o+1) = 1;
end

for o = 1:(Md+N+1)*l-1
    Aw(o,o+1) = 1;
end

for o = 1:N
    A{o} = AA(size(AA,1)/N*(o-1)+1:size(AA,1)/N*o,:);
    B{o} = BB(size(BB,1)/N*(o-1)+1:size(BB,1)/N*o,:);
    C{o} = CC(size(CC,1)/N*(o-1)+1:size(CC,1)/N*o,:);
    D{o} = DD(size(DD,1)/N*(o-1)+1:size(DD,1)/N*o,:);
    E{o} = EE(size(EE,1)/N*(o-1)+1:size(EE,1)/N*o,:);


    Gxw_k = [D{o} zeros(size(D{o},1),Md+N)];
    An_k(1:size(A{o},1),1:size(A{o},2)) = A{o};
    An_k(size(A{o},1)+1:size(A{o},1)+(Mr+N+1)*q,size(A{o},2)+1:size(A{o},2)+(Mr+N+1)*q) = Ar;
    An_k(1:size(Gxw_k,1),end-size(Gxw_k,2)+1:end) = Gxw_k;
    An_k(size(A{o},1)+(Mr+N+1)*q+1:end,size(A{o},2)+(Mr+N+1)*q+1:end) = Aw;
    Bn_k(1:size(B{o},1),1:q) = B{o};
    
    G1(1:size(CC,1)/N,1:size(CC,1)/N) = -eye(size(CC,1)/N);
    G2_k(1:size(E{o},1),1:size(E{o},2)) = E{o};
    Cn_k = [C{o} G1 G2_k];

    Am_k(1:size(An_k,1),1:size(An_k,2)) = An_k;
    Am_k(size(An_k,2)+1:end,1:size(Cn_k,2)) = Cn_k; 
    Am_k(end-size(Cn_k,1)+1:end, end-size(Cn_k,1)+1:end) = eye(size(Cn_k,1));
    
    Bm_k(1:size(Bn_k,1),1:size(Bn_k,2)) = Bn_k;
    Cz_k(1:size(C{o},1),1:size(C{o},2)) = C{o};
    Cz_k(size(C{o},1)+1:size(C{o},1)+(Mr+N+1)*q,size(C{o},2)+1:size(C{o},2)+(Mr+N+1)*q) = eye((Mr+N+1)*q);
    Cz_k(1:size(Gxw_k,1),size(C{o},2)+(Mr+N+1)*q+1:size(C{o},2)+(Mr+N+1)*q+size(Gxw_k,2)) = Gxw_k;
    Cz_k(size(C{o},1)+(Mr+N+1)*q+1:size(C{o},1)+(Mr+N+1)*q+(Md+N+1)*l,size(C{o},2)+(Mr+N+1)*q+1:size(C{o},2)+(Mr+N+1)*q+(Md+N+1)*l) = eye((Md+N+1)*l);
    Cz_k(size(C{o},1)+(Mr+N+1)*q+(Md+N+1)*l+1:size(C{o},1)+(Mr+N+1)*q+(Md+N+1)*l+q,size(C{o},2)+(Mr+N+1)*q+(Md+N+1)*l+1:size(C{o},2)+(Mr+N+1)*q+(Md+N+1)*l+q) = eye(q);

    Cn{o} = Cn_k;
    Am{o} = Am_k;
    Bm{o} = Bm_k;
    Cz{o} = Cz_k;
end

%% Mo phong tin hieu dat va nhieu
for i = 1:length(t)
    % Tin hieu dat
        if i >= 30
            r(i) = 2;
        else
            r(i) = 0;
        end
    % Tin hieu nhieu
        if i >= 30 && i<= 60
            w(i) = sin(1/(30*pi*i));
        else
            w(i) = 0;
        end
    % Tinh delta r, delta w
        if i <= N
            dr = 0;
            dw = 0;
        else
            dr(i) = r(i) - r(i-N);
            dw(i) = w(i) - w(i-N);
        end
end


%% Mo phong LMI Mr = 0, Md = 0
%% Mo phong LMI 
% Bo dieu khien LMI
K_00 = [1.00471 0.33084; 1.20883 0.36974];
for i = 1:N
    K_LMI_00{i} = K_00(i,:);
    Ky(i) = K_LMI_00{i}(1:size(y_k,1));
    Kv(i) = K_LMI_00{i}(size(y_k,1)+1:end);
end

% He thong chay trong thoi gian
for i = 1:length(t)
    % Thu tu trong chu ky
    i_N = mod(i,N)+1;

    % Tinh tin hieu dieu khien
    uMr = 0; uMd = 0;
    
    du_k = Ky(i_N)*dy_k + uMr + uMd + Kv(i_N)*dv_k;
    
        if i-N < 1
            u_k = du_k;
        else
            u_k = du_k + uu_LMI_00(:,i-N);   
        end

    % Trang thai
    x_k = A{i_N}*x_k + B{i_N}*u_k + D{i_N}*w(i);
    
        if i-N < 1
            dx_k = x_k;
        else
            dx_k = x_k - xx_LMI_00(:,i-N);   
        end

    % Cap nhat
    dy_k = C{i_N}*dx_k + E{i_N}*dw(i);

        if i-N < 1
            y_k = dy_k;
        else
            y_k = dy_k + yy_LMI_00(:,i-N);  
        end

    e_k = e_ts*(y_k - r(i));
    v_k = v_k + e_k ;
    
        if i-N < 1
            dv_k = v_k;
        else
            dv_k = v_k - vv_LMI_00(:,i-N);  
        end
    
    % Luu trang thai theo thoi gian
    xx_LMI_00(:,i) = x_k;
    uu_LMI_00(:,i)  = u_k;
    yy_LMI_00(:,i)  = y_k;
    vv_LMI_00(:,i)  = v_k;
    dv_LMI_00(:,i)  = dv_k;
end

%% Mo phong LQR
% Cap nhat lai cac bien
x_k = zeros(n,1);
u_k = zeros(q,1);
y_k = zeros(q,1);
v_k = zeros(q,1);
dx_k = zeros(n,1);
du_k = zeros(q,1);
dy_k = zeros(q,1);
dv_k = zeros(q,1);

xx_LQR_00 = [];
uu_LQR_00 = [];
yy_LQR_00 = [];
vv_LQR_00 = [];

% Chuyen ve he thong LTI
[m1,m2] = size(Bm{1});
Q0 = eye(m1);
R0 = eye(m2);
M = eye(m1);
F = []; Q = []; R = [];

    for o = 1:N
        Q = blkdiag(Q,Q0);
        R = blkdiag(R,R0);
    end
    
    for o = 1:N
        M=Am{o}*M;
        F=[F;M];
    end
    
AM = [zeros(N*m1,(N-1)*m1),F];
CM = [];BM = [];
    
    %Xu ly tung cot mot roi ghep lai
    for H = 1:N
        M = eye(m1);
        CM=[];
        for V = 1:N
            if V < H
                CM = [CM;zeros(m1,m2)];
            elseif V == H
                C0 = Bm{V};
                CM = [CM;C0];
            else
                M = Am{V}*M;
                CM = [CM;M*C0];
            end
        end
        BM=[BM,CM];
    end

% Bo dieu khien
[L,P,~]=dlqr(AM,BM,Q,R);
for i = 1:N
    K_LQR_00{i} = -L(i,size(L,2)/N*(N-1)+1:end);
    Kx{i} = K_LQR_00{i}(1:size(x_k,1));
    Kv(i) = K_LQR_00{i}(size(x_k,1)+Mr+N+Md+N+3:end);
end

% He thong chay trong thoi gian
for i = 1:length(t)
    % Thu tu trong chu ky
    i_N = mod(i,N)+1;

    % Tinh toan tin hieu dieu khien
    uMr = 0; uMd = 0;
    
    du_k = Kx{i_N}*dx_k + uMr + uMd + Kv(i_N)*dv_k;
    
        if i-N <= 1
            u_k = du_k;
        else
            u_k = du_k + uu_LQR_00(:,i-N);
        end
    
    % Trang thai
    x_k = A{i_N}*x_k + B{i_N}*u_k + D{i_N}*w(i);
    
        if i-N < 1
            dx_k = x_k;
        else
            dx_k = x_k - xx_LQR_00(:,i-N);
        end
    
    % Cap nhat
    dy_k = C{i_N}*dx_k + E{i_N}*dw(i);

        if i-N < 1
            y_k = dy_k;
        else
            y_k = dy_k + yy_LQR_00(:,i-N);  
        end

    % sai so
    e_k = e_ts*(y_k - r(i));
    v_k = v_k + e_k;
    
        if i-N < 1
            dv_k = v_k;
        else
            dv_k = v_k - vv_LQR_00(:,i-N);
        end
    
    % Luu trang thai theo thoi gian
    xx_LQR_00(:,i)  = x_k;
    uu_LQR_00(:,i)  = u_k;
    yy_LQR_00(:,i)  = y_k;
    vv_LQR_00(:,i)  = v_k;
    dv_LQR_00(:,i)  = dv_k;
end

%% Mo phong LQR co Review Control
% Cap nhat lai cac bien
x_k = zeros(n,1);
u_k = zeros(q,1);
y_k = zeros(q,1);
v_k = zeros(q,1);
dx_k = zeros(n,1);
du_k = zeros(q,1);
dy_k = zeros(q,1);
dv_k = zeros(q,1);

xx_LQR_cRC = [];
uu_LQR_cRC = [];
yy_LQR_cRC = [];
vv_LQR_cRC = [];

% Chuyen ve he thong LTI
[m1,m2] = size(Bm{1});
Q0 = eye(m1);
R0 = eye(m2);
M = eye(m1);
F = []; Q = []; R = [];

    for o = 1:N
        Q = blkdiag(Q,Q0);
        R = blkdiag(R,R0);
    end
    
    for o = 1:N
        M=Am{o}*M;
        F=[F;M];
    end
    
AM = [zeros(N*m1,(N-1)*m1),F];
CM = [];BM = [];
    
    %Xu ly tung cot mot roi ghep lai
    for H = 1:N
        M = eye(m1);
        CM=[];
        for V = 1:N
            if V < H
                CM = [CM;zeros(m1,m2)];
            elseif V == H
                C0 = Bm{V};
                CM = [CM;C0];
            else
                M = Am{V}*M;
                CM = [CM;M*C0];
            end
        end
        BM=[BM,CM];
    end

% Bo dieu khien
[L1,P1,~] = dlqr(AM,BM,Q,R);
for i = 1:N
    K_LQR_cRC{i} = -L1(i,size(L1,2)/N*(N-1)+1:end);
    Kx{i} = K_LQR_cRC{i}(1:size(x_k,1));
    Kv(i) = K_LQR_cRC{i}(size(x_k,1)+Mr+N+Md+N+3:end);
end

% He thong chay trong thoi gian
for i = 1:length(t)
    % Thu tu trong chu ky
    i_N = mod(i,N)+1;

    % Tinh toan tin hieu dieu khien
    uMr = 0; uMd = 0;
    
    du_k = Kx{i_N}*dx_k + uMr + uMd + Kv(i_N)*dv_k;
    
        if i-N < 1
            u_k = du_k;
        else
            u_k = du_k + uu_LQR_cRC(:,i-N);
        end
    
    % Trang thai
    x_k = A{i_N}*x_k + B{i_N}*u_k + D{i_N}*w(i);
    
        if i-N < 1
            dx_k = x_k;
        else
            dx_k = x_k - xx_LQR_cRC(:,i-N);
        end
    
    % Cap nhat
    dy_k = C{i_N}*dx_k + E{i_N}*dw(i);

        if i-N < 1
            y_k = dy_k;
        else
            y_k = dy_k + yy_LQR_cRC(:,i-N);  
        end

    % sai so
    e_k = e_ts*(y_k - r(i));
    v_k = v_k + e_k;
    
        if i-N < 1
            dv_k = v_k;
        else
            dv_k = v_k - vv_LQR_cRC(:,i-N);
        end
    
    % Luu trang thai theo thoi gian
    xx_LQR_cRC(:,i)  = x_k;
    uu_LQR_cRC(:,i)  = u_k;
    yy_LQR_cRC(:,i)  = y_k;
    vv_LQR_cRC(:,i)  = v_k;
    dv_LQR_cRC(:,i)  = dv_k;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mo phong LQR khong co Review Control
% Cap nhat lai cac bien
Mr = 0;
Md = 0;

x_k = zeros(n,1);
u_k = zeros(q,1);
y_k = zeros(q,1);
v_k = zeros(q,1);
dx_k = zeros(n,1);
du_k = zeros(q,1);
dy_k = zeros(q,1);
dv_k = zeros(q,1);

xx_LQR_kRC = [];
uu_LQR_kRC = [];
yy_LQR_kRC = [];
vv_LQR_kRC = [];

% Chuyen ve he thong LTI
[m1,m2] = size(B{1});
Q0 = eye(m1);
R0 = eye(m2);
M = eye(m1);
F = []; Q = []; R = [];

    for o = 1:N
        Q = blkdiag(Q,Q0);
        R = blkdiag(R,R0);
    end
    
    for o = 1:N
        M=A{o}*M;
        F=[F;M];
    end
    
AM = [zeros(N*m1,(N-1)*m1),F];
CM = [];BM = [];
    
    %Xu ly tung cot mot roi ghep lai
    for H = 1:N
        M = eye(m1);
        CM=[];
        for V = 1:N
            if V < H
                CM = [CM;zeros(m1,m2)];
            elseif V == H
                C0 = B{V};
                CM = [CM;C0];
            else
                M = A{V}*M;
                CM = [CM;M*C0];
            end
        end
        BM=[BM,CM];
    end

% Bo dieu khien
[L2,P,~] = dlqr(AM,BM,Q,R);
K_LQR_cRC = -L2(:,size(x_k,1)+1:end);

% He thong chay trong thoi gian
for i = 1:length(t)
    % Thu tu trong chu ky
    i_N = mod(i,N)+1;

    % Tinh toan tin hieu dieu khien
    u_k = K_LQR_cRC(i_N,:)*x_k + v_k;
    
    % Cap nhat he thong
    x_k = A{i_N}*x_k + B{i_N}*u_k + D{i_N}*w(i);
    y_k = C{i_N}*x_k + E{i_N}*u_k;
    e_k = e_ts_kRC*(y_k - r(i));
    v_k = v_k + e_k;
    
    % Luu trang thai theo thoi gian
    xx_LQR_kRC(:,i)  = x_k;
    uu_LQR_kRC(:,i)  = u_k;
    yy_LQR_kRC(:,i)  = y_k;
    vv_LQR_kRC(:,i)  = v_k;
end

%% Ve do thi
figure


subplot(2, 1, 1);
hold on
plot(t,r,'Color', [0, 0, 1],'LineStyle', '-', 'LineWidth', 1);
plot(t,yy_LMI_00,'Color', [1, 0, 0],'LineStyle', '-', 'LineWidth', 1);
plot(t,yy_LQR_cRC,'Color', [0, 0.8, 0],'LineStyle', '-', 'LineWidth', 1);
plot(t,yy_LQR_kRC,'Color', [1, 0.7, 0],'LineStyle', '-', 'LineWidth', 1);
legend({'$r$','$LMIPS$','$OCPS$','$OCNPS$'}, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 8);
ylabel('y(k), r(k)', 'FontSize', 10, 'Interpreter', 'latex');
hold off

subplot(2, 1, 2);
hold on
plot(t,uu_LMI_00,'Color', [1, 0, 0],'LineStyle', '-', 'LineWidth', 1);
plot(t,uu_LQR_cRC,'Color', [0, 0.8, 0],'LineStyle', '-', 'LineWidth', 1);
plot(t,uu_LQR_kRC,'Color', [1, 0.7, 0],'LineStyle', '-', 'LineWidth', 1);
hold off

%{
hold on
plot(t,r,'Color', [0, 0, 1],'LineStyle', '-', 'LineWidth', 1);
plot(t,yy_LMI_00,'Color', [1, 0, 0],'LineStyle', '-', 'LineWidth', 1);
plot(t,yy_LQR_cRC,'Color', [0, 0.8, 0],'LineStyle', '-', 'LineWidth', 1);
plot(t,yy_LQR_kRC,'Color', [1, 0.7, 0],'LineStyle', '-', 'LineWidth', 1);
legend({'$r$','$LMIPS$','$OCPS$','$OCNPS$'}, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 8);
ylabel('y(k), r(k)', 'FontSize', 10, 'Interpreter', 'latex');
hold off
%}