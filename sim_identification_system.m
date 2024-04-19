clc;
clear all;
close all;
%% create sampling time
dt = 0.01; %second
T = 10; %second
n = T/dt; %number of samples
u(1) = 0;
for i = 2:n
    t(i) = i*dt;
    u(i) = 1;
end
%% Motor parameters
R = 5.2;    %Ohm 
L = 1.1;    %H
Kt = 0.75;  %N.m/A
Ke = 0.75;  %Vs/rad
J = 0.01;   %Kg.m^2
b = 0.02;   %N.m.s

x = zeros(2,n);
A = [-R/L -Ke/L; Kt/J -b/J];
B = [ 1/L; 0];
%% Dynamic system
for i = 2:n
    x(:,i) = (A*x(:,i-1) + B*u(i-1))*dt + x(:, i-1);
end
%% Estimated parameter using Model reference adaptive
% Solve for P using MATLAB's care function
A_bar = [-7.5 2.5; 2.5 -2.5];
Q = [10 0; 0 10]; 

P = lyap(A_bar,Q)
Gamma_a1 = diag([1 1]);
Gamma_b1 = diag([1 1]);
A_hat1(:,:,1) = [0    0;0  0];
B_hat1(:,:,1) = [0; 0];
x_m1 = zeros(2,n);
for i = 2:n
    x_m1(:,i) = (A_hat1(:,:,i-1)*x_m1(:,i-1) + B_hat1(:,i-1)*u(i-1))*dt + x_m1(:,i-1);
    e(:,i) = x(:,i) - x_m1(:,i);
    A_hat1(:,:,i) = (Gamma_a1*P*e(:,i)*x_m1(:,i)')*dt + A_hat1(:,:,i-1);
    B_hat1(:,i) = (Gamma_b1*P*e(:,i)*u(i-1))*dt + B_hat1(:,i-1);
end
A_hat_ss1 = A_hat1(:,:,end);
B_hat_ss1 = B_hat1(:,end);

Gamma_a2 = diag([100 100]);
Gamma_b2 = diag([1 1]);
A_hat2(:,:,1) = [0    0;0  0];
B_hat2(:,:,1) = [0; 0];
x_m2 = zeros(2,n);
for i = 2:n
    x_m2(:,i) = (A_hat2(:,:,i-1)*x_m2(:,i-1) + B_hat2(:,i-1)*u(i-1))*dt + x_m2(:,i-1);
    e(:,i) = x(:,i) - x_m2(:,i);
    A_hat2(:,:,i) = (Gamma_a2*P*e(:,i)*x_m2(:,i)')*dt + A_hat2(:,:,i-1);
    B_hat2(:,i) = (Gamma_b2*P*e(:,i)*u(i-1))*dt + B_hat2(:,i-1);
end
A_hat_ss2 = A_hat2(:,:,end);
B_hat_ss2 = B_hat2(:,end);

Gamma_a3 = diag([1 1]);
Gamma_b3 = diag([100 100]);
A_hat3(:,:,1) = [0    0;0  0];
B_hat3(:,:,1) = [0; 0];
x_m3 = zeros(2,n);
for i = 2:n
    x_m3(:,i) = (A_hat3(:,:,i-1)*x_m3(:,i-1) + B_hat3(:,i-1)*u(i-1))*dt + x_m3(:,i-1);
    e(:,i) = x(:,i) - x_m3(:,i);
    A_hat3(:,:,i) = (Gamma_a3*P*e(:,i)*x_m3(:,i)')*dt + A_hat3(:,:,i-1);
    B_hat3(:,i) = (Gamma_b3*P*e(:,i)*u(i-1))*dt + B_hat3(:,i-1);
end
A_hat_ss3 = A_hat3(:,:,end);
B_hat_ss3 = B_hat3(:,end);

Gamma_a4 = diag([100 100]);
Gamma_b4 = diag([100 100]);
A_hat4(:,:,1) = [0    0;0  0];
B_hat4(:,:,1) = [0; 0];
x_m4 = zeros(2,n);
for i = 2:n
    x_m4(:,i) = (A_hat4(:,:,i-1)*x_m4(:,i-1) + B_hat4(:,i-1)*u(i-1))*dt + x_m4(:,i-1);
    e(:,i) = x(:,i) - x_m4(:,i);
    A_hat4(:,:,i) = (Gamma_a4*P*e(:,i)*x_m4(:,i)')*dt + A_hat4(:,:,i-1);
    B_hat4(:,i) = (Gamma_b4*P*e(:,i)*u(i-1))*dt + B_hat4(:,i-1);
end
A_hat_ss4 = A_hat4(:,:,end);
B_hat_ss4 = B_hat4(:,end);


figure('Position',[100 100 1200 500]);
plot(t,x(1,:),'-g', t, x_m1(1,:),':k', t, x_m2(1,:),'-r', t, x_m3(1,:),'-.m', t, x_m4(1,:),'--b','LineWidth', 1.25);
xlabel('Time - ($s$)', 'interpreter', 'latex', 'fontname','Times New Roman', 'fontsize', 16);
ylabel('Armature Current $i(t)$ - (A)', 'interpreter', 'latex','fontname','Times New Roman', 'fontsize', 16);
legend('Plant', 'Model 1', 'Model 2', 'Model 3', 'Model 4', 'Orientation','Horizontal', 'fontname','Times New Roman', 'fontsize', 16);

figure('Position',[100 100 1200 500]);
plot(t,x(2,:),'-g', t, x_m1(2,:),':k', t, x_m2(2,:),'-r', t, x_m3(2,:),'-.m', t, x_m4(2,:),'--b','LineWidth', 1.25);
xlabel('Time - ($s$)', 'interpreter', 'latex', 'fontname','Times New Roman', 'fontsize', 16);
ylabel('Angular Velocity ${\theta}\dot(t)$ - (rad/s)', 'interpreter', 'latex','fontname','Times New Roman', 'fontsize', 16);
legend('Plant', 'Model 1', 'Model 2', 'Model 3', 'Model 4', 'Orientation','Horizontal', 'fontname','Times New Roman', 'fontsize', 16);

u_t(1) = 0;
for i = 2:n
    t(i) = i*dt;
    if t(i) < 2
        u_t(i) = 1.5; 
    elseif t(i) < 4
        u_t(i) = 5; 
    elseif t(i) < 6
        u_t(i) = 2;
    elseif t(i) < 8
        u_t(i) = 4;
    else
        u_t(i) = 1;
    end
end
x_0 = zeros(2, n);
x_1 = zeros(2, n);
x_2 = zeros(2, n);
x_3 = zeros(2, n);
x_4 = zeros(2, n);
for i = 2:n
    x_0(:,i) = (A*x_0(:,i-1) + B*u_t(i-1))*dt + x_0(:, i-1);
    x_1(:,i) = (A_hat_ss1*x_1(:,i-1) + B_hat_ss1*u_t(i-1))*dt + x_1(:, i-1);
    x_2(:,i) = (A_hat_ss2*x_2(:,i-1) + B_hat_ss2*u_t(i-1))*dt + x_2(:, i-1);
    x_3(:,i) = (A_hat_ss3*x_3(:,i-1) + B_hat_ss3*u_t(i-1))*dt + x_3(:, i-1);
    x_4(:,i) = (A_hat_ss4*x_4(:,i-1) + B_hat_ss4*u_t(i-1))*dt + x_4(:, i-1);
end

figure('Position',[100 100 1200 500]);
plot(t, x_0(1,:), '-g', t, x_1(1,:),':k', t, x_2(1,:),'-r', t, x_3(1,:),'-.m', t, x_4(1,:),'--b','LineWidth', 1.25);
xlabel('Time - ($s$)', 'interpreter', 'latex', 'fontname','Times New Roman', 'fontsize', 16);
ylabel('Armature Current $i(t)$ - (A)', 'interpreter', 'latex','fontname','Times New Roman', 'fontsize', 16);
legend('Plant','Model 1', 'Model 2', 'Model 3', 'Model 4', 'Orientation','Horizontal', 'fontname','Times New Roman', 'fontsize', 16);

figure('Position',[100 100 1200 500]);
plot(t,x_0(2,:),'-g', t, x_1(2,:),':k', t, x_2(2,:),'-r', t, x_3(2,:),'-.m', t, x_4(2,:),'--b','LineWidth', 1.25);
xlabel('Time - ($s$)', 'interpreter', 'latex', 'fontname','Times New Roman', 'fontsize', 16);
ylabel('Angular Velocity ${\theta}\dot(t)$ - (rad/s)', 'interpreter', 'latex','fontname','Times New Roman', 'fontsize', 16);
legend('Plant', 'Model 1', 'Model 2', 'Model 3', 'Model 4', 'Orientation','Horizontal', 'fontname','Times New Roman', 'fontsize', 16);


