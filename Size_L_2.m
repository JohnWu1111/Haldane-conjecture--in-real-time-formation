clear;
tic

L = 2;
len1 = 2^L;
len2 = 3^L;
h0 = 0.3;
sigma_z_1 = [ 1 0; 0 -1];
sigma_x_1 = [ 0 1; 1 0];
sigma_z_2 = [ 1 0 0; 0 0 0; 0 0 -1];
sigma_x_2 = [ 0 1 0; 1 0 1; 0 1 0].*sqrt(2);
I1 = eye(2);
I2 = eye(3);
h = h0.*[1 -1];
dt = 0.01;
t = 0:dt:200;

H1 = -2*kron(sigma_z_1,sigma_z_1) + h(1).*kron(sigma_x_1,I1) + h(2).*kron(I1,sigma_x_1);
H2 = -2*kron(sigma_z_2,sigma_z_2) + h(1).*kron(sigma_x_2,I2) + h(2).*kron(I2,sigma_x_2);
[V1,D1] = eig(H1);
[V2,D2] = eig(H2);

phi0_1 = zeros(len1,1);
phi1_1 = zeros(len1,1);
phi0_2 = zeros(len2,1);
phi1_2 = zeros(len2,1);
phi0_1(1) = 1;
phi1_1(end) = 1;
phi0_2(1) = 1;
phi1_2(end) = 1;
G1 = zeros(length(t),1);
G2 = zeros(length(t),1);

G1(1) = phi1_1'*phi0_1;
G2(1) = phi1_2'*phi0_2;
a1 = V1'*phi0_1;
b1 = V1'*phi1_1;
e1 = diag(D1);
a2 = V2'*phi0_2;
b2 = V2'*phi1_2;
e2 = diag(D2);
temp1 = exp(-1i*e1*dt);
tran1 = diag(temp1);
temp2 = exp(-1i*e2*dt);
tran2 = diag(temp2);
for i = 2:length(t)
    a1 = tran1*a1;
    G1(i) = b1'*a1;
%     G(i) = -log(abs(G(i)))/L;
    G1(i) = norm(G1(i));
    
    a2 = tran2*a2;
    G2(i) = b2'*a2;
%     G(i) = -log(abs(G(i)))/L;
    G2(i) = norm(G2(i));
end

figure;
plot(t,G1);
xlabel('time')
ylabel('transition probility')
str = strcat('transition probility of spin half');
title(str)
% fname = ['transition probility_spin half','.png '];
% saveas(gcf, fname, 'png')

figure;
plot(t,G2);
xlabel('time')
ylabel('transition probility')
str = strcat('transition probility of spin int');
title(str)

toc