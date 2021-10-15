clear;
tic

L = 8;
len = 2^L;
h0 = 0.3;
sigma_z = [ 1 0; 0 -1];
sigma_x = [ 0 1; 1 0];
I = eye(2);
H = zeros(len,len);
% h = h0.*[1 0.8 0.6 0.4 -0.4 -0.6 -0.8 -1];
% h = h0.*[1 0.75 0.5 0.25 -0.25 -0.5 -0.75 -1];
h = h0.*[1 1 1 1 -1 -1 -1 -1];
% h = h0.*[7 5 3 1 -1 -3 -5 -7]./8;
% h = h0.*[1 -1];
dt = 0.01;
t = 0:dt:100;

% pos=1单独赋值
H1 = sigma_z;
H1 = kron(H1,sigma_z);
for j = 3:L
    H1 = kron(H1,I);
end
H = H - H1;

for i = 2:L-1
    H1 = I;
    for j = 2:i-1
        H1 = kron(H1,I);
    end
    H1 = kron(H1,sigma_z);
    H1 = kron(H1,sigma_z);
    for j = i+2:L
        H1 = kron(H1,I);
    end
    H = H - H1;
end

% pos=L单独赋值
H1 = sigma_z;
for j = 2:L-1
    H1 = kron(H1,I);
end
H1 = kron(H1,sigma_z);
H = H - H1;

% pos=1单独赋值
H2 = sigma_x;
for j = 2:L
    H2 = kron(H2,I);
end
H = H + h(1).*H2;

for i = 2:L
    H2 = I;
    for j = 2:i-1
        H2 = kron(H2,I);
    end
    H2 = kron(H2,sigma_x);
    for j = i+1:L
        H2 = kron(H2,I);
    end
    H = H + h(i).*H2;
end

[V,D] = eig(H);

phi0 = zeros(len,1);
phi1 = zeros(len,1);
phi0(1) = 1;
phi1(end) = 1;
G0 = zeros(length(t),1);
G = zeros(length(t),1);

G(1) = phi1'*phi0;
a = V'*phi0;
b = V'*phi1;
e = diag(D);
temp = exp(-1i*e*dt);
tran = diag(temp);
for i = 2:length(t)
    a = tran*a;
    G0(i) = b'*a;
%     G(i) = -log(abs(G(i)))/L;
    G(i) = norm(G0(i));
end

% for i = 2:length(t)
%     temp = exp(-1i*e*t(i));
%     tran = diag(temp);
%     phi = tran*a;
%     G(i) = b'*phi;
% %     G(i) = G(i)^2;
% end

figure;
plot(t,G);
xlabel('time')
ylabel('transition probility')
str = strcat('transition probility of spin half');
title(str)
fname = ['transition probility_spin half','.png '];
saveas(gcf, fname, 'png')

toc