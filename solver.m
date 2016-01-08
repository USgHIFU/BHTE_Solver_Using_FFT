clear;

x_min = -0.05;
x_max = 0.05;
y_min = -0.05;
y_max = 0.05;
z_min = -0.03;
z_max = 0.03;

rou_t = 1140;
C_t = 3570;
A = rou_t * C_t;
k_t = 0.498;
B = k_t;
omega_t = 2 * 3.97e-4;
rou_b = 1050;
C_b = 3850;
C = omega_t * rou_b * C_b;

delta_space = 0.005;

T_init = 37;
delta_t = 0.1;

load Q3DNew
Q = Q3D1;
[N, M, P] = size(Q);
Q_star = fft2(Q);

freq_x = 1 / (x_max - x_min);
freq_y = 1 / (y_max - y_min);
freq_z = 1 / (z_max - z_min);

[K_x, K_y, K_z] = meshgrid(0:N-1, 0:M-1, 0:P-1);
% nju_square = (K_x * freq_x / N).^2 +...
%              (K_y * freq_y / M).^2 +...
%              (K_z * freq_z / P).^2;

nju_square = (K_x).^2 +...
             (K_y).^2 +...
             (K_z).^2;

k = exp(-(4*pi^2 * nju_square * B) * delta_t / A);
b = (1 - k) ./ (4*pi^2 * nju_square * B);
b(1,1,1) = 0;

time = 20;
T_init_star = fft2(T_init * ones(size(nju_square)));


T = T_init * ones(N, M, P, floor(time/delta_t));

t_index = 1;
tic
while (time >= 0)
    T_star = k .* T_init_star + b .* Q_star;
    T_boundary = abs(ifft2(T_star));
%     T_boundary = T_boundary(:);
%     index = find(T_boundary < T_init);
%     T_boundary(index) = T_init;
%     T_boundary = reshape(T_boundary, N, M, P);
    T(:,:,:,t_index) = T_boundary;
    T_init_star = fft2(T_boundary);
    time = time - delta_t;
    t_index = t_index + 1;
end
toc