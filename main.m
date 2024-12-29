clear
clc
close all

tic;

%% Given
Vinf = 1;  % Freestream velocity
AoA = 5*pi/180;  % Angle of attack [rad]

%% Load Airfoil
airfoil_name = "naca2412";
fprintf("====================== Potential Flow ======================\n");
figure; hold on; grid on;
for N_B = 100:20:240 % Number of boundary points
[XB, YB] = LOAD_AIRFOIL(airfoil_name, rad2deg(AoA), N_B, true, "xfoil");
num_panels = N_B - 1;

%% check cw or ccw
edges = zeros(num_panels, 1);
edges(1:num_panels) = (XB(2:N_B)-XB(1:N_B-1)) .* (YB(2:N_B)+YB(1:N_B-1));
if(sum(edges)<0)
    XB = flipud(XB);
    YB = flipud(YB);
end
%% Centre of Panels
XC = (XB(2:N_B)+XB(1:N_B-1))/2;
YC = (YB(2:N_B)+YB(1:N_B-1))/2;

[CL, CD, CM, Vt, Vx, Vy, Vxy] = SVPM(XB, YB, XC, YC, Vinf, AoA, N_B, [], []);

plot(N_B, CL, "x", "LineWidth", 2);
plot(N_B, CM, "o", "LineWidth", 2);

end
% Print the results to the Command Window
fprintf('Lift Coefficient (CL) =  %2.4f\n',CL);
fprintf('Drag Coefficient (CD) =  %2.4f\n',CD);
fprintf('Moment Coefficient @c/4 (CM) = %2.4f\n',CM);

xlabel("N");
legend("CL", "CM");
title("Coefficient of Lift and Moment vs. no. of panels nodes");

fprintf("\n\n====================== Boundary Layer ======================\n");
c = 1; % m
L = c; % m
rho = 1.225; % kg.m^-3
nu = 1.4607e-5; % m^2/s
mu = 1.7894e-5; % kg.m^-1.s^-1
ReL = rho * Vinf * L / mu;
figure; hold on; grid on;
% N shouldn't be too large because of instabilities, guessing: for fast separating flow (usually high AoA) -> N=80, for low AoA and stable solution you might reach up to N=500 without any problem
% anyway you should conduct a convergence analysis
% ONLY after that you can choose a single value for N and see the results at this value
% otherwise the results at the end of the code might be totally rubish

for N = 10:10:70
x = linspace(0, c, N).'; %% x should be uniform

%% Upper Surface
XC_u = XC(YC>=0); % x points for upper surface
V_u = Vt(YC>=0); % external velocity on upper surface
[U_u, U_d_u, U_dd_u] = velocity(XC_u, V_u, x);
[ANS_upper, x_trans_u, x_sep_u] = pohlhausen(x, U_u, U_d_u, U_dd_u, ReL);

%% Lower Surface
XC_l = XC(YC<0);
V_l = -Vt(YC<0);
[U_l, U_d_l, U_dd_l] = velocity(XC_l, V_l, x, "CD4");
[ANS_lower, x_trans_l, x_sep_l] = pohlhausen(x, U_l, U_d_l, U_dd_l, ReL);

%% Calculate Drag
cf_u = ANS_upper(:, 16);
cf_l = ANS_lower(:, 16);

CD = trapz(x(~isnan(cf_u)), cf_u(~isnan(cf_u))) + trapz(x(~isnan(cf_l)), cf_l(~isnan(cf_l)));

plot(N, CD, "x", "LineWidth", 2);
end

xlim([0, N]);
xlabel("N");
ylabel("CD");
title("Drag Coefficient vs. no. of grid points");

fprintf("    The solution for the last value of N (no. of grid points)\n");
fprintf("Upper Surface:\n");
fprintf("X_tr = %d\n", x_trans_u);
fprintf("Separation occured just after x = %d\n", x_sep_u);
fprintf("x     U     U_d     U_dd     Rex     Z     K     Lambda     f1     f2     F     delta     delta_1     delta_2     tau     cf\n");
fprintf("%.2f   %.2f   %.2f   %.2f   %.0f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.0f   %.2f\n", ANS_upper.');

fprintf("\nLower Surface:\n");
fprintf("X_tr = %f\n", x_trans_l);
fprintf("Separation occured just after x = %d\n", x_sep_l);
fprintf("x      U      U_d      U_dd      Rex     Z     K     Lambda    f1    f2   F     delta     delta_1     delta_2    tau   cf\n");
fprintf("%.2f   %.2f   %.2f   %.2f   %.0f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.0f   %.2f\n", ANS_lower.');

fprintf("\n\nTotal Drag Coefficient CD = %f\n", CD);

toc;

