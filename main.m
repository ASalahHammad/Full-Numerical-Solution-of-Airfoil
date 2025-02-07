clear
clc
close all

tic;

%% Given
Vinf = 1;  % Freestream velocity
AoA = 0*pi/180;  % Angle of attack [rad]

%% Load Airfoil
airfoil_name = "naca2412";
fprintf("====================== Potential Flow ======================\n");
figure; hold on; grid on;
for N_B = [100,200] % Number of boundary points
    [XB, YB] = LOAD_AIRFOIL(airfoil_name, rad2deg(AoA), N_B, true, "xfoil");
    num_panels = N_B - 1;
    
    % check cw or ccw
    edges = zeros(num_panels, 1);
    edges(1:num_panels) = (XB(2:N_B)-XB(1:N_B-1)) .* (YB(2:N_B)+YB(1:N_B-1));
    if(sum(edges)<0)
        XB = flipud(XB);
        YB = flipud(YB);
    end
    
    % Centre of Panels
    XC = (XB(2:N_B)+XB(1:N_B-1))/2;
    YC = (YB(2:N_B)+YB(1:N_B-1))/2;
    
    [CL, CD, CM, Vt, Vx, Vy, Vxy] = SVPM(XB, YB, XC, YC, Vinf, AoA, N_B, [], []);
    
    plot(N_B, CL, "x", "LineWidth", 2);
    plot(N_B, CM, "o", "LineWidth", 2);

    % Print the results of the last solution to the Command Window
    fprintf("N = %i\n", N_B);
    fprintf('Lift Coefficient (CL) =  %2.4f\n',CL);
    fprintf('Drag Coefficient (CD) =  %2.4f\n',CD);
    fprintf('Moment Coefficient @c/4 (CM) = %2.4f\n',CM);
end

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
% N shouldn't be too large because of instabilities, guessing: for fast
% separating flow (usually high AoA) -> N=70, for low AoA and stable
% solution you might reach up to N=500 without any problem, anyway you
% should conduct a convergence analysis, ONLY after that you can choose a
% single value for N and see the results at this value otherwise the
% results at your popping at your terminal might be totally rubish

for N = 10:10:100
    x = linspace(0, 1, N).'; % x is from 0 -> 1, because it should be normalized
    % x = (1-cos(linspace(0, pi, N))).'/2;
    % Upper Surface
    XC_u = XC(YC>=0); % x points for upper surface
    V_u = Vt(YC>=0); % external velocity on upper surface
    U_u = interp1(XC_u/c, V_u/Vinf, x, "spline", "extrap"); % Velocity should be normalized
    U_d_u = gradient(U_u, x);
    ans_upper = boundary_layer(x, U_u, U_d_u, ReL, nu);

    % Lower Surface
    XC_l = XC(YC<0);
    V_l = -Vt(YC<0);
    U_l = interp1(XC_l, V_l/Vinf, x, "spline", "extrap"); % Velocity should be normalized
    U_d_l = gradient(U_l, x);
    ans_lower = boundary_layer(x, U_l, U_d_l, ReL, nu);
    
    % Calculate Drag
    % This needs fixing because you might need to check normalization
    cf_u = ans_upper.data(:, find(cellfun(@(entry) isequal(entry, "cf"), ans_upper.names)));
    cf_l = ans_lower.data(:, find(cellfun(@(entry) isequal(entry, "cf"), ans_lower.names)));
    CD = trapz(x(~isnan(cf_u))*c, cf_u(~isnan(cf_u))) + trapz(x(~isnan(cf_l)), cf_l(~isnan(cf_l)));

    plot(N, CD, "x", "LineWidth", 2);
end

xlim([0, N]);
xlabel("N");
ylabel("$C_d$", 'Interpreter', 'latex', 'FontSize', 15);
title("Drag Coefficient vs. no. of grid points");

fprintf("The solution for the last value of N (no. of grid points)\n");
fprintf("Upper Surface:\n");
fprintf("X_tr = %f\n", ans_upper.x_tr);
fprintf("x_sep = %f\n", ans_upper.x_sep);
fprintf("%-8s %-6s %-7s %-10s %-10s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n", ans_upper.names{:});
fprintf("%-8.3f %-6.2f %-7.2f %-10.2f %-10.2e %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.0f %-8.2f\n", ans_upper.data.');

fprintf("\nLower Surface:\n");
fprintf("X_tr = %f\n", ans_lower.x_tr);
fprintf("x_sep = %f\n", ans_lower.x_sep);
fprintf("%-8s %-6s %-7s %-10s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n", ans_lower.names{:});
fprintf("%-8.3f %-6.2f %-7.2f %-10.2e %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.0f %-8.2f\n", ans_lower.data.');

fprintf("\n\nTotal Drag Coefficient CD = %f\n", CD);

toc;

%% Figures
figure; hold on; grid on; xlabel("$x$", 'interpreter', 'latex', 'FontSize', 15); ylabel("$U_e$", 'Interpreter', 'latex', 'FontSize', 15);
plot(x, U_u); plot(x, U_l);
legend("Upper", "Lower");

figure; hold on; grid on; xlabel("$x$", 'interpreter', 'latex', 'FontSize', 15); ylabel("$\delta$", 'Interpreter', 'latex', 'FontSize', 15);
plot(x, ans_upper.data(:, find(cellfun(@(entry) isequal(entry, "delta"), ans_upper.names)))); plot(x, ans_lower.data(:, find(cellfun(@(entry) isequal(entry, "delta"), ans_lower.names))));
legend("Upper", "Lower");

figure; hold on; grid on; xlabel("$x$", 'interpreter', 'latex', 'FontSize', 15); ylabel("$\delta_1$", 'Interpreter', 'latex', 'FontSize', 15);
plot(x, ans_upper.data(:, find(cellfun(@(entry) isequal(entry, "delta_1"), ans_upper.names)))); plot(x, ans_lower.data(:, find(cellfun(@(entry) isequal(entry, "delta_1"), ans_lower.names))));

figure; hold on; grid on; xlabel("$x$", 'interpreter', 'latex', 'FontSize', 15); ylabel("$\delta_2$", 'Interpreter', 'latex', 'FontSize', 15);
plot(x, ans_upper.data(:, find(cellfun(@(entry) isequal(entry, "delta_2"), ans_upper.names)))); plot(x, ans_lower.data(:, find(cellfun(@(entry) isequal(entry, "delta_2"), ans_lower.names))));

figure; hold on; grid on; xlabel("$x$", 'interpreter', 'latex', 'FontSize', 15); ylabel("$C_f$", 'Interpreter', 'latex', 'FontSize', 15);
plot(x, ans_upper.data(:, find(cellfun(@(entry) isequal(entry, "cf"), ans_upper.names)))); plot(x, ans_lower.data(:, find(cellfun(@(entry) isequal(entry, "cf"), ans_lower.names))));

figure; hold on; grid on; xlabel("$x$", 'interpreter', 'latex', 'FontSize', 15); ylabel("$\tau$", 'Interpreter', 'latex', 'FontSize', 15);
plot(x, ans_upper.data(:, find(cellfun(@(entry) isequal(entry, "tau"), ans_upper.names)))); plot(x, ans_lower.data(:, find(cellfun(@(entry) isequal(entry, "tau"), ans_lower.names))));
