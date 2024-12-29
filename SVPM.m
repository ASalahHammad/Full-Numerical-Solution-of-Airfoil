function [CL, CD, CM, Vt, Vx, Vy, Vxy] = SVPM(XB, YB, XC, YC, Vinf, AoA, N, x_domain, y_domain)

if nargin == 8
    error("x_domain and y_domain must coexist");
end
if nargin < 8
    x_domain = [];
    y_domain = [];
end
if nargin < 7
    N = 200;
end

num_panels = N - 1;

S = sqrt((XB(2:N)-XB(1:N-1)).^2 + (YB(2:N)-YB(1:N-1)).^2);
phi = atan2(YB(2:N)-YB(1:N-1), XB(2:N)-XB(1:N-1));
phi(phi<0) = phi(phi<0) + 2*pi;
delta = phi + pi/2;
beta = delta - AoA;
beta(beta>2*pi) = beta(beta>2*pi) - 2*pi;

%% SVPM
%% Solve for I & J matrices
I = zeros(num_panels, num_panels);
J = zeros(num_panels, num_panels);
for i = 1:num_panels
  for j = 1:num_panels
    if i==j
      continue;
    end
    A = -(XC(i) - XB(j))*cos(phi(j)) - (YC(i) - YB(j))*sin(phi(j));
    B = (XC(i) - XB(j))^2 + (YC(i) - YB(j))^2;
    Cn = sin(phi(i) - phi(j));
    Dn = -(XC(i) - XB(j))*sin(phi(i)) + (YC(i) - YB(j))*cos(phi(i));
    Ct = -cos(phi(i) - phi(j));
    Dt = (XC(i) - XB(j))*cos(phi(i)) + (YC(i) - YB(j))*sin(phi(i));
    E = sqrt(B - A^2);
    E = E*(isreal(E));
    I(i, j) = Cn/2*log((S(j)^2+2*A*S(j)+B) / B)  +  (Dn-A*Cn)/E*(atan2((S(j)+A), E) - atan2(A, E));
    J(i, j) = Ct/2*log((S(j)^2+2*A*S(j)+B) / B)  +  (Dt-A*Ct)/E*(atan2((S(j)+A), E) - atan2(A, E));
  end
end
I(isnan(I) | ~isreal(I)) = 0;
J(isnan(J) | ~isreal(J)) = 0;

%% Solve for K & L matrices
K = zeros(num_panels, num_panels);
L = zeros(num_panels, num_panels);
for i = 1:num_panels
  for j = 1:num_panels
    if i==j
      continue;
    end
    A = -(XC(i) - XB(j))*cos(phi(j)) - (YC(i) - YB(j))*sin(phi(j));
    B = (XC(i) - XB(j))^2 + (YC(i) - YB(j))^2;
    Cn = -cos(phi(i) - phi(j));
    Dn = (XC(i) - XB(j))*cos(phi(i)) + (YC(i) - YB(j))*sin(phi(i));
    Ct = sin(phi(j) - phi(i));
    Dt = (XC(i) - XB(j))*sin(phi(i)) - (YC(i) - YB(j))*cos(phi(i));
    E = sqrt(B - A^2);
    E = E*(isreal(E));
    K(i, j) = Cn/2*log((S(j)^2+2*A*S(j)+B) / B)  +  (Dn-A*Cn)/E*(atan2((S(j)+A), E) - atan2(A, E));
    L(i, j) = Ct/2*log((S(j)^2+2*A*S(j)+B) / B)  +  (Dt-A*Ct)/E*(atan2((S(j)+A), E) - atan2(A, E));
  end
end
K(isnan(K) | ~isreal(K)) = 0;
L(isnan(L) | ~isreal(L)) = 0;

%% Solve linear system of equations
a = nan(N, N);
b = nan(N, 1);
a(1:num_panels, 1:num_panels) = I + pi*eye(size(I));
for i=1:size(I,1)
  a(i, end) = -sum(K(i, :));
end
for j=1:size(I,2)
  a(end, j) = J(1, j) + J(end, j);
end
a(end, end) = -sum(L(1,:) + L(end, :)) + 2*pi;

b(1:num_panels) = -Vinf*2*pi*cos(beta);
b(end) = -Vinf*2*pi*(sin(beta(1)) + sin(beta(end)));

result = a\b;
lambda = result(1:num_panels);
gamma = result(N);

%% Calculate Velocity Distribution on Airfoil
Vn = Vinf*cos(beta) + I * lambda / 2 / pi - K * gamma*ones(size(lambda)) / 2 / pi + 0.5*lambda;
Vt = Vinf*sin(beta) + J * lambda / 2 / pi - L * gamma*ones(size(lambda)) / 2 / pi + gamma/2*ones(size(lambda));

fprintf("Sanity Check: sum(Vn) = %f\n\n", sum(Vn));

Cp = 1 - (Vt / Vinf).^2;

%% Calculate Lift, Drag & Moment
CN = -Cp .* S .* sin(beta);
CA = -Cp .* S .* cos(beta);

CL = sum(CN .* cos(AoA) - CA .* sin(AoA));
CD = sum(CN .* sin(AoA) + CA .* cos(AoA));
CM = sum(Cp .* (XC-0.25).*S.*cos(phi));

Vx = zeros(size(x_domain));
Vy = zeros(size(y_domain));

%% Compute N matrices
for i=1:length(Vx)
  XP = x_domain(i);
  if y_domain(i)>=0
    YP = y_domain(i) + 0.000001;
  else
    YP = y_domain(i) - 0.000001;
  end

  A = -(XP - XB(1:num_panels)).*cos(phi) - (YP - YB(1:num_panels)).*sin(phi);
  B = (XP - XB(1:num_panels)).^2 + (YP - YB(1:num_panels)).^2;
  Cx = sin(phi);
  Dx = -(YP - YB(1:num_panels));
  Cy = -cos(phi);
  Dy = (XP - XB(1:num_panels));
  E = sqrt(B-A.^2);
  E = E*(isreal(E));
  Nx = Cx/2.*log((S.^2+2*A.*S+B)./B) + (Dx-A.*Cx)./E.*(atan2((S+A), E) - atan2(A, E));
  Ny = Cy/2.*log((S.^2+2*A.*S+B)./B) + (Dy-A.*Cy)./E.*(atan2((S+A), E) - atan2(A, E));

  Vx(i) = Vinf*cos(AoA) - sum(gamma.*Nx)/2/pi;
  Vy(i) = Vinf*sin(AoA) - sum(gamma.*Ny)/2/pi;
end

Vxy = sqrt(Vx.^2 + Vy.^2);
CpXY = 1 - (Vxy/Vinf).^2;

end % endfunction

