function [sol, x_trans, x_sep] = pohlhausen(x, U, U_d, U_dd, ReL, Re_tr)

if nargin < 6
    Re_tr = 5e5;
end

N = length(x);

x_trans = nan;
x_sep = nan;
found_transition = false;
Z = nan(size(x));
K = nan(size(x));
Lambda = nan(size(x));
f1 = nan(size(x));
f2 = nan(size(x));
F = nan(size(x));
Rex = nan(size(x));

Lambda(1) = 7.0523231;
K(1) = 0.0770356;
Z(1) = K(1) / U_d(1);
F(1) = 0;
f1(1) = (.3 - Lambda(1)/120) / (37/315 - Lambda(1)/945 - Lambda(1)^2/9072);
f2(1) = (2 + Lambda(1)/6) * (37/315 - Lambda(1)/945 - Lambda(1)^2/9072);

for i = 2:length(x)
  if i==2
    Z(i) = Z(1);
  else
    Z(i) = Z(i-1) + F(i-1) / U(i-1) * (x(i) - x(i-1));
  end
  K(i) = U_d(i) .* Z(i);
  Lambda(i) = find_lambda(K(i), Lambda(i-1), N);
  f1(i) = (.3 - Lambda(i)/120) / (37/315 - Lambda(i)/945 - Lambda(i)^2/9072);
  f2(i) = (2 + Lambda(i)/6) * (37/315 - Lambda(i)/945 - Lambda(i)^2/9072);
  F(i) = 2*f2(i) - 4*K(i) - 2*K(i)*f1(i);
  Rex = ReL * x .* U;
  if (Rex(i) >= Re_tr && ~found_transition)
    x_trans = x(i);
    found_transition = true;
  end
  if Lambda(i) <= -12
    x_sep = x(i);
    Lambda(i) = nan;
    Z(i) = nan;
    break;
  end
end

delta = sqrt(Lambda/ReL./U_d);
delta_2 = sqrt(Z/ReL);
delta_1 = f1 .* delta_2;
tau = sqrt(ReL) * f2 .* U ./ sqrt(Z);
cf = 2/ReL * tau;

%sol = [x', U', U_d', U_dd', Rex', Z', K', Lambda', f1', f2', F', delta', delta_1', delta_2', tau', cf'];
sol = [x, U, U_d, U_dd, Rex, Z, K, Lambda, f1, f2, F, delta, delta_1, delta_2, tau, cf];

end % endfunction

