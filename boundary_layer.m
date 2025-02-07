function [sol] = boundary_layer(x, U, U_d, ReL, nu, Lambda_0)

if(~exist("Lambda_0", "var"))
    Lambda_0 = 7.0523231;
end

N = length(x);

x_tr = nan;
x_sep = nan;
found_transition = false;
Z = nan(size(x));
K = nan(size(x));
Lambda = nan(size(x));
f1 = nan(size(x));
f2 = nan(size(x));
F = nan(size(x));
Rex = nan(size(x));
delta = nan(size(x));
delta_1 = nan(size(x));
delta_2 = nan(size(x));
tau = nan(size(x));
cf = nan(size(x));
Re_tr = nan(size(x));

Lambda(1) = Lambda_0;
K(1) = 0.0770356;
Z(1) = K(1) / U_d(1);
F(1) = 0;
f1(1) = (.3 - Lambda(1)/120) / (37/315 - Lambda(1)/945 - Lambda(1)^2/9072);
f2(1) = (2 + Lambda(1)/6) * (37/315 - Lambda(1)/945 - Lambda(1)^2/9072);
delta(1) = 0;
delta_1(1) = 0;
delta_2(1) = 0;
cf(1) = 0;
tau(1) = 0;
Re_tr(1) = 10^(-40.4557 + 64.8066*f1(1) - 26.7538*f1(1).^2 + 3.3819*f1(1).^3);
Rex(1) = 0;

%% pohlhausen's method
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
    Rex(i) = U(i)*x(i)*ReL;
    if(Lambda(i) <= -12) % check laminar separation
        x_sep = x(i);
        Lambda(i) = nan;
        Z(i) = nan;
        break;
        return;
    end
    delta(i) = sqrt(Lambda(i)/ReL/U_d(i));
    delta_2(i) = sqrt(Z(i)/ReL); % normalized
    delta_1(i) = f1(i) * delta_2(i);
    Re_tr(i) = 10^(-40.4557 + 64.8066*f1(i) - 26.7538*f1(i).^2 + 3.3819*f1(i).^3); % H = f1
    if (Rex(i) >= Re_tr(i) && ~found_transition) % check transition
        x_tr = x(i);
        found_transition = true;
        break;
    end
    tau(i) = sqrt(ReL) * f2(i) / U(i) / sqrt(Z(i)); % normalized, this needs fixing
    cf(i) = 2/ReL * tau(i);
end


%% head's method
H = nan(length(x)-i+1, 1);
H1 = nan(length(x)-i+1, 1);
H(1) = 1.4; % an approximate value corresponding to flat-plate flow
H1(1) = G(H(1));

[states] = RK4(@(x, states, nu, U, U_d) f_dot(x, states, nu, U, U_d), x(i-1:end), [delta_2(i-1); U(i-1)*delta_2(i-1)*H1(1)], nu, U(i-1:end), U_d(i-1:end));
states = states.';
delta_2(i-1:end) = states(:, 1);
H1 = states(2:end, 2) ./ U(i:end) ./ delta_2(i:end);
Rtheta = U(i:end).*delta_2(i:end)/nu;
H = arrayfun(@HofH1, H1);
cf(i:end) = 0.246*10.^(-0.678*H).*Rtheta.^(-0.268);
tau(i:end) = cf(i:end)*ReL/2;
delta_1(i:end) = H.*delta_2(i:end);
delta(i:end) = H1.*delta_2(i:end) + delta_1(i:end);
x_sep = find(isnan(states(:, 1))==1, 1);

sol.names = {"x", "U", "U_d", "Rex", "Z", "K", "Lambda", "f1", "f2", "F", "delta", "delta_1", "delta_2", "tau", "cf"};
sol.data = [x, U, U_d, Rex, Z, K, Lambda, f1, f2, F, delta, delta_1, delta_2, tau, cf];
sol.x_tr = x_tr;
sol.x_sep = x_sep;

end % endfunction

