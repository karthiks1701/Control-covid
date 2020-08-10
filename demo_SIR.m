function demo_SIR()
    dt = 1;                 
    D = 30;                 
    N_t = floor(D);  
    T = dt*N_t;             
    U_0 = [700 165 90];
    
    f_handle = @f;
    
    [u, t] = ode_FE(f_handle, U_0, dt, T);
    
    S = u(:,1);
    I = u(:,2);
    R = u(:,3);
    plot(t, S, 'b-', t, I, 'r-', t, R, 'g-');
    legend('S', 'I', 'R');
    xlabel('hours');
    % Consistency check:
    N = S(1) + I(1) + R(1);
    eps = 1E-12;           % Tolerance for comparing real numbers
    for n = 1:length(S)
        err = abs(S(n) + I(n) + R(n) - N);
        if (err > eps)
            error('demo_SIR: error=%g', err);                        
        end
    end
end

function result = f(u,t)
    beta = 0.002;   
    gamma = 0.25; 
    lambda = 0.56;
    N = 955;

    S = u(1); I = u(2); R = u(3);
    result = [-(beta)*S*I+gamma*N-gamma*S  (beta*S*I)-(gamma+lambda)*I (lambda*I-gamma*R)];
end



function [u, t] = ode_FE(f, U_0, dt, T)
    N_t = floor(T/dt);
    u = zeros(N_t+1, length(U_0));
    t = linspace(0, N_t*dt, length(u));
    u(1,:) = U_0;      % Initial values
    t(1) = 0;
    for n = 1:N_t
        u(n+1,:)  = u(n,:) + dt*f(u(n,:), t(n));
    end
end