
function sir_without_control()
D = 50;              
t = linspace(0, D, D+1);
S = zeros(D+1, 1);
I = zeros(D+1, 1);
R = zeros(D+1, 1);
I_w = zeros(D+1, 1);
S_w = zeros(D+1, 1);
R_w = zeros(D+1, 1);

S(1) = 1000;
I(1) = 100;
S_w(1) = 1000;
I_w(1) = 100;
R(1) = 50;
R_w(1) = 50;
dt = 1;

gamma = 0.25;
beta = 0.002;
lambda = 0.56;



N = S(1) + I(1) + R(1);
    for n=1:D
        ret = delay(I(n),n);
        S(n+1) = S(n) + gamma*N*dt-gamma*S(n)*dt -((ret)*beta*I(n)*S(n))*dt;
        I(n+1) = I(n) + (ret)*beta*I(n)*S(n)*dt - (lambda+gamma)*I(n)*dt;
        R(n+1) = R(n) + lambda*I(n)*dt - gamma*R(n)*dt;  
        
        S_w(n+1) = S_w(n) + gamma*N*dt -gamma*S_w(n)*dt -(beta*I_w(n)*S_w(n))*dt;
        I_w(n+1) = I_w(n) + (beta*I_w(n)*S_w(n))*dt - (lambda+gamma)*I_w(n)*dt;
        R_w(n+1) = R_w(n) + lambda*I_w(n)*dt - gamma*R_w(n)*dt;
    end

gcf;
plot(t,I,t,I_w);
title('SIR Model');
legend('i(t)','i_w');
xlabel('time, i(t),i_w(t)');
ylabel('Population');

end

function ret = delay(I,n)
  ret =1; 
 if (I>200)
      ret = 0.7;
  end
  
  
end