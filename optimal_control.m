function optimal_control()
  D = 30;
  t = linspace(0, D, D+1);
  S = zeros(D+1, 1);
  I = zeros(D+1, 1);
  R = zeros(D+1, 1);
  
  S_w = zeros(D+1, 1);
  I_w = zeros(D+1, 1);
  R_w = zeros(D+1, 1);
   
  lambda1 = zeros(D+1,1);
  lambda2 = zeros(D+1,1);
  lambda3 = zeros(D+1,1);
  
  
  u1 = zeros(D+1,1);
  u2 = zeros(D+1,1);
  
  
  u1(1) = 0;
  u2(1) = 0;
  
  S_w(1) = 800;
  I_w(1) = 110;
  R_w(1) = 70;
  S(1) = 800;
  I(1) = 110;
  R(1) = 70;
  dt = 1;
  gammaS = 0.25;  % death rates in suspetibles population 
  gammaI = 0.25;  % death rates in infected population
  gammaR = 0.25;  % death rates in recovered population
  beta = 0.002;
  lambda = 0.56;
  c1 = 0.06;
  c2 = 0.09;
  
  D = 30;
  
  N = S(1) + I(1) + R(1);
  for n=1:D-1
        S_w(n+1) = S_w(n) + gammaS*N*dt -gammaS*S_w(n)*dt-(beta*I_w(n)*S_w(n))*dt;
        I_w(n+1) = I_w(n) + (beta*I_w(n)*S_w(n))*dt -(lambda+gammaI)*I_w(n)*dt;
        R_w(n+1) = R_w(n) + lambda*I_w(n)*dt - gammaR*R_w(n)*dt  ; 
      
      
        S(n+1) = S(n) + gammaS*N*dt -gammaS*S(n)*dt -(beta*I(n)*S(n))*dt - u1(n)*dt;
        I(n+1) = I(n) + (beta*I(n)*S(n))*dt - (lambda+gammaI)*I(n)*dt - u2(n)*dt;
        R(n+1) = R(n) + lambda*I(n)*dt - gammaR*R(n)*dt + u1(n)*dt + u2(n)*dt ;  
        lambda1(D-n)= lambda1(D-n+1) - dt*( -1 +  (lambda1(D-n+1) -lambda2(D-n+1))*beta*I(n+1) + (gammaS)*lambda1(D-n+1));     
        lambda2(D-n)= lambda2(D-n+1) - dt*( -1 +  (lambda1(D-n+1) -lambda2(D-n+1))*beta*S(n+1) + (gammaI + lambda)*lambda2(D-n+1) - lambda*lambda3(D-n+1)); 
        lambda3(D-n)= lambda3(D-n+1) - dt*(gammaR*lambda3(D - n+1));
        u1(n+1) = (lambda1(D-n+1)-lambda3(D-n+1))/c1;
        u2(n+1) = (lambda2(D-n+1)-lambda3(D-n+1))/c2;
  
        if(u1(n+1)<0)
            u1(n+1) = 0;
        end
        if(u1(n+1)>500)
            u1(n+1) = 500;
        end
  
        if(u2(n+1)<0)
            u2(n+1) = 0;
        end
        
        if(u2(n+1)>500)
            u2(n+1) = 500;
        end
  end
plot(t,lambda1,t,lambda2,t,lambda3);
title('SIR Model');
legend('lambda1','lambda2','lambda3');
xlabel('time,l1,l2,l3');
ylabel('Population');

  
%plot(t,S,t,I,t,R,t,S_w,t,I_w,t,R_w);
%title('SIR Model');
%legend('S(t)','I(t)','R(t)','S_w(t)','I_w(t)','R_w(t)');
%xlabel('time, S(t),I(t),R(t),S_w(t),I_w(t),R_w(t)');
%ylabel('Population');
end
  

 