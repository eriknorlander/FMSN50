function [accepted, t] = drawt(lambda, t, tau, rho)
d = length(t) - 1;
accepted = zeros(1,d-1);
%Metropolis-Hasting algorithm
for i = 2:d
    R = rho(1)*(t(i+1)-t(i-1));
    Xstar =  t(i) - R + 2*R*rand;
    %     Xstar = sort(Xstar)';
    %
    while(Xstar < t(i-1) || Xstar >  t(i+1))
        Xstar =  t(i) - R + 2*R*rand;
        % Xstar = sort(Xstar)';
    end
    
    num = f(lambda, [t(1:i-1) Xstar t(i+1:end)], tau);
    den = f(lambda, t, tau);
    alpha = min(1, num/den);
    U = rand(1);
    if U <= alpha
        t(i) = Xstar;
        accepted(i-1) = accepted(i-1)+1;
    end
end
end


