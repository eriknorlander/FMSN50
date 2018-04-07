function theta = drawTheta(lambda, psi)
%Draw a new theta for the gibbs sampler
    d = length(lambda);
    theta = gamrnd(2*d + 2, 1./(psi + sum(lambda)));
end