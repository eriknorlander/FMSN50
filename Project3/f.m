function ft = f(lambda, t, tau)
d = length(t)-1;
tdiff = zeros(d, 1);

n = zeros(d, 1);
for i = 1:d
    bla = find(tau >= t(i) & tau<t(i+1));
    n(i) = length(bla);
end

for l = 1:d
    tdiff(l) = t(l+1) - t(l);
end
%-----EXP? Trodde man skulle ha på logaritmisk form!!
ft = exp(sum(log(lambda).*n + log(tdiff) - lambda.*tdiff));
end
