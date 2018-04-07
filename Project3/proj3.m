%% This section plots the trajectory and the breakpoints markings in the 
% coal-mine-disasters-plot.
close all
clear

%Number of disasters
load('coal_mine_disasters.mat')

%Defining start and end points
t1 = 1658;
td1 = 1980;
d = 6;

%Defining the initial breakpoints, this interval is chosen by empirical
%observations
tt1 = 1960;
tt2 = 1690;
step = (tt1-tt2)/d;
tmid = tt2:step:tt1;
t =[t1 tmid(2:end-1) td1];

samples = 1e5;
burn_in =  1e3;

rho = 0.005*ones(d, 1);
psi = 20;
tau = T;

burned = zeros(burn_in, length(t));
totalt = zeros(samples, length(t));

theta = gamrnd(2, 1/psi);
lambda = gamrnd(2, 1/theta, 1, d);

for i = burn_in
    theta = drawTheta(lambda, psi);
    lambda = drawLambda(theta, t, tau);
    [~, t] = drawt(lambda, t, tau, rho);
    burned(i, :) = t;
end

for i = 1:samples
    theta = drawTheta(lambda, psi);
    lambda = drawLambda(theta, t, tau);
    [~,t] = drawt(lambda, t, tau, rho);
    totalt(i, :) = t;
end

figure
plot(totalt)
title(['The behavior of the chain for d = ' num2str(d-1) ' of breakpoints'])
xlabel('Number of samples')
ylabel('Breakpoints')
set(gca, 'Fontsize', 16);

figure
plot(T, 1:length(T))
hold on
for i = 2:d
    line([mean(totalt(:, i)) mean(totalt(:, i))], [0 length(T)], 'Color', [rand rand rand])
end

title('Total number of accidents during 1658-1980')
xlabel('Breakpoints')
ylabel('Number of accidents')
set(gca, 'Fontsize', 16);
axis([t1 T(end) 0 length(T)])

%% 1c) Plotting how the number of breakpoints affect the results
close all
clear

%Number of disasters
load('coal_mine_disasters.mat')

samples = 1e4;
burn_in = 1e3;

%Defining start and end points
t1 = 1658;
td1 = 1980;

d=5; %This variable regulates the number of lambda values

for count = 2:d
    d = count;
    tt1 = 1960;
    tt2 = 1690;
    step = (tt1-tt2)/d;
    tmid = tt2:step:tt1;
    t =[t1 tmid(2:end-1) td1];
    
    rho = 0.01*ones(d, 1);
    psi = 20;
    tau = T;
    
    totalt = zeros(samples, length(t));
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);
    
    for i = 1:burn_in
        theta = drawTheta(lambda, psi);
        lambda = drawLambda(theta, t, tau);
        [~, t] = drawt(lambda, t, tau, rho);
    end
    
    for i = 1:samples
        theta = drawTheta(lambda, psi);
        lambda = drawLambda(theta, t, tau);
        [~, t] = drawt(lambda, t, tau, rho);
        totalt(i, :) = t;
    end
    
    figure
    plot(totalt)
    title(['The behavior of the chain for d = ' num2str(d-1) ' of breakpoints'])
    xlabel('Number of samples')
    ylabel('Breakpoints')
    set(gca, 'Fontsize', 16);
    
    figure
    for i = 2:d
        histogram(totalt(:,i), 20)
        hold on
    end
    title(['The histogram of the chain for d = ' num2str(d-1) ' of breakpoints'])
    xlabel('Breakpoints')
    ylabel('Intensity')
    set(gca, 'Fontsize', 16);
end

%% 1d) how does psi affect lambda?
clear
close all

load('coal_mine_disasters.mat')
tau = T;

%Defining start and end points
t1 = 1658;
td1 = 1980;

tt1 = 1960;
tt2 = 1690;
d = 5;

step = (tt1-tt2)/d;
tmid = tt2:step:tt1;
t =[t1 tmid(2:end-1) td1];

psi = 50;
rho = 0.01*ones(d,1);

n = 1e4;
burn_in = 1e3;

lambdaMean = zeros(psi, d);
lambdaVar = zeros(psi, d);

for psi = 1:psi
    lambdatemp = zeros(n, d);
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);

    for i = 1:burn_in
    theta = drawTheta(lambda, psi);
    lambda = drawLambda(theta, t, tau);
    [~, t] = drawt(lambda, t, tau, rho);
    end
    
    for i = 1:n
        theta = drawTheta(lambda, psi);
        lambda = drawLambda(theta, t, tau);
        [~, t] = drawt(lambda, t, tau, rho);
      
        lambdatemp(i, :) = lambda';
    end
    lambdaMean(psi, :) = mean(lambdatemp);
    lambdaVar(psi, :) = var(lambdatemp);
end

figure
plot(lambdaMean, '*')
title('The mean of the intensities \lambda''s dependency on \psi')
xlabel('\psi')
ylabel('\lambda')
legend('\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4', '\lambda_5')
set(gca, 'Fontsize', 16);

figure
plot(lambdaVar, '*')
title('The variance of intensities \lambda''s dependency on \psi')
xlabel('\psi')
ylabel('\lambda')
legend('\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4', '\lambda_5')
set(gca, 'Fontsize', 16);

%% 1d) how does psi affect theta?
clear
close all

load('coal_mine_disasters.mat')
tau = T;

%Defining start and end points
t1 = 1658;
td1 = 1980;

d = 5;
step = (td1-t1)/d;
t =t1:step:td1;

psi = 50;
rho = 0.025*ones(d,1);

n = 1e4;
burn_in = 1e3;

thetaMean = zeros(psi, 1);
thetaVar = zeros(psi, 1);
thetatemp = zeros(n);

for psi = 1:psi
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);
    thetatemp = zeros(n, 1);
    
    for i = 1:burn_in
        theta = drawTheta(lambda, psi);
        lambda = drawLambda(theta, t, tau);
        [~, t] = drawt(lambda, t, tau, rho);
    end
    
    for i = 1:n
        theta = drawTheta(lambda, psi);
        lambda = drawLambda(theta, t, tau);
        [~, t] = drawt(lambda, t, tau, rho);
        
        thetatemp(i) = theta;
    end
    
    thetaMean(psi) = mean(thetatemp);
    thetaVar(psi) = var(thetatemp);
end

figure
plot(thetaMean, '*')
title('The mean of the parameter \theta''s dependency on \psi')
xlabel('\psi')
ylabel('\theta')
set(gca, 'Fontsize', 16);

figure
plot(thetaVar, '*')
title('The variance of the parameter \theta''s dependency on \psi')
xlabel('\psi')
ylabel('\theta')
set(gca, 'Fontsize', 16);
%% 1d) how does psi affect t?
clear
close all

load('coal_mine_disasters.mat')
tau = T;

%Defining start and end points
t1 = 1658;
td1 = 1980;

tt1 = 1960;
tt2 = 1690;
d = 5;
step = (tt1-tt2)/d;
tmid = tt2:step:tt1;
t =[t1 tmid(2:end-1) td1];

%psi = [1 25 50 1e5];
psi = 50;
rho = 0.01*ones(d,1);

n = 1e4;
burn_in = 1e3;

tVec = zeros(n, length(t));
tMean = zeros(psi, length(t));
tVar = zeros(psi, length(t));

for psi = 1:psi
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);
    for i = 1:burn_in
        theta = drawTheta(lambda, psi);
        lambda = drawLambda(theta, t, tau);
        [~, t] = drawt(lambda, t, tau, rho);
    end
    
    for i = 1:n
        theta = drawTheta(lambda, psi);
        lambda = drawLambda(theta, t, tau);
        [~, t] = drawt(lambda, t, tau, rho);
        tVec(i, :) = t;
    end
    tMean(psi, :) = mean(tVec);
    tVar(psi, :) = var(tVec);
end

figure
plot(tMean(:, 2:end-1), '*')
title('The mean of the breakpoints t dependency on \psi')
xlabel('\psi')
ylabel('Breakpoints')
legend('t_1', 't_2', 't_3', 't_4', 'Orientation', 'horizontal')
set(gca, 'Fontsize', 16); 

figure
plot(tVar(:, 2:end-1), '*')
title('The variance of the breakpoints dependency of \psi')
xlabel('\psi')
ylabel('Variance of breakpoints')
legend('t_1', 't_2', 't_3', 't_4', 'Orientation', 'horizontal')
set(gca, 'Fontsize', 16); 

%% 1d) how does rho affect lambda and theta?
clear
close all

load('coal_mine_disasters.mat')
tau = T;

%Defining start and end points
t1 = 1658;
td1 = 1980;

tt1 = 1960;
tt2 = 1690;
d = 5;
step = (tt1-tt2)/d;
tmid = tt2:step:tt1;
t =[t1 tmid(2:end-1) td1];

nor = 30;
psi = 20;
rho1 = (1:nor)*1e-2;
rho = zeros(d,nor);
for i = 1:d
   rho(i,:) = rho1; 
end
n = 1e4;
burn_in = 1e3;

thetaMean = zeros(nor, 1);
lambdaMean = zeros(nor, d);

thetatemp = zeros(n, 1);
lambdaTemp = zeros(n, d);

for nor = 1:nor
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);
    thetatemp = zeros(nor, 1);
    lambdaTemp = zeros(nor, d);
    
    for i = 1:burn_in
        theta = drawTheta(lambda, psi);
        lambda = drawLambda(theta, t, tau);
        [~, t] = drawt(lambda, t, tau, rho(:, nor));
    end
    
    for i = 1:n
        theta = drawTheta(lambda, psi);
        lambda = drawLambda(theta, t, tau);
        [~, t] = drawt(lambda, t, tau, rho(:, nor));
        
        thetatemp(i) = theta;
        lambdaTemp(i, :) = lambda';
    end
    thetaMean(nor) = mean(thetatemp);
    lambdaMean(nor, :) = mean(lambdaTemp);
end

figure(1)
plot(rho(1,:),lambdaMean, '*')
title('The intensities \lambda dependency on \rho')
xlabel('\rho')
ylabel('\lambda')
set(gca, 'Fontsize', 16);

figure(2)
plot(rho(1,:),thetaMean, '*')
title('The parameter \theta dependency on \rho')
xlabel('\rho')
ylabel('\theta')
set(gca, 'Fontsize', 16);

%% 1d) how does rho affect t?
clear
close all

load('coal_mine_disasters.mat')
tau = T;

%Defining start and end points
t1 = 1658;
td1 = 1980;

tt1 = 1960;
tt2 = 1690;
d = 5;
step = (tt1-tt2)/d;
tmid = tt2:step:tt1;
t =[t1 tmid(2:end-1) td1];

psi = 20;
rho1 = [0.01 0.02 0.03 0.04];
nor = length(rho1);
rho = zeros(d,nor);
for i = 1:d
   rho(i,:) = rho1; 
end

n = 1e4;
burn_in = 1e3;

tVec = zeros(n, length(t));

for nor = 1:nor
    tVec = zeros(n, length(t));
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);
    
    for i = 1:burn_in
        theta = drawTheta(lambda, psi);
        lambda = drawLambda(theta, t, tau);
        [~, t] = drawt(lambda, t, tau, rho(:, nor));
    end
    
    for i = 1:n
        theta = drawTheta(lambda, psi);
        lambda = drawLambda(theta, t, tau);
        [~, t] = drawt(lambda, t, tau, rho(:, nor));
        tVec(i, :) = t;
    end
    
    figure
    subplot(2,2,1)
    autocorr(tVec(:, 2), 500)
    val = rho(1, nor);
    title(['The correlation function for t_1 with \rho = ' num2str(val)])
    xlabel('Time lag')
    ylabel('Dependency')
    set(gca, 'Fontsize', 16);
    
    subplot(2,2,2)
    autocorr(tVec(:, 3), 500)
    val = rho(1, nor);
    title(['The correlation function for t_2 with \rho = ' num2str(val)])
    xlabel('Time lag')
    ylabel('Dependency')
    set(gca, 'Fontsize', 16);
    
    subplot(2,2,3)
    autocorr(tVec(:, 4), 500)
    val = rho(1, nor);
    title(['The correlation function for t_3 with \rho = ' num2str(val)])
    xlabel('Time lag')
    ylabel('Dependency')
    set(gca, 'Fontsize', 16);
    
    subplot(2,2,4)
    autocorr(tVec(:, 4), 500)
    val = rho(1, nor);
    title(['The correlation function for t_4 with \rho = ' num2str(val)])
    xlabel('Time lag')
    ylabel('Dependency')
    set(gca, 'Fontsize', 16);
    
    figure
    plot(tVec)
    title(['The behavior of the chain for d = ' num2str(d-1) ' of breakpoints and \rho =' num2str(rho(nor))])
    xlabel('Number of samples')
    ylabel('Breakpoints')
    set(gca, 'Fontsize', 16);
end
%% 1d) Acceptance ratio plot one of the previous section has 
%      to be run before this one.
close all
rho = linspace(0.001,0.1);
accepted = zeros(length(rho),d-1);
accSamp = 100;
for r = 1:length(rho)
   for i = 1:accSamp
       theta = drawTheta(lambda, psi);
       lambda = drawLambda(theta, t, tau);
       [acc, t] = drawt(lambda, t, tau, rho(r));
       accepted(r,:) = accepted(r,:)+acc;
   end
end

ratio = sum(accepted,2)/(accSamp*(d-1));
plot(rho,ratio,'*');
hold on
plot([0,0.1],[0.3,0.3]);
title('Acceptance ratio as a function of \rho');
xlabel('Acceptance ratio');
ylabel('Value of \rho');
set(gca, 'Fontsize', 16);
%% Exercise 2b) 95% confidence interval for beta and mu using bootstrap
clear
close all
rng(56)
load('atlantic.txt')
[beta, mu] = est_gumbel(atlantic);
Finv = @(u, mu, beta) beta.*log(1./log(1./u)) + mu;
nbr = length(atlantic);
samples = 1000;
alpha1 = 0.025*samples;
alpha2 = 0.975*samples;
betaVec = zeros(samples,1);
muVec = zeros(samples,1);
delta = zeros(samples,2);

for i = 1:samples
    u = rand(582,1);
    drawWave = Finv(u, mu, beta);
    [betaDraw, muDraw]  = est_gumbel(drawWave);
    betaVec(i) = betaDraw;
    muVec(i) = muDraw;
    delta(i,:) = [beta-betaDraw, mu-muDraw];
end

delta = sort(delta);
% Confidence interval for beta and mu with bootstrap

CI = zeros(3,2);
CI(1,:) = beta + [delta(floor(alpha1)), delta(ceil(alpha2))];      % Conf Beta
CI(2,:) = mu + [delta(floor(alpha1)), delta(ceil(alpha2))];        % Conf mu
%% Exercise 2c) 95% confidence interval for the 100-year wave using bootsrap
T = 3*14*100;
alpha5 = 0.95*samples;
drawWaveVec = zeros(samples,1);
for i = 1:samples
    drawWaveVec(i) = Finv(1-1/T, muVec(i), betaVec(i));
end
avWave = Finv(1-1/T, mu, beta);
deltaWave = sort(avWave-drawWaveVec);
CI(3,:) =[0, avWave + deltaWave(ceil(alpha5))]; % Conf wave
%%
x = linspace(0,1,582);
plot(x,drawWave);
figure(2)
plot(1:582,atlantic);