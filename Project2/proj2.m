%% 3 A Naive approach - Random Walk 
clear
close all

d = 2; %Dimension
N = 5000; %Number of particles
walk = zeros(N, 2); %Initiating the walk-vector which remember the walk for a particle
dir = [eye(d); -1*eye(d)]; %Initiating possible direction for dimension d

for i = 3 : N + 1
    newDir = datasample(dir,1);
    next = newDir + walk(i-2, :);
    
    if ismember(next, walk, 'rows') == ones(1, 2)
        walk(i-1,:) = next;
        break
    end
    
    walk(i - 1, :) = next;
end
plot(walk(1:i-2, 1), walk(1:i-2, 2),'-*')
hold on
plot(walk(i-1,1),walk(i-1,2),'-*r')

%% Estimating Cn(2) using random walk
N = 5000; %Number of particles
walk = zeros(N, 2); 
n = [-1,0; 0,1; 1,0; 0,-1];
%c = zeros(10, 1);
nsa = zeros(10,1);

for NbrOfSteps = 1:10
    for NbrOfTries = 1:5000
        for i = 3:NbrOfSteps+2
            newDir = datasample(n,1);
            if i == NbrOfSteps+2
                nsa(NbrOfSteps) = nsa(NbrOfSteps) + 1;
                break
            end
            next = newDir + walk(i-2, :);
            if ismember(next, walk, 'rows') == ones(1, 2)
                walk(i-1,:) = next;
                break
            end
            walk(i - 1, :) = next;
        end
    end
    %c(NbrOfSteps)=(nsa(NbrOfSteps)./N).*(4.^NbrOfSteps);
end
k = 1:10;
k1 = 4.^k;
format long;
c1 = (nsa./N).*k1';

%% 4 Improving - actual self avoiding walk
close all
clear
d = 2; %Dimensions
N = 1e3; %Number of particles
walk = zeros(N, 2); %Initiating the walk-vector
g0 = zeros(N,1);
dir = [eye(d); -1*eye(d)];
neigh = zeros(2*d, d);

for i = 3 : N + 1
    for counter = 1:2*d
        neigh(counter, :) = walk(i-2, :) + dir(counter, :);
    end

    if ismember(neigh,walk,'rows') == ones(4,1)
        break
    end
    
    g0(i-2) = 1-(sum(ismember(neigh,walk,'rows'))/4);
    newDir = datasample(dir,1);
    
    while ismember(newDir + walk(i-2,:),walk,'rows') == ones(1,2)
        newDir = datasample(dir,1);
    end
    
    next = newDir + walk(i-2, :);
    walk(i-1, :) = next;
end

figure(1)
plot(walk(1:i-2, 1), walk(1:i-2, 2),'-*')
hold on
plot(walk(i-2, 1), walk(i-2,2), '*r')

%% 4 Estimating Cn(2) using SIS drawing from a gn = SAW
d = 2;
N = 1e3; %Number of particles
nbrOfSteps = 10;
nbrOfCycles = 10;
omega = zeros(N, nbrOfSteps);
g = zeros(N,1);
n = [eye(d); -1*eye(d)];
omega_0 = 1;

c2 = zeros(nbrOfCycles, nbrOfSteps);
for k = 1:nbrOfCycles
for nbrOfSteps = 1:nbrOfSteps
    for nbrOfTries = 1:N
        neigh = zeros(4, 2);
        walk = zeros(N, 2);
        for j = 3:nbrOfSteps+2
            for i = 1:4
                neigh(i, :) = walk(j-2, :) + n(i, :);
            end
            
            if ismember(neigh,walk,'rows') == ones(4,1)
                break
            end
            
            newDir = datasample(n,1);
            
            while ismember(newDir + walk(j-2,:),walk,'rows') == ones(1,2)
                newDir = datasample(n,1);
            end
            
            g(nbrOfTries) = 4-sum(ismember(neigh, walk, 'rows'));
            next = newDir + walk(j-2, :);
            walk(j-1, :) = next;
            
            if nbrOfSteps == 1
                omega(nbrOfTries, nbrOfSteps) = g(nbrOfTries)*omega_0;
            else
                omega(nbrOfTries, nbrOfSteps) = g(nbrOfTries).*omega(nbrOfTries, nbrOfSteps-1);
            end
        end
    end
end
c2(k,:) = mean(omega);
end
mean(c2);

%% Task 5. Estimating Cn(2) using SISR drawing from a gn = SAW
% This code works for arbitrary d and is used to calculate values for task
% 9
clear
close all

d = 2; %Dimension
d2 = 2*d;
N = 1e3; %Number of particles
nbrOfCycles = 10; %Number of cycles
stepLength = 10; %The length of the total step, i.e n
g = zeros(N,1); %The instrumental distribution function
n = [eye(d) ; -1*eye(d)]; %The possible directions in each step for a given d
neigh = zeros(d2, d); %The neighbours of a given point
cd = zeros(10, stepLength); %The estimate for each cycle

for k = 1:nbrOfCycles
    walk = zeros(N, d);
    omega = zeros(N, stepLength);
    %The rows are the particle and the cols are the variables in that dimension
    histOfCord = zeros(N, d*stepLength + 2);
    for nbrOfSteps = 1:stepLength
        for particle = 1:N
            %For a given particle, we gather its history of walk.
            historyOfWalk = zeros(nbrOfSteps, d);
            for s = 1:d
                historyOfWalk(:, s) = histOfCord(particle, s:d:nbrOfSteps*d)';
            end
            
            for i = 1:d2
                neigh(i, :) = walk(particle, :) + n(i, :);
            end
            
            if ismember(neigh, historyOfWalk,'rows') == ones(d2,1)
                continue
            end
            
            newDir = datasample(n,1);
            next = walk(particle, :) + newDir;
            
            while ismember(next, historyOfWalk, 'rows') == ones(1,d)
                newDir = datasample(n,1);
                next = walk(particle, :) + newDir;
            end
            
            g(particle) = d2-sum(ismember(neigh, historyOfWalk, 'rows'));
            walk(particle, :) = next;
            histOfCord(particle, nbrOfSteps*d+1:nbrOfSteps*d+d) = next;
            omega(particle, nbrOfSteps) = g(particle);
        end
        %Selection
        ind = randsample(N, N, true, omega(:, nbrOfSteps));
        %Mutation
        walk = walk(ind, :);
        histOfCord = histOfCord(ind, :);
    end
    cd(k, :) = mean(omega);
end

for k = 1:10
    for i = 2:nbrOfSteps
        cd(k,i) = cd(k,i-1).*cd(k, i);
    end
end

cd_var = var(cd);
cd_sisr = mean(cd);

m = size(cd);
m = m(1,2);
col1 = ones(m,1);
col2 = (1:m)';
col3 = log(col2);
X = [col1, col2, col3];
y = log(cd);
beta = ((X'*X)\X')*y';

A_d = exp(beta(1,:));
mu_d = exp(beta(2,:));
gam_d = beta(3,:) + 1;

A_d_var = var(A_d);
mu_d_var = var(mu_d);
gam_d_var = var(gam_d);
vars = [A_d_var, mu_d_var, gam_d_var]

A_d_mean = mean(A_d);
mu_d_mean = mean(mu_d);
gam_d_mean = mean(gam_d);
means = [A_d_mean, mu_d_mean, gam_d_mean]

ctest = @(n) A_d_mean.*mu_d_mean.^n.*n.^(gam_d_mean-1);
cEst = ctest(1:10);

%% Plot comparing the asymptotic behaviour of cn^(1/n) with our estimate of mu_2 as n grows large
clear
close all

d = 2; %Dimension
d2 = 2*d;
N = 1e3; %Number of particles
nbrOfCycles = 10; %Number of cycles
stepLength = 50; %The length of the total step, i.e n
g = zeros(N,1); %The instrumental distribution function
n = [eye(d) ; -1*eye(d)]; %The possible directions in each step for a given d
neigh = zeros(d2, d); %The neighbours of a given point
cd = zeros(10, stepLength); %The estimate for each cycle

for k = 1:nbrOfCycles
    walk = zeros(N, d);
    omega = zeros(N, stepLength);
    %The rows are the particle and the cols are the variables in that dimension
    histOfCord = zeros(N, d*stepLength + 2);
    for nbrOfSteps = 1:stepLength
        for particle = 1:N
            %For a given particle, we gather its history of walk.
            historyOfWalk = zeros(nbrOfSteps, d);
            for s = 1:d
                historyOfWalk(:, s) = histOfCord(particle, s:d:nbrOfSteps*d)';
            end
            
            for i = 1:d2
                neigh(i, :) = walk(particle, :) + n(i, :);
            end
            
            if ismember(neigh, historyOfWalk,'rows') == ones(d2,1)
                continue
            end
            
            newDir = datasample(n,1);
            next = walk(particle, :) + newDir;
            
            while ismember(next, historyOfWalk, 'rows') == ones(1,d)
                newDir = datasample(n,1);
                next = walk(particle, :) + newDir;
            end
            
            g(particle) = d2-sum(ismember(neigh, historyOfWalk, 'rows'));
            walk(particle, :) = next;
            histOfCord(particle, nbrOfSteps*d+1:nbrOfSteps*d+d) = next;
            omega(particle, nbrOfSteps) = g(particle);
        end
        %Selection
        ind = randsample(N, N, true, omega(:, nbrOfSteps));
        %Mutation
        walk = walk(ind, :);
        histOfCord = histOfCord(ind, :);
    end
    cd(k, :) = mean(omega);
end

for k = 1:nbrOfCycles
    for i = 2:nbrOfSteps
        cd(k,i) = cd(k,i-1).*cd(k, i);
    end
end

cd = mean(cd);

stepLength = 50; % n = 50
close all
plot(1:stepLength,cd)
plot(1:stepLength,cd.^(1./(1:stepLength)))
hold on
plot(1:50,2.6691*ones(stepLength)) %Our estimate of mu_2
xlabel('Number of steps (n)')
ylabel('\mu_2(n)')
set(gca,'fontsize',13)
legend('Convergence of (c_2)^(1/n) towards \mu_2','Computed \mu_2')

%% Boundary analysis using Graham 
close all
mu_d_bound = @(d) (2.*d-1) - 1./(2.*d) - 3./(2.*d).^2 - 16./(2.*d).^3;
mu_d = [mu_d_bound(3:10)];
mu_d_mine = [4.7172, 8.8321, 18.9219];
figure(2)
plot(3:10,mu_d)
hold on
plot([3,5,10], mu_d_mine,'*')
xlabel('Dimension (d)')
ylabel('\mu(d)')
set(gca,'fontsize', 13)
legend('Graham','Our estimates')