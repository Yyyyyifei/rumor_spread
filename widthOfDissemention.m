num_intervals = 100;
alpha_list = linspace(0.05, 0.25, num_intervals);
gamma_list = linspace(0.5, 1, num_intervals);
mu_list = linspace(0.02, 0.5, num_intervals);
r_mu_1 = zeros(num_intervals);
r_mu_2 = zeros(num_intervals);
r_mu_3 = zeros(num_intervals);
r_alpha_1 = zeros(num_intervals);
r_alpha_2 = zeros(num_intervals);
r_alpha_3 = zeros(num_intervals);
r_gamma_1 = zeros(num_intervals);
r_gamma_2 = zeros(num_intervals);
r_gamma_3 = zeros(num_intervals);

for n=1:num_intervals
% parameters
    param.N     = 100;         % total population size
    param.alpha  = 25/param.N; % infection rate, 1/beta = typical time between contacts
    param.gamma = 0.5;         % recovery rate, 1/gamma = typical time until recovery
    param.k = 25;
    % final time
    tf = 5;
    
    % initial condition
    x0 = param.N * [0.990;0.010;0.000];
    
    % simulate under different m
    param.mu = mu_list(n);

    [t,x] = ode45( @(t,x) SIR_rhs(t,x,param), [0,tf], x0 );
    R = x(:,3); 
    r_mu_1(n) = R(end);
    
    param.gamma = 1; param.alpha=0.25;
    param.mu = mu_list(n);

    [t,x] = ode45( @(t,x) SIR_rhs(t,x,param), [0,tf], x0 );
    R = x(:,3); 

    r_mu_2(n) = R(end);

    param.gamma = 0.75; param.alpha = 0.1;
    param.mu = mu_list(n);

    [t,x] = ode45( @(t,x) SIR_rhs(t,x,param), [0,tf], x0 );
    R = x(:,3); 

    r_mu_3(n) = R(end);
    % simulate under different alpha
    param.mu=0.1;
    param.gamma = 0.5;
    param.alpha = alpha_list(n);
    [t,x] = ode45( @(t,x) SIR_rhs(t,x,param), [0,tf], x0 );
    S = x(:,1); I = x(:,2); R = x(:,3);
    r_alpha_1(n) = R(end);
    
    param.mu=0.2; param.gamma=1;param.alpha=alpha_list(n);
    [t,x] = ode45( @(t,x) SIR_rhs(t,x,param), [0,tf], x0 );
    S = x(:,1); I = x(:,2); R = x(:,3);
    r_alpha_2(n) = R(end);

    param.mu=0.1; param.gamma=1;param.alpha=alpha_list(n);
    [t,x] = ode45( @(t,x) SIR_rhs(t,x,param), [0,tf], x0 );
    S = x(:,1); I = x(:,2); R = x(:,3);
    r_alpha_3(n) = R(end);

    %simulate under different gamma
    param.alpha=0.1; param.mu=0.1; param.gamma = gamma_list(n);

    [t,x] = ode45( @(t,x) SIR_rhs(t,x,param), [0,tf], x0 );
    S = x(:,1); I = x(:,2); R = x(:,3);
    r_gamma_1(n) = R(end);

    param.alpha=0.1; param.mu=0.2; param.gamma = gamma_list(n);

    [t,x] = ode45( @(t,x) SIR_rhs(t,x,param), [0,tf], x0 );
    S = x(:,1); I = x(:,2); R = x(:,3);
    r_gamma_2(n) = R(end);

    param.alpha=0.3; param.mu=0.2; param.gamma = gamma_list(n);

    [t,x] = ode45( @(t,x) SIR_rhs(t,x,param), [0,tf], x0 );
    S = x(:,1); I = x(:,2); R = x(:,3);
    r_gamma_3(n) = R(end);
end

figure();

title("Width of Dissemenation with respect to parameters")
subplot(3,1,1) 
hold on
P1=plot(alpha_list, r_alpha_1,'color', "r",  LineWidth=0.8);
P2=plot(alpha_list,r_alpha_2, 'color', 'g', LineWidth=0.8);
P3=plot(alpha_list,r_alpha_3, 'color', 'b',  LineWidth=0.8);
legend([P1(1) P2(1) P3(1)],{'mu=0.1, gamma=0.5', 'mu=0.2, gamma=1', 'mu=0.1, gamma=1'})
xlabel('Infectious Rate')
ylabel('Width of Dissemenation')

subplot(3,1,2)
hold on
P1=plot(gamma_list, r_gamma_1, 'r', LineWidth=0.8);
P2=plot(gamma_list, r_gamma_2, 'g', LineWidth=0.8);
P3=plot(gamma_list, r_gamma_3, 'b', LineWidth=0.8);
legend([P1(1) P2(1) P3(1)], {'mu=0.1, alpha=0.1', 'mu=0.2, alpha=0.1', 'mu=0.2, alpha=0.2'})
xlabel('Recover rate')
ylabel('Width of Dissemenation')

subplot(3,1,3)
hold on
P1=plot(mu_list, r_mu_1,  "r", LineWidth=0.8);
P2=plot(mu_list, r_mu_2,  "g", LineWidth=0.8);
P3=plot(mu_list, r_mu_3,  "b", LineWidth=0.8);
legend([P1(1) P2(1) P3(1)], {'alpha=0.25, gamma=0.5', 'alpha=0.25, gamma=1', 'alpha=0.1, gamma=0.75'})
xlabel("Threshold of convincement")
ylabel('Width of Dissemenation')


function dxdt = SIR_rhs(t,x,param)
% right-hand side of the ODE of the
% Kermack-McKendrick model 

    S = x(1); I = x(2); R = x(3); % extract the variables
    
    P = 1/(1+exp(-param.k*(I/param.N-param.mu)));
    
    dxdt = [ -param.alpha * S * I;                   ... dS/dt
              P*param.alpha * S * I - param.gamma * I; ... dI/dt
                                   (1-P)*param.alpha * S * I + param.gamma * I  ... dR/dt
           ];
end