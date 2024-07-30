  % parameters
    param.N     = 100;         % total population size
    param.alpha  = 25/param.N; % infection rate, 1/beta = typical time between contacts
    param.gamma = 0.5;         % recovery rate, 1/gamma = typical time until recovery
    param.k = 10;
    param.mu = 0.04;

    % final time
    tf = 5;
    
    % initial condition
    x0 = param.N * [0.990;0.010;0.000];
    
    R0 = param.alpha * x0(1) / (param.gamma*(1+exp(param.k*param.mu)));
    
    % simulate
    [t,x] = ode45( @(t,x) SIR_rhs(t,x,param), [0,tf], x0 );
    S = x(:,1); I = x(:,2); R = x(:,3);
    
    % plot
    figure();
    clrs = [0 0 240; 220 20 60; 60 60 60] / 256;
    
    subplot(2,1,1); hold on; % time courses
    plot(t,S,'Color',clrs(1,:),'LineWidth',3)
    plot(t,I,'Color',clrs(2,:),'LineWidth',3)
    plot(t,R,'Color',clrs(3,:),'LineWidth',3)
    legend({'S','I','R'})
    xlabel('time')
    ylabel('compartments')
    title(['R_0 = ' num2str(R0)])
    
    subplot(2,1,2); hold on; % ternary plot
    plot( 1/2 * (2*R+I) / param.N, sqrt(3)/2 * I / param.N, '-k', 'LineWidth', 2 )
    plot( [0 0.5 1 0], [0 sqrt(3)/2 0 0], '-k', 'LineWidth', 1 )
    set(gca,'XTick',[],'YTick',[])
    axis([-0.1 1.1 -0.1 sqrt(3)/2+0.1])
    axis equal
    box on
    
    text(0.5,sqrt(3)/2,'I=N','VerticalAlignment','bottom');
    text(1,0,'R=N','HorizontalAlignment','left');
    text(0,0,'S=N','HorizontalAlignment','right');

    drawnow

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