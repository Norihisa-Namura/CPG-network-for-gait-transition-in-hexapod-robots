%%
clear
close all
mystyle
rng('default')

%% Settings
class_name = @fitzhugh_nagumo;
cla = class_name();
    
% 1: wave, 2: tetrapod, 3: tripod
gait_name = ["wave","tetrapod","tripod"];
gaits = [1,2,3];
%gaits = [3,2,1];

speed = [1/6,1/3,1/2];

save_ind = 1;
savename = gait_name(gaits(1)) + "_" + gait_name(gaits(2)) + "_" + gait_name(gaits(3));

%% Limit cycle and PSF
dt = cla.dt;
[T,omega,initial_lc] = funcs.period(cla.dt,cla);
[X_0,theta_lc] = funcs.phase_map(T,dt,initial_lc,cla);

L = length(theta_lc); 

Z = PSF(T,dt,X_0,cla);
[threshold,threshold_state]= obtain_threshold(X_0(1,:),X_0);

Z_norm = sqrt(mean(sum(Z.^2,1)));
dtheta = omega*dt;

%% Parameters
alpha_star = [1/2,1/3,1/2]*2*pi;
beta_star = [1/6,1/3,1/2]*2*pi;

epsilon = 0.1;
c = [4,8];
s = speed(gaits(1));

init_time = 5;
trans_time = 10;
over_time = 5;

%% Time evolution
for trans = 1:2
    func_type = gait(gaits(trans:trans+1));
    b = binary(func_type);
     
    alpha = alpha_star(gaits(trans:trans+1)); % [start,end]
    beta = beta_star(gaits(trans:trans+1)); % [start,end]

    
    phi = linspace(0,2*pi,L+1)';
    M = L/2;
    
    % Even and odd function
    Gamma_odd = 0;
    for k = 1:10
        Gamma_odd = Gamma_odd + k * exp(-k^2/2) * sin(k*phi);
    end
    Gamma_odd = 10*Gamma_odd;
    Gamma_even = 2*cos(phi)+1;

    phi_star = [pi,2*pi/3];
    phi_bar = [0,-2*pi/3];

    % PCF for alpha and beta
    if func_type(1) == 2
        Gamma1 = Gamma_even;
    elseif func_type(1) == 1
        Gamma1 = Gamma_odd;
    end
    
    if func_type(2) == 2
        if beta(2) == 2*pi/3
            Gamma2 = Gamma_even;
        elseif beta(2) == pi/3
            Gamma2 = -[Gamma_even(M+2:end);Gamma_even(1:M+1)];
        end
    elseif func_type(2) == 1
        Gamma2 = Gamma_odd;
    end
    
    Gamma_alpha = c(1)*(Gamma1 - b(1)*Gamma1(end:-1:1));
    Gamma_beta = c(2)*Gamma2;
    Gamma1(end) = [];
    Gamma2(end) = [];
    Gamma_alpha(end) = [];
    Gamma_beta(end) = [];
    phi = linspace(-pi,pi,L+1)';
    phi(end) = [];

    y_labels = ["$\dot{\alpha}$","$\dot{\beta}$"];
    show_pcf_transition(phi,Gamma_alpha,Gamma_beta,M,alpha,beta,y_labels)
    
    % MCF
    H1 = zeros(2,L,L);
    H2 = zeros(2,L,L);
    
    for k = 1:L
        for l = 1:L
            H1(:,k,mod(k-l,L)+1) = Z(:,k)*Gamma1(l)/(Z_norm^2);
            H2(:,k,mod(k-l,L)+1) = Z(:,k)*Gamma2(l)/(Z_norm^2);
        end
    end

    % Time evolution for gait transition
    if trans == 1
        disp("initial gait")
        
        % LF,LM,LH,RF,RM,RH
        shift_init = funcs.state_to_phase(threshold_state(:,gaits(1)),X_0,theta_lc);
        theta_init = mod([0;beta(1);2*beta(1);alpha(1);alpha(1)+beta(1);alpha(1)+2*beta(1)]+shift_init+speed(gaits(1))*2*pi,2*pi);
        index_init = phase2index(theta_init,dtheta);
        init_state = [X_0(:,index_init(1));X_0(:,index_init(2));X_0(:,index_init(3));X_0(:,index_init(4));X_0(:,index_init(5));X_0(:,index_init(6))];
        
        i_end = round(init_time/dt) + 1;
        
        [alpha_init_sim,beta_init_sim,~,x_init,theta_init_direct] =...
            funcs.synchronization(init_state,X_0,i_end,alpha,beta,theta_init,theta_lc,0,H1,H2,Gamma_alpha,Gamma_beta,s,b,c,dt,omega,cla);

    end

    disp("gait transition" + num2str(trans))
    s = speed(gaits(trans+1));
    
    if trans == 1
        i_end = round(trans_time/dt) + 1;
        [alpha_trans_sim1,beta_trans_sim1,T_conv1,x_trans1,theta_trans_direct1] =...
            funcs.synchronization(x_init(:,end),X_0,i_end,alpha,beta,theta_init_direct(:,end),theta_lc,epsilon,H1,H2,Gamma_alpha,Gamma_beta,s,b,c,dt,omega,cla);
    elseif trans == 2
        i_end = round((trans_time + over_time)/dt) + 1;
        [alpha_trans_sim2,beta_trans_sim2,T_conv2,x_trans2,theta_trans_direct2] =...
            funcs.synchronization(x_trans1(:,end),X_0,i_end,alpha,beta,theta_trans_direct1(:,end),theta_lc,epsilon,H1,H2,Gamma_alpha,Gamma_beta,s,b,c,dt,omega,cla);
    end
end

%% Concatenate (t, alpha, beta, theta)
sim_index = 1:(init_time + 2*trans_time)/dt+1;

t_init_sim = 0:dt:init_time;
t_trans_sim1 = init_time + (0:dt:trans_time);
t_trans_sim2 = init_time + trans_time + (0:dt:trans_time+over_time);
t_sim = [t_init_sim(1:end-1),t_trans_sim1(1:end-1),t_trans_sim2];
t_sim = t_sim(sim_index);

x = [x_init(:,1:end-1),x_trans1(:,1:end-1),x_trans2];
x_out = x([1,3,5,7,9,11],:);

theta_direct = [theta_init_direct(:,1:end-1),theta_trans_direct1(:,1:end-1),theta_trans_direct2];
alpha_sim = [alpha_init_sim(:,1:end-1),alpha_trans_sim1(:,1:end-1),alpha_trans_sim2];
beta_sim = [beta_init_sim(:,1:end-1),beta_trans_sim1(:,1:end-1),beta_trans_sim2];
theta_direct = theta_direct(:,sim_index);
alpha_sim = alpha_sim(sim_index);
beta_sim = beta_sim(sim_index);

alpha_direct = [theta_direct(4,:) - theta_direct(1,:);theta_direct(5,:) - theta_direct(2,:);theta_direct(6,:) - theta_direct(3,:)];
beta_direct = [theta_direct(2,:) - theta_direct(1,:);theta_direct(3,:) - theta_direct(2,:);...
    theta_direct(5,:) - theta_direct(4,:);theta_direct(6,:) - theta_direct(5,:)];

clearvars alpha_init_sim alpha_trans_sim1 alpha_trans_sim2
clearvars beta_init_sim beta_trans_sim1 beta_trans_sim2
clearvars theta_init_direct theta_trans_direct1 theta_trans_direct2
clearvars x_init x_trans1 x_trans2

%% Threshold for swing/stance
kappa = 1;
sigma1 = ones(size(t_init_sim))*threshold(gaits(1));
sigma2 = threshold(gaits(1)) + (threshold(gaits(2)) - threshold(gaits(1))) * (1 - exp(-kappa*(t_trans_sim1(2:end) - t_trans_sim1(1))));
sigma3 = threshold(gaits(2)) + (threshold(gaits(3)) - threshold(gaits(2))) * (1 - exp(-kappa*(t_trans_sim2(2:end) - t_trans_sim2(1))));
sigma = [sigma1,sigma2,sigma3];
clearvars sigma1 sigma2 sigma3

%% Reference joint angles
ref_angle1 = reference_trajectory(x_out(1,:),sigma);
ref_angle2 = reference_trajectory(x_out(2,:),sigma);
ref_angle3 = reference_trajectory(x_out(3,:),sigma);
ref_angle4 = reference_trajectory(x_out(4,:),sigma);
ref_angle5 = reference_trajectory(x_out(5,:),sigma);
ref_angle6 = reference_trajectory(x_out(6,:),sigma);

ref_angle = [ref_angle4(:,sim_index);
    ref_angle5(:,sim_index);
    ref_angle6(:,sim_index);
    -ref_angle1(:,sim_index);
    -ref_angle2(:,sim_index);
    -ref_angle3(:,sim_index)];

save("data/" + cla.name + "_ref_angle","ref_angle");
clearvars ref_angle1 ref_angle2 ref_angle3 ref_angle4 ref_angle5 ref_angle6

%%
x = x(:,sim_index);
x_out = x_out(:,sim_index);
sigma = sigma(sim_index);
t_trans_sim2 = init_time + trans_time + (0:dt:trans_time);

%% Save data for matlab
if save_ind == 1
    close all
    save("data/" + cla.name + "_" + savename);
end

%% Load data
%load("data/" + cla.name + "_" + program_name); 

%% show figures
%% LC and PSF
fig = utils.show_lc_psf(theta_lc,X_0,Z,cla);
utils.save_fig(save_ind,fig,"lc_psf",cla);

%% PCF
fig = utils.show_pcf2(linspace(-pi,pi,L+1),[Gamma_odd(M+1:end-1);Gamma_odd(1:M+1)],[Gamma_even(M+1:end-1);Gamma_even(1:M+1)],phi_star,phi_bar);
utils.save_fig(save_ind,fig,"odd_even");

%% Gait results
show_results;
utils.save_fig(save_ind,fig,gait_name(gaits(1))+"_"+gait_name(gaits(2))+"_"+gait_name(gaits(3)),cla);

%%
function index = phase2index(theta,dtheta)
    index = round(mod(theta,2*pi)/dtheta)+1;
    L = round(2*pi/dtheta);
    index(index > L) = 1;
end

function [sigma,state] = obtain_threshold(x_out,X_0)
    L = length(x_out);
    M = round(L/2);
    index = 1:L;
    [~,i_min] = min(x_out);

    % min -> max -> min
    x_out = x_out([i_min:L,1:i_min-1]);
    X_0 = X_0(:,[i_min:L,1:i_min-1]);
    [~,i_max] = max(x_out);

    if i_max - 1 < M
        i_query = 2:i_max-1;
    else
        i_query = i_max+1:L;
    end

    i_others = setdiff(index,i_query);
    x_query = x_out(i_query);
    x_others = x_out(i_others);

    i_nn = dsearchn(x_others',x_query');
    i_nn_orig = i_others(i_nn);
    i_diff = abs(i_nn_orig - i_query);

    sigma = [0,0,0];
    i_wave = dsearchn(i_diff',round(L/6));
    i_tetrapod = dsearchn(i_diff',round(L/3));
    i_tripod = dsearchn(i_diff',round(L/2));
    
    sigma(1) = x_out(i_nn_orig(i_wave));
    sigma(2) = x_out(i_nn_orig(i_tetrapod));
    sigma(3) = x_out(i_nn_orig(i_tripod));

    state = zeros(2,3);
    state(:,1) = X_0(:,i_nn_orig(i_wave));
    state(:,2) = X_0(:,i_nn_orig(i_tetrapod));
    state(:,3) = X_0(:,i_nn_orig(i_tripod));
end

function func_type = gait(gaits)
    switch(gaits(2))
        case 1
          func_type = [1,2]; 
        case 2
          func_type = [2,2];
        case 3
          func_type = [1,1]; 
    end
end

function b = binary(func_type)
    b = zeros(1,2);
    
    if func_type(1) == 2
        b(1) = -1;
    elseif func_type(1) == 1
        b(1) = 1;
    end
    if func_type(2) == 2
        b(2) = -1;
    elseif func_type(2) == 1
        b(2) = 1;
    end
end

function show_pcf_transition(phi,Gamma1,Gamma2,M,alpha,beta,y_labels)
    figure()
    subplot(2,1,1)
    plot(phi,[Gamma1(M+1:end);Gamma1(1:M)])
    hold on
    scatter(alpha,0,1000,'.')
    xlim([-pi,pi])
    xticks([0,1,2,3,4,5,6]*pi/3-pi)
    xticklabels(["$-\pi$","$-\frac{2}{3}\pi$","$-\frac{1}{3}\pi$","$0$","$\frac{1}{3}\pi$","$\frac{2}{3}\pi$","$\pi$"])
    ylabel(y_labels(1))

    subplot(2,1,2)
    plot(phi,[Gamma2(M+1:end);Gamma2(1:M)])
    hold on
    scatter(beta,0,1000,'.')
    xlim([-pi,pi])
    xticks([0,1,2,3,4,5,6]*pi/3-pi)
    xticklabels(["$-\pi$","$-\frac{2}{3}\pi$","$-\frac{1}{3}\pi$","$0$","$\frac{1}{3}\pi$","$\frac{2}{3}\pi$","$\pi$"])
    ylabel(y_labels(2))
end
