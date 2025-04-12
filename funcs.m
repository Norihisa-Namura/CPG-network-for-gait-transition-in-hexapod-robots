classdef funcs
    methods
    end
    methods (Static)
        function [T,omega,initial_lc] = period(dt,cla)
            % Calculate period from data

            M = 1e8; % Max iteration
            xx = zeros(length(cla.initial),3);
            xx(:,2) = cla.initial;
            xx(:,3) = funcs.rk4(xx(:,2),dt,cla);
            count = 0;
            icount = 0;
            for i = 1:M-1
                xx(:,1:2) = xx(:,2:3);
                xx(:,3) = funcs.rk4(xx(:,2),dt,cla);
                cond_var1 = (xx(2,1) - cla.y_basis) * (xx(2,2) - cla.y_basis);
                if i > icount + 10 && cond_var1 < 0
                    count = count + 1;
                    icount = i;
                    disp(xx(:,3).')

                    if count == 28
                        i_start = i;
                    elseif count == 30 
                        T = (i - i_start) * dt;
                    	break;
                    end
                else
                    T = 0;
                end
            end
            if T == 0
                omega = Inf;
            else
                omega = 2*pi/T;
            end
            initial_lc = xx(:,3);
        end

        function [x_lc,theta_lc] = phase_map(T,dt,initial,cla)
            % Phase map on the limit cycle
            
            omega = 2*pi/T;
            n = round(T/dt);
            x_lc = zeros(length(initial),n);
            theta_lc = zeros(n,1);
            
            x_lc(:,1) = initial;
            
            for i = 1:n-1
                x_lc(:,i+1) = funcs.rk4(x_lc(:,i),dt,cla);
                theta_lc(i+1) = omega*dt*i;
            end
        end
        
        function theta = state_to_phase(x,x_lc,theta_lc,Z)
            % Calculation of asymptotic phase near the limit cycle

            index = dsearchn(x_lc.',x.');
            if exist("Z",'var')
                theta = theta_lc(index) + sum(Z(:,index).*(x - x_lc(:,index)),2);
            else
                theta = theta_lc(index);
            end
        end
        
        function fx = rk4(x,dt,cla)
            % Fourth-order Runge–Kutta method for oscillator dynamics

            ka = cla.func(x);
            kb = cla.func(x + dt*ka/2);
            kc = cla.func(x + dt*kb/2);
            kd = cla.func(x + dt*kc);
            fx = x + dt*(ka + 2*kb + 2*kc + kd)/6;
        end

        function fx = rk4_pcf(phi,dt,Gamma)
            % Fourth-order Runge–Kutta method for the dynamics of phase differences

            L = length(Gamma);
            phi_grid = linspace(0,2*pi,L+1)';
            phi_grid(end) = [];

            index = dsearchn(phi_grid,mod(phi,2*pi));
            ka = Gamma(index);
            index_a = dsearchn(phi_grid,mod(phi + dt*ka/2,2*pi));
            kb = Gamma(index_a);
            index_b = dsearchn(phi_grid,mod(phi + dt*kb/2,2*pi));
            kc = Gamma(index_b);
            index_c = dsearchn(phi_grid,mod(phi + dt*kc,2*pi));
            kd = Gamma(index_c);
            fx = phi + dt*(ka + 2*kb + 2*kc + kd)/6;
        end

        function x_out = rk4_network(x,X_0,H1,H2,dt,s,b,c,cla)
            % Fourth-order Runge–Kutta method for CPG network

            % x: 12 * 1
            x = reshape(x,[2,6]);

            index = dsearchn(X_0.',x.');
            ka = s*cla.func(x) + funcs.network(index,H1,H2,b,c);
            
            x_a = x + dt*ka/2;
            index_a = dsearchn(X_0.',x_a.');
            kb = s*cla.func(x_a) + funcs.network(index_a,H1,H2,b,c);

            x_b = x + dt*kb/2;
            index_b = dsearchn(X_0.',x_b.');
            kc = s*cla.func(x_b) + funcs.network(index_b,H1,H2,b,c);

            x_c = x + dt*kc;
            index_c = dsearchn(X_0.',x_c.');
            kd = s*cla.func(x_c) + funcs.network(index_c,H1,H2,b,c);

            x_out = x + dt*(ka + 2*kb + 2*kc + kd)/6;
            x_out = x_out(:);
        end

        function k = network(index,H1,H2,b,c)
            % Definition of CPG network structure

            k = [b(1)*c(1)*H1(:,index(1),index(4)) + b(2)*c(2)*H2(:,index(1),index(2)),...
                b(1)*c(1)*H1(:,index(2),index(5)) + c(2)*H2(:,index(2),index(1)) + b(2)*c(2)*H2(:,index(2),index(3)),...
                b(1)*c(1)*H1(:,index(3),index(6)) + c(2)*H2(:,index(3),index(2)),...
                c(1)*H1(:,index(4),index(1)) + b(2)*c(2)*H2(:,index(4),index(5)),...
                c(1)*H1(:,index(5),index(2)) + c(2)*H2(:,index(5),index(4)) + b(2)*c(2)*H2(:,index(5),index(6)),...
                c(1)*H1(:,index(6),index(3)) + c(2)*H2(:,index(6),index(5))];
        end

        function [alpha_sim,beta_sim,T_conv,x,theta_direct] = synchronization(x_init,X_0,i_end,alpha,beta,theta,theta_lc,epsilon,H1,H2,Gamma_alpha,Gamma_beta,s,b,c,dt,omega,cla)
            % Oscillator dynamics and phase difference dynamics
            
            alpha_sim = zeros(1,i_end);
            beta_sim = zeros(1,i_end);
            theta_direct = zeros(6,i_end);
            x = zeros(size(x_init,1),i_end);
            alpha_sim(:,1) = alpha(1);
            beta_sim(:,1) = beta(1);
            theta_direct(:,1) = theta;
            x(:,1) = x_init;

            convergence_time = zeros(2,1);
            phi_epsilon = dt*omega; % Threshold for convergence

            alpha_count = 0;
            beta_count = 0;

            H1_in = epsilon*H1;
            H2_in = epsilon*H2;
            Gamma_alpha_in = epsilon*Gamma_alpha;
            Gamma_beta_in = epsilon*Gamma_beta;

            for i = 1:i_end-1
                alpha_sim(i+1) = funcs.rk4_pcf(alpha_sim(i),dt,Gamma_alpha_in);
                beta_sim(i+1) = funcs.rk4_pcf(beta_sim(i),dt,Gamma_beta_in);
                x(:,i+1) = funcs.rk4_network(x(:,i),X_0,H1_in,H2_in,dt,s,b,c,cla);
                theta_direct(:,i+1) = funcs.state_to_phase(reshape(x(:,i+1),[2,6]),X_0,theta_lc);
        
                if alpha_count == 0 && abs(alpha_sim(i) - alpha(2)) < phi_epsilon 
                    convergence_time(1) = i*dt;
                    alpha_count = 1;
                end
                if beta_count == 0 && abs(beta_sim(i) - beta(2)) < phi_epsilon 
                    convergence_time(2) = i*dt;
                    beta_count = 1;
                end
            end

            disp(convergence_time)
            T_conv = max(convergence_time);
        end
    end
end