function Z = PSF(T,dt,X_0,cla)
    % This function calculates the PSF of two-dimensional limit-cycle oscillators.
    %
    %
    % input: 
    %
    % - T: Period of oscillation
    %
    % - dt: Time step for numerical simualtions
    %
    % - X_0: Limit-cycle trajectory
    %
    % - cla: Class of a limit-cycle oscillator
    %
    % output: 
    %
    % - Z: PSF

    Tnum = round(T/dt);
    omega = 2*pi/T;

    disp("calc PSF")
    Z_ = Calc_Z(X_0, dt, Tnum, cla); 
    Z = omega*Z_;
end

%% functions
function Z = Jacobi_evol_adj(X, dt, X_half_next, X_next, Z, cla)
    % 4-th order Runge-Kutta method for transposed Jacobian matrix

    k1 = dt * (cla.func_jacobi(X(1,:),X(2,:))'*Z);
    k2 = dt * (cla.func_jacobi(X_half_next(1,:),X_half_next(2,:)))'*(Z+k1/2);
    k3 = dt * (cla.func_jacobi(X_half_next(1,:),X_half_next(2,:)))'*(Z+k2/2);
    k4 = dt * (cla.func_jacobi(X_next(1,:),X_next(2,:)))'*(Z+k3);
    Z = Z + (k1+2*k2+2*k3+k4)/6;
end

function Z_ = Calc_Z(X_0, dt, Tnum,cla)
    Z_ = zeros(size(X_0, 1),Tnum);
    Z = ones(size(X_0, 1),1); 

    for rep = 1:3
        for tt = 1:Tnum
            Z = evol_Z(Z,X_0,Tnum,tt,dt,cla);
        end
    end

    % Logging
    for tt = 1:Tnum
        Z_(:,Tnum-tt+1) = Z;
        Z = evol_Z(Z,X_0,Tnum,tt,dt,cla);
    end
end

function Z = evol_Z(Z, X_0, Tnum, tt, dt, cla)
    X = X_0(:,Tnum-tt+1);
    h = -1/2*dt;
    X_half_next = funcs.rk4(X,h,cla);
    X_next = X_0(:,mod((Tnum-tt-1),Tnum)+1);
    Z = Jacobi_evol_adj(X, dt, X_half_next, X_next, Z, cla);
    f = cla.func(X_next);
    prob = Z'* f; 
    Z = Z / prob; % Normalization condition
end
