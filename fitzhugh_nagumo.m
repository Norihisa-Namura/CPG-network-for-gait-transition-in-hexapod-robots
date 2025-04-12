classdef fitzhugh_nagumo
    % This class contains information on FitzHugh–Nagumo (FHN) oscillators.
    %
    %
    % properties (principals only):
    %
    % - a: Parameter "a"
    %
    % - b: Parameter "b"
    %
    % - c: Parameter "c"
    %
    % - c: Parameter "d"
    %
    % - dt: Time step for numerical simulations
    %
    %
    % functions: 
    %
    % - fitzhugh-nagumo: constructor for parameter setting
    %
    % - func: Dynamics of FHN oscillator
    %
    % - func_jacobi: Jacobian matrix of FHN oscillator

    properties
        % Parameters
        a = 1/3;
        b = 0.25;
        c = 0.15;
        d = 40;

        dt = 0.001;

        initial = [2;0];
        x_tick = [-2,0,2];
        y_tick = [-2,0,2];
        x_lim = [-2.5,2.5];
        y_lim = [-2,2];
        output_lim = [-4,4];
        name = "fhn";
        evol_basin = 3;
        y_basis = 0;
        dim = 2;
    end
    
    methods
        function obj = fitzhugh_nagumo(a,b,c,d)
            if exist('a','var')
                obj.a = a;
            end
            if exist('b','var')
                obj.b = b;
            end
            if exist('c','var')
                obj.c = c;
            end
            if exist('d','var')
                obj.d = d;
            end
        end
        
        function dxdt = func(obj,x)
            dx1dt = obj.d*(x(1,:) - obj.a * x(1,:).^3 - x(2,:));
            dx2dt = obj.d*((x(1,:) + obj.b) * obj.c);
            dxdt = [dx1dt;dx2dt];
        end
        
        function J = func_jacobi(obj,x1_lc,x2_lc)
            J = obj.d*[1 - 3*obj.a*x1_lc.^2,-1;obj.c,0];
        end
    end
end