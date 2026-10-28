classdef Pipe
    properties
        d
        r
        L
        h
        beta
        rho
        V
        A
        mu
        B0
        B1
    end

    methods
        function obj = Pipe(d, L, h, beta, rho, mu, B0, B1)
            % Define pipe parameters
            %
            % d: pipe's axial section diameter [m]
            % L: pipe's length [m]
            % h: pipe's height [m]
            % beta: fluid's dP/drho [Pa m^3 /kg]
            % rho: fluid density [kg/m^3]
            % mu: fluid's viscosity [Pa s]
            % B0: friction parameter 1
            % B1: friction parameter 2

            obj.d = d;
            obj.r = d / 2;                 % [m] pipe's axial section radius
            obj.L = L;
            obj.h = h;
            obj.beta = beta;
            obj.rho = rho;
            obj.V = pi * obj.r^2 * L;      % [m^3] pipe's volume
            obj.A = pi * obj.r^2;          % [m^2] pipe's axial section area
            obj.mu = mu;
            obj.B0 = B0;
            obj.B1 = B1;
        end
    end
end
