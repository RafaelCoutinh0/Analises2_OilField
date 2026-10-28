classdef Well
    properties
        pipe1
        pipe2
        bcs
        choke
        PI
        P_reservoir
    end

    methods
        function obj = Well(pipe1, pipe2, bcs, choke, PI, P_reservoir)
            % Well setup
            % pipe1: first section of well (before BCS)
            % pipe2: second section of well (after BCS)
            % bcs: pump object
            % choke: choke object
            % PI: production index [m^3/s/Pa]
            % P_reservoir: reservoir pressure [Pa]
            obj.pipe1 = pipe1;
            obj.pipe2 = pipe2;
            obj.bcs = bcs;
            obj.choke = choke;
            obj.PI = PI;
            obj.P_reservoir = P_reservoir;
        end

        function [dpfbdt, dpchokedt, dqaveragedt, g_pin, g_bcs] = model(obj, t, x, z, u, bound_conditions) %#ok<INUSD>
            % Differential-algebraic system to describe Well dynamics
            %
            % t: time [s]
            % x: [P_fbhp, P_choke, q_average]
            % z: [P_intake_BCS, dP_BCS]
            % u: [f_BCS, choke_opening]
            % bound_conditions: P_manifold
            % return: 3 ODE, 2 algebraic equations

            P_fbhp = x(1);
            P_choke = x(2);
            q_average = x(3);
            P_intake = z(1);
            bcs_dP = z(2);

            frequency = u(1);
            opening = u(2);

            P_manifold = bound_conditions;

            q_reservoir = obj.PI * (obj.P_reservoir - P_fbhp * 1e5);

            friction_1 = obj.friction(q_average / 3600, obj.pipe1);
            friction_2 = obj.friction(q_average / 3600, obj.pipe2);

            q_choke = obj.flow_rate_choke(P_choke * 1e5, P_manifold * 1e5, opening);

            [l_bar, r_bar, A_bar, rho_bar] = obj.average_properties(); %#ok<ASGLU>

            dp = obj.bcs.head_fun(frequency, q_average / 3600) * 9.81 * rho_bar;

            height_1 = obj.pipe1.h * obj.pipe1.rho * 9.81;
            height_2 = obj.pipe2.h * obj.pipe2.rho * 9.81;

            % ODE's
            dpfbdt = obj.pipe1.beta / obj.pipe1.V * ...
                (q_reservoir - q_average / 3600) / 1e5;

            dpchokedt = obj.pipe2.beta / obj.pipe2.V * ...
                (q_average / 3600 - q_choke) / 1e5;

            dqaveragedt = A_bar / rho_bar / l_bar * ...
                (P_fbhp * 1e5 - P_choke * 1e5 - friction_1 - friction_2 ...
                - height_1 - height_2 + dp) * 3600;

            % ALGEBRAIC
            P_discharge = P_choke * 1e5 + friction_2 + height_2;

            g_pin = P_fbhp - (friction_1 + height_1) / 1e5 - P_intake;

            g_bcs = bcs_dP - (P_discharge / 1e5 - P_intake);
        end

        function fric = friction(~, q, pipe)
            % Evaluation of pipe friction
            % q: flow [m^3/s]
            % pipe: Pipe object
            % return: pressure variation due to friction in pipe
            Re = (2 * pipe.rho * q) / (pipe.r * pi * pipe.mu);
            fric_factor = 0.36 * Re^(-0.25);
            fric = (pipe.B0 + pipe.B1 * fric_factor) * q^2 * pipe.rho / 2;
        end

        function [l_bar, r_bar, A_bar, rho_bar] = average_properties(obj)
            % Average properties of Well, evaluated from pipes
            l_bar = (obj.pipe1.L + obj.pipe2.L) / 2;
            r_bar = (obj.pipe1.r * obj.pipe1.L + obj.pipe2.r * obj.pipe2.L) * 2 / l_bar;
            A_bar = pi * r_bar^2;
            rho_bar = (obj.pipe1.rho * obj.pipe1.V + obj.pipe2.rho * obj.pipe2.V) / ...
                (obj.pipe1.V + obj.pipe2.V);
        end

        function q = flow_rate_choke(obj, p_in, p_out, opening)
            % Flow rate through choke
            % p_in: pressure in valve intake [Pa]
            % p_out: pressure in valve discharge [Pa]
            % opening: choke opening (0 to 1)
            % return: flow [m^3/s]
            import casadi.*
            choke_factor = obj.choke.model(opening);
            delta_p = p_in - p_out;
            q = sign(delta_p) * choke_factor * abs(delta_p)^0.5;
        end
    end
end
