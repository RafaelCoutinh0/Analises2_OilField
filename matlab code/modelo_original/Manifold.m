classdef Manifold
    % Class to describe the manifold model
    % Adapted from:
    % Rasmus (2011) - Automatic Start-up and Control of Artificially Lifted Wells

    properties
        wells
        pipe
        booster
        q_wv
        avrWC
    end

    methods
        function obj = Manifold(pipe, booster, q_wv, avrWC, wells)
            % Manifold setup
            % pipe: Pipe object of manifold
            % booster: Pump object of booster pump
            % q_wv: [%] Initial water cut
            % avrWC: [%] Initial, average water cut
            % wells: cell array with Well objects connected to manifold

            obj.wells = wells;
            obj.pipe = pipe;
            obj.booster = booster;

            % not in use
            obj.q_wv = q_wv;
            obj.avrWC = avrWC;
        end

        function out = model(obj, t, x, z, u)
            % Manifold complete model, including associated wells
            % x: [P_manifold, q_tr, P_fbhp_1, P_choke_1, q_average_1, ...]
            % z: [P_intake_1, dP_BCS_1, ...]
            % u: [f_BP, P_topside, f_ESP_1, alpha_1, ...]
            % return: 2 manifold ODE, n*3 well ODE, n*2 algebraic equations

            % Inputs
            f_BP = u(1);
            p_topside = u(2);

            f_ESP = cell(1, numel(obj.wells));
            alpha = cell(1, numel(obj.wells));
            cont = 2;
            for i = 1:numel(obj.wells)
                f_ESP{i} = u(cont + 1);
                alpha{i} = u(cont + 2);
                cont = cont + 2;
            end

            % States
            p_man = x(1);
            q_tr = x(2);

            P_fbhp = cell(1, numel(obj.wells));
            P_choke = cell(1, numel(obj.wells));
            P_intake = cell(1, numel(obj.wells));
            P_discharge = cell(1, numel(obj.wells)); %#ok<NASGU>
            q_mean = cell(1, numel(obj.wells));
            q_choke = cell(1, numel(obj.wells));
            dot_x = {};
            g_z = {};
            qin_man = 0;

            contx = 2;
            contz = 0;
            for i = 1:numel(obj.wells)
                well = obj.wells{i};

                P_fbhp{i} = x(contx + 1);
                P_choke{i} = x(contx + 2);
                q_mean{i} = x(contx + 3);
                P_intake{i} = z(contz + 1);
                P_discharge{i} = z(contz + 2);

                % for each well
                q_choke{i} = well.flow_rate_choke(P_choke{i} * 1e5, ...
                    p_man * 1e5, alpha{i});

                qin_man = qin_man + q_choke{i};

                [dot_p_fbhp, dot_p_choke, dot_q_ave, g_p_in, g_dp_bcs] = ...
                    well.model(t, ...
                    [x(contx + 1); x(contx + 2); x(contx + 3)], ...
                    [z(contz + 1); z(contz + 2)], ...
                    [f_ESP{i}; alpha{i}], p_man);

                contx = contx + 3;
                contz = contz + 2;

                dot_x = [dot_x, {dot_p_fbhp, dot_p_choke, dot_q_ave}]; %#ok<AGROW>
                g_z = [g_z, {g_p_in, g_dp_bcs}]; %#ok<AGROW>
            end

            fric_manifold = obj.friction(q_tr / 3600, obj.pipe);

            dp_boost = obj.booster.head_fun(f_BP, q_tr / 3600) * obj.pipe.rho * 9.81;

            height = obj.pipe.h * obj.pipe.rho * 9.81;

            dot_p_man = (obj.pipe.beta / obj.pipe.V) * ...
                (qin_man - q_tr / 3600) / 1e5;

            dot_qtr = (obj.pipe.A / obj.pipe.rho / obj.pipe.L) * ...
                (p_man * 1e5 - p_topside - height - fric_manifold + dp_boost) * 3600;

            out = [{dot_p_man, dot_qtr}, dot_x, g_z];
        end

        function fric = friction(~, q, pipe)
            % Friction model of manifold
            Re = (2 * pipe.rho * q) / (pipe.r * pi * pipe.mu);
            fric_factor = 0.36 * Re^(-0.25);
            fric = (pipe.B0 + pipe.B1 * fric_factor) * q^2 * pipe.rho / 2;
        end
    end
end
