classdef Choke
    properties
        k_choke
        choke_fun
    end

    methods
        function obj = Choke(k_choke, choke_fun)
            % Choke setup
            % k_choke: valve constant [m^3/s/Pa^0.5]
            % choke_fun: [adm -> adm] valve opening function in CasADi Function format
            obj.k_choke = k_choke;
            obj.choke_fun = choke_fun;
        end

        function value = model(obj, opening)
            % Multiplication of k_choke and choke characteristic as function
            % of valve opening. Equation needed to evaluate the flow rate
            % through the valve.
            % opening: valve opening (0 to 1)
            % return: k_choke*choke characteristic
            value = obj.k_choke * obj.choke_fun(opening);
        end
    end
end
