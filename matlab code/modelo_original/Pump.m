classdef Pump
    properties
        head_fun
        eff_fun
        pot_fun
    end

    methods
        function obj = Pump(head_fun, eff_fun, pot_fun)
            % Pump setup
            % head_fun: head of BCS in function of frequency and volumetric
            %           flow in CasADi Function format [Hz, m^3/s -> m]
            % eff_fun: efficiency of BCS [Hz, m^3/s -> adm]
            % pot_fun: power of BCS [Hz, m^3/s -> W]
            obj.head_fun = head_fun;
            obj.eff_fun = eff_fun;
            obj.pot_fun = pot_fun;
        end
    end
end
