function n = nomes_modelo()
%NOMES_MODELO Convenções fixadas pelo sym_dynamics.m do usuário.
    n.x = {'P_man', 'q_tr', 'P_fbhp_1', 'P_choke_1', ...
        'q_mean_1', 'P_fbhp_2', 'P_choke_2', 'q_mean_2'};
    n.x_labels = {'P_{man}', 'q_{tr}', 'P_{FBHP,1}', 'P_{choke,1}', ...
        'q_{mean,1}', 'P_{FBHP,2}', 'P_{choke,2}', 'q_{mean,2}'};
    n.x_units = {'bar', 'm^3/h', 'bar', 'bar', 'm^3/h', ...
        'bar', 'bar', 'm^3/h'};
    n.z = {'P_intake_1','dP_BCS_1','P_intake_2','dP_BCS_2'};
    n.z_labels = {'P_{intake,1}', '\Delta P_{BCS,1}', ...
        'P_{intake,2}', '\Delta P_{BCS,2}'};
    n.z_units = {'bar','bar','bar','bar'};
    n.u_full = {'f_BP', 'p_topside', 'f_BCS_1', 'alpha_1', ...
        'f_BCS_2', 'alpha_2'};
    n.idx_u = [1 3 4 5 6];  % p_topside não é entrada manipulada
    n.u = n.u_full(n.idx_u);
    n.u_labels = {'f_{BP}', 'f_{BCS,1}', '\alpha_1', ...
        'f_{BCS,2}', '\alpha_2'};
    n.u_units = {'Hz','Hz','-','Hz','-'};
end
