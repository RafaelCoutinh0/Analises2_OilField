function lin = linearizar_dae(mdl)
%LINEARIZAR_DAE Elimina algebricas e preserva a matriz de entrada original.
%  Modo principal: DESVIO DIMENSIONAL, como nos notebooks de referência.
%  Modo alternativo: desvio percentual, como em linearizacao_ss_ft.m.
    import casadi.*
    n = nomes_modelo();
    f = mdl.f; g = mdl.g; x = mdl.x; z = mdl.z; u = mdl.u;
    Jfun = Function('jac_bcs_completa', {x,z,u}, { ...
        jacobian(f,x),jacobian(f,z),jacobian(f,u), ...
        jacobian(g,x),jacobian(g,z),jacobian(g,u)});
    [fx, fz, fu, gx, gz, gu] = Jfun(mdl.x_ss,mdl.z_ss,mdl.u0);
    fx = full(fx); fz = full(fz); fu = full(fu);
    gx = full(gx); gz = full(gz); gu = full(gu);
    if rcond(gz)<1e-12
        error('BCS:DAE', 'A Jacobiana G_z é singular ou mal condicionada.');
    end

    Zx = -(gz\gx); Zu_total = -(gz\gu);
    A_dim = fx+fz*Zx;
    B_total = fu+fz*Zu_total;
    idx = n.idx_u;
    B_dim = B_total(:,idx);
    C_dim = eye(mdl.nx);
    D_dim = zeros(mdl.nx,numel(idx));

    % Também possibilita saídas algébricas, linearizadas pela restrição g=0.
    % dz = Zx*dx + Zu*du; não devem ser confundidas com estados dinâmicos.
    C_total = [C_dim;Zx];
    D_total = [D_dim; Zu_total(:,idx)];

    % De desvio dimensional para desvio relativo; multiplicar por 100, se quiser %.
    xs = full(mdl.x_ss(:));
    us = full(mdl.u0(idx));
    Tx = diag(xs); Tu = diag(us);
    if any(abs(xs)<eps) || any(abs(us)<eps)
        warning('Ponto de operação tem componente nula: escala percentual indefinida.');
        A_pct = []; B_pct = []; C_pct = []; D_pct = [];
    else
        A_pct = Tx\(A_dim*Tx);
        B_pct = Tx\(B_dim*Tu);
        C_pct = eye(mdl.nx);
        D_pct = zeros(mdl.nx,numel(idx));
    end
    lin = struct('A',A_dim,'B',B_dim,'C',C_dim,'D',D_dim, ...
        'A_pct',A_pct,'B_pct',B_pct,'C_pct',C_pct,'D_pct',D_pct, ...
        'Zx',Zx,'Zu',Zu_total(:,idx),'Zu_total',Zu_total, ...
        'C_total',C_total,'D_total',D_total,'Gz',gz, ...
        'x_ss',xs,'z_ss',full(mdl.z_ss(:)),'u0',full(mdl.u0(:)), ...
        'u_ss',us,'idx_u',idx,'nomes',n, ...
        'f_x',fx,'f_z',fz,'f_u',fu,'g_x',gx,'g_u',gu);
end
