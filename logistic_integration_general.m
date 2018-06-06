function [r_out, t_out] = logistic_integration_general( stim, p, q, tonic_activation )

    if nargin<3
        tonic_activation = zeros(p.N,1);
    end

    stim_vec = zeros(p.N,1);
    if length(stim) == 1
        stim_vec(1) = stim;
    else
        stim_vec = stim;
    end
       
    options = odeset('NonNegative',1:p.N);
    [t_out, r_out] = ode45(@(t,y)evolution_eq(t,y,stim_vec,p,q,tonic_activation),[0 p.tmax],zeros(p.N,1),options);

end

function dydt = evolution_eq( t, y, stim, p, q, tonic_activation)

    dydt = diag(p.tauinv) * ( -p.V0 - y + stim * (t>p.t0) * (t<p.t1) +   logistic( p.A*y + tonic_activation, q ));

end

function g = logistic( y, q )
    g = q.rel_strength ./ (1 + exp( -q.k .* (y - q.yh) ) ) - 1 ./ (1 + exp( -q.k .* -q.yh) );
end