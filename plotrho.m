
% syms s z w Q 
% 
% eqn1 = -(1+w)*(1+s)*Q+w/z*Q - s*(1-s*(1+1/w))+w*s*z*Q+Q^2;
% eqn2 = -(1+w)*(1+s)*Q+w/z*Q - s*(1-s*(1+1/w))+w*s*z*Q+Q^2+(1-s)*Q;
% sol = solve(eqn1 == 0, Q);
% sol2 = solve(eqn2 == 0, Q);
N = 100000;
w = 1.*ones(N,1);
s  = 0.5;
s = linspace(0,1,N)';
% s = (w + 2.*ones(N,1) - 2.*sqrt(w+1))./(w.*s);
    z1 = (w + 2.*ones(N,1) - 2.*sqrt(w+1))./(w.*s);
    z2 = (w + 2.*ones(N,1) + 2.*sqrt(w+1))./(w.*s);
    rho  = (w.*(ones(N,1)-s)-w.*s.*sqrt((z2-ones(N,1)).*(z1-ones(N,1))));
%%
tiledlayout(2,2)
nexttile
    plot(s,z1, 'k', 'LineWidth',3) 
    xlabel('s')
    ylabel('z_1')
    nexttile
    plot(s,z2, 'k', 'LineWidth',3) 
    xlabel('s')
    ylabel('z_2')
nexttile(3,[1 2])
    plot(s,real(rho), 'k', 'LineWidth',3) 
    xlabel('s')
    ylabel('rho')


    %%
tiledlayout(2,2)
nexttile
    plot(s,z1, 'k', 'LineWidth',3) 
    xlabel('s')
    ylabel('z_1')
    nexttile
    plot(s,z2, 'k', 'LineWidth',3) 
    xlabel('s')
    ylabel('z_2')
nexttile(3,[1 2])
    plot(s,real(rho), 'k', 'LineWidth',3) 
    xlabel('s')
    ylabel('rho')


    %%
syms x a b m
taylor(sqrt((z-a)*(z-b)),z)
%%


fun  = @(y) (1-y)./y.*(sqrt((y-a).*(y-b)));
C = [1+1i -1+1i -1-1i 1-1i];
q2 = integral(fun,1,1,'Waypoints',C)