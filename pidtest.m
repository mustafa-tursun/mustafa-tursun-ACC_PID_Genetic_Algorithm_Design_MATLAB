function J = pidtest(G,dt,parms)
    s = tf('s');
    K = parms(1) + parms(2)/s + parms(3)*s/(1+.001*s);
    Loop = series(K,G);
    ClosedLoop = feedback(Loop,1+2*s);
    t = 0:dt:20;
    [y,t] = step(ClosedLoop,t);
    %H=1+2*s
    CTRLtf = K/(1+K*G*(1+2*s));
    u = lsim(CTRLtf,1-y,t);
     
    Q = 100;
    R = 0.001;
    J = dt*sum(Q*(1-y(:)).^2+R*u(:).^2)
    

    [y,t] = step(ClosedLoop,t);
    plot(t,y,'LineWidth',2,'color','b')
    % h = findobj(gcf,'type','line');
    set(gca, 'color','k','color', 'w', 'color','w')
    set(gcf, 'color','w')
    grid on
    %set(h, 'linewidth',2, 'color', 'r');
    drawnow
end
