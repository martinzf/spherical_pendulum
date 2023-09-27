function animate(dt,tf,alpha,phi,l)
% Function animates pendulum
    % Time frames
    t = 0:dt:tf;
    % Cartesian coords
    xyz = l*[sin(alpha).*cos(phi);
        sin(alpha).*sin(phi);
        -cos(alpha)];
    % pendulum and wake line objects
    hold on
    pend = animatedline('MaximumNumPoints',2,'Marker','.','MarkerSize',12, ...
        'Color','k');
    wake = animatedline('MaximumNumPoints',5,'LineWidth',2, ...
        'Color',[.8 .8 .8 .5]);
    title('Spherical pendulum')
    % Spherical grid
    [X,Y,Z] = sphere;
    X = l*X;
    Y = l*Y;
    Z = l*Z;
    mesh(X,Y,Z,'FaceAlpha','0')
    hold off
    % Force 3D view
    view(3) 
    % Axis limits
    lim = 1.2*l; 
    xlim("manual")
    ylim("manual")
    zlim("manual")
    xlim([-lim lim])
    ylim([-lim lim])
    zlim([-lim lim])
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    grid on
    % Time text
    txt = annotation('textbox',[.1,.875,.1,.1],'String','t=0');
    for i = 1:length(t)
        addpoints(pend,[0 xyz(1,i)],[0 xyz(2,i)],[0 xyz(3,i)])
        addpoints(wake,xyz(1,i),xyz(2,i),xyz(3,i))
        set(txt,'String',['t = ',num2str(t(i)),'s'])
        pause(dt)
    end
end