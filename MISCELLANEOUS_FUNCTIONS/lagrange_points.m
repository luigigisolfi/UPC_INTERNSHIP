%----------------------------------------------------------------------------------
% When called, this function plots the 5 lagrange points for a system with given mu
%----------------------------------------------------------------------------------
function lagrange_points(mu)

    function y = collinear_lagrange(xstar)
    r_s1= (xstar- mu); %distance to primary 
    r_s2= (xstar-mu+1); %distance to secondary 
    rho_s13= norm(r_s1)^3;
    rho_s23=norm(r_s2)^3;
    y = xstar- ((1-mu)*r_s1/rho_s13 + mu*r_s2/rho_s23);
    end

    options = optimset('Display','iter');
    L_2 = fzero(@collinear_lagrange,[-1, -1.5],options);
    L_1 = fzero(@collinear_lagrange,[-0.5, -0.97],options);
    L_3 = fzero(@collinear_lagrange,[0.5, 3],options);
    fprintf('L_1=%f, L_2=%f, L_3=%f\n', L_1, L_2, L_3)
    
    function output
        figure()
        title(sprintf('Lagrange Points for mu = %d', mu))
        hold on
        ylim([-1.2,1.2])
        xlim([-1.2,1.2])
        axis equal
        xlabel('x')
        ylabel('y')
        plot(mu ,0,'bo','MarkerFaceColor','b')
        plot(mu-1,0,'go','MarkerFaceColor','g')
        plot(L_1,0,'rv','MarkerFaceColor','r')
        plot(L_2,0,'r^','MarkerFaceColor','r')
        plot(L_3,0,'rp','MarkerFaceColor','r')
        plot(0.5 + mu- 1,sqrt(3)/2,'rX','MarkerFaceColor','r')
        plot(0.5 + mu- 1,-sqrt(3)/2,'rs','MarkerFaceColor','r')
    end
    output
end