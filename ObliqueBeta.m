function [Beta, err]=ObliqueBeta(Theta,M, gamma)
%%OBLIQUEBETA Oblique Beta Function by Joshua McMillin
% Elements of Gas Turbine Propulsion, J.D. Mattingly, 1996
% Function uses a quadratic fit on 4 points iteratively to converge on a
% solution.  Calls the ObliqueTheta Function for initial assumed Beta
% values at given Mach and gamma.
% Beta and Theta are expressed in Radians
% err(1) is calculated theta - input theta
% err(2) is estimated error in Beta due to err(1)
    err_tol = 1e-10;
    y(1) = pi()/3; % 60 deg
    y(2) = pi()/4; % 45 deg
    y(3) = pi()/5; % 36 deg
    
    x(1) = ObliqueTheta(y(1),M,gamma);
    x(2) = ObliqueTheta(y(2),M,gamma);
    x(3) = ObliqueTheta(y(3),M,gamma);
    
    f=polyfit(x(:),y(:),2);
    y(4) = polyval(f,Theta);
    x(4) = ObliqueTheta(y(4),M,gamma);
    err(1) = x(4)-Theta;
    
    while abs(err(1)) > err_tol
       [S, r] = sortrows([y(:) x(:)],'descend');
       i=find(r==4);
       switch i
           case {1, 2}
               y=S(1:3,1)';
               x=S(1:3,2)';
           case {3,4}
               y=S(2:4,1)';
               x=S(2:4,2)';
       end
       f=polyfit(x(:),y(:),2);
       y(4) = polyval(f,Theta);
       x(4) = ObliqueTheta(y(4),M,gamma);
       err(1) = x(4)-Theta;
    end
    Beta=y(4);
    err(2) = (x(4)-Theta)*(2*f(1)*x(4) + f(2));
end
