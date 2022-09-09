%% Linear 1D Advection equation
% Routine for solving the following PDE
%
%    du/dt + v*du/dx = 0
%
% Using the Method of Lines
clc; clear; close all; 

%% Parameters
v=1.25;
d=1.0;
zl=1.0;
n=21;
dz=zl/(n-1);

%% Select three point, finite differencing (FD) of spatial derivative
% ifd = 1: Centered approximations
% ifd = 2: Two point upwind approximation
% ifd = 3: Five point, biased upwind approximation
% ifd = 4: van Leer flux limiter
ifd=2;

%% Initial Condition
% Initial, final times, integration interval, number of Euler
% steps for each output
t=0.0; tf=100.0; h=0.0001; nout=100;
t0=0;
t_in_Start=0.25;       % t_in_Start=pipeLength/m_flow_Start*rho*dh^2/4*pi;
x0=0;                      % x��ʼλ����0
pipe_Length=39;     

% Initial conditions
for i=1:n
    Tin(i)=t0;
end
x=zl;
%% x��v�Ĺ�ϵ
% dx/dt=f(t,x)=v
% x(t0)=x0;
f=v;

%% Integrate until t = tf
    while t < tf*1.0001
    % Monitor solution by displaying t
        current_t=t;
        current_Length=x+x0;
    % Write selected output, including analytical solution and difference
    % between numerical and analytical solutions
        
        if current_Length-pipe_Length >1e-6
            time_out=current_t;
            break;
        end
    
        % Take nout Euler steps
        for iout=1:nout
            % Temporal derivatives
            
            % Centered approximations
            if ifd==1
            % Boundary condition at z = 0
                Tin(1)=t0+t_in_Start;
            % Spatial derivative
                [Tin_dz]=dss002(0.0,zl,n,Tin);
            % Temporal derivative
                Tin_dt(1)=t0;
                for i=2:n
                    Tin_dt(i)=-v*Tin_dz(i);
                end
            % End of three point centered approximation
            end
            
            % Two point upwind approximation
            if(ifd==2)
                Tin(1)=t0+t_in_Start;
                [Tin_dz]=dss012(0.0,zl,n,Tin,v);
                Tin_dt(1)=t0;
                for i=2:n
                    Tin_dt(i)=-v*Tin_dz(i);
                end
            end
            
        % Take Euler step   //���currentlength��currenttime����������ѭ������nout�йأ�������h�й�
            for i=1:n
                Tin(i)=Tin(i)+Tin_dt(i)*h;      % h��С����
            end
            x=x+f*h;         % ���ô󲽳�
            
            t=t+h;
            zl=x;               % ���������Ƕ�̬��
        % Next Euler step    
        end
    end
    
%% Centered approximations

function [ux]=dss002(xl,xu,n,u)
%...
%... FUNCTION DSS002 COMPUTES THE FIRST DERIVATIVE, U , OF A
%...                                                                                        X
%... VARIABLE U OVER THE SPATIAL DOMAIN XL LE X LE XU
%...
%... ARGUMENT LIST
%...
%... XL LOWER BOUNDARY VALUE OF X (INPUT)
%...
%... XU UPPER BOUNDARY VALUE OF X (INPUT)
%...
%... N NUMBER OF GRID POINTS IN THE X DOMAIN INCLUDING THE
%... BOUNDARY POINTS (INPUT)
%...
%... U ONE-DIMENSIONAL ARRAY CONTAINING THE VALUES OF U AT
%... THE N GRID POINT POINTS FOR WHICH THE DERIVATIVE IS
%... TO BE COMPUTED (INPUT)
%...
%... UX ONE-DIMENSIONAL ARRAY CONTAINING THE NUMERICAL
%... VALUES OF THE DERIVATIVES OF U AT THE N GRID POINTS
%... (OUTPUT)
%...
%... THE WEIGHTING COEFFICIENTS CAN BE SUMMARIZED AS
%...
%...        -3 4 -1
%...
%... 1/2  -1 0  1
%...
%...         1 -4 3
%...
%... WHICH ARE THE COEFFICIENTS REPORTED BY BICKLEY FOR N = 2, M =
%... 1, P = 0, 1, 2 (BICKLEY, W. G., FORMULAE FOR NUMERICAL DIFFER-
%... ENTIATION, MATH. GAZ., VOL. 25, 1941).
%...
%... COMPUTE THE SPATIAL INCREMENT
        dx=(xu-xl)/(n-1);
        r2fdx=1./(2.*dx);
        nm1=n-1;
%...
%... LEFT END POINT. THE FOLLOWING CODING HAS BEEN FORMATTED
% SO THAT THE NUMERICAL WEIGHTING COEFFICIENTS CAN BE MORE
%... EASILY ASSOCIATED WITH THE BICKLEY MATRIX LISTED ABOVE)
        ux(1)=r2fdx*( -3.*u( 1) +4.*u( 2) -1.*u( 3));
%...
%... INTERIOR POINTS
        for i=2:nm1
            ux(i)=r2fdx*( -1.*u(i-1) +0.*u( i) +1.*u(i+1));
        end
%...
%... RIGHT END POINT
        ux(n)=r2fdx*( 1.*u(n-2) -4.*u(n-1) +3.*u(n));
end

%% upwind FD
function [ux]=dss012(xl,xu,n,u,v)
%...
%... FUNCTION DSS012 IS AN APPLICATION OF FIRST-ORDER DIRECTIONAL
%... DIFFERENCING IN THE NUMERICAL METHOD OF LINES. IT IS INTENDED
%... SPECIFICALLY FOR THE ANALYSIS OF CONVECTIVE SYSTEMS MODELLED BY
%... FIRST-ORDER HYPERBOLIC PARTIAL DIFFERENTIAL EQUATIONS WITH THE
%... SIMPLEST FORM
%...
%... U + V*U = 0                                                    (1)
%...    T        X
%...
%... THE FIRST FOUR PARAMETERS, XL, XU, N AND U, ARE THE SAME AS
%... FOR FUNCTION DSS002. THE FIFTH PARAMETER, V, MUST BE PROVIDED
%... TO DSS012 SO THAT THE DIRECTION OF FLOW IN EQUATION (1) CAN BE
%... USED TO SELECT THE APPROPRIATE FINITE DIFFERENCE APPROXIMATION
%... FOR THE FIRST-ORDER SPATIAL DERIVATIVE IN EQUATION (1), U .
%... THE CONVENTION FOR THE SIGN OF V IS X
%...
%... FLOW LEFT TO RIGHT V GT 0
%... (I.E., IN THE DIRECTION (I.E., THE SIXTH ARGUMENT IS
%... OF INCREASING X) POSITIVE IN CALLING DSS012)
%...
%... FLOW RIGHT TO LEFT V LT 0
%... (I.E., IN THE DIRECTION (I.E., THE SIXTH ARGUMENT IS
%... OF DECREASING X) NEGATIVE IN CALLING DSS012)
%...
%... COMPUTE THE SPATIAL INCREMENT, THEN SELECT THE FINITE DIFFERENCE
%... APPROXIMATION DEPENDING ON THE SIGN OF V IN EQUATION (1).
        dx=(xu-xl)/(n-1);
        if v > 0
%...
%... (1) FINITE DIFFERENCE APPROXIMATION FOR POSITIVE V
            ux(1)=(u(2)-u(1))/dx;
                for i=2:n
                    ux(i)=(u(i)-u(i-1))/dx;
                end
        end
%...
%... (2) FINITE DIFFERENCE APPROXIMATION FOR NEGATIVE V
        if v < 0
            nm1=n-1;
                for i=1:nm1
                    ux(i)=(u(i+1)-u(i))/dx;
                end
            ux(n)=(u(n)-u(n-1))/dx;
        end
end

%%The van Leer limiter
function [ux]=vanl2(xl,xu,n,u,v)
    dx=(xu-xl)/(n-1);
    delta=1.0e-05;
    if v >= 0.0
        for i=3:n-1
            if(abs(u(i)-u(i-1))<delta)
                phi(i)=0.0;
            else
                r(i)=(u(i+1)-u(i))/(u(i)-u(i-1));
                phi(i)=(r(i)+abs(r(i)))/(1.0+abs(r(i)));
            end
            if(abs(u(i-1)-u(i-2))<delta)
                phi(i-1)=0.0;
            else
                r(i-1)=(u(i)-u(i-1))/(u(i-1)-u(i-2));
                phi(i-1)=(r(i-1)+abs(r(i-1)))/(1.0+abs(r(i-1)));
            end
            flux2=u(i )+(u(i )-u(i-1))*phi(i )/2.0;
            flux1=u(i-1)+(u(i-1)-u(i-2))*phi(i-1)/2.0;
            ux(i)=(flux2-flux1)/dx;
        end
        ux(1)=(-u(1)+u(2))/dx;
        ux(2)=(-u(1)+u(2))/dx;
        ux(n)=(u(n)-u(n-1))/dx;
    end
    if v < 0.0
        for i=2:n-2
            if(abs(u(i)-u(i+1))<delta)
                phi(i)=0.0;
            else
                r(i)=(u(i-1)-u(i))/(u(i)-u(i+1));
                phi(i)=(r(i)+abs(r(i)))/(1.0+abs(r(i)));
            end
            if(abs(u(i+1)-u(i+2))<delta)
                phi(i+1)=0.0;
            else
                r(i+1)=(u(i)-u(i+1))/(u(i+1)-u(i+2));
                phi(i+1)=(r(i+1)+abs(r(i+1)))/(1.0+abs(r(i+1)));
            end
            flux2=u(i )+(u(i )-u(i+1))*phi(i )/2.0;
            flux1=u(i+1)+(u(i+1)-u(i+2))*phi(i+1)/2.0;
            ux(i)=(flux2-flux1)/dx;
        end
        ux(1)=(u(2)-u(1))/dx;
        ux(n-1)=(-u(n-1)+u(n))/dx;
        ux(n )=(-u(n-1)+u(n))/dx;
    end
end