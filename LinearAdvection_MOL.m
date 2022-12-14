%% Linear 1D Advection equation
% Routine for solving the following PDE
%
%    du/dt + v*du/dx = 0
%
% Using the Method of Lines
clc; clear; close all; 

%% Parameters
v=0.589/1000/pi/0.05248^2*4;
zl=39;
n=21;
dz=zl/(n-1);

R=2.1641;                   % (m`K)/W,  Thermal resistance per unit length from fluid to boundary temperature
Cp=4200;            % J/(kg.K), 水的比热容
A=0.05248^2*4;          %管道横截面积
C=A*Cp*1000;
tauRC=R*C;

ambTemperature=291-273.15;
initialTemperature=18.2;

% 加载管道入口及出口温度数据
% 第一列：时刻/s；第二列：质量流量kg/s；
% 第三列：出口管道温度℃；第四列：出口水温；
% 第五列：入口管道温度；第六列：入口水温
PipeData=load('E:\\programs\\matlab\\LinearAdvection_MOL-main\\DataULg2.txt');
inletTemperature=PipeData(:,6);

%% Select three point, finite differencing (FD) of spatial derivative
% ifd = 1: Centered approximations
% ifd = 2: Two point upwind approximation
% ifd = 3: Five point, biased upwind approximation
% ifd = 4: van Leer flux limiter
ifd=1;

%% Initial Condition
% Initial, final times, integration interval, number of Euler
% steps for each output
t=0.0; tf=875.0; 
nout=500;
h=(tf-t)/nout;

t_in_Start=min(0,39/0.589*1000*0.05248^2/4*pi);       % t_in_Start=pipeLength/m_flow_Start*rho*dh^2/4*pi;   
t_out_Start=min(-0,39/0.589*1000*0.05248^2/4*pi);    
t0=t+t_in_Start;

PipeLength=39;     

% Initial conditions
for i=1:n
    Timein(i)=t0;
    Temperaturein(i)=initialTemperature;
end
x=0;

TimeSpan=PipeData(:,1);
inletT=inletTemperature(1);
delayInletTemperature=[];
outletTemperature=[];

%% Integrate until t = tf
% calculate the time delay,v(x,t)=inletTime(x,t)
tic
for j=1:length(TimeSpan)
    % Take nout Euler steps
    for iout=1:nout
            % Monitor solution by displaying x
            currentTime=t;
            currentLength=x;
            if currentLength>=PipeLength
                TimeDelay=Timein(1)-Timein(n);
                delayInletTemperature(end+1)=Temperaturein(n);               
                break;
            end
            
            if ifd==1       % Centered approximations
            % Boundary condition at z = 0
                Timein(1)=currentTime;
                Temperaturein(1)=inletTemperature(j);  
            % Spatial derivative
                [Timein_dz]=dss002(0.0,zl,n,Timein);
                [Temperaturein_dz]=dss002(0.0,zl,n,Temperaturein);
            % End of three point centered approximation
            end
            
            if ifd==2      % Two point upwind approximation
                Timein(1)=t;
                [Timein_dz]=dss012(0.0,zl,n,Timein,v);
            end
            
            % Temporal derivative
            Timein_dt(1)=0;
            Temperaturein_dt(1)=0;
            for i=2:n
                Timein_dt(i)=-func(x,t)*Timein_dz(i);
                Temperaturein_dt(i)=-func(x,t)*Temperaturein_dz(i);
            end            

            % Take Euler step                
                x=x+h*func(x,t);

            for i=1:n
                Timein(i)=Timein(i)+Timein_dt(i)*h;
                Temperaturein(i)=Temperaturein(i)+Temperaturein_dt(i)*h;            %这是一组ODE，后面的变化速率要远远慢于前面
            end
                  
        t=t+h;
    end
    t=0;
    x=0;
    for i=1:n
        Timein(i)=t0;
        Temperaturein(i)=Temperaturein(n);
    end
end  

toc
%% Calculate the outlet temperature
delayInletTemperature=delayInletTemperature';
h2=0.01;
outTime=[];
outletTem=[initialTemperature,initialTemperature];
t=0;

for j=1:length(TimeSpan)          
count=0;
    while t<tf
                
        if  t>=TimeDelay       
            tau=max(0,t-TimeDelay-TimeSpan(j));
            if tau >0
                outletTem(end+1)=ambTemperature+(delayInletTemperature(j)-ambTemperature)*exp(-tau/tauRC);
                
                if abs(outletTem(end)-outletTem(end-1:end))<1e-2
                        count=count+h2;
                end
                if count>7
                    outletTemperature(end+1)=outletTem(end);
                    outTime(end+1)=t;
                    t=t+h2;
                    break;
                end
            end

        end
        t=t+h2;
    end

end    

%% 绘图
measureTemperature=PipeData(:,4);
outletTem(1:2)=[];
delayInletTemperature=[initialTemperature;delayInletTemperature];
delayInletTimeSpan=[0;TimeSpan+TimeDelay];
for k =1:length(delayInletTimeSpan)
    if delayInletTimeSpan(k)>600
        delayInletTimeSpan(k:end)=[];
        delayInletTemperature(k:end)=[];
        break;
    end
end
outletTemperature=[initialTemperature,outletTemperature];
outTime=[0,outTime];
for k=1:length(outTime)
    if outTime(k)>600
        outTime(k:end)=[];
        outletTemperature(k:end)=[];
        break;
    end
end
plot(TimeSpan,inletTemperature,'black-.',delayInletTimeSpan,delayInletTemperature,'red--', ...
        outTime,outletTemperature,'green',TimeSpan,measureTemperature,'blue');

%% v(t)
function dx=func(x,t)
v=0.589/1000/pi/0.05248^2*4;      
dx=v;
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