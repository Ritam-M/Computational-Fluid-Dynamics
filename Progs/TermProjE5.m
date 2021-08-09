% Problem Paramters
RaT=10000.0; RaS=1000.0;
% Grid Paramters
n=51;                                               % imax=jmax=n
h=1/(n-1);                                          % Grid step size
% Time March Paramters
maxit=6000;                                         % Maximum number of time steps
dt=.0001;                                           % Time step size
t=0;                                                % Current time step
% SOR Paramters
b=1.5;                                              % Over-relaxation factor
mk=100;                                             % Maximum number of SOR iterations
% Declare and Initialize
[w,si]=deal(zeros(n,n));                            % Vorticity and stream function
[T,C]=deal(zeros(n,n));                             % Temperature and Concentration
% Time March Loop
while t<maxit
    % Compute the boundary nodes values
    for i=2:n-1
        w(i,1)=2*(si(i,1)-si(i,2)+0.5*si(i,3))/(h*h);             % vorticity on bottom wall
        T(i,1)=1.0; C(i,1)=1.0;
        w(i,n)=2*(si(i,n)-si(i,n-1)+0.5*si(i,n-2))/(h*h);           % vorticity on top wall
        T(i,n)=0.0; C(i,n)=0.0;
        w(1,i)=2*(si(1,i)-si(2,i)+0.5*si(3,i))/(h*h);             % vorticity on left wall
        T(1,i)=T(2,i); C(1,i)=C(2,i);
        w(n,i)=2*(si(n,i)-si(n-1,i)+0.5*si(n-2,i))/(h*h);           % vorticity on right wall
        T(n,i)=T(n-1,i); C(n,i)=C(n-1,i);
    end
    % Compute the interior nodes initial vorticity (predictor)
%     for i=2:n-1
%         for j=2:n-1
%             w(i,j)=(4*si(i,j)-si(i+1,j)-si(i-1,j)-si(i,j+1)-si(i,j-1))/(h*h);
%         end
%     end
    % Compute the interior nodes new time step vorticity, temp and concn.
    for i=2:n-1
        for j=2:n-1
            w(i,j)=w(i,j)+dt*((((si(i+1,j)-si(i-1,j))*(w(i,j+1)-w(i,j-1)))/(4*h*h))...
                -(((w(i+1,j)-w(i-1,j))*(si(i,j+1)-si(i,j-1)))/(4*h*h))...
                +(w(i+1,j)+w(i-1,j)+w(i,j+1)+w(i,j-1)-4*w(i,j))*(1/(h*h))+RaT*(T(i+1,j)-T(i-1,j))/(2*h)+RaS*(C(i+1,j)-C(i-1,j))/(2*h));
            T(i,j)=T(i,j)+dt*((((si(i+1,j)-si(i-1,j))*(T(i,j+1)-T(i,j-1)))/(4*h*h))...
                -(((T(i+1,j)-T(i-1,j))*(si(i,j+1)-si(i,j-1)))/(4*h*h))...
                +(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)-4*T(i,j))*(1/(h*h)));
            C(i,j)=C(i,j)+dt*((((si(i+1,j)-si(i-1,j))*(C(i,j+1)-C(i,j-1)))/(4*h*h))...
                -(((C(i+1,j)-C(i-1,j))*(si(i,j+1)-si(i,j-1)))/(4*h*h))...
                +(10^-0.5)*(C(i+1,j)+C(i-1,j)+C(i,j+1)+C(i,j-1)-4*C(i,j))*(1/(h*h)));
        end
    end
    % Compute the interior nodes new time step stream function
    ck=0;                                           % Current iteration number
    while ck<mk
        for i=2:n-1
            for j=2:n-1
                si(i,j)=(1-b)*si(i,j)+(si(i+1,j)+si(i-1,j)+si(i,j+1)+si(i,j-1)+w(i,j)*(h*h))*(b/4);
            end
        end
        ck=ck+1;
    end
    t=t+1;                                      % Increment current time
end
[u,v]=deal(zeros(n,n));
for i=2:n-1
    for j=2:n-1
        u(i,j)=(si(i,j+1)-si(i,j-1))/(2*h);
        v(i,j)=-(si(i+1,j)-si(i-1,j))/(2*h);
    end
end

x=linspace(0,1,n); y=linspace(0,1,n);
data=zeros(n*n,8); k=1;
for i=1:n
    for j=1:n
        data(k,:)=[x(1,i), y(1,i), u(i,j), v(i,j), w(i,j), si(i,j), T(i,j), C(i,j)];
        k=k+1;
    end
end

F=figure(1);
subplot(1,3,1);
[c1,h1]=contour(linspace(0,1,n),linspace(0,1,n),si.',15);
set(h1,'edgecolor','k');
title('Streamline');
subplot(1,3,2);
cMap=jet(256);
[c2,h2]=contourf(linspace(0,1,n),linspace(0,1,n),T.',100);
set(h2,'edgecolor','none');
colormap(cMap);
title('Isotherm');
subplot(1,3,3);
cMap=jet(256);
[c3,h3]=contourf(linspace(0,1,n),linspace(0,1,n),C.',100);
set(h3,'edgecolor','none');
colormap(cMap);
title('Iso-concentration');
set(F,'WindowStyle','docked');