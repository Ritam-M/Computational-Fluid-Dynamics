%  CL-613 Term Project, Group-E
%  Numerical solution of double diffusive convection
%  Numerical details : No. of grid/cells = 50*50
%                      Methodology : Finite volume discretisation, method of underrelaxation for solving discretised equations
%                      Data generated : Plot of streamlines, isotherms and iso-concentration lines
nx=50; ny=50; theta=60; Dy=1.0; Dx=tand(theta);
dx=Dx/nx; dy=Dy/ny;
Pr=0.7; Ra_T=2000; Ra_S=1000;
alpha_V=0.3; alpha_P=0.7;
maxit=120000; maxdiv=1e-6;

u_k=zeros(nx+1,ny+2); u_star=zeros(nx+1,ny+2); u_dash=zeros(nx+1,ny+2); u_cell=zeros(nx+1,ny+1);
aUE=zeros(nx+1,ny+2); aUW=zeros(nx+1,ny+2); aUN=zeros(nx+1,ny+2); aUS=zeros(nx+1,ny+2); aUP=zeros(nx+1,ny+2); dU=zeros(nx+1,ny+2);
v_k=zeros(nx+2,ny+1); v_star=zeros(nx+2,ny+1); v_dash=zeros(nx+2,ny+1); v_cell=zeros(nx+1,ny+1);
aVE=zeros(nx+2,ny+1); aVW=zeros(nx+2,ny+1); aVN=zeros(nx+2,ny+1); aVS=zeros(nx+2,ny+1); aVP=zeros(nx+2,ny+1); dV=zeros(nx+2,ny+1);
p=zeros(nx+2,ny+2); p_dash=zeros(nx+2,ny+2); p_cell=zeros(nx+1,ny+1);
aPE=zeros(nx+2,ny+2); aPW=zeros(nx+2,ny+2); aPN=zeros(nx+2,ny+2); aPS=zeros(nx+2,ny+2); aPP=zeros(nx+2,ny+2); bPP=zeros(nx+2,ny+2);
%The problem lies here : If we initialise temperatures and concenterations
%to 0, we are getting reverse circulation. On the other hand, if we
%initialise them as 1, we get almost similar results.
T_k=0.0*ones(nx+2,ny+2); T_k_minus_1=0.0*ones(nx+2,ny+2); T_cell=zeros(nx+1,ny+1);
aTE=zeros(nx+2,ny+2); aTW=zeros(nx+2,ny+2); aTN=zeros(nx+2,ny+2); aTS=zeros(nx+2,ny+2); aTP=zeros(nx+2,ny+2);
C_k=0.0*ones(nx+2,ny+2); C_k_minus_1=0.0*ones(nx+2,ny+2); C_cell=zeros(nx+1,ny+1);
aCE=zeros(nx+2,ny+2); aCW=zeros(nx+2,ny+2); aCN=zeros(nx+2,ny+2); aCS=zeros(nx+2,ny+2); aCP=zeros(nx+2,ny+2);

u_k_minus_1=u_k; v_k_minus_1=v_k; p_star=p;

for t=1:maxit
    %x-momentum eqn :
    for j=2:ny
        for i=2:nx+1-j
            %convection coefficients
            aCe=dy*0.5*(u_k_minus_1(i+1,j)+u_k_minus_1(i,j))/Pr;
            aCw=dy*0.5*(u_k_minus_1(i,j)+u_k_minus_1(i-1,j))/Pr;
            aCn=dx*0.5*(v_k_minus_1(i+1,j)+v_k_minus_1(i,j))/Pr;
            aCs=dx*0.5*(v_k_minus_1(i+1,j-1)+v_k_minus_1(i,j-1))/Pr;
            %diffusive coefficients
            aDx=dy/dx; aDy=dx/dy;
            %hybrid scheme
            aUE(i,j)=max([-aCe,(aDx-0.5*aCe),0]);
            aUW(i,j)=max([aCw,(aDx+0.5*aCw),0]);
            aUN(i,j)=max([-aCn,(aDy-0.5*aCn),0]);
            aUS(i,j)=max([aCs,(aDy+0.5*aCs),0]);
            aUP(i,j)=(aUE(i,j)+aUW(i,j)+aUN(i,j)+aUS(i,j))/alpha_V;
            %coefficient of pressure correction equation
            dU(i,j)=dy/aUP(i,j);
        end
    end
    u_star=u_k_minus_1;
    for j=2:ny
        for i=2:nx+1-j
            u_star(i,j)=(aUE(i,j)*u_star(i+1,j)+aUW(i,j)*u_star(i-1,j)+aUN(i,j)*u_star(i,j+1)+aUS(i,j)*u_star(i,j-1)-dy*(p_star(i+1,j)-p_star(i,j))+(1-alpha_V)*aUP(i,j)*u_k_minus_1(i,j))/aUP(i,j);
        end
    end
    %y-momentum eqn
    for j=2:ny
        for i=2:nx+1-j
            %convection coefficients
            aCe=dy*0.5*(u_k_minus_1(i,j+1)+u_k_minus_1(i,j))/Pr;
            aCw=dy*0.5*(u_k_minus_1(i-1,j+1)+u_k_minus_1(i-1,j))/Pr;
            aCn=dx*0.5*(v_k_minus_1(i,j+1)+v_k_minus_1(i,j))/Pr;
            aCs=dx*0.5*(v_k_minus_1(i,j)+v_k_minus_1(i,j-1))/Pr;
            %diffusive coefficients
            aDx=dy/dx; aDy=dx/dy;
            %hybrid scheme
            aVE(i,j)=max([-aCe,(aDx-0.5*aCe),0]);
            aVW(i,j)=max([aCw,(aDx+0.5*aCw),0]);
            aVN(i,j)=max([-aCn,(aDy-0.5*aCn),0]);
            aVS(i,j)=max([aCs,(aDy+0.5*aCs),0]);
            aVP(i,j)=(aVE(i,j)+aVW(i,j)+aVN(i,j)+aVS(i,j))/alpha_V;
            %coefficient of pressure correction equation
            dV(i,j)=dx/aVP(i,j);
        end
    end
    v_star=v_k_minus_1;
    for j=2:ny
        for i=2:nx+1-j
            v_star(i,j)=(aVE(i,j)*v_star(i+1,j)+aVW(i,j)*v_star(i-1,j)+aVN(i,j)*v_star(i,j+1)+aVS(i,j)*v_star(i,j-1)-dy*(p_star(i,j+1)-p_star(i,j))+(1-alpha_V)*aVP(i,j)*v_k_minus_1(i,j)+0.5*Ra_T*(T_k_minus_1(i,j+1)+T_k_minus_1(i,j))*dx*dy+0.5*Ra_S*(C_k_minus_1(i,j+1)+C_k_minus_1(i,j))*dx*dy)/aVP(i,j);
        end
    end
    %pressure correction equation (poisson eqn) :
    for j=2:ny
        for i=2:nx+1-j+1
            aPE(i,j)=dU(i,j)*dy;
            aPW(i,j)=dU(i-1,j)*dy;
            aPN(i,j)=dV(i,j)*dx;
            aPS(i,j)=dV(i,j-1)*dx;
            aPP(i,j)=aPE(i,j)+aPW(i,j)+aPN(i,j)+aPS(i,j);
            bPP(i,j)=(u_star(i-1,j)-u_star(i,j))*dy+(v_star(i,j-1)-v_star(i,j))*dx;
        end
    end
    for j=1:ny+2
        for i=1:nx+2
            p_dash(i,j)=0.0;
        end
    end
    for j=2:ny
        for i=2:nx+1-j+1
            p_dash(i,j)=(aPE(i,j)*p_dash(i+1,j)+aPW(i,j)*p_dash(i-1,j)+aPN(i,j)*p_dash(i,j+1)+aPS(i,j)*p_dash(i,j-1)+bPP(i,j))/aPP(i,j);
        end
    end
    %velocity correction terms :
    for j=2:ny
        for i=2:nx+1-j+1
            u_dash(i,j)=dU(i,j)*(p_dash(i,j)-p_dash(i+1,j));
        end
    end
    for j=2:ny
        for i=2:nx+1-j+1
            v_dash(i,j)=dV(i,j)*(p_dash(i,j)-p_dash(i,j+1));
        end
    end
    %pressure updation (p_new=p* + alpha_P.p_dash, to reduce variables p_new and p_dash are taken as the same).
    for j=2:ny
        for i=2:nx+1-j+1
            p_dash(i,j)=p_star(i,j)+p_dash(i,j)*alpha_P;
        end
    end
    for j=2:ny
        for i=2:nx+1-j+1
            u_k(i,j)=u_star(i,j)+u_dash(i,j);
        end
    end
    for j=2:ny
        for i=2:nx+1-j+1
            v_k(i,j)=v_star(i,j)+v_dash(i,j);
        end
    end
    %Temperqature and concentration equations :
    for j=2:ny
        for i=2:nx+1-j+1
            aCe=u_k(i,j)*dy; aCw=u_k(i-1,j)*dy; aCn=v_k(i,j)*dx; aCs=v_k(i,j-1)*dx; aDx=dy/dx; aDy=dx/dy;
            %hybrid scheme
            aTE(i,j)=max([-aCe,(aDx-0.5*aCe),0]);
            aTW(i,j)=max([aCw,(aDx+0.5*aCw),0]);
            aTN(i,j)=max([-aCn,(aDy-0.5*aCn),0]);
            aTS(i,j)=max([aCs,(aDy+0.5*aCs),0]);
            aTP(i,j)=(aTE(i,j)+aTW(i,j)+aTN(i,j)+aTS(i,j))/alpha_P;
        end
    end
    T_k=T_k_minus_1;
    for j=2:ny
        for i=2:nx+1-j+1
            T_k(i,j)=(aTE(i,j)*T_k(i+1,j)+aTW(i,j)*T_k(i-1,j)+aTN(i,j)*T_k(i,j+1)+aTS(i,j)*T_k(i,j-1)+(1-alpha_P)*T_k_minus_1(i,j)*aTP(i,j))/aTP(i,j);
        end
    end
    %alpha_C=(10^-0.5)*alpha_V;
    for j=2:ny
        for i=2:nx+1-j+1
            aCe=u_k(i,j)*dy; aCw=u_k(i-1,j)*dy; aCn=v_k(i,j)*dx; aCs=v_k(i,j-1)*dx; aDx=(10^-0.5)*dy/dx; aDy=(10^-0.5)*dx/dy;
            aCE(i,j)=max([-aCe,(aDx-0.5*aCe),0]);
            aCW(i,j)=max([aCw,(aDx+0.5*aCw),0]);
            aCN(i,j)=max([-aCn,(aDy-0.5*aCn),0]);
            aCS(i,j)=max([aCs,(aDy+0.5*aCs),0]);
            aCP(i,j)=(aCE(i,j)+aCW(i,j)+aCN(i,j)+aCS(i,j))/alpha_P;
        end
    end
    C_k=C_k_minus_1;
    for j=2:ny
        for i=2:nx+1-j+1
            C_k(i,j)=(aCE(i,j)*C_k(i+1,j)+aCW(i,j)*C_k(i-1,j)+aCN(i,j)*C_k(i,j+1)+aCS(i,j)*C_k(i,j-1)+(1-alpha_P)*C_k_minus_1(i,j)*aCP(i,j))/aCP(i,j);
        end
    end
    error1=max(max(abs(u_k_minus_1(2:nx,2:ny+1)-u_k(2:nx,2:ny+1)))); error2=max(max(abs(v_k_minus_1(2:nx+1,2:ny)-v_k(2:nx+1,2:ny))));
    error3=max(max(abs(T_k_minus_1(2:nx+1,2:ny+1)-T_k(2:nx+1,2:ny+1)))); error4=max(max(abs(C_k_minus_1(2:nx+1,2:ny+1)-C_k(2:nx+1,2:ny+1))));
    if(t>100)
        error=max([error1,error2,error3,error4]);
        if(error<maxdiv)
            break;
        end
    end
    u_k_minus_1=u_k; v_k_minus_1=v_k; p_star=p_dash; T_k_minus_1=T_k; C_k_minus_1=C_k;
    %boundary conditions :
    for i=2:nx
        u_k_minus_1(i,1)=u_k_minus_1(i,2);
    end
    for j=1:ny+2
        u_k_minus_1(1,j)=0.0;
    end
    for j=1:ny+2
        for i=1:nx+1
            if(i+j==nx+3)
                u_k_minus_1(i,j)=u_k_minus_1(i,j-1);
            end
        end
    end
    for j=1:ny+1
        v_k_minus_1(1,j)=v_k_minus_1(2,j);
    end
    for j=1:ny+1
        for i=1:nx+2
            if(i+j==nx+3)
                v_k_minus_1(i,j)=v_k_minus_1(i-1,j);
            end
        end
    end
    for i=1:nx+2
        v_k_minus_1(i,1)=0.0; v_k_minus_1(i,ny+1)=0.0;
    end
    for i=1:nx+2
        p_dash(i,1)=-p_dash(i,2);
    end
    for j=1:ny+2
        for i=1:nx+2
            if(i+j==nx+2)
               p_dash(i+1,j+1)=-p_dash(i,j);
            end
        end
    end
    for j=2:ny+1
        p_dash(1,j)=p_dash(2,j);
    end
    for i=2:nx+1
        T_k_minus_1(i,1)=T_k_minus_1(i,2);
    end
    for i=2:nx+1
        for j=2:ny+1
            if(i+j==nx+2)
                T_k_minus_1(i+1,j+1)=0-T_k_minus_1(i,j);
            end
        end
    end
    for j=1:ny+2
        T_k_minus_1(1,j)=2-T_k_minus_1(2,j);
    end
    for i=2:nx+1
        C_k_minus_1(i,1)=C_k_minus_1(i,2);
    end
    for i=2:nx+1
        for j=2:ny+1
            if(i+j==nx+2)
                C_k_minus_1(i+1,j+1)=0-C_k_minus_1(i,j);
            end
        end
    end
    for j=1:ny+2
        C_k_minus_1(1,j)=2-C_k_minus_1(2,j);
    end
end

%calculation of grid point parameters by taking average of staggered grid parameters :
for i=1:nx+1
    for j=1:ny+1
        u_cell(i,j)=(u_k(i,j+1)+u_k(i,j))/2;
        v_cell(i,j)=(v_k(i+1,j)+v_k(i,j))/2;
        p_cell(i,j)=(p_dash(i,j)+p_dash(i+1,j)+p_dash(i+1,j+1)+p_dash(i,j+1))/4;
        T_cell(i,j)=(T_k(i,j)+T_k(i+1,j)+T_k(i+1,j+1)+T_k(i,j+1))/4;
        C_cell(i,j)=(C_k(i,j)+C_k(i+1,j)+C_k(i+1,j+1)+C_k(i,j+1))/4;
    end
end

%vorticity :
w=zeros(nx+1,ny+1);
for i=2:nx
    for j=2:ny
        w(i,j)=(-(u_cell(i,j+1)-u_cell(i,j-1))/(2*dy))+((v_cell(i+1,j)-v_cell(i-1,j))/(2*dx));
    end
end

%generation of stream function for streamlines :
psi=zeros(nx+1,ny+1); psi_new=zeros(nx+1,ny+1); errorS=10.0; maxS=0.0;
while(errorS>0.0001)
    for i=2:nx
        for j=2:ny
            psi_new(i,j)=0.25*0.5*(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1)+(dx*dx*w(i,j)))+((1-0.5)*psi(i,j));
        end
    end
    maxS=0.0;
    for i=2:nx-1
        for j=2:ny-1
            errorS=abs(psi_new(i,j)-psi(i,j));
            if(errorS>maxS)
                maxS=errorS;
            end
        end
    end
    errorS=maxS;
    psi=psi_new;
end

F=figure(1);
subplot(1,3,1);
[c1,h1]=contour(linspace(0,Dx,nx+1),linspace(0,Dy,ny+1),psi_new.',15);
set(h1,'edgecolor','k');
title('Streamline');
patch([0,tand(theta),tand(theta),0],[1,0,1,1],'k');
subplot(1,3,2);
cMap=jet(256);
[c2,h2]=contourf(linspace(0,Dx,nx+1),linspace(0,Dy,ny+1),T_cell.',100);
set(h2,'edgecolor','none');
colormap(cMap);
title('Isotherm');
patch([0,tand(theta),tand(theta),0],[1,0,1,1],'k');
subplot(1,3,3);
cMap=jet(256);
[c3,h3]=contourf(linspace(0,Dx,nx+1),linspace(0,Dy,ny+1),C_cell.',100);
set(h3,'edgecolor','none');
colormap(cMap);
title('Iso-concentration');
patch([0,tand(theta),tand(theta),0],[1,0,1,1],'k');
set(F,'WindowStyle','docked');

x=linspace(0,Dx,nx+1); y=linspace(0,Dy,ny+1); data=zeros((nx+1)*(ny+1),9); k=1;
for i=1:nx+1
    for j=1:ny+1
        data(k,:)=[x(1,i),y(1,j),u_cell(i,j),v_cell(i,j),p_cell(i,j),w(i,j),psi_new(i,j),T_cell(i,j),C_cell(i,j)];
        k=k+1;
    end
end