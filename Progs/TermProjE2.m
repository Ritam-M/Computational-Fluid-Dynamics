%  CL-613 Term Project, Group-E
%  Numerical solution of double diffusive convection
%  Numerical details : No. of grid/cells = 50*50
%                      Methodology : Finite volume discretisation, method of underrelaxation for solving discretised equations
%                      Data generated : Plot of streamlines, isotherms and iso-concentration lines

theta=30; ny=50; Dy=1.0*sind(theta);
nx=ny+ceil(ny/cosd(theta)); Dx=1.0+1.0*cosd(theta);
dx=Dx/nx; dy=Dy/ny; b=floor(1.0/dx);
Pr=1.0; Ra_T=10000; Ra_S=1000;
alpha_V=0.3; alpha_P=0.7;
maxit=120000; maxdiv=1e-6;

u_k=zeros(nx+1,ny+2); u_star=zeros(nx+1,ny+2); u_dash=zeros(nx+1,ny+2); u_cell=zeros(nx+1,ny+1);
aUE=zeros(nx+1,ny+2); aUW=zeros(nx+1,ny+2); aUN=zeros(nx+1,ny+2); aUS=zeros(nx+1,ny+2); aUP=zeros(nx+1,ny+2); dU=zeros(nx+1,ny+2);
v_k=zeros(nx+2,ny+1); v_star=zeros(nx+2,ny+1); v_dash=zeros(nx+2,ny+1); v_cell=zeros(nx+1,ny+1);
aVE=zeros(nx+2,ny+1); aVW=zeros(nx+2,ny+1); aVN=zeros(nx+2,ny+1); aVS=zeros(nx+2,ny+1); aVP=zeros(nx+2,ny+1); dV=zeros(nx+2,ny+1);
p=zeros(nx+2,ny+2); p_dash=zeros(nx+2,ny+2); p_cell=zeros(nx+1,ny+1);
aPE=zeros(nx+2,ny+2); aPW=zeros(nx+2,ny+2); aPN=zeros(nx+2,ny+2); aPS=zeros(nx+2,ny+2); aPP=zeros(nx+2,ny+2); bPP=zeros(nx+2,ny+2);
T_k=0.5*ones(nx+2,ny+2); T_k_minus_1=0.5*ones(nx+2,ny+2); T_cell=zeros(nx+1,ny+1);
aTE=zeros(nx+2,ny+2); aTW=zeros(nx+2,ny+2); aTN=zeros(nx+2,ny+2); aTS=zeros(nx+2,ny+2); aTP=zeros(nx+2,ny+2);
C_k=1.0*ones(nx+2,ny+2); C_k_minus_1=1.0*ones(nx+2,ny+2); C_cell=zeros(nx+1,ny+1);
aCE=zeros(nx+2,ny+2); aCW=zeros(nx+2,ny+2); aCN=zeros(nx+2,ny+2); aCS=zeros(nx+2,ny+2); aCP=zeros(nx+2,ny+2);

u_k_minus_1=u_k; v_k_minus_1=v_k; p_star=p;

error1=zeros(size(u_k)); error2=zeros(size(v_k)); error3=zeros(size(T_k)); error4=zeros(size(C_k));

for t=1:maxit
    %x-momentum eqn :
    for j=2:ny+1
        for i=j:j+b-1
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
            aUP(i,j)=(aUE(i,j)+aUW(i,j)+aUN(i,j)+aUS(i,j)+aCe-aCw+aCn-aCs)/alpha_V;
            %coefficient of pressure correction equation
            dU(i,j)=dy/aUP(i,j);
        end
    end
    u_star=u_k_minus_1;
    for j=2:ny+1
        for i=j:j+b-1
            u_star(i,j)=(aUE(i,j)*u_star(i+1,j)+aUW(i,j)*u_star(i-1,j)+aUN(i,j)*u_star(i,j+1)+aUS(i,j)*u_star(i,j-1)-dy*(p_star(i+1,j)-p_star(i,j))+(1-alpha_V)*aUP(i,j)*u_k_minus_1(i,j))/aUP(i,j);
        end
    end
    %y-momentum eqn
    for j=2:ny
        for i=j+1:j+b
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
            aVP(i,j)=(aVE(i,j)+aVW(i,j)+aVN(i,j)+aVS(i,j)+aCe-aCw+aCn-aCs)/alpha_V;
            %coefficient of pressure correction equation
            dV(i,j)=dx/aVP(i,j);
        end
    end
    v_star=v_k_minus_1;
    for j=2:ny
        for i=j+1:j+b
            v_star(i,j)=(aVE(i,j)*v_star(i+1,j)+aVW(i,j)*v_star(i-1,j)+aVN(i,j)*v_star(i,j+1)+aVS(i,j)*v_star(i,j-1)-dx*(p_star(i,j+1)-p_star(i,j))+(1-alpha_V)*aVP(i,j)*v_k_minus_1(i,j)+0.5*Ra_T*(T_k_minus_1(i,j+1)+T_k_minus_1(i,j))*dx*dy+0.5*Ra_S*(C_k_minus_1(i,j+1)+C_k_minus_1(i,j))*dx*dy)/aVP(i,j);
        end
    end
    %pressure correction equation (poisson eqn) :
    for j=2:ny+1
        for i=j+1:j+b
            aPE(i,j)=dU(i,j)*dy;
            aPW(i,j)=dU(i-1,j)*dy;
            aPN(i,j)=dV(i,j)*dx;
            aPS(i,j)=dV(i,j-1)*dx;
            aPP(i,j)=aPE(i,j)+aPW(i,j)+aPN(i,j)+aPS(i,j);
            bPP(i,j)=(u_star(i-1,j)-u_star(i,j))*dy+(v_star(i,j-1)-v_star(i,j))*dx;
        end
    end
    for i=1:nx+2
        for j=1:ny+2
            p_dash(i,j)=0.0;
        end
    end
    for j=2:ny+1
        for i=j+1:j+b
            p_dash(i,j)=(aPE(i,j)*p_dash(i+1,j)+aPW(i,j)*p_dash(i-1,j)+aPN(i,j)*p_dash(i,j+1)+aPS(i,j)*p_dash(i,j-1)+bPP(i,j))/aPP(i,j);
        end
    end
    %velocity correction terms :
    for j=2:ny+1
        for i=j:j+b-1
            u_dash(i,j)=dU(i,j)*(p_dash(i,j)-p_dash(i+1,j));
        end
    end
    for j=2:ny
        for i=j+1:j+b
            v_dash(i,j)=dV(i,j)*(p_dash(i,j)-p_dash(i,j+1));
        end
    end
    %pressure updation (p_new=p* + alpha_P.p_dash, to reduce variables p_new and p_dash are taken as the same).
    for j=2:ny+1
        for i=j+1:j+b
            p_dash(i,j)=p_star(i,j)+p_dash(i,j)*alpha_P;
        end
    end
    for j=2:ny+1
        for i=j:j+b-1
            u_k(i,j)=u_star(i,j)+u_dash(i,j);
        end
    end
    for j=2:ny
        for i=j+1:j+b
            v_k(i,j)=v_star(i,j)+v_dash(i,j);
        end
    end
    %Temperqature and concentration equations :
    for j=2:ny+1
        for i=j+1:j+b
            aCe=u_k(i,j)*dy; aCw=u_k(i-1,j)*dy; aCn=v_k(i,j)*dx; aCs=v_k(i,j-1)*dx; aDx=dy/dx; aDy=dx/dy;
            %hybrid scheme
            aTE(i,j)=max([-aCe,(aDx-0.5*aCe),0]);
            aTW(i,j)=max([aCw,(aDx+0.5*aCw),0]);
            aTN(i,j)=max([-aCn,(aDy-0.5*aCn),0]);
            aTS(i,j)=max([aCs,(aDy+0.5*aCs),0]);
            aTP(i,j)=(aTE(i,j)+aTW(i,j)+aTN(i,j)+aTS(i,j)+aCe-aCw+aCn-aCs)/alpha_P;
        end
    end
    T_k=T_k_minus_1;
    for j=2:ny+1
        for i=j+1:j+b
            T_k(i,j)=(aTE(i,j)*T_k(i+1,j)+aTW(i,j)*T_k(i-1,j)+aTN(i,j)*T_k(i,j+1)+aTS(i,j)*T_k(i,j-1)+(1-alpha_P)*T_k_minus_1(i,j)*aTP(i,j))/aTP(i,j);
        end
    end
    for j=2:ny+1
        for i=j+1:j+b
            aCe=u_k(i,j)*dy; aCw=u_k(i-1,j)*dy; aCn=v_k(i,j)*dx; aCs=v_k(i,j-1)*dx; aDx=(10^-0.5)*dy/dx; aDy=(10^-0.5)*dx/dy;
            aCE(i,j)=max([-aCe,(aDx-0.5*aCe),0]);
            aCW(i,j)=max([aCw,(aDx+0.5*aCw),0]);
            aCN(i,j)=max([-aCn,(aDy-0.5*aCn),0]);
            aCS(i,j)=max([aCs,(aDy+0.5*aCs),0]);
            aCP(i,j)=(aCE(i,j)+aCW(i,j)+aCN(i,j)+aCS(i,j)+aCe-aCw+aCn-aCs)/alpha_P;
        end
    end
    C_k=C_k_minus_1;
    for j=2:ny+1
        for i=j+1:j+b
            C_k(i,j)=(aCE(i,j)*C_k(i+1,j)+aCW(i,j)*C_k(i-1,j)+aCN(i,j)*C_k(i,j+1)+aCS(i,j)*C_k(i,j-1)+(1-alpha_P)*C_k_minus_1(i,j)*aCP(i,j))/aCP(i,j);
        end
    end
%     error1=max(max(abs(u_k_minus_1(2:nx,2:ny+1)-u_k(2:nx,2:ny+1)))); error2=max(max(abs(v_k_minus_1(2:nx+1,2:ny)-v_k(2:nx+1,2:ny))));
%     error3=max(max(abs(T_k_minus_1(2:nx+1,2:ny+1)-T_k(2:nx+1,2:ny+1)))); error4=max(max(abs(C_k_minus_1(2:nx+1,2:ny+1)-C_k(2:nx+1,2:ny+1))));
    for j=2:ny+1
        for i=j+2:j+b-1
            error1(i,j)=abs(u_k_minus_1(i,j)-u_k(i,j));
        end
    end
    for j=2:ny
        for i=j+1:j+b-1-1
            error2(i,j)=abs(v_k_minus_1(i,j)-v_k(i,j));
        end
    end
    for j=2:ny+1
        for i=j+1+1:j+b-1
            error3(i,j)=abs(T_k_minus_1(i,j)-T_k(i,j)); error4(i,j)=abs(T_k_minus_1(i,j)-T_k(i,j));
        end
    end
    error11=max(max(error1)); error22=max(max(error2)); error33=max(max(error3)); error44=max(max(error4));
    if(t>100)
        error=max([error11,error22,error33,error44]);
        if(error<maxdiv)
            break;
        end
    end
    u_k_minus_1=u_k; v_k_minus_1=v_k; p_star=p_dash; T_k_minus_1=T_k; C_k_minus_1=C_k;
    %boundary conditions :
    %traction free u
    for i=2:1:b-1
        u_k_minus_1(i,1)=u_k_minus_1(i,2);%bottom wall
        u_k_minus_1(i,ny+2)=u_k_minus_1(i,ny+1);%top wall
    end
    %symmetry u
    for j=1:ny+2
        for i=1:nx+1
            if(i==j)
                u_k_minus_1(i,j)=0.0; %left wall
            end
            if(i==j+b)
                u_k_minus_1(i,j)=0.0; %right wall
            end
        end
    end
    %symmetry v
    for j=2:ny
        for i=2:nx+1
            if(i==j)
                v_k_minus_1(i,j)=v_k_minus_1(i+1,j);%left wall
            end
            if(i==j+b)
                v_k_minus_1(i+1,j)=v_k_minus_1(i,j);%right wall
            end
        end
    end
    %traction free v
    for i=1:nx+2
        v_k_minus_1(i,1)=0.0; %bottom wall
        v_k_minus_1(i,ny+1)=0.0;%top wall
    end
    %traction free p
    for i=1:1+b
        p_dash(i,1)=-p_dash(i,2);%bottom wall
    end
    for i=ny:ny+b
        p_dash(i,ny+2)=-p_dash(i,ny+1);%top wall
    end
    %symmetry p
    for j=2:ny+1
        for i=2:nx+1
            if(i==j)
                p_dash(i,j+1)=p_dash(i+1,j);%left wall
            end
            if(i==j+b)
                p_dash(i+1,j)=p_dash(i,j+1);%right wall
            end
        end
    end
    %constant wall temperature
    for i=2:1+b
        T_k_minus_1(i,1)=2-T_k_minus_1(i,2);%bottom wall
    end
    for i=ny+1:ny+1+b
        T_k_minus_1(i,ny+2)=-T_k_minus_1(i,ny+1);%top wall
    end
    %constant wall temperature flux
    for j=2:ny+1
        for i=2:nx+1
            if(i==j)
                T_k_minus_1(i,j+1)=T_k_minus_1(i+1,j);%left wall
            end
            if(i==j+b)
                T_k_minus_1(i+1,j)=T_k_minus_1(i,j+1);%right wall
            end
        end
    end
    %constant wall concentration
    for i=2:1+b
        C_k_minus_1(i,1)=2-C_k_minus_1(i,2);%bottom wall
    end
    for i=ny+1:ny+1+b
        C_k_minus_1(i,ny+2)=-C_k_minus_1(i,ny+1);%top wall
    end
    %constant wall concentration flux
    for j=2:ny+1
        for i=2:nx+1
            if(i==j)
                C_k_minus_1(i,j+1)=C_k_minus_1(i+1,j);%left wall
            end
            if(i==j+b)
                C_k_minus_1(i+1,j)=C_k_minus_1(i,j+1);%right wall
            end
        end
    end
end

%calculation of grid point parameters by taking average of staggered grid parameters :
for j=1:ny+1
    for i=j:j+b
        u_cell(i,j)=(u_k(i,j+1)+u_k(i,j))/2;
        v_cell(i,j)=(v_k(i+1,j)+v_k(i,j))/2;
        p_cell(i,j)=(p_dash(i,j)+p_dash(i+1,j)+p_dash(i+1,j+1)+p_dash(i,j+1))/4;
        T_cell(i,j)=(T_k(i,j)+T_k(i+1,j)+T_k(i+1,j+1)+T_k(i,j+1))/4;
        C_cell(i,j)=(C_k(i,j)+C_k(i+1,j)+C_k(i+1,j+1)+C_k(i,j+1))/4;
    end
end

%vorticity :
w=zeros(nx+1,ny+1);
for j=2:ny
    for i=j:j+b-1
        w(i,j)=(-(u_cell(i,j+1)-u_cell(i,j-1))/(2*dy))+((v_cell(i+1,j)-v_cell(i-1,j))/(2*dx));
    end
end

%generation of stream function for streamlines :
psi=zeros(nx+1,ny+1); psi_new=zeros(nx+1,ny+1); errorS=10.0; maxS=0.0;
while(errorS>0.0001)
    for j=2:ny
        for i=j+1:j+b-1
            psi_new(i,j)=0.25*0.5*(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1)+(dx*dx*w(i,j)))+((1-0.5)*psi(i,j));
        end
    end
    maxS=0.0;
    for j=2:ny-1
        for i=j+1+1:j+b-1-1
            errorS=abs(psi_new(i,j)-psi(i,j));
            if(errorS>maxS)
                maxS=errorS;
            end
        end
    end
    errorS=maxS;
    psi=psi_new;
end

%plot data :
F=figure(1);
subplot(1,3,1);
[c1,h1]=contour(linspace(0,Dx,nx+1),linspace(0,Dy,ny+1),psi_new.',15);
patch([0,0,cosd(theta),0],[0,sind(theta),sind(theta),0],[1,1,1]);
hold on;
patch([1,1+cosd(theta),1+cosd(theta),1],[0,0,sind(theta),0],[1,1,1]);
set(h1,'edgecolor','k');
title('Streamline');
subplot(1,3,2);
cMap=jet(256);
[c2,h2]=contourf(linspace(0,Dx,nx+1),linspace(0,Dy,ny+1),T_cell.',100);
patch([0,0,cosd(theta),0],[0,sind(theta),sind(theta),0],[1,1,1]);
hold on;
patch([1,1+cosd(theta),1+cosd(theta),1],[0,0,sind(theta),0],[1,1,1]);
set(h2,'edgecolor','none');
colormap(cMap);
title('Isotherm');
subplot(1,3,3);
cMap=jet(256);
[c3,h3]=contourf(linspace(0,Dx,nx+1),linspace(0,Dy,ny+1),C_cell.',100);
patch([0,0,cosd(theta),0],[0,sind(theta),sind(theta),0],[1,1,1]);
hold on;
patch([1,1+cosd(theta),1+cosd(theta),1],[0,0,sind(theta),0],[1,1,1]);
set(h3,'edgecolor','none');
colormap(cMap);
title('Iso-concentration');
set(F,'WindowStyle','docked');

%form data file :
x=linspace(0,Dx,nx+1); y=linspace(0,Dy,ny+1); data=zeros((nx+1)*(ny+1),9); k=1;
for i=1:nx+1
    for j=1:ny+1
        data(k,:)=[x(1,i),y(1,j),u_cell(i,j),v_cell(i,j),p_cell(i,j),w(i,j),psi_new(i,j),T_cell(i,j),C_cell(i,j)];
        k=k+1;
    end
end