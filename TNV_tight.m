function [u]=TNV_tight(f,xk,dt,lambda,ak,T,epsilon, maxIter)
f = max(double(f),1);
[Image_h,Image_w]=size(f);

int_f_Txk=sum(sum(f./(imfilter(xk,T,'symmetric','same'))));
norm_f_Txk = norm(f./(imfilter(xk,T,'symmetric','same')),'fro')^2;
TVxk = TV(xk);
F=@(p) p^3*ak*TVxk + 2*p*int_f_Txk - 2*norm_f_Txk;
Fp = @(p) 3*p^2*ak*TVxk+2*int_f_Txk; 

init = newton(256 ,F,Fp,10^-6)

%Txk = imfilter(xk,T,'symmetric','same'); 
%init = norm(f./Txk,'fro')^2/(sum(sum(f./Txk)))

%initialization: u^0=f/Txk Here, I'm epsilon regularizing to prevent /0
%init=sum(sum(f./(imfilter(xk,T,'symmetric','same'))))/(Image_h*Image_w); 
u=init.*ones(size(f)); %f./max(1,xk);
for n=1:maxIter
    Z=u.*xk;
    %These aren't square matrices
    DXF=dxf(u); %forward diff in x, for all y indexes
    DXB=dxb(u); % back diff in x, for all y indexes
    DYF=dyf(u); %forward diff in y, for all x indexes
    DYB=dyb(u); %back diff in y, for all x indexes
    DXFZ=dxf(Z); %forward diff in x, for all y indexes
    DXBZ=dxb(Z); % back diff in x, for all y indexes
    DYFZ=dyf(Z); %forward diff in y, for all x indexes
    DYBZ=dyb(Z); %
% EPSILON REGUALRIZATION of div terms
%     %u terms
    c1= 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    c2= 1./sqrt(epsilon^2 + DXB(:,2:end-1).^2 + DYF(1:end-2,:).^2);
    c3= c1; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    c4= 1./sqrt(epsilon^2 + DXF(:,1:end-2).^2 + DYB(2:end-1,:).^2);
    %z terms
    dz= 1./sqrt(epsilon^2 + DXFZ(:,2:end-1).^2 + DYFZ(2:end-1,:).^2);
    dz2= 1./sqrt(epsilon^2 + DXBZ(:,2:end-1).^2 + DYFZ(1:end-2,:).^2);
    dz3= dz; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    dz4= 1./sqrt(epsilon^2 + DXFZ(:,1:end-2).^2 + DYBZ(2:end-1,:).^2);
% fidelity term
    Tuxk=imfilter(u.*xk,T,'symmetric','same');
    fidelity = (f(2:end-1,2:end-1)./Tuxk(2:end-1,2:end-1) - 1).*(-f(2:end-1,2:end-1)./(Tuxk(2:end-1,2:end-1).^2));
    fidelity = imfilter(fidelity,T,'symmetric','same');
    fidelity = fidelity.*xk(2:end-1,2:end-1);
%update interior grid points            
    u(2:end-1,2:end-1)=1./(1+dt/2/lambda.*(c1+c2+c3+c4)+dt/2*ak.*(dz+dz2+dz3+dz4).*xk(2:end-1,2:end-1).^2).*...
        (u(2:end-1,2:end-1) - dt.*fidelity +...
        dt/2/lambda.*(c1.*u(3:end,2:end-1)+c2.*u(1:end-2,2:end-1)+c3.*u(2:end-1,3:end)+c4.*u(2:end-1,1:end-2))+...
        dt/2*ak.*xk(2:end-1,2:end-1).*(dz.*Z(3:end,2:end-1)+dz2.*Z(1:end-2,2:end-1)+dz3.*Z(2:end-1,3:end)+dz4.*Z(2:end-1,1:end-2)) );
            
%Update boundary conditions
    %Left and Right Boundary Conditions
    u(1,2:Image_w-1)=u(2,2:Image_w-1); 
    u(Image_h,2:Image_w-1)=u(Image_h-1,2:Image_w-1); 
    
    %Top and Bottom boundary conditions
    u(2:Image_h-1,1)=u(2:Image_h-1,2); 
    u(2:Image_h-1,Image_w)=u(2:Image_h-1,Image_w-1);
    
    %Corner Boundary Conditions:
    u(1,1)=u(2,2); u(1,Image_w) = u(2,Image_w-1); 
    u(Image_h,1) =u(Image_h-1,2); 
    u(Image_h,Image_w)=u(Image_h-1,Image_w-1);
    
end
end

%function dxf performs the forward difference in x, the first index-position
%of matrix M, and returns DXF, a matrix of size size(M,1)-2 x size(M,2)
%representing the forward-difference approximations to the x-derviative of
%the matrix M at interior grid points (2:end-1,1:end)
function DXF=dxf(M)
    DXF =M(3:end,:) -M(2:end-1,:);
end

%backwards difference x-derivative approximation at interior grid points
% (2:end-1,1:end)
function DXB=dxb(M) 
    DXB=M(2:end-1,:) - M(1:end-2,:);
end

%forwards difference y-derivative approximation at interior grid points
% (1:end, 2:end-1)
function DYF=dyf(M) 
    DYF=M(:,3:end) - M(:,2:end-1);
end

%backward difference y-derivative approximation at interior grid points
% (1:end, 2:end-1)
function DYB=dyb(M)
    DYB=M(:,2:end-1) - M(:,1:end-2);
end

%Cetnered difference x and y derivatives at interior grid points
function DXC=dxc(M)
    DXC = (M(3:end,:) - M(1:end-2,:))./2;
end
function DYC= dyc(M)
    DYC = (M(:,3:end) - M(:,1:end-2))./2;
end

function p = newton(p0,F,Fp,tol)
%function newton performs newton's iteration to solve F(p)+c=0.
    while abs(F(p0)) > tol
        p=p0 - F(p0)/Fp(p0);
        p0=p;
    end
end

%function to compute the non-regularized total variation of image A
function tvA=TV(A)
    A=padarray(A,[1,1],'replicate');
    DXF=A(3:end,2:end-1)-A(2:end-1,2:end-1);
    DYF=A(2:end-1,3:end)-A(2:end-1,2:end-1);
    tvA=sum(sum(sqrt(DXF.^2+DYF.^2)));
end


