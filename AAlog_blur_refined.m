%Multiplicative denoising and deblurring  (Semi-implicit)
%Author: Joel Barnett

function [u]=AAlog_blur_refined(f,xk,dt,lambda,ak,T,epsilon, maxIter)
f = max(double(f),1.1);
[Image_h,Image_w]=size(f);

%initialization: u^0=average(f/Txk )
%init=sum(sum(f./(max(1.1,imfilter(xk,T,'symmetric','same')))))/(Image_h*Image_w); 
%u=init.*ones(size(f));

%optimal init minimizing TV(log(u^0x_k)) choosing u=c/x_k, with c =
%average(f);
%u = mean(f(:))./max(1.1,xk);

% %non-constant initialization
%u = f./max(1.1,xk);

u=ones(size(f));

phi=zeros(size(f));
for n=1:maxIter
    Z=u.*xk;
    %These aren't square matrices   
    DXFphi=dxf(phi); %forward diff in x, for all y indexes
    DXBphi=dxb(phi); % back diff in x, for all y indexes
    DYFphi=dyf(phi); %forward diff in y, for all x indexes
    DYBphi=dyb(phi); %back diff in y, for all x indexes
    
    DXFZ=dxf(Z); %forward diff in x, for all y indexes
    DXBZ=dxb(Z); % back diff in x, for all y indexes
    DYFZ=dyf(Z); %forward diff in y, for all x indexes
    DYBZ=dyb(Z); %
% EPSILON REGUALRIZATION of div terms
    %phi terms
    cPhi1= 1./sqrt(epsilon^2 + DXFphi(:,2:end-1).^2 + DYFphi(2:end-1,:).^2);
    cPhi2= 1./sqrt(epsilon^2 + DXBphi(:,2:end-1).^2 + DYFphi(1:end-2,:).^2);
    cPhi3= cPhi1; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    cPhi4= 1./sqrt(epsilon^2 + DXFphi(:,1:end-2).^2 + DYBphi(2:end-1,:).^2);
     
% Compute phi update
    if n==1 %after first iteration, gradPhi is already updated 
        DXCphi=dxc(phi); DYCphi=dyc(phi);
        gradPhi = sqrt(epsilon^2+DXCphi(:,2:end-1).^2+DYCphi(2:end-1,:).^2);
        TV_phi= sum(gradPhi(:));
    end
    %S_n=integral( log(u)*phi)/integral(|gradPhi|)
    S_n=sum(sum(log(u(2:end-1,2:end-1)).*phi(2:end-1,2:end-1)))./TV_phi;
    
    phi(2:end-1,2:end-1) = 1./(1+S_n*dt*(cPhi1+cPhi2+cPhi3+cPhi4)).*(phi(2:end-1,2:end-1)+...
        dt.*(log(u(2:end-1,2:end-1)) + S_n*(cPhi1.*phi(3:end,2:end-1)+...
        cPhi2.*phi(1:end-2,2:end-1)+cPhi3.*phi(2:end-1,3:end)+cPhi4.*phi(2:end-1,1:end-2))) );
    %now, phi=phi^{n+1}, so we can recompute the gradient for u update
    DXCphi=dxc(phi); DYCphi=dyc(phi);
    gradPhi = sqrt(epsilon^2+DXCphi(:,2:end-1).^2+DYCphi(2:end-1,:).^2);
    TV_phi = sum(gradPhi(:));
% for div(grad(u)/|u||grad u|), we need to divide the ci terms by |u(i,j)|
% appropriately regularized
    %Z terms
    dz1= 1./(sqrt(epsilon^2 + DXFZ(:,2:end-1).^2 + DYFZ(2:end-1,:).^2).*abs(Z(2:end-1,2:end-1)));
    dz2= 1./(sqrt(epsilon^2 + DXBZ(:,2:end-1).^2 + DYFZ(1:end-2,:).^2).*abs(Z(1:end-2,2:end-1)));
    dz3= dz1; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    dz4= 1./(sqrt(epsilon^2 + DXFZ(:,1:end-2).^2 + DYBZ(2:end-1,:).^2).*abs(Z(2:end-1,1:end-2)));
% form |grad u|/u|u| term. Here, I tried centered differences for |grad(u)|
    %Z terms
    %DXFZ=dxc(Z); DYFZ=dyc(Z);
    ZabsZ= Z(2:end-1,2:end-1).*abs(Z(2:end-1,2:end-1));
    gradZ_ZabsZ = sqrt(DXFZ(:,2:end-1).^2 + DYFZ(2:end-1,:).^2)./ZabsZ;

% fidelity term
    Tuxk=imfilter(u.*xk,T,'symmetric','same');
    fidelity = (1./(Tuxk(2:end-1,2:end-1)) - f(2:end-1,2:end-1)./(Tuxk(2:end-1,2:end-1).^2));
    fidelity = imfilter(fidelity,T,'symmetric','same');
    fidelity = lambda.*fidelity.*xk(2:end-1,2:end-1);
%update interior grid points of u            
    u(2:end-1,2:end-1)=1./(1+dt*lambda*ak.*(dz1+dz2+dz3+dz4).*xk(2:end-1,2:end-1).^2).*...
        (u(2:end-1,2:end-1) - dt.*fidelity + dt*ak*lambda.*gradZ_ZabsZ.*xk(2:end-1,2:end-1)+...
        dt*ak*lambda.*xk(2:end-1,2:end-1).*(dz1.*Z(3:end,2:end-1)+dz2.*Z(1:end-2,2:end-1)+dz3.*Z(2:end-1,3:end)+dz4.*Z(2:end-1,1:end-2)) -dt.*phi(2:end-1,2:end-1)./(u(2:end-1,2:end-1).*TV_phi) );
    %chop values over 256 or under 1
    %     overflow_coords= (u(2:end-1,2:end-1).*xk(2:end-1,2:end-1))>256;
    %     underflow_coords = (u(2:end-1,2:end-1).*xk(2:end-1,2:end-1))<1;
    %     u(overflow_coords) = 256./xk(overflow_coords);
    %     u(underflow_coords)= 1./xk(underflow_coords);        
%Update boundary conditions for u and phi
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
    
    %Left and Right Boundary Conditions
    phi(1,2:Image_w-1)=phi(2,2:Image_w-1); 
    phi(Image_h,2:Image_w-1)=phi(Image_h-1,2:Image_w-1); 
    
    %Top and Bottom boundary conditions
    phi(2:Image_h-1,1)=phi(2:Image_h-1,2); 
    u(2:Image_h-1,Image_w)=u(2:Image_h-1,Image_w-1);
    
    %Corner Boundary Conditions:
    phi(1,1)=phi(2,2); phi(1,Image_w) = phi(2,Image_w-1); 
    phi(Image_h,1) =phi(Image_h-1,2); 
    phi(Image_h,Image_w)=phi(Image_h-1,Image_w-1);
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

%Centered difference x and y derivatives at interior grid points
function DXC=dxc(M)
    DXC = (M(3:end,:) - M(1:end-2,:))./2;
end
function DYC= dyc(M)
    DYC = (M(:,3:end) - M(:,1:end-2))./2;
end






