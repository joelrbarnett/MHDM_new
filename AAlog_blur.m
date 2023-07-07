%Multiplicative denoising and deblurring  (Semi-implicit)
%Author: Joel Barnett

function [u,energy]=AAlog_blur(f,xk,dt,lambda,T,epsilon, maxIter)

f = double(f);
[Image_h,Image_w]=size(f);

%constant initialization: u^0=average(f/Txk )
int_f_Txk=sum(sum(f./(imfilter(xk,T,'symmetric','same'))));
u=int_f_Txk/(Image_h*Image_w)*ones(size(f));

%non-constant initialization: u^0 = f./Txk
%u = f./max(1.1,imfilter(xk,T,'symmetric','same')); %poor results

energy=zeros(maxIter,1); 
for n=1:maxIter    
        %These aren't square matrices
        DXF=dxf(u); %forward diff in x, for all y indexes
        DXB=dxb(u); % back diff in x, for all y indexes
        DYF=dyf(u); %forward diff in y, for all x indexes
        DYB=dyb(u); %back diff in y, for all x indexes
% EPSILON REGUALRIZATION,
        d1= 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
        d2= 1./sqrt(epsilon^2 + DXB(:,2:end-1).^2 + DYF(1:end-2,:).^2);
        d3= d1; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
        d4= 1./sqrt(epsilon^2 + DXF(:,1:end-2).^2 + DYB(2:end-1,:).^2);
% for div(grad(u)/|u||grad u|), we need to divide the ci terms by |u(i,j)|
% appropriately regularized
        d1=d1./abs(u(2:end-1,2:end-1)); %divide by |u(i,j)|
        d2=d2./abs(u(1:end-2,2:end-1)); %divide by |u(i-1,j)|
        d3= d1; %d3./abs(u(2:end-1,2:end-1)); %divide by |u(i,j)|
        d4=d4./abs(u(2:end-1,1:end-2)); %divide by |u(i,j-1)|
% form |grad u|/u|u| term. Here, I tried centered differences for |grad(u)|
    DXC = dxc(u); DYC=dyc(u);
    UabsU= u(2:end-1,2:end-1).*abs(u(2:end-1,2:end-1));
    gradU_UabsU = sqrt(DXC(:,2:end-1).^2 + DYC(2:end-1,:).^2)./UabsU;

% fidelity term
    Tuxk=imfilter(u.*xk,T,'symmetric','same');
    fidelity = (1./Tuxk(2:end-1,2:end-1) - f(2:end-1,2:end-1)./(Tuxk(2:end-1,2:end-1).^2));
    fidelity = imfilter(fidelity,T,'symmetric','same');
    fidelity = lambda.*fidelity.*xk(2:end-1,2:end-1);

%store energy of nth iterate
    energy(n) = TV(log(u)) + lambda*sum(log(Tuxk(:)) + f(:)./Tuxk(:));
%     if n>=200 && ((abs(energy(n) -energy(n-5))/energy(n-5)<10^(-7)))
%         energy=energy(1:n);
%         break
%     end
%update interior grid points            
        u(2:end-1,2:end-1)=1./(1+dt.*(d1+d2+d3+d4)).*...
            (u(2:end-1,2:end-1) - dt.*fidelity + dt.* gradU_UabsU+...
            dt.*(d1.*u(3:end,2:end-1)+d2.*u(1:end-2,2:end-1)+d3.*u(2:end-1,3:end)+d4.*u(2:end-1,1:end-2)));
    %chop values over 256 or under 1
    %     overflow_coords= (u(2:end-1,2:end-1).*xk(2:end-1,2:end-1))>256;
    %     underflow_coords = (u(2:end-1,2:end-1).*xk(2:end-1,2:end-1))<1;
    %     u(overflow_coords) = 256./xk(overflow_coords);
    %     u(underflow_coords)= 1./xk(underflow_coords);
        
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

%Centered difference x and y derivatives at interior grid points
function DXC=dxc(M)
    DXC = (M(3:end,:) - M(1:end-2,:))./2;
end
function DYC= dyc(M)
    DYC = (M(:,3:end) - M(:,1:end-2))./2;
end

%function to compute the non-regularized total variation of image A
function tvA=TV(A)
    A=padarray(A,[1,1],'replicate');
    DXF=A(3:end,2:end-1)-A(2:end-1,2:end-1);
    DYF=A(2:end-1,3:end)-A(2:end-1,2:end-1);
    tvA=sum(sum(sqrt(DXF.^2+DYF.^2)));
end






