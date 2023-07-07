%ROF Denoising Algorithm (Semi-implicit)
%Author: Joel Barnett

function [w]=Osher(f,xk,dt,lambda, epsilon, maxIter)
f = double(f);
[Image_h,Image_w]=size(f);
Image_h=Image_h+2; Image_w=Image_w+2;
%Initial condition


%initialize picking w to be nearest residual remaining. 
%w=zeros(Image_h,Image_w); %w is padded by 1 cell around boundary
%w(2:end-1,2:end-1)=log(f)-xk; %u=f./xk; %log(f/(u0*u1*u2...*uk)) = log(f) - log(u0)-...-log(uk) = log(f) - xk

%constInit optimal constant initialization
int_f_xk = sum(sum(f.*exp(-xk)));
w = log(int_f_xk/(size(f,1)*size(f,2)))*ones(Image_h,Image_w);

% u=ones(size(u0))*sum(sum(log(u0./exp(xk))))/(size(u0,1)*size(u0,2)); %initialization v2
% u = ones(size(u0))*log(sum(sum(u0./exp(xk)))/(size(u0,1)*size(u0,2))); %initialization v3

%set intial boundaries values to satisfy bdy condition
%Left and Right Boundary Conditions
    w(1,2:Image_w-1)=w(2,2:Image_w-1); 
    w(Image_h,2:Image_w-1)=w(Image_h-1,2:Image_w-1); 
    
    %Top and Bottom boundary conditions
    w(2:Image_h-1,1)=w(2:Image_h-1,2); 
    w(2:Image_h-1,Image_w)=w(2:Image_h-1,Image_w-1);
    
    %Corner Boundary Conditions:
    w(1,1)=w(2,2); w(1,Image_w) = w(2,Image_w-1); 
    w(Image_h,1) =w(Image_h-1,2); 
    w(Image_h,Image_w)=w(Image_h-1,Image_w-1);

for n=1:maxIter
        %These aren't square matrices
        DXF=dxf(w); %forward diff in x, for all y indexes
        DXB=dxb(w); % back diff in x, for all y indexes
        DYF=dyf(w); %forward diff in y, for all x indexes
        DYB=dyb(w); %back diff in y, for all x indexes
% EPSILON REGUALRIZATION,
        c1= 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
        c2= 1./sqrt(epsilon^2 + DXB(:,2:end-1).^2 + DYF(1:end-2,:).^2);
        c3= c1; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
        c4= 1./sqrt(epsilon^2 + DXF(:,1:end-2).^2 + DYB(2:end-1,:).^2);
%update interior grid points            
        w(2:end-1,2:end-1)=1./(1+dt.*(c1+c2+c3+c4)).*...
            (w(2:end-1,2:end-1) - dt.*lambda.*(1-f.*exp(-w(2:end-1,2:end-1)-xk))+...
            dt.*(c1.*w(3:end,2:end-1)+c2.*w(1:end-2,2:end-1)+c3.*w(2:end-1,3:end)+c4.*w(2:end-1,1:end-2)));


%Update boundary conditions
    %Left and Right Boundary Conditions
    w(1,2:Image_w-1)=w(2,2:Image_w-1); 
    w(Image_h,2:Image_w-1)=w(Image_h-1,2:Image_w-1); 
    
    %Top and Bottom boundary conditions
    w(2:Image_h-1,1)=w(2:Image_h-1,2); 
    w(2:Image_h-1,Image_w)=w(2:Image_h-1,Image_w-1);
    
    %Corner Boundary Conditions:
    w(1,1)=w(2,2); w(1,Image_w) = w(2,Image_w-1); 
    w(Image_h,1) =w(Image_h-1,2); 
    w(Image_h,Image_w)=w(Image_h-1,Image_w-1);
end
%return resized image: 
    w=w(2:end-1,2:end-1);
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





