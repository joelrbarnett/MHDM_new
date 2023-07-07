%Denoising Algorithm (Semi-implicit)
%Author: Joel Barnett
% returns w, the argmin of lambda*integral{ f*exp(-w-xk)+(w+xk)} +
% alp*lambda*J(w+xk)+J(w), where J(u) = TV(u). 
% Here, xk is the log of the partially restored multiscale image
% u=u1*u2...uk. I.e., xk=log(u1)+log(u2)+...+log(uk), and w will
% log(u_{k+1}).
function [w]=OsherTight(f,xk,dt,lambda, alp, epsilon, maxIter)
f = double(f);
[Image_h,Image_w]=size(f);
Image_h=Image_h+2; Image_w=Image_w+2;

%Initial condition: choose to be closest to residual f/xk 
%w=zeros(Image_h,Image_w); %w is padded by 1 cell around boundary
%w(2:end-1,2:end-1)=log(f)-xk; %u=f./xk; %log(f/(u0*u1*u2...*uk)) = log(f) - log(u0)-...-log(uk) = log(f) - xk

%constInit optimal constant initialization
int_f_xk = sum(sum(f.*exp(-xk)));
w = log(int_f_xk/(size(f,1)*size(f,2)))*ones(Image_h,Image_w);

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
%extend xk with pad
    xk=padarray(xk,[1,1],'replicate');
for n=1:maxIter
    %For tighter term J(w+xk),  wx=w+xk;
    z=w+xk;
    
    %These aren't square matrices
    DXF=dxf(w); %forward diff in x, for all y indexes
    DXB=dxb(w); % back diff in x, for all y indexes
    DYF=dyf(w); %forward diff in y, for all x indexes
    DYB=dyb(w); %back diff in y, for all x indexes

    DXFz=dxf(z); %forward diff in x, for all y indexes
    DXBz=dxb(z); % back diff in x, for all y indexes
    DYFz=dyf(z); %forward diff in y, for all x indexes
    DYBz=dyb(z); %
% EPSILON REGUALRIZATION,
    %w terms
    c1= 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    c2= 1./sqrt(epsilon^2 + DXB(:,2:end-1).^2 + DYF(1:end-2,:).^2);
    c3= c1; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    c4= 1./sqrt(epsilon^2 + DXF(:,1:end-2).^2 + DYB(2:end-1,:).^2);
    
    %w+xk terms
    cz1= 1./sqrt(epsilon^2 + DXFz(:,2:end-1).^2 + DYFz(2:end-1,:).^2);
    cz2= 1./sqrt(epsilon^2 + DXBz(:,2:end-1).^2 + DYFz(1:end-2,:).^2);
    cz3= cz1; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    cz4= 1./sqrt(epsilon^2 + DXFz(:,1:end-2).^2 + DYBz(2:end-1,:).^2);
    
    %div(xk/|grad(w+xk)|)
    DXFxk=dxf(xk); DXBxk=dxb(xk); DYFxk=dyf(xk); DYBxk=dyb(xk);
    div_xk_grad_z=cz1.*DXFxk(:,2:end-1)-cz2.*DXBxk(:,2:end-1)+cz3.*DYFxk(2:end-1,:)-cz4.*DYBxk(2:end-1,:);
%update interior grid points            
    w(2:end-1,2:end-1)=1./(1+dt.*((c1+c2+c3+c4)+lambda*alp.*(cz1+cz2+cz3+cz4))).*...
        (w(2:end-1,2:end-1) - dt.*lambda.*(1-f.*exp(-w(2:end-1,2:end-1)-xk(2:end-1,2:end-1)))+...
        dt.*(c1.*w(3:end,2:end-1)+c2.*w(1:end-2,2:end-1)+c3.*w(2:end-1,3:end)+c4.*w(2:end-1,1:end-2))+...
        dt*lambda*alp.*(cz1.*w(3:end,2:end-1)+cz2.*w(1:end-2,2:end-1)+cz3.*w(2:end-1,3:end)+cz4.*w(2:end-1,1:end-2))+...
        dt*lambda*alp.*div_xk_grad_z);


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





