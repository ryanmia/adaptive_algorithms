function [y,e,w]=NLMS_equalizer(x,d,mu,Nw)

Nx=length(x);
centerTap=fix(Nw/2)+1;
Ny=Nx-Nw+1;
y=zeros(Ny,1);
e=zeros(Ny,1);
w=zeros(Nw,1);
w(centerTap)=1;

x=x(:);
X=zeros(Nw,Ny);
for i=1:Ny
    X(:,i)=x(i+Nw-1:-1:i);
end

for i=1:Ny
    y(i)=w'*X(:,i);
    e(i)=d(i+centerTap-1)-y(i);

    w=w+mu*e(i)'*X(:,i)/norm(X(:,i))^2;
end
end