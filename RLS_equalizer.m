function [y,e,w]=RLS_equalizer(x,d,c_reg,c_exp,Nw)
    Nx=length(x);
    centerTap=fix(Nw/2)+1;
    Ny=Nx-Nw+1;
    y=zeros(Ny,1);
    e=zeros(Ny,1);
    w=zeros(Nw,1);
    
    c_inv=1/c_exp;
    w(centerTap)=1;
    P=c_reg*eye(Nw,Nw);
    
    x=x(:);
    X=zeros(Nw,Ny);
    for i=1:Ny
        X(:,i)=x(i+Nw-1:-1:i);
    end
    
    for i=1:Ny
        y(i)=w'*X(:,i);
        K=(c_inv*P*X(:,i))/(1+c_inv*X(:,i)'*P*X(:,i));
        e(i)=d(i+centerTap-1)-w'*X(:,i);
        w=w+K*e(i)';
        P=c_inv*P-c_inv*K*X(:,i)'*P;    
    end
end