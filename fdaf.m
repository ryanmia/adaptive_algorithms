function [y,e]=fdaf(d,x,mu,Nw)
    Nx=length(x);
    Nd=length(d);
    fftSize=2*Nw;
    nBlk=round(Nx/Nw);
    tLen=(nBlk+2)*Nw;
    y=zeros(tLen,1);
    e=zeros(tLen,1);
    x=[zeros(Nw,1); x(:); zeros(tLen-(Nx+Nw),1)];
    d=[zeros(Nw,1); d(:); zeros(tLen-(Nx+Nw),1)];
    W=ones(fftSize,1);
    estPwr=eps*W;
    c=.1;%fget fact
    
    for i=1:nBlk
        bInd=(i-1)*Nw+(1:Nw);
        xtemp=x((i-1)*Nw+(1:fftSize),1);
        X=fft(xtemp,fftSize);
        
        Y=W.*X;
        ytemp=ifft(Y,fftSize);
        y(bInd)=ytemp(Nw+1:fftSize);
        
        e(bInd)=d(bInd)-y(bInd);
        etemp=zeros(fftSize,1);
        etemp(Nw+1:fftSize,1)=e(bInd,1);
        E=fft(etemp,fftSize);
        
        estPwr=c*estPwr+(1-c)*abs(X).^2;
        P=diag(estPwr);
        
        V=P\conj(X).*E;
        v=ifft(V,fftSize);
        v(Nw+1:fftSize,1)=zeros(Nw,1);
        V=fft(v,fftSize);
        W=W+mu*V;
    end
    w=ifft(W,fftSize);
    w=w(1:Nw,1);
    y=y(1:Nx,1);
    y=y((Nw+1)/2:end-(Nw-1)/2);
    e=e(1:Nx,1);
end