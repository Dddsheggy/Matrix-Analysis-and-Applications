clear all;clc;close all;
N=100;     %signal length
sigma=0;
m = N/2;
L=N-m;
s_amp=[1.31*exp(1i*pi/4),2.07*exp(1i*pi/3),1.88*exp(1i*pi/5)];
s_omega=[0.12*pi,0.37*pi,0.72*pi];
MSE=[];

for t=1:100
    sigma=sigma+0.1;
    tmp=0;
    for K=1:200
        x=zeros(1,N);%initialize
        w = sqrt(sigma)*randn(1,N);
        n = [1:N];
        for slen=1:length(s_omega)
            x = x+s_amp(slen)*exp(1j*s_omega(slen)*n)  ;
        end
        x=x+w;
        for n = 1:L
            X(:,n) = x(n:(n+m-1));
        end
        for n = 1:L
            Y(:,n) = x((n+1):(n+m));
        end
        %Rxx\Rxy
        Rxx = 0;
        for i = 1:L
            Rxx = Rxx+X(:,i)*X(:,i)';
        end
        Rxx = Rxx/L;
        Rxy = 0;
        for i = 1:L
            Rxy = Rxy+X(:,i)*Y(:,i)';
        end
        Rxy = Rxy/L;
        [A,B] = eig(Rxx);
        var = min(diag(B));
        I = eye(m);
        Z = diag(ones(1,m-1),-1);
        Cxx = Rxx - I*var;
        Cxy = Rxy - Z*var;
        [~,B] = eig(Cxx,Cxy);
        f=angle(diag(B));
        [~,fpos]=sort(abs(abs(diag(B))-1));
        f=f(fpos);
        fval=f(f>0);
        omega_est=sort(fval(1:length(s_amp)));
        tmp=tmp+sum((abs(omega_est-s_omega')).^2);
    end
    MSE=[MSE,tmp/K];
end

xax=1:length(MSE);
xax=xax/10;
plot(xax,MSE)