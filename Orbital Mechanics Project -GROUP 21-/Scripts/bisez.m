function [xvect,iter]=bisez(a,b,toll,fun)
xvect=[];
iter=-1;
err=toll+1;
nmax=ceil(log2((b-a)/toll)-1);
while (iter<nmax && err>toll)
    iter=iter+1;
    x=(b+a)/2;
    fc=fun(x);
    if (fc==0)
        err=0;
    else
        err=abs(fc);
    end
    xvect=[xvect;x];
    if (fc*fun(a)>0)
        a=x;
    else
        b=x;
    end
end