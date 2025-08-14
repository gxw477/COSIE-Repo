function[a,b,corr,err,res]=lsqn(X,Y,dep)
% least squares method
% [a,b,corr,err,res] = lsqn(X,Y,dep)
% dep(1) is the minimum depth value in mm
% dep(2) specifies the depth step resolution in mm, if set to 0 does not perform partial fit analysis
% dep(3) specifies what type of partial fit routines to try:
%                           0 means start from shallow end only
%                           1 means start from deep end only
%                           2 means try analysis from both ends and pick
%                           overall best result

s=size(X);
s1=size(Y);

n=min(s(1)*s(2),s1(1)*s1(2));
X=reshape(X,n,1);
Y=reshape(Y,n,1);

dec='y';
if dep(2)==0
    dec='n';
end

l=length(X);

if dec=='y'

    dat=[X(1:n),Y(1:n)];
dat=sortrows(dat);
X=dat(:,1);
Y=dat(:,2);
clear dat

    a2=inf;b2=inf;res2=inf;corr2=inf;err2=inf;
    
    if (dep(3)==0)|(dep(3)==2)
        % From bottom (shallow depth) up
        if (X(end)-X(1))<dep(1)
            f_i=l;
        else
            f_i=find((X-min(X))>=dep(1));
            f_i=min(f_i);
        end

        [Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,n]=assign_var(X(1:f_i),Y(1:f_i));
        [a,b,corr,err,res]=les(Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,n);
        
        check=0;
        Xe=X(f_i);
        
        f_i=[];
        for f_x=Xe:dep(2):max(X)
            f_i=[f_i,min(find(X>=f_x))];
        end
        f_i=[f_i,length(X)];
        
        for k=1:length(f_i)
            [Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,li]=assign_var(X(1:f_i(k)),Y(1:f_i(k)));
            [a1,b1,corr1,err1,res1]=les(Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,li);

            if (res1<res)
                a=a1;b=b1;corr=corr1;err=err1;res=res1;
            end
        end
        a2=a;b2=b;corr2=corr;err2=err;res2=res;
        
    end

    if (dep(3)==1)|(dep(3)==2)
        %-----from the top end
        
        X=flipud(X);
        Y=flipud(Y);
        
        if (X(1)-X(end))<dep(1)
            f_i=l;
        else
            f_i=find((max(X)-X)>=dep(1));
            f_i=min(f_i);
        end

        [Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,n]=assign_var(X(1:f_i),Y(1:f_i));
        [a,b,corr,err,res]=les(Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,n);
        
        check=0;
        Xe=X(f_i);
        
        f_i=[];
        for f_x=Xe:-dep(2):min(X)
            f_i=[f_i,min(find(X<=f_x))];
        end
        f_i=[f_i,length(X)];
        
        for k=1:length(f_i)
            [Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,li]=assign_var(X(1:f_i(k)),Y(1:f_i(k)));
            [a1,b1,corr1,err1,res1]=les(Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,li);

            if (res1<res)
                a=a1;b=b1;corr=corr1;err=err1;res=res1;
            end
        end
    end

    if (res>res2)
        a=a2;b=b2;corr=corr2;err=err2;res=res2;
    end
        
else
    [Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,n]=assign_var(X,Y);
    [a,b,corr,err,res]=les(Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,n);
end

% ----------------------------------------------------------------
function[Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,n]=assign_var(x,y)

Ex=sum(x(:),'omitnan');
Ey=sum(y(:),'omitnan');
y2=y.*y;
Ey2=sum(y2(:),'omitnan');
xy=x.*y;
Exy=sum(xy(:),'omitnan');
d=size(x);
n=d(1).*d(2),'omitnan';
x2=x.*x;
Ex2=sum(x2(:),'omitnan');

% ----------------------------------------------------------------
function[a,b,corr,err,res]=les(Ex,Ey,Ex2,Ey2,Exy,x,y,xy,x2,y2,n);

a=(Ex*Exy-Ey*Ex2)/(Ex*Ex-n*Ex2);
b=(n*Exy-Ex*Ey)/(n*Ex2-Ex*Ex);
corr=((n*Exy-Ex*Ey)^2)/((n*Ex2-Ex^2)*(n*Ey2-Ey^2));
err(1)=(((Ey2-(Ey^2)/n-b*(Exy-Ex*Ey/n))/(n-2))^0.5)*sqrt((1/n)+((Ex/n)^2)/(Ex2-(Ex^2)/n));
err(2)=((Ey2-(Ey^2)/n-b*(Exy-Ex*Ey/n))/((n-2)*(Ex2-(Ex^2)/n)))^0.5;
res=err(2)*(Ex2-(Ex^2)/n)^(0.5);