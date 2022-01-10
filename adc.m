clc
clear
close all

%m=[1 1 0 1 1; 0 1 1 1 1; 0 0 1 1 1; 0 1 1 0 1; 1 1 1 1 0];
%m=[1 0 1;0 1 -1];
m=[2 -4 3;0 -4 3;0 0 3];

y=(m)';
[x z]=size(y);
k=0.01;
T=10;
t= 0:k:T;

for j=1:x
    up=0; 
    for i=1:z
        if(i==1)
           s(j,:) = y(j,i)*rectangularPulse((i-1)*T/z, i*T/z,t);        
        else
           s(j,:) = s(j,:) + y(j,i)*rectangularPulse((i-1)*T/z, i*T/z,t);
        end
    end
end

phi(1,:) =  s(1,:) ./sqrt(trapz(t,s(1,:).^2 ));
figure;plot(t,phi(1,:))


for j = 2: x  
    sum=0;
    for i = 1: j-1
        sum = sum + phi(i,:)*trapz(t,phi(i,:).*s(j,:));
    end 

    g(j,:)= s(j,:) - sum;
    if abs(g(j,:))< 0.00001 
        break;
    end
    
phi(j,:) = g(j,:)/sqrt(trapz(t,g(j,:).^2));


figure;plot(t,phi(j,:));

end
%%
[wp lp]=size(phi);
for j=1:wp
    phi(j,1)=phi(j,2);
    phi(j,lp)=phi(j,lp -1);
    
    for i=2:lp-1
        if ((phi(j,i)~=phi(j,i-1) && phi(j,i)~=phi(j,i+1)))
            phi(j,i)=phi(j,i-1);
        end
    end
end

%%

a=sqrt(trapz(t,s(1,:).^2 )) ;  

for i=2:x
    for  j=1:x
        try
            b(i,j)=trapz(t,phi(j,:).*s(i,:));
        end
    end
end
[wb lb]=size(b);
if lb<3
    figure;
    plot([0 a],[0 0]);
    for j=2:wb
        hold on
        plot([0 b(j,1)],[0 b(j,2)]);
    end
    grid on;
else
    for i=1:(ceil((lb)/3))
    r=(i-1)*round((5/3 - floor(5/3))*3);
    figure;
    plot3([0 a],[0 0],[0 0]);
    for j=2:wb
        hold on
        plot3([0 b(j,1+r)],[0 b(j,2+r)] ,[0 b(j,3+r)]);
    end
    grid on; 
    end
end
