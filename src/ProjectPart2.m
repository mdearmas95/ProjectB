%%Miguel De Armas
%1351046
%%Diffusion Equation Project(Explicit Method)
clc;clear;close all;
N=40;
h=2*pi/(N-1);
%given lengths
a_x= -pi;
a_y= -pi;
b_x= pi;
b_y= pi;

%arrays of lengths
deltat= h.^2/4;
t= 0:deltat:40;
x= a_x:h:b_x;
y= a_y:h:b_y;
%coefficient shortcuts
a= deltat/(h.^2);
c=1-4*a;
%length shortcuts
ly=length(y);
lt=length(t);

%Boundary conditions imposed on boundaries and initial condition imposed on
%inner points
f_a= y.*(y-a_y).^2;
g_a= (y-a_y).^2.*cos(pi*y/a_y);
U= [f_a;zeros(N-2,N);g_a];
Uzeros= zeros(ly,ly);
center=zeros(1,lt);
[X,Y]=meshgrid(x,y);


%time for loop computations
for q=1:lt
    for j=1:ly
        for i=2:ly-1
            if j==1
                Uzeros(i,1)=a*U(i-1,1)+c*U(i,1)+a*U(i+1,1)+2*a*U(i,2);%Ghost node doubles the final term
            elseif j==ly
                Uzeros(i,j)=a*U(i,j-1)+a*U(i-1,j)+c*U(i,j)+a*U(i+1,j); %the general explicit form without the j+1 value
            else 
                Uzeros(i,j)=a*U(i,j-1)+a*U(i-1,j)+c*U(i,j)+a*U(i+1,j)+a*U(i,j+1);%the general explicit form
            end
        end
        center(q)=Uzeros(N/2,N/2);
    end
    U=[U(1,:); Uzeros(2:end-1,:);U(end,:)]; %the final U vector with all known displacements
    
    if mod(q,10)==0 % Throughout the time steps plots are generated to confirm diffusion
        figure(1)
        surf(X,Y,U)
        drawnow
    end
end
figure(2) %The centrally located grid node is tracked for use of steady state determination
plot(t,center)
                
                