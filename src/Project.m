%Miguel De Armas 
%1351046
%Diffusion Equation Project (ADI Method)
clc;clear; close all;
N=40; %number of partitions in the x and y direction
h= (2*pi)/(N-1); %The increment value
lambda= -pi; %given

%given lengths
a_x= -pi;
a_y= -pi;
b_x= pi;
b_y= pi;

%arrays of lengths
x= a_x:h:b_x;
y= a_y:h:b_y;
deltat= h.^2/2;
t= 0:deltat:40; %abritrary time vector made
lx= length(x);
lt=length(t);

%mesh required for plotting
[X,Y]=meshgrid(x,y);

%Dirichlet boundary conditions set in vector
f_a= y.*(y-a_y).^2;
g_a= (y-a_y).^2.*cos(pi*y/a_y);
U= [f_a; zeros(lx-2,lx); g_a];

%%For solving the implicit part, use the tridiagonal function created.
%The Boundary conditions are reflected in the diagonals of the matrix
%We will use coefficients to make computation faster
a= deltat/(2*h^2);
b= (2*a+1);
c= (1-2*a);

%Diagonals for implicit in x
e= [-a*ones(1,N-2) 0];
f= [1 b*ones(1,N-2) 1];
g= [0 -a*ones(1, N-2)];
r= zeros(1,N);

%Diagonals for implicit in y
e2= [-a*ones(1, N-2) -2*a];
f2= b*ones(N);
g2= [0 -a*ones(1,N-2)];
r2= zeros(1,N);

%Vectors have been set for tridiagonal matrix function to run.

%%Now we have to take explicit in y in the first step and explicit in the
%%second step

for q=1:lt %time loop
    for j=1:lx %goes through loop and takes explicit in y with proper boundary conditions 
        for i= 2:lx-1
            if j==1
                r(i)=c*U(i,1)+2*a*U(i,2); %the ghost node doubles the j+1 value
            elseif j==lx 
                 r(i)= a*U(i,j-1) +c*U(i,j); 
            else
                r(i)=a*U(i,j-1)+c*U(i,j)+a*U(i,j+1);
            end
        end
        r= [U(1,j) r(2:lx-1) U(lx,j)]; %The Dirichlet boundary conditions are imposed
        x= tridiag(e,f,g,r);
        U(:,j)=x; %U(n+1/2) is dtermined by combining explicit in y and implicit in x
    end
    
    %%Now comes the second half step, explicit in x and implicit in y
    for i=2:lx-1
        for j=1:lx
            r2(j)=a*U(i-1,j)+c*U(i,j)+a*U(i+1,j);%explicit solution in x
        end
    x2= tridiag(e2,f2,g2,r2);
    U(i,:)=x2; %combining explicit in x and implicit in y
    end
    center(q)= U(N/2,N/2); %centrally located grid node is tracked for steady state state determination and diffusion confirmation
    figure(1)
    if mod(q,5)==0;
        surf(X,Y,U)
        drawnow
    end
end
figure(2) %the centrally located grid node plotted against time to confirm steady state 
plot(t,center) %MAY TAKE A FEW SECONDS TO COMPUTE
gi=norm(U,2);
