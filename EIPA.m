%Yousef Nour
%Elec 4700-LAB5-Harmonic Wave Equation in 2D FD and Modes-EIPA
close all
clear all
clc

set(0,'DefaultFigureWindowStyle', 'docked')
set(0,'defaultaxesfontsize', 20)
set(0, 'defaultaxesfontname', 'Times New Roman')
set(0,'DefaultLineLineWidth', 2);

nx = 50;
ny = 50;
V = zeros(nx,ny);
G = sparse(nx*ny, nx*ny);

Inculsion = 0;

map = @(i,j) j + (i - 1)*ny;
for i=1:nx
    for j=1:ny
        mm = map(i,j);
        nxNeg = map(i-1,j);
        nxPos = map(i+1,j);
        nyNeg = map(i,j-1);
        nyPos = map(i,j+1);
        
        if (i == 1 || i == nx)
            G(mm,mm) = 1 / 1^2;
        elseif (j == 1 || j == ny)
            G(mm,mm) = 1 / 1^2;
            
        elseif (i > 10 & i <20 & j > 10 & j <20)
            G(mm,mm) = -1 / 1^2 + -1 / 1^2;
            G(mm,nxNeg) = 1 / 1^2;
            G(mm,nxPos) = 1 / 1^2;
            G(mm,nyNeg) = 1 / 1^2;
            G(mm,nyPos) = 1 / 1^2;
        else
            G(mm,mm) = -2 / 1^2 + -2 / 1^2;
            G(mm,nxNeg) = 1 / 1^2;
            G(mm,nxPos) = 1 / 1^2;
            G(mm,nyNeg) = 1 / 1^2;
            G(mm,nyPos) = 1 / 1^2;
        end
    end 
end 

figure('name', 'Matrix')
spy(G)

nmodes = 20; %eigen values
[E,D] = eigs(G,nmodes,'SM');

figure('name', 'EigenValues')
plot(diag(D),'*');

np = ceil(sqrt(nmodes));
figure('name','Modes')
for k=1:nmodes
    M = E(:,k);
    for i=1:nx
        for j=1:ny
            n = i+(j-1)*nx;
            V(i,j) = M(n);
        end
        subplot(np,np,k), surf(V,'linestyle', 'none')
        title(['EV = ' num2str(D(k,k))])
    end
end