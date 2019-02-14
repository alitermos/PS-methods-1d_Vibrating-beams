% ------------------------------------------------------------------
% > Non-Ideal Boundary Condition Simulation for a FIXED BEAM
% ------------------------------------------------------------------
% > By Ali Termos & Salim Laaguel
% > Contributors: Alfa Heryudono & Jinhee Lee
% > University of Massachusetts Dartmouth, Mathematics Department 
% > Date: November 21, 2018
% ------------------------------------------------------------------

% Fresh Start
clear; close all; clc;

% Initiate Variables
L=2;
n =100;
k_f = 0.95;
k_step = 0.01;

for k = 0.00:k_step:k_f
    
    % Create Chebyshev Points
    x = - cos((0:n).*(pi/n)).';
    % Create matrix T
    T = ones(n+1,n+1);
    T(:,2) = x;
    % Create matrix Tp
    Tp = zeros(n+1,n+1);
    Tp(:,2) = 1;
    % Create matrix Tp2
    Tp2 = zeros(n+1,n+1);
    % Create matrix Tp3
    Tp3 = zeros(n+1,n+1);
    % Create matrix Tp4
    Tp4 = zeros(n+1,n+1);
    % Loop to get T,Tp,Tpp Matrices
    for i=3:n+1
        T(:,i) = 2*x.*T(:,i-1)-T(:,i-2);
        Tp(:,i) = 2*T(:,i-1) + 2*x.*Tp(:,i-1) - Tp(:,i-2);
        Tp2(:,i) = 4*Tp(:,i-1) + 2*x.*Tp2(:,i-1) - Tp2(:,i-2);
        Tp3(:,i) = 6*Tp2(:,i-1) + 2*x.*Tp3(:,i-1) - Tp3(:,i-2);
        Tp4(:,i) = 8*Tp3(:,i-1) + 2*x.*Tp4(:,i-1) - Tp4(:,i-2);
    end
    % Solve for A.x = lambda.B.x
    % Form matrix A
    A = Tp4; % Form matrix A
    A(1,:) = T(1,:);
    A(2,:) = k*L*Tp2(1,:)-(1.0-k)*Tp(1,:);
    A(n,:) = T(n+1,:);
    A(n+1,:) = k*L*Tp2(n+1,:)+(1.0-k)*Tp(n+1,:);
    % Form matrix B
    B = T;
    B(1,:) = 0;
    B(2,:) = 0;
    B(n,:) = 0;
    B(n+1,:) = 0;
    % Get Eigen Values
    [V,D]=eig(A\B);
    l_s = diag(D);
    l = (16/(L.^4)).*(1./l_s);
    [l,ind]=sort(l);
    % Display Eigenfunctions
    U = T*V(:,ind);
    xstar = (L/2)*(x+1);
    
    % Generating Plots for Modes
    plot(xstar,U(:,1),'-r.');hold on;
    plot(xstar,U(:,2),'-b.');hold on;
    plot(xstar,U(:,3),'-m.');hold on;
    plot(0.0,0,'ko','MarkerSize',55,'MarkerEdgeColor',...
        [0.85,0.65,0.5],'LineWidth',2.5)
    plot(L,0,'ko','MarkerSize',55,'MarkerEdgeColor',...
        [0.85,0.65,0.5],'LineWidth',2.5)
    legend({'1^{st} mode','2^{nd} mode','3^{rd} mode'});
    xlabel('{\bf x (Beam Length)}','Interpreter','latex')
    ylabel('{\bf Eigen Modes}','Interpreter','latex')
    title({'First Three Modes of','${\bf Y^{4}(x)=\lambda.Y(x)}$',...
        'for a Fixed Beam with Non-Ideal B.C.s',...
        strcat('{\bf kl} =',num2str(k)),...
        strcat('{\bf kr} =',num2str(k))},...
        'Interpreter','latex')
    hold off; 
    
    pause(0.0002)
    
end

% Validating Data
fprintf('-------------------------------------\n')
fprintf('    kl = %.2f      kr = %.2f\n',k,k)
fprintf('-------------------------------------\n')
fprintf('Eigen Values      Natural Frequencies\n')
fprintf('---------------   -------------------\n')
fprintf('l1 = %.4e | BL1 = %.4f\n', l(1), L*nthroot(l(1),4))
fprintf('l2 = %.4e | BL2 = %.4f\n', l(2), L*nthroot(l(2),4))
fprintf('l3 = %.4e | BL3 = %.4f\n', l(3), L*nthroot(l(3),4))
fprintf('l4 = %.4e | BL4 = %.4f\n', l(4), L*nthroot(l(4),4))
fprintf('l5 = %.4e | BL5 = %.4f\n', l(5), L*nthroot(l(5),4))
fprintf('-------------------------------------\n')