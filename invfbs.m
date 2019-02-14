% -------------------------------------------------------------------------------------
% > Function Name: invfbs (INVERSE FIXED BEAM SOLVER)
% -------------------------------------------------------------------------------------
%
% > Purpose: Inverse Solver that Computes the Assesement Parameters kl & kr of
%            a Fixed Beam given First Three Natural Frequency Measurments.
%            ( Accurately Detects up to 97.7% damage )
%
% > Function Call: [kl,kr] = invfbs(bl_1, bl_2, bl_3)   
%
% > Input:
%     First 3 Non-Dimentional Natural Frequencies
%     bl_1, bl_2, bl_3
%
% > Outputs:
%     kl & kr: Assessement Parameter with range: [0,1)
%             - 0 Corresponds to Ideal Cantilever (no Damage)
%             - (0,1) Domain Corresponds to Severity of Damage
%             - 1 Corresponds to a Non-Functional Cantilever (Simply-Supported Model),
%               where "forfbs" no Longer Applies 
% 
% -------------------------------------------------------------------------------------
% > By Ali Termos & Salim Laaguel
% > Contributors: Alfa Heryudono & Jinhee Lee
% > University of Massachusetts Dartmouth, Mathematics Department 
% > Date: November 21, 2018
% -------------------------------------------------------------------------------------

function [kl,kr] = invfbs(wm_1,wm_2,wm_3)    
    
    % Initial Guesses
    kl= 0.50;
    kr= 0.50;
    delta= 0.02202; % Step Size
    iter= 90; % Iterations
    r= 0.25; % Relaxation Factor
    
    for i=1:iter % Loop to get Appropriate 'kl & kr' Values
        
        % Calling Forward Solver with Initial Guesses
        [w_1, w_2, w_3]= forfbs(kl,kr);
        % Calling Forward Solver with Initial Guesses plus Step Size
        [wkl_1, wkl_2, wkl_3]= forfbs(kl+delta,kr);
        [wkr_1, wkr_2, wkr_3]= forfbs(kl,kr+delta);
        
        % Forming the Jacobian Matrix
        J(1,1)= (wkl_1-w_1)/delta;
        J(2,1)= (wkl_2-w_2)/delta;
        J(3,1)= (wkl_3-w_3)/delta;
        J(1,2)= (wkr_1-w_1)/delta;
        J(2,2)= (wkr_2-w_2)/delta;
        J(3,2)= (wkr_3-w_3)/delta;
        
        % Finding the Residual
        R(1,1)= w_1-wm_1;
        R(2,1)= w_2-wm_2;
        R(3,1)= w_3-wm_3;
        
        % Applying Singular Value Decomposition to Jacobian Matrix
        [U,S,V]= svd(J,0);
        Ut= U';
        Dk= -(V/S*Ut*R);
        
        % Controling 'kl' Propagation with Relaxation Factor
        if kl+Dk(1)*r <0
            kl= 0;
        else
            kl= kl+Dk(1)*r;
        end
    
        % Controling 'kr' Propagation with Relaxation Factor
        if kr+Dk(2)*r <0
            kr= 0;
        else
            kr= kr+Dk(2)*r;
        end
        
        % Printing 'kl & kr' Convergence
        fprintf('%.4f | %.4f\n',kl,kr);
        disp('---------------');
        
    end
    
end