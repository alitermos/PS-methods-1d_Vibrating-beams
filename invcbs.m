% ---------------------------------------------------------------------------------
% > Function Name: invcbs (INVERSE CANTILEVER BEAM SOLVER)
% ---------------------------------------------------------------------------------
%
% > Purpose: Inverse Solver that Computes the Assesement Parameter kl of
%            a Catilever Beam given First Three Natural Frequency Measurments.
%            ( Accurately Detects up to 97.7% damage )
%
% > Function Call: [kl] = invcbs(bl_1, bl_2, bl_3)   
%
% > Input:
%     First 3 Non-Dimentional Natural Frequencies
%     bl_1, bl_2, bl_3
%
% > Outputs:
%     kl: Assessement Parameter with range: [0,1)
%         - 0 Corresponds to Ideal Cantilever (no Damage)
%         - (0,1) Domain Corresponds to Severity of Damage
%         - 1 Corresponds to a Non-Functional Cantilever (Simply-Supported Model),
%           where "forcbs" no Longer Applies 
% ---------------------------------------------------------------------------------
% > By Ali Termos & Salim Laaguel
% > Contributors: Alfa Heryudono & Jinhee Lee
% > University of Massachusetts Dartmouth, Mathematics Department 
% > Date: November 21, 2018
% ---------------------------------------------------------------------------------

function [kl] = invcbs(wm_1,wm_2,wm_3)    
 
    kl = 0.5; % Initial Guess
    delta = 0.02202; % Step Size
    iter = 30; % Iterations
    r = 0.25; % Relaxation Factor
    
    for i = 1:iter % Loop to get Appropriate 'kl' Value
        
        % Calling Forward Solver with Initial Guess
        [w_1, w_2, w_3] = forcbs(kl); 
        % Calling Forward Solver with Initial Guess plus Step Size
        [wkl_1, wkl_2, wkl_3] = forcbs(kl+delta); 
        
        % Forming the Jacobian Matrix
        J(1,1) = (wkl_1-w_1) / delta;
        J(2,1) = (wkl_2-w_2) / delta;
        J(3,1) = (wkl_3-w_3) / delta;
        
        % Finding the Residual
        R(1,1) = w_1-wm_1;
        R(2,1) = w_2-wm_2;
        R(3,1) = w_3-wm_3;
        
        % Calculating the 'new kl' Approximation
        kl_new = -J \ R;
        
        % Controling 'kl' Propagation with Relaxation Factor
        if kl + (kl_new)*r < 0
            kl = 0;
        else
            kl = kl + (kl_new)*r;
        end
        
        % Printing 'kl' Convergence
        fprintf('%.5f\n', kl);   
        
    end
    
end