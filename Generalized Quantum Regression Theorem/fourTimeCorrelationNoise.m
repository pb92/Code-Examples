function returnValue = fourTimeCorrelationNoise(rho_ss, M, N, dt, A_mat, B_mat, C_mat, D_mat, t_vec)
% Philip Daniel Blocher
% Quantum Optics Group, Aarhus University
% pblocher@phys.au.dk - AU458007
% Last updated: March 2018.
%
% Description: Generalized first order function, which calculates the 
% second order correlation function for input operators A(tA), B(tB), C(tC)
% and D(tD).
% rho_ss must be the steadystate density matrix for the system, M the 
% superoperator acting upon rho_ss (vectorized) to yield the EoM for 
% rho_ss. t_vec must be a vector with the times tA, tB, tC, and tD.
% A_mat, B_mat, C_mat and D_mat must be the matrix representation of the 
% operators A, B, C, and D.
%
% The function returns corrFun, the correlation function 
% <A(tA)B(tB)C(tC)D(tD)>.

%% Rename elements for convenience:
tA = t_vec(1);
tB = t_vec(2);
tC = t_vec(3);
tD = t_vec(4);

[n, ~] = size(A_mat);

%% Calculate the correlation function (no noise terms)

if tA == min([tA,tB,tC,tD])
    
    if tB == min([tB,tC,tD])
        
        if tC <= tD
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tB <= tC <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho_A_gd at time tA.
            rho_A_tA = rho_ss*A_mat;
            
            % Evolve rho_A_gd from time tA to time tB.
            rho_A_tB = timeEvolveSingleEXPM(M,rho_A_tA,tA,tB);         
            
            % Calculate rho_AB_gf at time tB.
            rho_AB_tB = NaN(n,n);
            
            for g = 1:1:n
                for f = 1:1:n
                    
                    value = 0;
                    for d = 1:1:n
                        value = value + B_mat(d,f) * rho_A_tB(g,d);
                    end
                    rho_AB_tB(g,f) = value;
                    
                end
            end
            
            % Evolve rho_AB_gf from time tB to time tC.
            rho_AB_tC = timeEvolveSingleEXPM(M,rho_AB_tB,tB,tC);
            
            % Calculate rho_ABC_gh at time tC.
            rho_ABC_tC = NaN(n,n);
            
            for g = 1:1:n
                for h = 1:1:n
                    
                    value = 0;
                    for f = 1:1:n
                        value = value + C_mat(f,h) * rho_AB_tC(g,f);                     
                    end
                    rho_ABC_tC(g,h) = value;
                    
                end
            end
            
            % Evolve rho_ABC_gh from time tC to time tD.
            rho_ABC_tD = timeEvolveSingleEXPM(M,rho_ABC_tC,tC,tD);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tD.
            returnValue = 0;
            
            for g = 1:1:n
                for h = 1:1:n
                    returnValue = returnValue + D_mat(h,g) * rho_ABC_tD(g,h);
                end
            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tB <= tC <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tB <= tD <= tC %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Calculate rho_A_gd at time tA.
            rho_A_tA = rho_ss*A_mat;
            
            % Evolve rho_A_gd from time tA to time tB.
            rho_A_tB = timeEvolveSingleEXPM(M,rho_A_tA,tA,tB);

            % Calculate rho_AB_gf at time tB.
            rho_AB_tB = NaN(n,n);
            
            for g = 1:1:n
                for f = 1:1:n
                    
                    value = 0;
                    for d = 1:1:n
                        value = value + B_mat(d,f) * rho_A_tB(g,d);                        
                    end
                    rho_AB_tB(g,f) = value;
                    
                end
            end
            
            % Evolve rho_AB_gf from time tB to time tD.
            rho_AB_tD = timeEvolveSingleEXPM(M,rho_AB_tB,tB,tD);
            
            % Calculate rho_AB_D_ef at time tD.
            rho_AB_D_tD = NaN(n,n);
            
            for e = 1:1:n
                for f = 1:1:n
                    
                    value = 0;
                    for g = 1:1:n
                        value = value + D_mat(e,g) * rho_AB_tD(g,f);
                    end
                    rho_AB_D_tD(e,f) = value;
                    
                end
            end
            
            % Evolve rho_AB_D_ef from time tD to time tC.
            rho_AB_D_tC = timeEvolveSingleEXPM(M,rho_AB_D_tD,tD,tC);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tC.
            returnValue = 0;
            
            for e = 1:1:n
                for f = 1:1:n
                    returnValue = returnValue + C_mat(f,e) * rho_AB_D_tC(e,f);
                end
            end            
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tB <= tD <= tC %%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
    elseif tC == min([tB,tC,tD])
        
        if tB <= tD
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tC <= tB <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%   

            % Calculate rho_A_gd at time tA.
            rho_A_tA = rho_ss*A_mat;
             
            % Evolve rho_A_gd from time tA to time tC.
            rho_A_tC = timeEvolveSingleEXPM(M,rho_A_tA,tA,tC);
            
            % Calculate rho_A_C__cdgh at time tC.
            rho_A_C__tC = NaN(n^4,1);
            
            for c = 1:1:n
                for d = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            
                            index = (c-1)*n^3 + (d-1)*n^2 + (g-1)*n + h;
                            rho_A_C__tC(index) = C_mat(c,h) * rho_A_tC(g,d);
                            
                        end
                    end
                end
            end

            % Evolve rho_A_C__cdgh from time tC to time tB.
            % rho_A_C__tB = timeEvolveDoubleEXPM(M,rho_A_C__tC,tC,tB);
            rho_A_C__tB = timeEvolveDoubleNoise(M,N,dt,rho_A_C__tC,tC,tB);
            
            % Calculate rho_ABC__gh at time tB.
            rho_ABC__tB = NaN(n,n);
            
            for g = 1:1:n
                for h = 1:1:n
                    
                    value = 0;                    
                    for c = 1:1:n
                        for d = 1:1:n
                            index = (c-1)*n^3 + (d-1)*n^2 + (g-1)*n + h;
                            value = value + B_mat(d,c) * rho_A_C__tB(index);
                        end
                    end
                    rho_ABC__tB(g,h) = value;
                    
                end
            end
            
            % Evolve rho_ABC__gh from time tB to time tD.
            rho_ABC__tD = timeEvolveSingleEXPM(M,rho_ABC__tB,tB,tD);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tD.
            returnValue = 0;
            
            for g = 1:1:n
                for h = 1:1:n
                    returnValue = returnValue + D_mat(h,g) * rho_ABC__tD(g,h);
                end
            end            
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tC <= tB <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%               
            
        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tC <= tD <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Calculate rho_A_gd at time tA.
            rho_A_tA = rho_ss*A_mat;
                   
            % Evolve rho_A_gd from time tA to time tC.
            rho_A_tC = timeEvolveSingleEXPM(M,rho_A_tA,tA,tC);          
            
            % Calculate rho_A_C__cdgh at time tC.
            rho_A_C__tC = NaN(n^4,1);
            
            for c = 1:1:n
                for d = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            
                            index = (c-1)*n^3 + (d-1)*n^2 + (g-1)*n + h;
                            rho_A_C__tC(index) = C_mat(c,h) * rho_A_tC(g,d);
                            
                        end
                    end
                end
            end           
            
            % Evolve rho_A_C__cdgh from time tC to time tD.
            % rho_A_C__tD = timeEvolveDoubleEXPM(M,rho_A_C__tC,tC,tD);
            rho_A_C__tD = timeEvolveDoubleNoise(M,N,dt,rho_A_C__tC,tC,tD);
            
            
            % Calculate rho_A_CD_cd at time tD.
            rho_A_CD_tD = NaN(n,n);
            
            for c = 1:1:n
                for d = 1:1:n
                    
                    value = 0;
                    for g = 1:1:n
                        for h = 1:1:n
                            
                            index = (c-1)*n^3 + (d-1)*n^2 + (g-1)*n + h;
                            value = value + D_mat(h,g) * rho_A_C__tD(index);
                            
                        end
                    end
                    rho_A_CD_tD(c,d) = value;
                    
                end
            end
            
            % Evolve rho_A_CD_cd from time tD to time tB.
            rho_A_CD_tB = timeEvolveSingleEXPM(M,rho_A_CD_tD,tD,tB);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tB. 
            returnValue = 0;
            
            for c = 1:1:n
                for d = 1:1:n
                    returnValue = returnValue + B_mat(d,c) * rho_A_CD_tB(c,d);
                end
            end            
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tC <= tD <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%           

        end
        
    else
        
        if tB <= tC
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tD <= tB <= tC %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Calculate rho_A_gd at time tA.
            rho_A_tA = rho_ss*A_mat;
            
            % Evolve rho_A_gd from time tA to time tD.
            rho_A_tD = timeEvolveSingleEXPM(M,rho_A_tA,tA,tD);            
            
            % Calculate rho_A__D_ed at time tD.
            rho_A__D_tD = NaN(n,n);
            
            for e = 1:1:n
                for d = 1:1:n
                    
                    value = 0;
                    for g = 1:1:n
                        value = value + D_mat(e,g) * rho_A_tD(g,d);                        
                    end
                    rho_A__D_tD(e,d) = value;
                    
                end
            end
            
            % Evolve rho_A__D_ed from time tD to time tB.
            rho_A__D_tB = timeEvolveSingleEXPM(M,rho_A__D_tD,tD,tB);         
            
            % Calculate rho_AB_D_ef at time tB.
            rho_AB_D_tB = NaN(n,n);
            
            for e = 1:1:n
                for f = 1:1:n
                    
                    value = 0;
                    for d = 1:1:n
                        value = value + B_mat(d,f) * rho_A__D_tB(e,d);
                    end
                    rho_AB_D_tB(e,f) = value;
                    
                end
            end
            
            % Evolve rho_AB_D_ef from time tB to time tC.
            rho_AB_D_tC = timeEvolveSingleEXPM(M,rho_AB_D_tB,tB,tC);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tC.
            returnValue = 0;
            
            for e = 1:1:n
                for f = 1:1:n
                    returnValue = returnValue + C_mat(f,e) * rho_AB_D_tC(e,f);
                end
            end                  
        
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tD <= tB <= tC %%%%%%%%%%%%%%%%%%%%%%%%%%%

        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tD <= tC <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho_A_gd at time tA.
            rho_A_tA = rho_ss*A_mat;
            
            % Evolve rho_A_gd from time tA to time tD.
            rho_A_tD = timeEvolveSingleEXPM(M,rho_A_tA,tA,tD);
            
            % Calculate rho_A__D_ed at time tD.
            rho_A__D_tD = NaN(n, n);
            
            for e = 1:1:n
                for d = 1:1:n
                    
                    value = 0;
                    for g = 1:1:n
                        value = value + D_mat(e,g) * rho_A_tD(g,d);                        
                    end
                    rho_A__D_tD(e,d) = value;
                    
                end
            end
            
            % Evolve rho_A__D_ed from time tD to time tC.
            rho_A__D_tC = timeEvolveSingleEXPM(M,rho_A__D_tD,tD,tC);       
            
            % Calculate rho_A_CD_cd at time tC.
            rho_A_CD_tC = zeros(n, n);
            
            for c = 1:1:n
                for d = 1:1:n
                    
                    value = 0;
                    for e = 1:1:n
                        value = value + C_mat(c,e) * rho_A__D_tC(e,d);
                    end
                    rho_A_CD_tC(c,d) = value;
                    
                end
            end
            
            % Evolve rho_A_CD_cd from time tC to time tB.
            rho_A_CD_tB = timeEvolveSingleEXPM(M,rho_A_CD_tC,tC,tB);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tB.
            returnValue = 0;
            
            for c = 1:1:n
                for d = 1:1:n
                    returnValue = returnValue + B_mat(d,c) * rho_A_CD_tB(c,d);
                end
            end                     
%%%%%%%%%%%%%%%%%%%%%%%%%% tA <= tD <= tC <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
        
    end
    
elseif tB == min([tA,tB,tC,tD])
    
    if tA == min([tA,tC,tD])
        
        if tC <= tD
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tA <= tC <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Calculate rho__B__abgf at time tB.
            rho__B__tB = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for f = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                            rho__B__tB(index) = B_mat(a,f) * rho_ss(g,b);
                        end
                    end
                end
            end
            
            % Evolve rho__B__abgf from time tB to time tA.
            % rho__B__tA = timeEvolveDoubleEXPM(M,rho__B__tB,tB,tA);
            rho__B__tA = timeEvolveDoubleNoise(M,N,dt,rho__B__tB,tB,tA);
            
            % Evaluate rho_AB__gf at time tA.
            rho_AB__tA = NaN(n,n);
            
            for g = 1:1:n
                for f = 1:1:n
                    
                    value = 0;
                    for a = 1:1:n
                        for b = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                            value = value + A_mat(b,a) * rho__B__tA(index);
                        end
                    end
                    rho_AB__tA(g,f) = value;
                    
                end
            end
            
            % Evolve rho_AB__gf from time tA to time tC.
            rho_AB__tC = timeEvolveSingleEXPM(M,rho_AB__tA,tA,tC);
            
            % Evaluate rho_ABC_gh at time tC.
            rho_ABC_tC = NaN(n,n);
            
            for g = 1:1:n
                for h = 1:1:n
                    
                    value = 0;
                    for f = 1:1:n
                        value = value + C_mat(f,h) * rho_AB__tC(g,f);
                    end
                    rho_ABC_tC(g,h) = value;
                    
                end
            end
            
            % Evolve rho_ABC_gh from time tC to time tD.
            rho_ABC_tD = timeEvolveSingleEXPM(M,rho_ABC_tC,tC,tD);
                    
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tD.
            returnValue = 0;
            
            for g = 1:1:n
                for h = 1:1:n
                    returnValue = returnValue + D_mat(h,g) * rho_ABC_tD(g,h);
                end
            end                     
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tA <= tC <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%           
           
        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tA <= tD <= tC %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Calculate rho__B__abgf at time tB.
            rho__B__tB = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for f = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                            rho__B__tB(index) = B_mat(a,f) * rho_ss(g,b);
                        end
                    end
                end
            end
            
            % Evolve rho__B__abgf from time tB to time tA.
            % rho__B__tA = timeEvolveDoubleEXPM(M,rho__B__tB,tB,tA);
            rho__B__tA = timeEvolveDoubleNoise(M,N,dt,rho__B__tB,tB,tA);
            
            % Evaluate rho_AB__gf at time tA.
            rho_AB__tA = NaN(n,n);
            
            for g = 1:1:n
                for f = 1:1:n
                    
                    value = 0;
                    for a = 1:1:n
                        for b = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                            value = value + A_mat(b,a) * rho__B__tA(index);
                        end
                    end
                    rho_AB__tA(g,f) = value;
                    
                end
            end
            
            % Evolve rho_AB__gf from time tA to time tD.
            rho_AB__tD = timeEvolveSingleEXPM(M,rho_AB__tA,tA,tD);           
            
            % Evaluate rho_AB_D_ef at time tD.
            rho_AB_D_tD = NaN(n,n);
            
            for e = 1:1:n
                for f = 1:1:n
                    
                    value = 0;
                    for g = 1:1:n
                        value = value + D_mat(e,g) * rho_AB__tD(g,f);
                    end
                    rho_AB_D_tD(e,f) = value;
                    
                end
            end
            
            % Evolve rho_AB_D_ef from time tD to time tC.
            rho_AB_D_tC = timeEvolveSingleEXPM(M,rho_AB_D_tD,tD,tC);
                    
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tC.
            returnValue = 0;
            
            for e = 1:1:n
                for f = 1:1:n
                    returnValue = returnValue + C_mat(f,e) * rho_AB_D_tC(e,f);
                end
            end                     
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tA <= tC <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%            
                 
        end
        
    elseif tC == min([tA,tC,tD])
        
        if tA <= tD
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tC <= tA <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__B__abgf at time tB.
            rho__B__tB = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for f = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                            rho__B__tB(index) = B_mat(a,f) * rho_ss(g,b);
                        end
                    end
                end
            end
            
            % Evolve rho__B__abgf from time tB to time tC.
            % rho__B__tC = timeEvolveDoubleEXPM(M,rho__B__tB,tB,tC);
            rho__B__tC = timeEvolveDoubleNoise(M,N,dt,rho__B__tB,tB,tC);

            % Evaluate rho__BC__abgh at time tC.
            rho__BC__tC = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            
                            index1 = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            
                            value = 0;
                            for f = 1:1:n
                                index2 = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                                value = value + C_mat(f,h) * rho__B__tC(index2);
                            end
                            rho__BC__tC(index1) = value;
                            
                        end
                    end
                end
            end
            
            % Evolve rho__BC__abgh from time tC to time tA.
            % rho__BC__tA = timeEvolveDoubleEXPM(M,rho__BC__tC,tC,tA);
            rho__BC__tA = timeEvolveDoubleNoise(M,N,dt,rho__BC__tC,tC,tA);
            
            % Evaluate rho_ABC__gh at time tA.
            rho_ABC__tA = NaN(n,n);
            
            for g = 1:1:n
                for h = 1:1:n
                    
                    value = 0;
                    for a = 1:1:n
                        for b = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            value = value + A_mat(b,a) * rho__BC__tA(index);
                        end
                    end
                    rho_ABC__tA(g,h) = value;
                    
                end
            end
            
            % Evolve rho_ABC__gh from time tA to time tD.
            rho_ABC__tD = timeEvolveSingleEXPM(M,rho_ABC__tA,tA,tD);
                            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tD.
            returnValue = 0;
            
            for g = 1:1:n
                for h = 1:1:n
                    returnValue = returnValue + D_mat(h,g) * rho_ABC__tD(g,h);
                end
            end                     
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tC <= tA <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tC <= tD <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__B__abgf at time tB.
            rho__B__tB = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for f = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                            rho__B__tB(index) = B_mat(a,f) * rho_ss(g,b);
                        end
                    end
                end
            end
            
            % Evolve rho__B__abgf from time tB to time tC.
            % rho__B__tC = timeEvolveDoubleEXPM(M,rho__B__tB,tB,tC);
            rho__B__tC = timeEvolveDoubleNoise(M,N,dt,rho__B__tB,tB,tC);
            
            % Evaluate rho__BC__abgh at time tC.
            rho__BC__tC = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            
                            index1 = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            
                            value = 0;
                            for f = 1:1:n
                                index2 = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                                value = value + C_mat(f,h) * rho__B__tC(index2);
                            end
                            rho__BC__tC(index1) = value;
                            
                        end
                    end
                end
            end
            
            % Evolve rho__BC__abgh from time tC to time tD.
            % rho__BC__tD = timeEvolveDoubleEXPM(M,rho__BC__tC,tC,tD);
            rho__BC__tD = timeEvolveDoubleNoise(M,N,dt,rho__BC__tC,tC,tD);
            
            % Evaluate rho__BCD_ab at time tD.
            rho__BCD_tD = NaN(n,n);
            
            for a = 1:1:n
                for b = 1:1:n
                    
                    value = 0;
                    for g = 1:1:n
                        for h = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            value = value + D_mat(h,g) * rho__BC__tD(index);
                        end
                    end
                    rho__BCD_tD(a,b) = value;
                    
                end
            end

            % Evolve rho__BCD_ab from time tD to time tA.
            rho__BCD_tA = timeEvolveSingleEXPM(M,rho__BCD_tD,tD,tA);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tA.
            returnValue = 0;
            
            for a = 1:1:n
                for b = 1:1:n
                    returnValue = returnValue + A_mat(b,a) * rho__BCD_tA(a,b);
                end
            end                     
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tC <= tD <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
        end
        
    else
        
        if tA <= tC
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tD <= tA <= tC %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__B__abgf at time tB.
            rho__B__tB = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for f = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                            rho__B__tB(index) = B_mat(a,f) * rho_ss(g,b);
                        end
                    end
                end
            end
            
            % Evolve rho__B__abgf from time tB to time tD.
            % rho__B__tD = timeEvolveDoubleEXPM(M,rho__B__tB,tB,tD);
            rho__B__tD = timeEvolveDoubleNoise(M,N,dt,rho__B__tB,tB,tD);
            
            % Evaluate rho__B_D_abef at time tD.
            rho__B_D_tD = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for e = 1:1:n
                        for f = 1:1:n
                            
                            index1 = (a-1)*n^3 + (b-1)*n^2 + (e-1)*n + f;
                            
                            value = 0;
                            for g = 1:1:n
                                index2 = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                                value = value + D_mat(e,g) * rho__B__tD(index2);
                            end
                            rho__B_D_tD(index1) = value;
                            
                        end
                    end
                end
            end
            
            % Evolve rho__B_D_abef from time tD to time tA.   
            % rho__B_D_tA = timeEvolveDoubleEXPM(M,rho__B_D_tD,tD,tA);
            rho__B_D_tA = timeEvolveDoubleNoise(M,N,dt,rho__B_D_tD,tD,tA);
            
            % Evaluate rho_AB_D_ef at time tA.
            rho_AB_D_tA = NaN(n,n);
            
            for e = 1:1:n
                for f = 1:1:n
                    
                    value = 0;
                    for a = 1:1:n
                        for b = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (e-1)*n + f;
                            value = value + A_mat(b,a) * rho__B_D_tA(index);
                        end
                    end
                    rho_AB_D_tA(e,f) = value;
                    
                end
            end
            
            % Evolve rho_AB_D_ef from time tA to time tC.
            rho_AB_D_tC = timeEvolveSingleEXPM(M,rho_AB_D_tA,tA,tC);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tC.
            returnValue = 0;
            
            for e = 1:1:n
                for f = 1:1:n
                    returnValue = returnValue + C_mat(f,e) * rho_AB_D_tC(e,f);
                end
            end
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tD <= tA <= tC %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tD <= tC <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__B__abgf at time tB.
            rho__B__tB = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for f = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                            rho__B__tB(index) = B_mat(a,f) * rho_ss(g,b);
                        end
                    end
                end
            end
            
            % Evolve rho__B__abgf from time tB to time tD.
            %  rho__B__tD = timeEvolveDoubleEXPM(M,rho__B__tB,tB,tD);
            rho__B__tD = timeEvolveDoubleNoise(M,N,dt,rho__B__tB,tB,tD);
            
            % Evaluate rho__B_D_abef at time tD.
            rho__B_D_tD = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for e = 1:1:n
                        for f = 1:1:n
                            
                            index1 = (a-1)*n^3 + (b-1)*n^2 + (e-1)*n + f;
                            
                            value = 0;
                            for g = 1:1:n
                                index2 = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + f;
                                value = value + D_mat(e,g) * rho__B__tD(index2);
                            end
                            rho__B_D_tD(index1) = value;
                            
                        end
                    end
                end
            end
            
            % Evolve rho__B_D_abef from time tD to time tC.   
            % rho__B_D_tC = timeEvolveDoubleEXPM(M,rho__B_D_tD,tD,tC);
            rho__B_D_tC = timeEvolveDoubleNoise(M,N,dt,rho__B_D_tD,tD,tC);
            
            % Evaluate rho__BCD_ab at time tC.
            rho__BCD_tC = NaN(n,n);
            
            for a = 1:1:n
                for b = 1:1:n
                    
                    value = 0;
                    for e = 1:1:n
                        for f = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (e-1)*n + f;
                            value = value + C_mat(f,e) * rho__B_D_tC(index);
                        end
                    end
                    rho__BCD_tC(a,b) = value;
                    
                end
            end
            
            % Evolve rho__BCD_ab from time tC to time tA.
            rho__BCD_tA = timeEvolveSingleEXPM(M,rho__BCD_tC,tC,tA);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tA.
            returnValue = 0;
            
            for a = 1:1:n
                for b = 1:1:n
                    returnValue = returnValue + A_mat(b,a) * rho__BCD_tA(a,b);
                end
            end                     
%%%%%%%%%%%%%%%%%%%%%%%%%% tB <= tD <= tC <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%    
            
        end
        
    end
    
elseif tC == min([tA,tB,tC,tD])
    
    if tA == min([tA,tB,tD])
        
        if tB <= tD
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tA <= tB <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__C__cbgh at time tC.
            rho__C__tC = NaN(n^4,1);
            
            for c = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            index = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            rho__C__tC(index) = C_mat(c,h) * rho_ss(g,b);
                        end
                    end
                end
            end

            % Evolve rho__C__cbgh from time tC to time tA.
            % rho__C__tA = timeEvolveDoubleEXPM(M,rho__C__tC,tC,tA);
            rho__C__tA = timeEvolveDoubleNoise(M,N,dt,rho__C__tC,tC,tA);
            
            % Evaluate rho_A_C__cdgh at time tA.
            rho_A_C__tA = NaN(n^4,1);
            
            for c = 1:1:n
                for d = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            
                            index1 = (c-1)*n^3 + (d-1)*n^2 + (g-1)*n + h;
                            
                            value = 0;
                            for b = 1:1:n
                                index2 = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                                value = value + A_mat(b,d) * rho__C__tA(index2);
                            end
                            rho_A_C__tA(index1) = value;
                            
                        end
                    end
                end
            end
                            
            % Evolve rho_A_C__cdgh from time tA to time tB.   
            % rho_A_C__tB = timeEvolveDoubleEXPM(M,rho_A_C__tA,tA,tB);
            rho_A_C__tB = timeEvolveDoubleNoise(M,N,dt,rho_A_C__tA,tA,tB);
            
            % Evaluate rho_ABC_gh at time tB.
            rho_ABC__tB = NaN(n,n);
            
            for g = 1:1:n
                for h = 1:1:n
                    
                    value = 0;
                    for c = 1:1:n
                        for d = 1:1:n
                            index = (c-1)*n^3 + (d-1)*n^2 + (g-1)*n + h;
                            value = value + B_mat(d,c) * rho_A_C__tB(index);
                        end
                    end
                    rho_ABC__tB(g,h) = value;
                end
            end

            % Evolve rho_ABC_gh from time tB to time tD.
            rho_ABC__tD = timeEvolveSingleEXPM(M,rho_ABC__tB,tB,tD);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tD.
            returnValue = 0;
            
            for g = 1:1:n
                for h = 1:1:n
                    returnValue = returnValue + D_mat(h,g) * rho_ABC__tD(g,h);
                end
            end                     
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tA <= tB <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tA <= tD <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
            % Calculate rho__C__cbgh at time tC.
            rho__C__tC = NaN(n^4,1);
            
            for c = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            index = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            rho__C__tC(index) = C_mat(c,h) * rho_ss(g,b);
                        end
                    end
                end
            end

            % Evolve rho__C__cbgh from time tC to time tA.
            % rho__C__tA = timeEvolveDoubleEXPM(M,rho__C__tC,tC,tA);
            rho__C__tA = timeEvolveDoubleNoise(M,N,dt,rho__C__tC,tC,tA);
            
            % Evaluate rho_A_C__cdgh at time tA.
            rho_A_C__tA = NaN(n^4,1);
            
            for c = 1:1:n
                for d = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            
                            index1 = (c-1)*n^3 + (d-1)*n^2 + (g-1)*n + h;
                            
                            value = 0;
                            for b = 1:1:n
                                index2 = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                                value = value + A_mat(b,d) * rho__C__tA(index2);
                            end
                            rho_A_C__tA(index1) = value;
                            
                        end
                    end
                end
            end
                            
            % Evolve rho_A_C__cdgh from time tA to time tD.
            % rho_A_C__tD = timeEvolveDoubleEXPM(M,rho_A_C__tA,tA,tD);
            rho_A_C__tD = timeEvolveDoubleNoise(M,N,dt,rho_A_C__tA,tA,tD);
            
            % Evaluate rho_A_CD_cd at time tD.
            rho_A_CD_tD = NaN(n,n);
            
            for c = 1:1:n
                for d = 1:1:n
                    
                    value = 0;
                    for g = 1:1:n
                        for h = 1:1:n
                            index = (c-1)*n^3 + (d-1)*n^2 + (g-1)*n + h;
                            value = value + D_mat(h,g) * rho_A_C__tD(index);
                        end
                    end
                    rho_A_CD_tD(c,d) = value;
                end
            end
            
            % Evolve rho_A_CD_cd from time tD to time tB.
            rho_A_CD_tB = timeEvolveSingleEXPM(M,rho_A_CD_tD,tD,tB);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tB.
            returnValue = 0;
            
            for c = 1:1:n
                for d = 1:1:n
                    returnValue = returnValue + B_mat(d,c) * rho_A_CD_tB(c,d);
                end
            end                     
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tA <= tD <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
    elseif tB == min([tA,tB,tD])
        
        if tA <= tD
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tB <= tA <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%            

            % Calculate rho__C__cbgh at time tC.
            rho__C__tC = NaN(n^4,1);
            
            for c = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            index = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            rho__C__tC(index) = C_mat(c,h) * rho_ss(g,b);
                        end
                    end
                end
            end

            % Evolve rho__C__cbgh from time tC to time tB.
            % rho__C__tB = timeEvolveDoubleEXPM(M,rho__C__tC,tC,tB);
            rho__C__tB = timeEvolveDoubleNoise(M,N,dt,rho__C__tC,tC,tB);
            
            % Evaluate rho__BC__abgh at time tB.
            rho__BC__tB = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            
                            index1 = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            
                            value = 0;
                            for c = 1:1:n
                                index2 = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                                value = value + B_mat(a,c) * rho__C__tB(index2);
                            end
                            rho__BC__tB(index1) = value;
                            
                        end
                    end
                end
            end
                                
            % Evolve rho__BC__abgh from time tB to time tA.
            % rho__BC__tA = timeEvolveDoubleEXPM(M,rho__BC__tB,tB,tA);
            rho__BC__tA = timeEvolveDoubleNoise(M,N,dt,rho__BC__tB,tB,tA);
            
            % Evaluate rho_ABC__gh at time tA.
            rho_ABC__tA = NaN(n,n);
            
            for g = 1:1:n
                for h = 1:1:n
                    
                    value = 0;
                    for a = 1:1:n
                        for b = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            value = value + A_mat(b,a) * rho__BC__tA(index);
                        end
                    end
                    rho_ABC__tA(g,h) = value;
                    
                end
            end          
            
            % Evolve rho_ABC__gh from time tA to time tD.
            rho_ABC__tD = timeEvolveSingleEXPM(M,rho_ABC__tA,tA,tD);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tD.
            returnValue = 0;
            
            for g = 1:1:n
                for h = 1:1:n
                    returnValue = returnValue + D_mat(h,g) * rho_ABC__tD(g,h);
                end
            end                     
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tB <= tA <= tD %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tB <= tD <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__C__cbgh at time tC.
            rho__C__tC = NaN(n^4,1);
            
            for c = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            index = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            rho__C__tC(index) = C_mat(c,h) * rho_ss(g,b);
                        end
                    end
                end
            end

            % Evolve rho__C__cbgh from time tC to time tB.
            % rho__C__tB = timeEvolveDoubleEXPM(M,rho__C__tC,tC,tB);
            rho__C__tB = timeEvolveDoubleNoise(M,N,dt,rho__C__tC,tC,tB);
            
            % Evaluate rho__BC__abgh at time tB.
            rho__BC__tB = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            
                            index1 = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            
                            value = 0;
                            for c = 1:1:n
                                index2 = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                                value = value + B_mat(a,c) * rho__C__tB(index2);
                            end
                            rho__BC__tB(index1) = value;
                            
                        end
                    end
                end
            end
                                
            % Evolve rho__BC__abgh from time tB to time tD.
            % rho__BC__tD = timeEvolveDoubleEXPM(M,rho__BC__tB,tB,tD);
            rho__BC__tD = timeEvolveDoubleNoise(M,N,dt,rho__BC__tB,tB,tD);
            
            % Evaluate rho__BCD_ab at time tD.
            rho__BCD_tD = NaN(n,n);
            
            for a = 1:1:n
                for b = 1:1:n
                    
                    value = 0;
                    for g = 1:1:n
                        for h = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            value = value + D_mat(h,g) * rho__BC__tD(index);
                        end
                    end
                    rho__BCD_tD(a,b) = value;
                end
            end
            
            % Evolve rho__BCD_ab from time tD to time tA.
            rho__BCD_tA = timeEvolveSingleEXPM(M,rho__BCD_tD,tD,tA);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tD.
            returnValue = 0;
            
            for a = 1:1:n
                for b = 1:1:n
                    returnValue = returnValue + A_mat(b,a) * rho__BCD_tA(a,b);
                end
            end          
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tB <= tD <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
        
    else
        
        if tA <= tB
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tD <= tA <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__C__cbgh at time tC.
            rho__C__tC = NaN(n^4,1);
            
            for c = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            index = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            rho__C__tC(index) = C_mat(c,h) * rho_ss(g,b);
                        end
                    end
                end
            end

            % Evolve rho__C__cbgh from time tC to time tD.
            % rho__C__tD = timeEvolveDoubleEXPM(M,rho__C__tC,tC,tD);
            rho__C__tD = timeEvolveDoubleNoise(M,N,dt,rho__C__tC,tC,tD);
            
            % Evaluate rho__CD_cb at time tD.
            rho__CD_tD = NaN(n,n);
            
            for c = 1:1:n
                for b = 1:1:n
                    
                    value = 0;
                    for g = 1:1:n
                        for h = 1:1:n
                            index = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            value = value + D_mat(h,g) * rho__C__tD(index);
                        end
                    end
                    rho__CD_tD(c,b) = value;
                    
                end
            end
            
            % Evolve rho__CD_cb from time tD to time tA.
            rho__CD_tA = timeEvolveSingleEXPM(M,rho__CD_tD,tD,tA);
            
            % Evaluate rho_A_CD_cd at time tA.
            rho_A_CD_tA = NaN(n,n);
            
            for c = 1:1:n
                for d = 1:1:n
                    
                    value = 0;
                    for b = 1:1:n
                        value = value + A_mat(b,d) * rho__CD_tA(c,b);
                    end
                    rho_A_CD_tA(c,d) = value;
                    
                end
            end

            % Evolve rho_A_CD_cd from time tA to time tB.
            rho_A_CD_tB = timeEvolveSingleEXPM(M,rho_A_CD_tA,tA,tB);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tB.
            returnValue = 0;
            
            for c = 1:1:n
                for d = 1:1:n
                    returnValue = returnValue + B_mat(d,c) * rho_A_CD_tB(c,d);
                end
            end          
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tD <= tA <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%

        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tD <= tB <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__C__cbgh at time tC.
            rho__C__tC = NaN(n^4,1);
            
            for c = 1:1:n
                for b = 1:1:n
                    for g = 1:1:n
                        for h = 1:1:n
                            index = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            rho__C__tC(index) = C_mat(c,h) * rho_ss(g,b);
                        end
                    end
                end
            end

            % Evolve rho__C__cbgh from time tC to time tD.
            % rho__C__tD = timeEvolveDoubleEXPM(M,rho__C__tC,tC,tD);
            rho__C__tD = timeEvolveDoubleNoise(M,N,dt,rho__C__tC,tC,tD);
            
            % Evaluate rho__CD_cb at time tD.
            rho__CD_tD = NaN(n,n);
            
            for c = 1:1:n
                for b = 1:1:n
                    
                    value = 0;
                    for g = 1:1:n
                        for h = 1:1:n
                            index = (c-1)*n^3 + (b-1)*n^2 + (g-1)*n + h;
                            value = value + D_mat(h,g) * rho__C__tD(index);
                        end
                    end
                    rho__CD_tD(c,b) = value;
                    
                end
            end
            
            % Evolve rho__CD_cb from time tD to time tB.
            rho__CD_tB = timeEvolveSingleEXPM(M,rho__CD_tD,tD,tB);
            
            % Evaluate rho__BCD_ab at time tB.
            rho__BCD_tB = NaN(n,n);
            
            for a = 1:1:n
                for b = 1:1:n
                    
                    value = 0;
                    for c = 1:1:n
                        value = value + B_mat(a,c) * rho__CD_tB(c,b);
                    end
                    rho__BCD_tB(a,b) = value;
                    
                end
            end
            
            % Evolve rho__BCD_ab from time tB to time tA.
            rho__BCD_tA = timeEvolveSingleEXPM(M,rho__BCD_tB,tB,tA);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tA.
            returnValue = 0;
            
            for a = 1:1:n
                for b = 1:1:n
                    returnValue = returnValue + A_mat(b,a) * rho__BCD_tA(a,b);
                end
            end          
%%%%%%%%%%%%%%%%%%%%%%%%%% tC <= tD <= tB <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%        
            
        end
        
    end
    
else
    
    if tA == min([tA,tB,tC])
        
        if tB <= tC
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tA <= tB <= tC %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Calculate rho__D_eb at time tD.
            rho__D_tD = D_mat*rho_ss;

            % Evolve rho__D_eb from time tD to time tA.
            rho__D_tA = timeEvolveSingleEXPM(M,rho__D_tD,tD,tA);
            
            % Evaluate rho_A__D_ed at time tA.
            rho_A__D_tA = NaN(n,n);
            
            for e = 1:1:n
                for d = 1:1:n
                    
                    value = 0;
                    for b = 1:1:n
                        value = value + A_mat(b,d) * rho__D_tA(e,b);
                    end
                    rho_A__D_tA(e,d) = value;
                    
                end
            end
            
            % Evolve rho_A__D_ed from time tA to time tB.
            rho_A__D_tB = timeEvolveSingleEXPM(M,rho_A__D_tA,tA,tB);

            % Evaluate rho_AB_D_ef at time tB.
            rho_AB_D_tB = NaN(n,n);
            
            for e = 1:1:n
                for f = 1:1:n
                    
                    value = 0;
                    for d = 1:1:n
                        value = value + B_mat(d,f) * rho_A__D_tB(e,d);
                    end
                    rho_AB_D_tB(e,f) = value;
                    
                end
            end
            
            % Evolve rho_AB_D_ef from time tB to time tC.
            rho_AB_D_tC = timeEvolveSingleEXPM(M,rho_AB_D_tB,tB,tC);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tC.
            returnValue = 0;
            
            for e = 1:1:n
                for f = 1:1:n
                    returnValue = returnValue + C_mat(f,e) * rho_AB_D_tC(e,f);
                end
            end          
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tA <= tB <= tC %%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tA <= tC <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__D_eb at time tD.
            rho__D_tD = D_mat*rho_ss;

            % Evolve rho__D_eb from time tD to time tA.
            rho__D_tA = timeEvolveSingleEXPM(M,rho__D_tD,tD,tA);
            
            % Evaluate rho_A__D_ed at time tA.
            rho_A__D_tA = NaN(n,n);
            
            for e = 1:1:n
                for d = 1:1:n
                    
                    value = 0;
                    for b = 1:1:n
                        value = value + A_mat(b,d) * rho__D_tA(e,b);
                    end
                    rho_A__D_tA(e,d) = value;
                    
                end
            end

            % Evolve rho_A__D_ed from time tA to time tC.
            rho_A__D_tC = timeEvolveSingleEXPM(M,rho_A__D_tA,tA,tC);

            % Evaluate rho_A_CD_cd at time tB.
            rho_A_CD_tC = NaN(n,n);
            
            for c = 1:1:n
                for d = 1:1:n
                    
                    value = 0;
                    for e = 1:1:n
                        value = value + C_mat(c,e) * rho_A__D_tC(e,d);
                    end
                    rho_A_CD_tC(c,d) = value;
                    
                end
            end

            % Evolve rho_A_CD_cd from time tC to time tB.
            rho_A_CD_tB = timeEvolveSingleEXPM(M,rho_A_CD_tC,tC,tB);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tB.
            returnValue = 0;
            
            for c = 1:1:n
                for d = 1:1:n
                    returnValue = returnValue + B_mat(d,c) * rho_A_CD_tB(c,d);
                end
            end          
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tA <= tC <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        end
        
    elseif tB == min([tA,tB,tC])
        
        if tA <= tC
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tB <= tA <= tC %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__D_eb at time tD.
            rho__D_tD = D_mat*rho_ss;

            % Evolve rho__D_eb from time tD to time tB.
            rho__D_tB = timeEvolveSingleEXPM(M,rho__D_tD,tD,tB);
            
            % Evaluate rho__B_D_abef at time tB.
            rho__B_D_tB = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for e = 1:1:n
                        for f = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (e-1)*n + f;
                            rho__B_D_tB(index) = B_mat(a,f) * rho__D_tB(e,b);
                        end
                    end
                end
            end
            
            % Evolve rho__B_D_abef from time tB to time tA.
            % rho__B_D_tA = timeEvolveDoubleEXPM(M,rho__B_D_tB,tB,tA);
            rho__B_D_tA = timeEvolveDoubleNoise(M,N,dt,rho__B_D_tB,tB,tA);
            
            % Evaluate rho_AB_D_ef at time tA.
            rho_AB_D_tA = NaN(n,n);
            
            for e = 1:1:n
                for f = 1:1:n
                    
                    value = 0;
                    for a = 1:1:n
                        for b = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (e-1)*n + f;
                            value = value + A_mat(b,a) * rho__B_D_tA(index);
                        end
                    end
                    rho_AB_D_tA(e,f) = value;
                    
                end
            end

            % Evolve rho_AB_D_ef from time tA to time tC.
            rho_AB_D_tC = timeEvolveSingleEXPM(M,rho_AB_D_tA,tA,tC);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tC.
            returnValue = 0;
            
            for e = 1:1:n
                for f = 1:1:n
                    returnValue = returnValue + C_mat(f,e) * rho_AB_D_tC(e,f);
                end
            end          
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tB <= tA <= tC %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tB <= tC <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__D_eb at time tD.
            rho__D_tD = D_mat*rho_ss;

            % Evolve rho__D_eb from time tD to time tB.
            rho__D_tB = timeEvolveSingleEXPM(M,rho__D_tD,tD,tB);
            
            % Evaluate rho__B_D_abef at time tB.
            rho__B_D_tB = NaN(n^4,1);
            
            for a = 1:1:n
                for b = 1:1:n
                    for e = 1:1:n
                        for f = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (e-1)*n + f;
                            rho__B_D_tB(index) = B_mat(a,f) * rho__D_tB(e,b);
                        end
                    end
                end
            end
            
            % Evolve rho__B_D_abef from time tB to time tC.   
            % rho__B_D_tC = timeEvolveDoubleEXPM(M,rho__B_D_tB,tB,tC);
            rho__B_D_tC = timeEvolveDoubleNoise(M,N,dt,rho__B_D_tB,tB,tC);
            
            % Evaluate rho__BCD_ab at time tC.
            rho__BCD_tC = NaN(n,n);
            
            for a = 1:1:n
                for b = 1:1:n
                    
                    value = 0;
                    for e = 1:1:n
                        for f = 1:1:n
                            index = (a-1)*n^3 + (b-1)*n^2 + (e-1)*n + f;
                            value = value + C_mat(f,e) * rho__B_D_tC(index);
                        end
                    end
                    rho__BCD_tC(a,b) = value;
                    
                end
            end

            % Evolve rho__BCD_ab from time tC to time tA.
            rho__BCD_tA = timeEvolveSingleEXPM(M,rho__BCD_tC,tC,tA);
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tA.
            returnValue = 0;
            
            for a = 1:1:n
                for b = 1:1:n
                    returnValue = returnValue + A_mat(b,a) * rho__BCD_tA(a,b);
                end
            end          
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tB <= tC <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
    else
        
        if tA <= tB
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tC <= tA <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__D_eb at time tD.
            rho__D_tD = D_mat*rho_ss;

            % Evolve rho__D_eb from time tD to time tC.
            rho__D_tC = timeEvolveSingleEXPM(M,rho__D_tD,tD,tC);
            
            % Evaluate rho__CD_cb at time tC.
            rho__CD_tC = NaN(n,n);
            
            for c = 1:1:n
                for b = 1:1:n
                    
                    value = 0;
                    for e = 1:1:n
                        value = value + C_mat(c,e) * rho__D_tC(e,b);
                    end
                    rho__CD_tC(c,b) = value;
                    
                end
            end

            % Evolve rho__CD_cb from time tC to time tA.
            rho__CD_tA = timeEvolveSingleEXPM(M,rho__CD_tC,tC,tA);            
            
            % Evaluate rho_A_CD_cd at time tA.
            rho_A_CD_tA = NaN(n,n);
            
            for c = 1:1:n
                for d = 1:1:n
                    
                    value = 0;
                    for b = 1:1:n
                        value = value + A_mat(b,d) * rho__CD_tA(c,b);
                    end
                    rho_A_CD_tA(c,d) = value;
                    
                end
            end

            % Evolve rho_A_CD_cd from time tA to time tB.
            rho_A_CD_tB = timeEvolveSingleEXPM(M,rho_A_CD_tA,tA,tB);           
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tB.
            returnValue = 0;
            
            for c = 1:1:n
                for d = 1:1:n
                    returnValue = returnValue + B_mat(d,c) * rho_A_CD_tB(c,d);
                end
            end          
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tC <= tA <= tB %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        else
            
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tC <= tB <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Calculate rho__D_eb at time tD.
            rho__D_tD = D_mat*rho_ss;

            % Evolve rho__D_eb from time tD to time tC.
            rho__D_tC = timeEvolveSingleEXPM(M,rho__D_tD,tD,tC);
            
            % Evaluate rho__CD_cb at time tC.
            rho__CD_tC = NaN(n,n);
            
            for c = 1:1:n
                for b = 1:1:n
                    
                    value = 0;
                    for e = 1:1:n
                        value = value + C_mat(c,e) * rho__D_tC(e,b);
                    end
                    rho__CD_tC(c,b) = value;
                    
                end
            end

            % Evolve rho__CD_cb from time tC to time tB.
            rho__CD_tB = timeEvolveSingleEXPM(M,rho__CD_tC,tC,tB);          
            
            % Evaluate rho__BCD_ab at time tB.
            rho__BCD_tB = NaN(n,n);
            
            for a = 1:1:n
                for b = 1:1:n
                    
                    value = 0;
                    for c = 1:1:n
                        value = value + B_mat(a,c) * rho__CD_tB(c,b);
                    end
                    rho__BCD_tB(a,b) = value;
                    
                end
            end

            % Evolve rho__BCD_ab from time tB to time tA.
            rho__BCD_tA = timeEvolveSingleEXPM(M,rho__BCD_tB,tB,tA);          
            
            % Evaluate the correlation function <A(tA)B(tB)C(tC)D(tD)> at
            % time tA.
            returnValue = 0;
            
            for a = 1:1:n
                for b = 1:1:n
                    returnValue = returnValue + A_mat(b,a) * rho__BCD_tA(a,b);
                end
            end          
%%%%%%%%%%%%%%%%%%%%%%%%%% tD <= tC <= tB <= tA %%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
        end
        
    end

end


end


% Time evolution function for a single dyadic product via matrix exponentiation.
function rho_final = timeEvolveSingleEXPM(M, rho_initial, t_initial, t_final)

    [n,~] = size(rho_initial);
    
    rho_final = rho_initial;
    rho_final = reshape(rho_final.',n^2,1);
    rho_final = expm(M*(t_final-t_initial))*rho_final;
    rho_final = reshape(rho_final,n,n).';
    
end


% Time evolution function for two dyadic products via matrix exponentiation (with no noise contributions)
function rho_final = timeEvolveDoubleEXPM(M, rho_initial, t_initial, t_final)

    [n,~] = size(rho_initial);
    n = sqrt(sqrt(n));

    Mbig = kron(M, eye(n^2)) + kron(eye(n^2), M);

    rho_final = rho_initial;
    rho_final = expm(Mbig*(t_final-t_initial))*rho_final;

end


% Time evolution function for two dyadic products via diff. eq. (with noise contributions)
function rho_final = timeEvolveDoubleNoise(M, N, dt, rho_initial, t_initial, t_final)

    % Initialization of support variables
    [n,~] = size(rho_initial);
    n = sqrt(sqrt(n));

    t = t_initial;

    % Initialization of return variable. We will update this at every
    % single time step.
    rho_final = rho_initial;

    % Time step loop: Continue until we reach t_final.
    while t < t_final
        
        % Update time variable.
        t = t + dt;
        
        % rho0 is the initial density matrix of this time step.
        rho0 = rho_final;
        
        % Structure: rho_{ijkl} = <(|j><i|)(t) B(0) (|l><k|)(t)>

        for i = 1:1:n
            for j = 1:1:n
                for k = 1:1:n
                    for l = 1:1:n

                        index1 = (i-1)*n^3 + (j-1)*n^2 + (k-1)*n + l;

                        drho = 0;

                        for iprime = 1:1:n
                            for jprime = 1:1:n

                                kprime = iprime; lprime = jprime;

                                index2 = (iprime-1)*n^3 + (jprime-1)*n^2 + (k-1)*n + l;
                                index3 = (i-1)*n^3 + (j-1)*n^2 + (kprime-1)*n + lprime;

                                mIndex1 = (i-1)*n + j;
                                mIndex2 = (iprime-1)*n + jprime;
                                mIndex3 = (k-1)*n + l;
                                mIndex4 = (kprime-1)*n + lprime;

                                drho = drho + M(mIndex1,mIndex2)*rho0(index2) + M(mIndex3,mIndex4)*rho0(index3);

                            end
                        end
                        
                        drho_noise = 0;

                        for noiseIndex = 1:1:n^4

                            drho = drho + rho0(noiseIndex)*N(index1,noiseIndex);
                            drho_noise = drho_noise + rho0(noiseIndex)*N(index1,noiseIndex);

                        end

                        rho_final(index1) = rho0(index1) + drho*dt;

                    end
                end
            end
        end

    end

end