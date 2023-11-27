% The SALT algorithm (Reinforcment learning and Trust mechanisms of the Selfish Algorithm)
% Korosh Mahmoodi 111618
% Cite: Mahmoodi, Korosh, Bruce J. West, and Cleotilde Gonzalez. "Selfish algorithm and emergence of collective intelligence." Journal of Complex Networks 8.3 (2020): cnaa019.
tic
clc ;
clear all ;
close all ;
% In SALT, there are two decision mechanisms for each agent:
% Defect or Cooperate (1 or 2)
% not-Trust or Trsut (1 or 2)
% For each decision, there is a moving threshold that divides the interval [0 1] into two parts, P and 1-P, where P is the propensity of the corresponding decision
Trial = 1e5 ; % Simulaton length
Size = 20 ; % Number of the agents (players)

% Constants relating payoff to the change of the propensity/probability of
Chi_RL = 0.1 ; %  SAL mechanism (Cooperation or Defection decision)
Chi_Trust = 0.1 ; % SAT mechanism (trust decision)

% Prisoner's Dilema game's payoffs
s = 0 ; % s is the payoff of the cooperator agent if another agent defected
p = 0 ; % p is the payoff of the agents if both defected
tt = 0.9 ; % (1+ tt) is the payoff of the defector agent if the other agent cooperated (tt < s + 1 and  tt < 1)

Ratio_CC = zeros(Trial, 1) ; CC = 0 ; % Ration of mutual cooperation
TracPofC = zeros(Trial, 1) ;

Out1 = zeros(Size, 1) ; % Previous payoff
Out2 = zeros(Size, 1) ; % Current payoff

D_C_decision = zeros(Size, 1) ; % D (defection) as 1 and C (cooperation) as 2
Trust_decision = zeros(Size, 1) ; % Not-Trust as 1 and Trust as 2

P_RL = zeros(Size, Size -1, 2) ; % P_RL and (1-P_RL) are the propensity of decison D and C, respectively
P_Trust = zeros(Size,Size -1, 2) ; % P_Trust and (1-P_Trust) are the propensity of decison "not to trust" and "trust" the decison of another agent, respectively

P_RL_Used = zeros(Size, 1) ; % Tracks if the RL decision mechanism has used (2) in the final decision making or not (1)
P_Trust_Used = zeros(Size, 1) ; % Tracks if the Trust decision mechanism has used (2) in the final decision making or not (1)

% Initial conditions
for jjj = 1 : Size
    r = rand ;
    if r < 0.5
        D_C_decision(jjj) = 1 ;
    else
        D_C_decision(jjj) = 2 ;
    end

    if r < 0.5
        Trust_decision(jjj) = 1 ;
    else
        Trust_decision(jjj) = 2 ;
    end
end
for i = 1 : Size
    for j = 1 : Size
        for k = 1 : 2
            if i ~= j
                P_RL(i , j, k) = 0.5 ;
                P_Trust(i , j, k) =  0.5 ;
            end
        end
    end
end


gh =  1 ;
for ti = 2 : Trial

% Selecting the two agents to play
    m = floor(1 + Size * rand) ; 
    n = floor(1 + Size * rand) ;

    while m == n
        n = floor(1 + Size * rand) ;
    end

    % Decision of C or D
    P_RL_m_cumsum = cumsum(P_RL(m, n, :) ) ;
    r = rand ;
    for i = 1 : 2
        if r < P_RL_m_cumsum(1,1,i)
            D_C_decision(m) = i ;
            break
        end
    end

    P_RL_n_cumsum = cumsum(P_RL(n, m, :) ) ;
    r = rand ;
    for i = 1 : 2
        if r < P_RL_n_cumsum(1,1,i)
            D_C_decision(n) = i ;
            break
        end
    end

    P_RL_Used(m) = 2 ; % RL decision mechanism used
    P_RL_Used(n) = 2 ;

    %  Decision to not-trust or trust the decision of the other agent
    P_Trust_m_cumsum = cumsum(P_Trust(m , n, :) ) ;
    r = rand ;
    for i = 1 : 2
        if r < P_Trust_m_cumsum(1, 1, i)
            Trust_decision(m) = i ;
            break
        end
    end

    P_Trust_n_cumsum = cumsum(P_Trust(n ,m , :) ) ;
    r = rand ;
    for i = 1 : 2
        if r < P_Trust_n_cumsum(1,1,i)
            Trust_decision(n) = i ;
            break
        end
    end


    D_C_m = D_C_decision(m) ;
    D_C_n = D_C_decision(n) ;

    if  Trust_decision(m) == 2
        D_C_decision(m) = D_C_n ;

        P_Trust_Used(m) = 2 ; % Trust decision mechanism used
        P_RL_Used(m) = 1 ; % RL decision mechanism eliminated
    end
    if  Trust_decision(n) == 2
        D_C_decision(n) = D_C_m ;

        P_Trust_Used(n) = 2 ;
        P_RL_Used(n) = 1 ;
    end


    % Payoffs from the Prisoner's Dilemma game
    if D_C_decision(m) == 2
        ggn = 1 ;
    else
        ggn = 0 ;
    end

    if D_C_decision(n) == 2
        ggm = 1 ;
    else
        ggm = 0 ;
    end

    if D_C_decision(m) == 2
        Out2(m) = ggm * (1) + (1 - ggm) * (-s) ;
    else
        Out2(m) = ggm * (1 + tt) + (1 - ggm) * (p) ;
    end

    if D_C_decision(n) == 2
        Out2(n) = ggn * (1) + (1 - ggn) * (-s) ;
    else
        Out2(n) = ggn * (1 + tt) + (1 - ggn) * (p) ;
    end


    %   Updating the thresholds

    PRLm = zeros(1, 2);
    PRLn = zeros(1, 2);
    PTrustm = zeros(1, 2);
    PTrustn = zeros(1, 2);
    for k = 1 : 2
        PRLm(k) = P_RL(m, n, k) ;
        PRLn(k) = P_RL(n, m, k) ;

        PTrustm(k) = P_Trust(m, n, k) ;
        PTrustn(k) = P_Trust(n, m, k) ;
    end

            TracPofC(ti) = PRLm(2) ;


    P_RL(m, n, :) = SA_Update(P_RL_Used(m), PRLm , D_C_decision(m), Out2(m), Out1(m), Chi_RL) ;
    P_RL(n, m, :) = SA_Update(P_RL_Used(n), PRLn , D_C_decision(n)  , Out2(n) , Out1(n), Chi_RL) ;


    P_Trust(m, n, :) = SA_Update(P_Trust_Used(m), PTrustm, Trust_decision(m) , Out2(m), Out1(m), Chi_Trust) ;
    P_Trust(n, m, :) = SA_Update(P_Trust_Used(n), PTrustn, Trust_decision(n) , Out2(n), Out1(n), Chi_Trust) ;



    Out1(m) = Out2(m) ;
    Out1(n) = Out2(n) ;

    if  D_C_decision(m) == 2  &&  D_C_decision(n) == 2
        CC = CC + 1 ;
    end
    Ratio_CC(ti) = CC/ ti ;


end   % end trial


plot(Ratio_CC)
xlabel('Trial (time)'), ylabel('Ratio of mutual cooperation');


toc



function P = SA_Update(MechanismUsed, pp, Decision, Pay, Paybefore, Chi)

% This function updates the propensities of the decision mechanisms which
% contributed to the final decision of the agent
% MechanismUsed: 2 if the decision mechanism is used and 1 otherwise
% pp: The input propensities
% Decision: Indicates the index of the propensity that resulted in the decision
% Pay: Payoff of the agent at the current trial
% Paybefore: Payoff of the agent at the previous trial
% Chi: Constants relating payoff to the change of the propensity/probability

if MechanismUsed == 2 % If the decision mechanism has used (2) in the final decision making or not (1)
    g = length(pp) ;
    P = zeros(1 , g) ;

    Delta =  Chi * (Pay - Paybefore) / (abs(Pay) + abs(Paybefore)) ;
    if   Pay == Paybefore
        P = pp ;
    else

        P(Decision) = pp(Decision) + Delta ;

        for i = 1 : g
            if i ~= Decision
                P(i) = pp(i) - Delta/(g-1) ;
            end
        end

        % Boundary condition
        for j = 1 : g
            if  P(j) < 0
                P(j) = 0 ;
            end
            if  P(j) > 1
                P(j) = 1 ;
            end
        end

        % Normalization
        R = cumsum(P);
        P = P ./ R(g) ;

    end

else
    P = pp ;
end

end




