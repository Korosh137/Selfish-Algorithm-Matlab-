% The SAL Algorithm (Reinforcment Learning mechanism of the Selfish Algorithm)
% Korosh Mahmoodi 111618
% Cite: Mahmoodi, Korosh, Bruce J. West, and Cleotilde Gonzalez. "Selfish algorithm and emergence of collective intelligence." Journal of Complex Networks 8.3 (2020): cnaa019.
tic
clc ;
clear all ;
close all ;
% In SAL, there is one decision mechanism for each agent:
% Defect or Cooperate (1 or 2)
% For this decision, there is a moving threshold that divides the interval [0 1] into two parts, P and 1-P, where P is the propensity of the corresponding decision
Trial = 1e5 ; % Simulaton length
Size = 20 ; % Number of the agents (players)

TracPofC = zeros(Trial, 1) ;
% Constant relating payoff to the change of the propensity/probability of
Chi_RL = 0.1 ; %  SAL mechanism (Cooperation or Defection decision)

% Prisoner's Dilema game's payoffs
s = 0 ; % s is the payoff of the cooperator agent if another agent defected
p = 0 ; % p is the payoff of the agents if both defected
tt = 0.9 ; % (1+ tt) is the payoff of the defector agent if the other agent cooperated (tt < s + 1 and  tt < 1)

Ratio_CC = zeros(Trial, 1) ; CC = 0 ; % Ration of mutual cooperation

Out1 = zeros(Size, 1) ; % Previous payoff
Out2 = zeros(Size, 1) ; % Current payoff

D_C_decision = zeros(Size, 1) ; % D (defection) as 1 and C (cooperation) as 2
P_RL = zeros(Size, Size, 2) ; % P_RL and (1-P_RL) are the propensity of decison D and C, respectively

P_RL_Used = zeros(Size, 1) ; % Tracks if the RL decision mechanism has used (2) in the final decision making or not (1)

% Initial conditions
for jjj = 1 : Size
    r = rand ;
    if r < 0.5
        D_C_decision(jjj) = 1 ;
    else
        D_C_decision(jjj) = 2 ;
    end
end
for i = 1 : Size
    for j = 1 : Size
        for k = 1 : 2
            if i ~= j
                P_RL(i, j, k) = 0.5 ;
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
        if r < P_RL_m_cumsum(1, 1, i)
            D_C_decision(m) = i ;
            break
        end
    end
    P_RL_Used(m) = 2 ; % The RL mechanism is used in decision making

    P_RL_n_cumsum = cumsum(P_RL(n, m, :) ) ;
    r = rand ;
    for i = 1 : 2
        if r < P_RL_n_cumsum(1, 1, i)
            D_C_decision(n) = i ;
            break
        end
    end
    P_RL_Used(n) = 2 ; % The RL mechanism is used in decision making


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
    for k = 1 : 2
        PRLm(k) = P_RL(m, n, k) ;
        PRLn(k) = P_RL(n, m, k) ;
    end

    TracPofC(ti) = PRLm(2) ;

    P_RL(m, n, :) = SA_Update(P_RL_Used(m), PRLm, D_C_decision(m), Out2(m), Out1(m), Chi_RL) ;
    P_RL(n, m, :) = SA_Update(P_RL_Used(n), PRLn, D_C_decision(n), Out2(n), Out1(n), Chi_RL) ;


    Out1(m) = Out2(m) ;
    Out1(n) = Out2(n) ;

    if  D_C_decision(m) == 2  &&  D_C_decision(n) == 2
        CC = CC + 1 ;
    end
    Ratio_CC(ti) = CC/ ti ;


end   % end trial


plot(Ratio_CC)
xlabel('Trial (time)'), ylabel('Ratio of mutual cooperation');



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



