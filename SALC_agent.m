% The SALC algorithm (Reinforcement Learning (SAL) and Connection (SAC) mechanisms of the Selfish Algorithm)
% Korosh Mahmoodi 111618
% Cite: Mahmoodi, Korosh, Bruce J. West, and Cleotilde Gonzalez. "Selfish algorithm and emergence of collective intelligence." Journal of Complex Networks 8.3 (2020): cnaa019.
tic
clc ;
clear all ;
close all ;
% In SA, there are three decisions for each agent to make:
% Connect or not-Connect ( 1 or 2)
% Defect or Cooperate ( 1 or 2)
% not-Trust or Trsut ( 1 or 2)
% For each decision, there is a moving threshold that divides the interval [0 1] into two parts, P and 1-P, where P is the propensity of the corresponding decision
Trial = 1e5 ; % Simulaton length
Size = 20 ; % Number of the agents

% Constants relating payoff to the change of the propensity/probability of
Chi_Connection = 0.1 ; % SAC mechanism (connect/play decision)
Chi_RL = 0.1 ; %  SAL mechanism (Cooperation or Defection decision)

DemoStart = Trial - 1 ; % Time for the start of the demo
Demoend = Trial ;

% Prisoner's Dilema game's payoffs
s = 0 ; % s is the payoff of the cooperator agent if another agent defected
p = 0 ; % p is the payoff of the agents if both defected
tt = 0.1 ; % (1+ tt) is the payoff of the defector agent if the other agent cooperated (tt < s + 1 and  tt < 1)

mov = VideoWriter('SA.avi') ;
open(mov)

Ratio_CC = zeros(Trial, 1) ; CC = 0 ; % Ration of mutual cooperation

Out1 = zeros(Size, 1) ; % Previous payoff
Out2 = zeros(Size, 1) ; % Current payoff

Connection_index = zeros(Size, 1) ; % Records the index of the connection propensity's matrix
D_C_decision = zeros(Size, 1) ; % D (defection) as 1 and C (cooperation) as 2

D_C_decisionColor = zeros(Trial, Size) ; % Assigns color for the nodes at each trial

P_Connection =  zeros(Size, Size -1) ; % Tendensy of the agents to connect/play
P_RL = zeros(Size, Size -1, 2) ; % P_RL and (1-P_RL) are the propensity of decison D and C, respectively

P_Connection_Used = zeros(Size, 1) ; % Tracks if the Connection decision mechanism has used (2) in the final decision making or not (1)
P_RL_Used = zeros(Size, 1) ; % Tracks if the RL decision mechanism has used (2) in the final decision making or not (1)

Conect2ty = zeros(Trial, Size, Size) ; % Collects the Connection matrix of each trial

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
                P_RL(i , j, k) = 0.5 ;
            end
        end
    end
end

for i = 1 : Size
    for j = 1 : Size - 1
        P_Connection(i, j) = 1/(Size - 1) ;
    end
end


gh =  1 ;
for ti = 2 : Trial

    % Selecting the two agents to play
xxm = 1 ; n = 0 ;
xxn = 1 ; m = 0 ;
     while  xxn  ~= m || xxm ~= n

        % Randomly selecting two agents to connect/play
        m = floor(1 + Size * rand) ;
        n = floor(1 + Size * rand) ;

        while m == n
            n = floor(1 + Size * rand) ;
        end

        P_Connection_m_cumsum = cumsum(P_Connection(m , :) ) ;
        r = rand ;
        for i = 1 : Size - 1
            if  i < m
                if r < P_Connection_m_cumsum(i)
                    xxm = i ;
                    Connection_index(m) = i ;
                    break
                end
            end
            if  i > m
                if r < P_Connection_m_cumsum(i)
                    xxm = i + 1 ;
                    Connection_index(m) = i ;
                    break
                end
            end
        end

        P_Connection_n_cumsum = cumsum(P_Connection(n , :) ) ;
        r = rand ;
        for i = 1 : Size - 1
            if  i < n
                if r < P_Connection_n_cumsum(i)
                    xxn = i ;
                    Connection_index(n) = i ;
                    break
                end
            end
            if  i > n
                if r < P_Connection_n_cumsum(i)
                    xxn = i + 1 ;
                    Connection_index(n) = i ;
                    break
                end
            end
        end

    end

    P_Connection_Used(m) = 2 ; % Connection decision mechanism used
    P_Connection_Used(n) = 2 ;


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

     % update the color of the nodes
    D_C_decisionColor(ti, :) =  D_C_decisionColor(ti-1, :) ;
    D_C_decisionColor(ti, m) =  D_C_decision(m) ;
    D_C_decisionColor(ti, n) =  D_C_decision(n) ;

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
    PConnectionm = zeros(1, Size-1);
    PConnectionn = zeros(1, Size-1);
    for k = 1 : Size-1
        PConnectionm(k) = P_Connection(m, k) ;
        PConnectionn(k) = P_Connection(n, k) ;
    end

    PRLm = zeros(1, 2);
    PRLn = zeros(1, 2);
    for k = 1 : 2
        PRLm(k) = P_RL(m, n, k) ;
        PRLn(k) = P_RL(n, m, k) ;
    end


    P_Connection(m, :) = SA_Update(P_Connection_Used(m), PConnectionm ,  Connection_index(m) , Out2(m) , Out1(m) , Chi_Connection) ;
    P_Connection(n, :) = SA_Update(P_Connection_Used(n), PConnectionn ,  Connection_index(n) , Out2(n) , Out1(n) , Chi_Connection) ;

    P_RL(m, n, :) = SA_Update(P_RL_Used(m), PRLm , D_C_decision(m), Out2(m), Out1(m), Chi_RL) ;
    P_RL(n, m, :) = SA_Update(P_RL_Used(n), PRLn , D_C_decision(n)  , Out2(n) , Out1(n), Chi_RL) ;


    Out1(m) = Out2(m) ;
    Out1(n) = Out2(n) ;

    if  D_C_decision(m) == 2  &&  D_C_decision(n) == 2
        CC = CC + 1 ;
    end
    Ratio_CC(ti) = CC/ ti ;


P_Connection2 = zeros(Size , Size) ;
for i = 1 : Size 
    for j = 1 : Size 
        if j < i
P_Connection2(i, j) = P_Connection(i, j) ;
        end
        if j > i
P_Connection2(i, j) = P_Connection(i, j-1) ;
        end
    end
end

    gh = gh + 1 ;
 Conect2ty(gh, :, :) = P_Connection2 ;
    


end   % end trial

%  Demo
for ty = DemoStart : Demoend

    weights(:, :) = Conect2ty(ty, :, :) ; % The intensity of the lines represents the propensity of the agents to connect/play with one another

    for tn  = 1 : Size
        Summ = cumsum(weights(tn, :)) ;
        weights(tn, :) = weights(tn, :) / Summ(Size) ;
    end

    G = digraph(weights) ;
    LWidths = 1*G.Edges.Weight ;

    plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths)
    intensityValue = LWidths  ;
    OOOO =   intensityValue  ;
    hhh = plot(G,'LineWidth',LWidths) ;
    set(gcf,'color','w') ;
    cccc = 0 ;
    for uuu = 1 : Size
        axis off
        for vvv = 1 : Size
            if uuu ~= vvv
                cccc =  cccc + 1 ;
                if   D_C_decisionColor(ty ,vvv) == 2
                    Col = 'g' ; % Cooperator
                else
                    Col = 'r' ; % Defector
                end
                highlight(hhh, uuu, vvv, 'EdgeColor', [ (1 -OOOO(cccc, 1))  (1-OOOO(cccc, 1))  (1-OOOO(cccc, 1)) ], 'LineWidth', 2, 'MarkerSize', 8)  ;
                highlight(hhh,  vvv, 'NodeColor', Col)  ;
            end
        end
    end
    frame = getframe(gcf) ;
    writeVideo(mov,frame)
end

close(mov)
figure ;
plot(Ratio_CC)
xlabel('Trial (time)'), ylabel('Ratio of mutual cooperation');
hold off


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
                P(j) = 1e-5 ;
            end
            if  P(j) > 1
                P(j) = 1-1e-5 ;
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




