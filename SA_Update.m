
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
                P(j) = 0 ; % 1e-6
            end
            if  P(j) > 1
                P(j) = 1 ; % 1-1e6
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



