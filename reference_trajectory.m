function ref_angle = reference_trajectory(x_out,sigma)
    % This function calculates reference joint angles from CPG outputs.
    %
    %
    % input: 
    %
    % - x_out: Output of the CPG network
    %
    % - sigma: Threshold for swing-stance determination
    %
    % output: 
    %
    % - ref_angle: Reference angles for 18 joints

    W = 0.2; % Distance from AEP to PEP
    H = 0.1; % Maximum height of the foot tip from the ground

    swing_stance = x_out > sigma; % 1: Swing, 0: Stance
    l = length(swing_stance);
    
    i_TD = -1;
    i_LO = -1;
    if swing_stance(1) == 1 % Swing
        i_TD = find(swing_stance == 0,1,"first"); % Index of t_TD
        l_swing = i_TD;
        ref_pos = zeros(3,l_swing);
    elseif swing_stance(1) == 0 % Stance
        i_LO = find(swing_stance == 1,1,"first");
        l_stance = i_LO;
        ref_pos = zeros(3,l_stance);
    end
    
    for i = 1:l
        if i == i_LO % Liftoff
            i_TD = find(swing_stance(i_LO:end) == 0,1,"first") + i_LO - 1; % Index of t_TD
            if i_TD
                l_swing = i_TD - i_LO + 1; % Swing duration
            else
                break
            end
            ref_pos_tmp = [zeros(1,l_swing);linspace(-W/2,W/2,l_swing);H*(x_out(i_LO:i_TD) - sigma(i_LO:i_TD))/(max(x_out(i_LO:i_TD) - sigma(i_LO:i_TD)))];
            ref_pos = [ref_pos(:,1:end-1),ref_pos_tmp];
        elseif i == i_TD % Touchdown
            i_LO = find(swing_stance(i_TD:end) == 1,1,"first") + i_TD - 1; % Index of t_LO
            if i_LO
                l_stance = i_LO - i_TD + 1; % Stance duration
            else
                break
            end
            ref_pos_stance_tmp = [zeros(1,l_stance);linspace(W/2,-W/2,l_stance);zeros(1,l_stance)];
            ref_pos = [ref_pos(:,1:end-1),ref_pos_stance_tmp];
        end
    end

    ref_angle = inverse_kinematics(ref_pos);
end

function angle = inverse_kinematics(pos)
    Px = pos(1,:);
    Py = pos(2,:);
    Pz = pos(3,:);

    L1 = 0.14;
    L2 = 0.24;
    L3 = 0.42;
    Px0 = L1 + L2/sqrt(2);
    Pz0 = L3 - L2/sqrt(2);
    L = sqrt((Pz0 - Pz).^2 + (sqrt(Px0^2 + Py.^2) - L1).^2);

    eta1 = atan(Py ./ (Px0 + Px));
    eta2 = acos((L2^2 + L.^2 - L3^2) ./ (2*L2*L)) + atan((sqrt(Px0^2 + Py.^2) - L1)./(Pz0 - Pz)) - pi/2;
    eta3 = acos((L2^2 + L3^2 - L.^2) ./ (2*L2*L3)) - pi;
    
    angle = [eta1;eta2;eta3];
end
