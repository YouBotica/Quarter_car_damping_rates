%% Homework 5: 

% 1.) Using the quarter car transfer functions derived in 
% this lecture compare the effect of relative and inertial damping rates.  
% Make a conclusion justified by your modeling results. 

%% Relative damping rates analysis:

Cs_array = [0.1, 0.2, 0.3, 0.4, 0.6]*10; 
Ct_array = [0.1, 0.2, 0.3, 0.4, 0.6]*10;

Kt = 1000; % lbs/in 
Ks = 100; % lbs/in
g = 32.174; % ft/sec^2
Ws = 1000; % lbs
Wu = 100; % lbs
ms = Ws / g;  % Sprung mass (slugs)
mu = Wu / g; % Unsprung mass (slugs)


% Generate transmissibility plots varying Cs (sprung mass damper rate):

for i=1:length(Cs_array)
    Ct = mean(Ct_array);
    Cs = Cs_array(i);
    fig_num = 1;
    twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu, fig_num);
end

for i=1:length(Ct_array)
    Ct = Ct_array(i);
    Cs = mean(Cs_array);
    fig_num = 2;
    twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu, fig_num);
end

%% Inertial damping rates analysis:

Cs_array = [0.1, 0.2, 0.3, 0.4, 0.6]*10; 
Ct_array = [0.1, 0.2, 0.3, 0.4, 0.6]*10;

Kt = 1000; % lbs/in 
Ks = 100; % lbs/in
g = 32.174; % ft/sec^2
Ws = 1000; % lbs
Wu = 100; % lbs
ms = Ws / g;  % Sprung mass (slugs)
mu = Wu / g; % Unsprung mass (slugs)

% Legend arrays:
legendEntries = cell(1, length(Cs_array) + length(Ct_array)); % Initialize the cell array
legendCounter = 1; % Initialize a counter for legend entries

figure(1); % Create the first figure outside the loop
hold on; % Hold on to plot multiple lines

for i=1:length(Cs_array)
    Ct = mean(Ct_array);
    Cs = Cs_array(i);
    % plot:
    [mgs, mgu, w] = twomass_inertial_damp(Ks, Kt, Cs, Ct, ms, mu);
    semilogx(w, mgs); % Plot the magnitude vs frequency for the sprung mass
    semilogx(w, mgu); % Plot the magnitude vs frequency for the unsprung mass
    % legend:
    legendEntries{legendCounter} = sprintf('Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounter = legendCounter + 1;
end

legend(legendEntries{1:length(Cs_array)}); % Create the legend for the first figure
title('Transmissibility plot varying Cs');
xlabel('frequency [rad/sec]');
ylabel('amplitude ratio');
hold off; % Release the hold after plotting


figure(2); % Create the second figure outside the loop
hold on; % Hold on to plot multiple lines


for i=1:length(Ct_array)
    Ct = Ct_array(i);
    Cs = mean(Cs_array);
    % plot:
    [mgs, mgu, w] = twomass_inertial_damp(Ks, Kt, Cs, Ct, ms, mu);
    semilogx(w, mgs); % Plot the magnitude vs frequency for the sprung mass
    semilogx(w, mgu); % Plot the magnitude vs frequency for the unsprung mass
    % legend:
    legendEntries{legendCounter} = sprintf('Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounter = legendCounter + 1;
end

legend(legendEntries{length(Cs_array)+1:end}); %legend(legendEntries{length(Cs_array)+1:end}); % Create the legend for the second figure
title('Transmissibility plot varying Ct');
xlabel('frequency [rad/sec]');
ylabel('amplitude ratio');
hold off; % Release the hold after plotting



function [mgs, mgu, w] = twomass_inertial_damp(Ks, Kt, Cs, Ct, ms, mu)

    % Mathematical model:
    
    A=[ 0,  1,  0,  0;
        -Ks/ms, -Cs/ms,  Ks/ms,  0;
         0,   0,   0,   1;
         Ks/mu,  Cs/mu,  -(Ks+Kt)/mu,  0];
     
    B=[0, 0;
        0, 0;
        0, 0;
        Kt/mu, Ct/mu];
     
    C=[1, 0, 0, 0;
        0, 0, 1, 0];
     
    D=[0, 0;0, 0];

    %twomass calculates the frequency response of a two-mass
    [mag, phase, w]=bode(ss(A,B(:,1),C(1,:),[0]),logspace(0,3)); % Outputs Xs
     mgs(1:50)=mag;

     
    [mag, phase, w]=bode(ss(A,B(:,1),C(2,:),[0]),logspace(0,3)); % Outputs Xu
    mgu(1:50)=mag;
end


% function [mgs, mgu, w] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu, fig_num)
% 
%     % Mathematical model:
% 
%     A=[ 0,  1,  0,  0;
%         -Ks/ms, -Cs/ms,  Ks/ms,  Cs/ms;
%          0,   0,   0,   1'
%          Ks/mu,  Cs/mu,  -(Ks+Kt)/mu,  -(Cs+Ct)/mu];
% 
%      B=[0, 0;
%         0, 0;
%         0, 0;
%         Kt/mu, Ct/mu];
% 
%     C=[1, 0, 0, 0;
%         0, 0, 1, 0];
% 
%     D=[0, 0;0, 0];
% 
%     %twomass calculates the frequency response of a two-mass
%     [mag,phase,w]=bode(ss(A,B(:,1),C(1,:),[0]),logspace(0,3)); % Outputs Xs
%      mgs(1:50)=mag;
% 
%     [mag,phase,w]=bode(ss(A,B(:,1),C(2,:),[0]),logspace(0,3)); % Outputs Xu
%     mgu(1:50)=mag;
% 
%     % semilogx(w,mgs,w,mgu),figure(fig_num)
%     %     ['Cs = ' num2str(Cs) ' lb-s/in']
% 
%     %semilogx(w, mgs, w, mgu); %, 'DisplayName', ['Cs = ' num2str(Cs)]);
%     xlabel('frequency [rad/sec]');
%     ylabel('amplitude ratio');
%     title(sprintf('Transmissibility plot with Ct = %2d and Cs = %2d', Ct, Cs))
%     legend(sprintf('Ct = %2d and Cs = %2d', Ct, Cs));
%     % legend('Ct = ' num2str(Ct) ' and Cs = ' num2str(Cs) '.');
% 
% end

