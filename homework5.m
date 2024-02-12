%% Cs damping rates analysis: 

Cs_array = logspace(log10(0.1), log10(1000), 10);
Ct_array = 0.1; %logspace(log10(0.1), log10(100), 7);


Kt = 1000*12; % lbs/ft A conversion factor of 12 is applied here to make units compatible with feet and slugs
Ks = 100*12; % lbs/ft A conversion factor of 12 is applied here to make units compatible with feet and slugs
g = 32.174; % ft/sec^2
Ws = 1000; % lbs
Wu = 100; % lbs
ms = Ws / g;  % Sprung mass (slugs)
mu = Wu / g; % Unsprung mass (slugs)

% Legend arrays:
legendTopEntries = cell(1, 2*length(Cs_array)); % Initialize the cell array
legendBottomEntries = cell(1, 2*length(Cs_array)); % Initialize the cell array
legendCounter = 1; % Initialize a counter for legend entries

figure; % Create the first figure outside the loop

% subplots:
subplot(2,1,1);
xscale log;
hold on;

subplot(2,1,2);
xscale log;
hold on;


for i=1:length(Cs_array)
    Ct = mean(Ct_array);
    Cs = Cs_array(i);

    % plot:
    subplot(2,1,1);
    [mgs, mgu, w] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu);

    semilogx(w', mgs, 'o-'); % Plot the magnitude vs frequency for the sprung mass
    legendTopEntries{legendCounter} = sprintf('Sprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    % legendCounter = legendCounter + 1;

    semilogx(w, mgu, '--'); % Plot the magnitude vs frequency for the unsprung mass
    legendBottomEntries{legendCounter} = sprintf('Unsprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounter = legendCounter + 1;

    
    % plot:
    subplot(2,1,2);
    [mgs, mgu, w] = twomass_inertial_damp(Ks, Kt, Cs, Ct, ms, mu);

    semilogx(w, mgs, 'o-'); % Plot the magnitude vs frequency for the sprung mass
    legendTopEntries{legendCounter} = sprintf('Sprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    % legendCounter = legendCounter + 1;

    semilogx(w, mgu, '--'); % Plot the magnitude vs frequency for the unsprung mass
    legendBottomEntries{legendCounter} = sprintf('Sprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounter = legendCounter + 1;

end

subplot(2,1,1);
title('Relative damping transmissibility plot varying Cs');
legend(legendTopEntries{1:2*length(Cs_array)}); % Create the legend for the first figure
xlabel('frequency [rad/sec]');
ylabel('amplitude ratio');
hold off;

subplot(2,1,2);
title('Inertial damping transmissibility plot varying Cs');
legend(legendBottomEntries{1:2*length(Cs_array)}); % Create the legend for the first figure
xlabel('frequency [rad/sec]');
ylabel('amplitude ratio');
hold off;


%% Ct damping rates analysis: 

Cs_array = logspace(log10(0.1), log10(100), 7);
Ct_array = logspace(log10(0.1), log10(100), 7);

Kt = 1000*12; % lbs/ft 
Ks = 100*12; % lbs/ft
g = 32.174; % ft/sec^2
Ws = 1000; % lbs
Wu = 100; % lbs
ms = Ws / g;  % Sprung mass (slugs)
mu = Wu / g; % Unsprung mass (slugs)


% Legend arrays:
legendTopEntries = cell(1, 2*length(Ct_array)); % Initialize the cell array
legendBottomEntries = cell(1, 2*length(Ct_array)); % Initialize the cell array
legendCounter = 1; % Initialize a counter for legend entries

figure; % Create the first figure outside the loop

% subplots:
subplot(2,1,1);
xscale log;
hold on;

subplot(2,1,2);
xscale log;
hold on;

for i=1:length(Ct_array)
    Ct = Ct_array(i);
    Cs = mean(Cs_array);

    % plot:
    subplot(2,1,1);
    [mgs, mgu, w] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu);

    semilogx(w', mgs, 'o-'); % Plot the magnitude vs frequency for the sprung mass
    legendTopEntries{legendCounter} = sprintf('Sprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    % legendCounter = legendCounter + 1;

    semilogx(w, mgu, '--'); % Plot the magnitude vs frequency for the unsprung mass
    legendBottomEntries{legendCounter} = sprintf('Unsprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounter = legendCounter + 1;

    
    % plot:
    subplot(2,1,2);
    [mgs, mgu, w] = twomass_inertial_damp(Ks, Kt, Cs, Ct, ms, mu);

    semilogx(w, mgs, 'o-'); % Plot the magnitude vs frequency for the sprung mass
    legendTopEntries{legendCounter} = sprintf('Sprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    % legendCounter = legendCounter + 1;

    semilogx(w, mgu, '--'); % Plot the magnitude vs frequency for the unsprung mass
    legendBottomEntries{legendCounter} = sprintf('Sprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounter = legendCounter + 1;
end


subplot(2,1,1);
title('Relative damping transmissibility plot varying Ct');
legend(legendTopEntries{1:2*length(Cs_array)}); % Create the legend for the first figure
xlabel('frequency [rad/sec]');
ylabel('amplitude ratio');
hold off;

subplot(2,1,2);
title('Inertial damping transmissibility plot varying Ct');
legend(legendBottomEntries{1:2*length(Cs_array)}); % Create the legend for the first figure
xlabel('frequency [rad/sec]');
ylabel('amplitude ratio');
hold off;




%% Functions:

function [mgs, mgu, w] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu)

    % Mathematical model:
    
    A=[ 0,  1,  0,  0;
        -Ks/ms, -Cs/ms,  Ks/ms,  Cs/ms;
         0,   0,   0,   1;
         Ks/mu,  Cs/mu,  -(Ks+Kt)/mu,  -(Cs+Ct)/mu];
     
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


function [mgs, mgu, w] = twomass_inertial_damp(Ks, Kt, Cs, Ct, ms, mu)

    % Mathematical model:
    
    A=[ 0,  1,  0,  0;
        -Ks/ms, -Cs/ms,  Ks/ms,  0;
         0,   0,   0,   1;
         Ks/mu,  0,  -(Ks+Kt)/mu,  -Ct/mu]; %Ks/mu,  Cs/mu,  -(Ks+Kt)/mu,  0];
     
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