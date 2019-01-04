close all
figure
hold on
str = {};
%%% Constantes
    
    c = 3*10^8;                   % Celerite de la lumiere (m*s^(-1))
    epsilon_0 = 8.85 * 10^(-12);  % Permittivite du vide (F*m^(-1))
    beta_2 = -2 * 10^(-26);       % Dispersion de vitesse de groupe (s^2*m^(-1))
    gamma = 2*10^(-3);            % Coefficient non lineaire (W^(-1)*m^(-1))
    d_eff = 10*10^(-6);           % Diametre effectif du mode guidé (m)
    A_eff = pi*(d_eff/2)^2;       % Aire effective(m^2)
    n = 1.47;                     % Indice de refraction
    lambda_p = 532 * 10^(-9);     % Longueur d'onde du faisceau (m)
    delta_0 = 0.02;               % Profondeur de modulation initiale
	
    Omega_p = 2*pi()*c/lambda_p;  % Frequence du faisceau (Hz)
    Omega_max = Omega_p/15000;    % Frequence de normalisation (Hz)
    Omega_min = 0;
    
%%% Evolution gain lineique (g)
    
        
    % Differentes valeurs de P_0
        P_0 = [10^(-2), 5*10^(-2), 10^(-1)];      % Puissance moyenne (W)
        
        % plage de frequence de modulation (abscisse)
            N = 1000; % Nombre de points par courbes
            Omega = [Omega_min : (Omega_max - Omega_min)/N : Omega_max];
        
        for p = P_0
            Omega_c = sqrt(4*gamma*p/abs(beta_2)); % Frequence de coupure
            % Calcul du gain linéique
                g = zeros(N,1); 
                k = 1;
                for omega = Omega
                    if Omega_c^2 - omega^2 < 0
                       g(k)=0; % Sinon negatif sous la racine
                    
                    else
                        g(k) = abs(beta_2)*omega*sqrt(Omega_c^2 - omega^2)/2;
                    end
                    k= k+1;
                    
                end
                % Ajout de la courbe 
                    plot(Omega/Omega_max, g)
                    str = [str, strcat('P_0 = ', num2str(p*10^3),' mW')]; % incrementation de la legende
        end
        title('Evolution gain lineique')
        set(gca, 'fontsize', 15);
        xlabel('Omega Normalise')
        ylabel('g (m^{-1})')
        legend(str{:})