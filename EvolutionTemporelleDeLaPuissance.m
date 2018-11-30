close all
figure
hold on
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
    Omega_p = 2*pi()*c/lambda_p;
    Omega_max = Omega_p/15000;
    Omega_min = 0;
    
%%% Constantes indicés fix
    P_0_fix = 100*10^(-3); 
    Omega_c_fix = sqrt(4*gamma*P_0_fix/abs(beta_2)); % Frequence de coupure pour P_0 fixe
    Omega_fix = 5*10^9 * 2*pi;
    g_fix = abs(beta_2)*Omega_fix*sqrt(Omega_c_fix^2 - Omega_fix^2)/2;
    
    L_0 = -log(delta_0)/(g_fix);
    
    Nb = 1000; % Nombre de points
    z_m = -log(delta_0)/(g_fix); % Fin du tronçon de fibre etudie 
    Z_fix = [0 : z_m/Nb : z_m];    % abscisse
    
    delta_fix = delta_0*exp(g_fix*L_0);
    
%%% Evolution temporelle de la puissance lumineuse en sortie d'une
%%% fibre L_0
        % Ici on fixe toutes les variables précédentes (on les indices fix)
            
            % Plage temporelle 
                Nb_t = 10000;
                t_max = 10^(-9);
                step = t_max/Nb_t;
                t = [step: step: t_max]; % On commence un peu apres zero pour 
                                         % ne pas recupere l indice 1 pour le maximum
            % Calcul P_out(t)
                P_out = (P_0_fix/(2*n*epsilon_0*c*A_eff)).*(1+delta_fix.*cos(Omega_fix.*t)).^2;
            % Tracer de P_out(t)
                plot(t*10^9,P_out)
                set(gca, 'fontsize', 15);
                title('Evolution temporelle de la puissance lumineuse en fin de fibre')
                xlabel('temps (ns)')
                ylabel('P_{out} (W)')
                
                [y,iy] = max(P_out); % Recuperation du maximum et du premier 
                                     % indice ou le trouver 
                
                sprintf('Le taux de repetition du train d impulsion est de % f GHz',10^(-9)/(iy*step))