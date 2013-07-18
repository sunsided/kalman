%% Beispiel für das Kalman-Filter im vierdimensionalen Fall
% In diesem Beispiel sind die Eingangswerte des Filters bekannt und
% werden nicht als Systemzustand modelliert. Kalman-Filter und 
% Systemmodell verwenden identische Matrizen A und B.

clc; home; close all;

% Zeitvektor
dt = 0.01;           % Abtastzeit in [s]
Tend = 1.5;          % Endzeit in [s]
t = 0:dt:Tend-dt;    % Zeitvektor [s]
N = length(t);
  
%% Durchführung der Simulation

% Systemparameter
ax = -1;            % Beschleunigung in x-Richtung [m/s^2]
vx0 = 2.1;          % Startgeschwindigkeit in x-Richtung [m/s]

ay = -9.81;         % Beschleunigung in y-Richtung [m/s^2]
vy0 = 7.28;         % Startgeschwindigkeit in y-Richtung [m/s]

px0 = 0;            % Startposition in x-Richtung [m]
py0 = 0;            % Startposition in y-Richtung [m]

% Initialer Zustandsvektor
x0 = [px0;           % Initiale x-Position
     py0;           % Initiale y-Position
     vx0;           % Initiale x-Geschwindigkeit
     vy0];          % Initiale y-Geschwindigkeit
 
% Zeitkontinuierliche Systemübergangsmatrix
A = [0, 0, dt,  0;    % Geschwindigkeit in x-Richtung
     0, 0,  0, dt;    % Geschwindigkeit in y-Richtung
     0, 0,  0,  0;    % Beschleunigung in x-Richtung
     0, 0,  0,  0];   % Beschleunigung in y-Richtung
    
% Zeitkontinuierliche Eingangsmatrix
B = [0.5*dt^2,        0;    % Einfluss von Eingangswert auf x-Geschwindigkeit
            0, 0.5*dt^2;    % Einfluss von Eingangswert auf y-Geschwindigkeit
           dt,        0;    % Einfluss von Eingangswert auf x-Beschleunigung
            0,      dt];    % Einfluss von Eingangswert auf y-Beschleunigung

% Ausgangsmatrix
C = [1, 0, 0, 0;    % Ausgang ist Identität des Zustandes
     0, 1, 0, 0;
     0, 0, 1, 0;
     0, 0, 0, 1];
 
% Durchgriffsmatrix
D = [0, 0;          % Einfluss von Eingangswert auf x-Position
     0, 0;          % Einfluss von Eingangswert auf y-Position
     0, 0;          % Einfluss von Eingangswert auf x-Geschwindigkeit
     0, 0];         % Einfluss von Eingangswert auf y-Geschwindigkeit
 
% Eingangsvektor (konstant)
u = [ax;            % Beschleunigung in x-Richtung
     ay];           % Beschleunigung in y-Richtung

% Ergebnisvektor
x_results = zeros(length(x0), length(t));

% Zustandsvektor initialisieren
x = x0;

for i=1:N
    % Simulationszeit ermitteln
    time = t(i);
    
    % Folgezustand ermitteln (Differential)
    dx = A*x + B*u;
    
    % Ausgang berechnen
    y = C*x + D*u;
    
    % "Integration" des Zustandes
    x = x + dx;

    % Vektor sichern
    x_results(:,i) = y;
end

%% Diskretisierung des zeitkontinuierlichen Systems

% Zeitvektor
ddt = 0.005;           % Abtastzeit in [s]
dt = 0:ddt:Tend-ddt;    % Zeitvektor [s]
dN = length(dt);

% Diskretisierung der Zustandsübergangsmatrix
Ad = [1, 0, ddt,   0;     % x = x0 + vx*dt + ...
      0, 1,   0, ddt;     % y = y0 + vy*dt + ...
      0, 0,   1,   0;     % vx = vx0
      0, 0,   0,   1];    % vy = vy0

% Diskretisierung der Eingangsmatrix
Bd = [0.5*ddt^2,         0;     % x = ... + 0.5*ax*dt^2
              0, 0.5*ddt^2;     % y = ... + 0.5*ay*dt^2
            ddt,         0;     % vx = ... ax*dt
              0,       ddt];    % vy = ... ay*dt
  
% Ergebnisvektor
dx_results = zeros(length(x0), length(dt));

% Zustandsvektor initialisieren
x = x0;

for i=1:dN
    % Simulationszeit ermitteln
    time = dt(i);
    
    % Folgezustand ermitteln (Differential)
    x_next = Ad*x + Bd*u;
    
    % Ausgang berechnen
    y = C*x + D*u;
    
    % "Integration" des Zustandes
    x = x_next;

    % Vektor sichern
    dx_results(:,i) = y;
end         
         
%% Approximation durch das Kalman-Filter

% Es wird definiert, dass nur die Position in x- und y-Richtung gemessen
% werden kann, nicht jedoch Geschwindigkeit oder Beschleunigung.

% Messrauschen erzeugen
vmess = 0.1;                       % Varianz in x/y-Richtung [m]
v = sqrt(vmess)*randn(2, dN);       % Messrauschen [m]

% Virtueller Messwertvektor
messungen = dx_results(1:2, :) + v;

% Messübergangsmatrix
H = [1, 0, 0, 0;            % Übergang von Zustandsvektor in x-Messung
     0, 1, 0, 0];           % Übergang von Zustandsvektor in y-Messung
 
% Messunsicherheit / Messkovarianzmatrix
R = [vmess, 0;
         0, vmess];

% Unsicherheit der Eingänge (Kovarianz des Prozessrauschens)
Q = [10E-5, 0;
     0, 10E-5];
     
% Konvergenzfaktor
lambda = 0.999;

% Initialisierung der Zustandsschätzung
x = [0;                     % Position in x-Richtung [m]
     0;                     % Position in y-Richtung [m]
     1;                     % Geschwindigkeit in x-Richtung [m/s]
     1];                    % Geschwindigkeit in y-Richtung [m/2]
 
% Initialisierung der Systemkovarianzen
P = [10, 0, 0, 0;
     0, 10, 0, 0;
     0, 0, 10, 0;
     0, 0, 0, 10];

% Identitätsmatrix
I = eye(size(P));
  
% Ergebnisvektor
kx_results = zeros(length(x), length(dt));
kx_results(:, 1) = x;
 
% Schleife beginnt mit einer a-priori-Schätzung,
% da der Initialzustand bereits eine a-posteriori-Schätzung ist.
x_pos = x;
P_pos = P;
clearvars x P A B D N;

% Iterieren
for i=1:dN        
    % Prädiktion des Systemzustandes
    x_pri = Ad*x_pos + Bd*u;
    
    % Prädiktion der Systemkovarianz
    P_pri = Ad*P_pos*Ad' * 1/(lambda^2) + Bd*Q*Bd';

    % Messung beziehen
    z = messungen(:, i);

    % Ermittlung von Innovation und Residualkovarianz
    w = z - H*x_pri;
    S = H*P_pri*H' + R;

    % Ermittlung des Kalman-Gains
    K = P_pri * H' / S;

    % Korrektur der Schätzung mittels Messwert
    x_pos = x_pri + K*w;

    % Korrektur der Kovarianzmatrix mittels Messwert
    % "Joseph-Form", insbesondere, wenn Kalman-Gain nicht optimal
    IKH = I-K*H;
    P_pos = IKH*P_pri*IKH' + K*R*K';

    % Vektor sichern
    kx_results(:,i) = x_pos;
end
     
%% Darstellung der Ergebnisse der Simulation
figure('Name', 'Kalman-Simulation: Ballistische Kurve', 'NumberTitle', 'Off');

% Plot der Position
subplot(2,2,1:2);
plot(x_results(1,:), x_results(2, :), 'k', 'LineWidth', 2);
title('Position');
xlabel('x [m]');
ylabel('h [m]');
grid on;
hold on;
stairs(dx_results(1,:), dx_results(2, :), 'b', 'LineWidth', 1);
plot(messungen(1,:), messungen(2, :), 'g+', 'LineWidth', 1);
plot(kx_results(1,:), kx_results(2, :), 'r-', 'LineWidth', 1);
legend('pos_{kont.}', 'pos_{disk.}', 'Messwerte', 'Kalman', 'Location', 'NorthWest');

% Plot der x-Geschwindigkeit
subplot(2,2,3);
plot(t, x_results(3, :), 'k', 'LineWidth', 2);
title('Geschwindigkeit');
xlabel('t [s]');
ylabel('v_x [m/s]');
grid on;
hold on;
plot(dt, kx_results(3, :), 'r-', 'LineWidth', 1);
legend('v_{x}', 'Kalman', 'Location', 'Northeast');

% Plot der y-Geschwindigkeit
subplot(2,2,4);
plot(t, x_results(4, :), 'k', 'LineWidth', 2);
title('Geschwindigkeit');
xlabel('t [s]');
ylabel('v_y [m/s]');
grid on;
hold on;
plot(dt, kx_results(4, :), 'r-', 'LineWidth', 1);
legend('v_{y}', 'Kalman', 'Location', 'Northeast');