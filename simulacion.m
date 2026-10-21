clear;
clc;
close all;

%% 1. PARÁMETROS DEL SISTEMA
% Unidades: kg, m, s y N/m

p.m1 = 2.0;       % Masa del soporte
p.m2 = 0.5;       % Masa del péndulo
p.k  = 15.0;      % Constante del resorte
p.l  = 1.0;       % Longitud del péndulo
p.g  = 9.81;      % Gravedad

l0 = 1.5;         % Longitud natural del resorte

%% 2. CONDICIONES INICIALES: MODIFICA ESTOS VALORES

x0 = 0.30;            % Alargamiento inicial del resorte [m]
vx0 = 0.0;            % Velocidad inicial del soporte [m/s]

theta0_grados = 25;    % Ángulo inicial desde la vertical hacia abajo [°]
omega0 = 0.0;         % Velocidad angular inicial [rad/s]

tiempoFinal = 20;      % Duración de la simulación [s]
fps = 40;             % Fotogramas por segundo

% Conversión del ángulo a radianes
theta0 = theta0_grados*pi/180;

% Vector de condiciones iniciales
u0 = [x0; vx0; theta0; omega0];

%% 3. RESOLUCIÓN DE LAS ECUACIONES

if calcularNormal(u0, p) <= 0
    error('Las condiciones iniciales no permiten mantener el contacto.');
end

opciones = odeset( ...
    'RelTol', 1e-9, ...
    'AbsTol', 1e-11, ...
    'MaxStep', 0.02, ...
    'Events', @(t,u) eventoContacto(t,u,p));

sol = ode45( ...
    @(t,u) ecuacionesMovimiento(t,u,p), ...
    [0 tiempoFinal], u0, opciones);

tFin = sol.x(end);
numeroFrames = max(2, ceil(tFin*fps) + 1);

t = linspace(0, tFin, numeroFrames);
U = deval(sol, t);

x = U(1,:);
theta = U(3,:);

perdidaContacto = isfield(sol, 'xe') && ~isempty(sol.xe);

%% 4. POSICIONES DE LAS MASAS

X1 = l0 + x;
Y1 = zeros(size(t));

X2 = X1 + p.l*sin(theta);
Y2 = -p.l*cos(theta);

%% 5. FIGURA DEL MOVIMIENTO

fig = figure( ...
    'Name', 'Movimiento del sistema masa-resorte-péndulo', ...
    'Color', 'w', ...
    'Position', [150 150 1000 650]);

ax = axes('Parent', fig);

hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');

xmin = min([0 X1 X2]) - 0.4;
xmax = max([0 X1 X2]) + 0.4;

xlim(ax, [xmin xmax]);
ylim(ax, [-1.3*p.l 1.3*p.l]);

xlabel(ax, 'Posición horizontal [m]');
ylabel(ax, 'Altura [m]');
title(ax, 'Sistema masa-resorte-péndulo sin rozamiento');

% Apoyo horizontal esquemático
plot(ax, [xmin xmax], [0 0], ...
    'Color', [0.65 0.65 0.65], ...
    'LineWidth', 2);

% Pared del extremo fijo
plot(ax, [0 0], [-0.25 0.25], ...
    'k', 'LineWidth', 5);

% Resorte
hResorte = plot(ax, NaN, NaN, ...
    'Color', [0 0.45 0.74], ...
    'LineWidth', 2);

% Varilla del péndulo
hBarra = plot(ax, NaN, NaN, ...
    'k', 'LineWidth', 2);

% Masa m1
hMasa1 = plot(ax, NaN, NaN, 's', ...
    'MarkerSize', 22, ...
    'MarkerFaceColor', [0 0.45 0.74], ...
    'MarkerEdgeColor', 'k');

% Masa m2
hMasa2 = plot(ax, NaN, NaN, 'o', ...
    'MarkerSize', 16, ...
    'MarkerFaceColor', [0.85 0.15 0.15], ...
    'MarkerEdgeColor', 'k');

% Información instantánea
hTexto = text(ax, 0.03, 0.95, '', ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 12);

legend(ax, [hMasa1 hMasa2], {'m_1', 'm_2'}, ...
    'Location', 'northeast');

%% 6. ANIMACIÓN

s = linspace(0, 1, 250);

for i = 1:numeroFrames

    if ~isgraphics(fig)
        break;
    end

    inicioFrame = tic;

    % Forma del resorte
    xs = X1(i)*s;
    ys = 0.07*sin(2*pi*10*s);

    set(hResorte, 'XData', xs, 'YData', ys);

    % Posiciones de las masas
    set(hMasa1, ...
        'XData', X1(i), ...
        'YData', Y1(i));

    set(hMasa2, ...
        'XData', X2(i), ...
        'YData', Y2(i));

    % Varilla que une ambas masas
    set(hBarra, ...
        'XData', [X1(i) X2(i)], ...
        'YData', [Y1(i) Y2(i)]);

    % Texto
    mensaje = sprintf( ...
        't = %.2f s\nx = %.3f m\nÁngulo = %.1f grados', ...
        t(i), x(i), theta(i)*180/pi);

    if perdidaContacto && i == numeroFrames
        mensaje = sprintf('%s\nPÉRDIDA DE CONTACTO', mensaje);
    end

    set(hTexto, 'String', mensaje);

    drawnow;

    % Reproducción aproximadamente en tiempo real
    if i < numeroFrames
        espera = (t(i+1) - t(i)) - toc(inicioFrame);

        if espera > 0
            pause(espera);
        end
    end
end

%% FUNCIONES AUXILIARES
% Mantener estas funciones al final del archivo.

function du = ecuacionesMovimiento(~, u, p)

    x     = u(1);
    vx    = u(2);
    theta = u(3);
    omega = u(4);

    [ax, alfa] = aceleraciones(x, theta, omega, p);

    du = [vx; ax; omega; alfa];
end


function [ax, alfa] = aceleraciones(x, theta, omega, p)

    A = [
        p.m1 + p.m2, p.m2*p.l*cos(theta);
        cos(theta),  p.l
    ];

    b = [
        p.m2*p.l*sin(theta)*omega^2 - p.k*x;
        -p.g*sin(theta)
    ];

    resultado = A\b;

    ax   = resultado(1);
    alfa = resultado(2);
end


function N = calcularNormal(u, p)

    x     = u(1);
    theta = u(3);
    omega = u(4);

    [~, alfa] = aceleraciones(x, theta, omega, p);

    N = (p.m1 + p.m2)*p.g ...
        + p.m2*p.l*( ...
            sin(theta)*alfa + cos(theta)*omega^2 ...
        );
end


function [valor, detener, direccion] = eventoContacto(~, u, p)

    valor = calcularNormal(u, p);
    detener = 1;
    direccion = -1;
end