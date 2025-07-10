%% --- Carga datos ---
filename = 'Fort Point.xlsx';            % Nombre del archivo Excel
sheetname = 'Fort Point';                % Nombre de la hoja dentro del archivo
datos = readtable(filename, 'Sheet', sheetname); % Carga la tabla con los datos

% agregado pcm
%graficaremos las magnitudes y sus alturas maximas, como mapa, para ver las
%dependencias
mags=unique(datos.Mw); %guarda magnitudes, pero tiene NaN
idx=~isnan(mags); %pregunta si NO ES NaN
mags=mags(idx); %salva solo los casos no NAN

figure(22)

for i=1:length(mags)
    idx=datos.Mw==mags(i);
    scatter(ones(size(datos.Lat_S(idx)))*mags(i)+datos.Dist_TC(idx)/(2*max(datos.Dist_TC)),(datos.Lat_S(idx)+datos.Lat_N(idx))*.5,datos.Slip_max(idx)*10,datos.H_max(idx),'filled',marker='s',MarkerEdgeColor='k')
    hold on;

end
hold off
colorbar
title('Tamano del marcador representan slip*10, Colores del marcador la altura maxima del tsunami, ')


%% --- Filtrado por magnitudes ---
datos_Mw8 = datos(datos.Mw == 8.0, :);   % Filtra eventos con magnitud Mw 8.0
datos_Mw9 = datos(datos.Mw == 9.0, :);   % Filtra eventos con magnitud Mw 9.0

% --- Calcular momentos sísmicos y factor de escala entre Mw9 y Mw8 ---
M0_8 = 10^((3/2)*(8 + 6.06));            % Momento sísmico para Mw 8.0, si esta dado en Nm
M0_9 = 10^((3/2)*(9 + 6.06));            % Momento sísmico para Mw 9.0
Factor = round(M0_9/M0_8);               % Factor aproximado de aumento de energía

% --- Casos Mw8: 32 unitarios de alto y bajo slip ---
Mw8_unitarios_AltoSlip = datos_Mw8(1:32, :);       % Primeros 32 eventos de alto slip
Mw8_unitarios_BajoSlip = datos_Mw8(36:67, :);      % Eventos 36 al 67 como bajo slip

% --- Construcción de grupo intercalado Mw8 (14 unitarios y 12 con factor 1.5) ---
n_unidades = 14;                         % Número de eventos unitarios (factor 1)
n_1_5 = 12;                              % Número de eventos con factor 1.5
total_events = n_unidades + n_1_5;       % Total de eventos = 26

% Selección de eventos para alto slip
indices_alto = 1:total_events;                                 
Mw8_intercalado_AltoSlip = datos_Mw8(indices_alto, :);         % Primeros 26 eventos

% Selección de eventos para bajo slip
indices_bajo = 36:(36+total_events-1);                         
Mw8_intercalado_BajoSlip = datos_Mw8(indices_bajo, :);         % Eventos 36 a 61

% Calcular factores intercalados y aplicar a H_max para alto slip
factors_AltoSlip = intercalar_indices(n_unidades, n_1_5);      % Vector de factores [1 1.5 1 ...]
H_max_intercalado_AltoSlip = Mw8_intercalado_AltoSlip.H_max .* factors_AltoSlip;

% Igual para bajo slip
factors_BajoSlip = intercalar_indices(n_unidades, n_1_5);
H_max_intercalado_BajoSlip = Mw8_intercalado_BajoSlip.H_max .* factors_BajoSlip;

% --- Calcular estadísticas para cada grupo ---

% Estadísticas de eventos unitarios
stats_Mw8_unit_AltoSlip = calc_stats(Mw8_unitarios_AltoSlip.H_max);
stats_Mw8_unit_BajoSlip = calc_stats(Mw8_unitarios_BajoSlip.H_max);

% Estadísticas de eventos intercalados
stats_Mw8_int_AltoSlip = calc_stats(H_max_intercalado_AltoSlip);
stats_Mw8_int_BajoSlip = calc_stats(H_max_intercalado_BajoSlip);

% Estadísticas para eventos Mw9 (alto y bajo slip)
Mw9_AltoSlip = datos_Mw9(1:32,:);                % Primeros 32 como alto slip
Mw9_BajoSlip = datos_Mw9(33:60,:);               % Eventos 33 al 60 como bajo slip
stats_Mw9_AltoSlip = calc_stats(Mw9_AltoSlip.H_max);
stats_Mw9_BajoSlip = calc_stats(Mw9_BajoSlip.H_max);

% --- Crear tabla resumen con resultados clave ---
Casos = { ...
    '32 unitarios Grupo Alto Slip';
    '32 unitarios Grupo Bajo Slip';
    '14x1 + 12x1.5 + 1x1 Grupo Alto Slip';
    '14x1 + 12x1.5 + 1x1 Grupo Bajo Slip';
    '32 Magnitud 9 Grupo Alto Slip';
    '28 Magnitud 9 Grupo Bajo Slip'};

% Número de terremotos por caso
N_Terremotos = [ ...
    stats_Mw8_unit_AltoSlip.N;
    stats_Mw8_unit_BajoSlip.N;
    stats_Mw8_int_AltoSlip.N;
    stats_Mw8_int_BajoSlip.N;
    stats_Mw9_AltoSlip.N;
    stats_Mw9_BajoSlip.N];

% Suma total de alturas máximas (H_max)
Total_Hmax = [ ...
    stats_Mw8_unit_AltoSlip.Sum;
    stats_Mw8_unit_BajoSlip.Sum;
    stats_Mw8_int_AltoSlip.Sum;
    stats_Mw8_int_BajoSlip.Sum;
    NaN;   % No se calcula suma para Mw9
    NaN];

% Promedio de alturas máximas
Promedio_Hmax = [ ...
    stats_Mw8_unit_AltoSlip.Mean;
    stats_Mw8_unit_BajoSlip.Mean;
    stats_Mw8_int_AltoSlip.Mean;
    stats_Mw8_int_BajoSlip.Mean;
    stats_Mw9_AltoSlip.Mean;
    stats_Mw9_BajoSlip.Mean];

% Máximos de H_max
Max_Hmax = [ ...
    stats_Mw8_unit_AltoSlip.Max;
    stats_Mw8_unit_BajoSlip.Max;
    stats_Mw8_int_AltoSlip.Max;
    stats_Mw8_int_BajoSlip.Max;
    stats_Mw9_AltoSlip.Max;
    stats_Mw9_BajoSlip.Max];

% Mínimos de H_max
Min_Hmax = [ ...
    stats_Mw8_unit_AltoSlip.Min;
    stats_Mw8_unit_BajoSlip.Min;
    stats_Mw8_int_AltoSlip.Min;
    stats_Mw8_int_BajoSlip.Min;
    stats_Mw9_AltoSlip.Min;
    stats_Mw9_BajoSlip.Min];

% Construcción de la tabla con todos los datos
tabla_resumen = table(Casos, N_Terremotos, Total_Hmax, Promedio_Hmax, Max_Hmax, Min_Hmax);
disp(tabla_resumen)      % Muestra la tabla en la consola

% --- Gráfico comparativo de alturas máximas ---

% Vector con valores a graficar
valores = [tabla_resumen.Total_Hmax(1);
           tabla_resumen.Total_Hmax(2);
           tabla_resumen.Total_Hmax(3);
           tabla_resumen.Total_Hmax(4);
           tabla_resumen.Promedio_Hmax(5);
           tabla_resumen.Promedio_Hmax(6);
           tabla_resumen.Max_Hmax(5);
           tabla_resumen.Max_Hmax(6)];

% Errores para barras de promedio Mw9 (Max - Promedio)
errores = [ ...
    tabla_resumen.Max_Hmax(5)-tabla_resumen.Promedio_Hmax(5);
    tabla_resumen.Max_Hmax(6)-tabla_resumen.Promedio_Hmax(6)];

% Etiquetas de cada barra
etiquetas = { ...
    'Mw8 Alto Slip 1', 'Mw8 Bajo Slip 1', ...
    'Mw8 Alto Slip 1.5', 'Mw8 Bajo Slip 1.5', ...
    'Mw9 Alto Slip prom', 'Mw9 Bajo Slip prom', ...
    'Mw9 Alto Slip max', 'Mw9 Bajo Slip max'};

% Crear gráfico de barras
figure
hold on

% Dibujar barras
b = bar(valores, 'FaceColor', [0.2 0.6 0.8]);  % Color base azul
b.FaceColor = 'flat';
b.CData(1:4,:) = repmat([0.2 0.6 0.8],4,1);    % Azul para Mw8
b.CData(5:6,:) = repmat([0.1 0.8 0.4],2,1);    % Verde para promedios Mw9
b.CData(7:8,:) = repmat([0.9 0.5 0.2],2,1);    % Naranjo para máximos Mw9

% Añadir barras de error (líneas verticales en promedio Mw9)
errorbar(5, valores(5), errores(1), 'k', 'LineWidth', 1.5)
errorbar(6, valores(6), errores(2), 'k', 'LineWidth', 1.5)

% Configurar ejes y etiquetas
xticks(1:8)
xticklabels(etiquetas)
xtickangle(45)
ylabel('Altura máxima (m)')
title('Comparación de alturas Mw8 vs Mw9')
grid on
box on
hold off

%% --- Funciones locales abajo del todo ---

% Función que genera una secuencia de factores intercalados 1 y 1.5
function idx = intercalar_indices(n_unidades, n_1_5)
    total = n_unidades + n_1_5;
    idx = zeros(total,1);
    i1 = 1; i15 = 1;
    for k=1:total
        if mod(k,2)==1 && i1 <= n_unidades
            idx(k) = 1;      % Evento unitario
            i1 = i1+1;
        elseif i15 <= n_1_5
            idx(k) = 1.5;    % Evento con factor 1.5
            i15 = i15+1;
        elseif i1 <= n_unidades
            idx(k) = 1;
            i1 = i1+1;
        end
    end
end

% Función que calcula estadísticas básicas de un vector H
function stats = calc_stats(H)
    stats.N = length(H);   % Número de eventos
    stats.Sum = sum(H);    % Suma de H_max
    stats.Mean = mean(H);  % Promedio de H_max
    stats.Max = max(H);    % Máximo H_max
    stats.Min = min(H);    % Mínimo H_max
end

%%hola