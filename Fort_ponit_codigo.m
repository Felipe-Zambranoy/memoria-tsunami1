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


% --- Filtrado por magnitudes ---
datos_Mw8 = datos(datos.Mw == 8.0, :);   % Filtra eventos con magnitud Mw 8.0
datos_Mw9 = datos(datos.Mw == 9.0, :);   % Filtra eventos con magnitud Mw 9.0
datos_Mw9_5 = datos(datos.Mw == 9.5, :);   % Filtra eventos con magnitud Mw 9.5


% --- Casos Mw8: 32 unitarios de alto y bajo slip ---
Mw8_unitarios_AltoSlip = datos_Mw8(1:32, :);       % Primeros 32 eventos de alto slip
Mw8_unitarios_BajoSlip = datos_Mw8(36:67, :);      % Eventos 36 al 67 como bajo slip

% --- Casos Mw9.5: unitarios de alto y bajo slip ---
Mw9_5_unitarios_AltoSlip = datos_Mw9_5(1:28, :);       % Primeros 28 eventos de alto slip
Mw9_5_unitarios_BajoSlip = datos_Mw9_5(29:48, :);      % Eventos 29 al 48 como bajo slip
alto = Mw9_5_unitarios_AltoSlip(14, :);                 % un evento de 9,5 de alto slip
bajo = Mw9_5_unitarios_BajoSlip(10, :);
h_max_alto_mw9_5=alto{1,3};
h_max_bajo_mw9_5=bajo{1,3};

% --- Calcular momento sismico de 9.5 Mw---
slip_mean_9_5_alto = alto{1,19};  % extrae el valor numérico, no la tabla
slip_mean_9_5_bajo = bajo{1,19};
Area_9_5_alto = alto{1,18};
Area_9_5_bajo = bajo{1,18};
Ms_9_5_alto = slip_mean_9_5_alto * Area_9_5_alto;
Ms_9_5_bajo = slip_mean_9_5_bajo * Area_9_5_bajo;

% --- latitudes norte y sur del sismo de  9.5 Mw---
ls_9_5_alto = round(alto{1,1}, 5);
ln_9_5_alto = round(alto{1,2}, 5);
ls_9_5_bajo = round(bajo{1,1}, 5);
ln_9_5_bajo = round(bajo{1,2}, 5);


% --- filtrar valores de 8mw que esten entre las latitudes N y S de 9.5 Mw---
filtro_alto = (Mw8_unitarios_AltoSlip.Lat_S >= ls_9_5_alto) & (Mw8_unitarios_AltoSlip.Lat_N <= ln_9_5_alto);
Mw8_alto = Mw8_unitarios_AltoSlip(filtro_alto, :);
filtro_bajo = (Mw8_unitarios_BajoSlip.Lat_S >= ls_9_5_bajo) & (Mw8_unitarios_BajoSlip.Lat_S <= ln_9_5_bajo);
Mw8_bajo = Mw8_unitarios_BajoSlip(filtro_alto, :);

n_alto = height(Mw8_alto);   % cantidad de filas de alto slip
n_bajo = height(Mw8_bajo);   % cantidad de filas de bajo slip

Ms_8_alto = zeros(n_alto,1);  % preasigno vector para momentos
Ms_8_bajo = zeros(n_bajo,1);

for i = 1:n_alto
    Ms_8_alto(i) = Mw8_alto{i,19} * Mw8_alto{i,18};
end

for i = 1:n_bajo
    Ms_8_bajo(i) = Mw8_bajo{i,19} * Mw8_bajo{i,18};
end

% Inicializar suma acumulada
Ms_acumulado_alto_mw8 = 0;
Ms_acumulado_bajo_mw8 = 0;
indice_alto = 0;
indice_bajo = 0;
h_max_acumulada_alto_mw8 = 0;
h_max_acumulada_bajo_mw8 = 0;

% Bandera para saber si se alcanzó el objetivo
alcanzado = false;

% Recorrer eventos y acumular
for i = 1:length(Ms_8_alto)
    Ms_acumulado_alto_mw8 = Ms_acumulado_alto_mw8 + Ms_8_alto(i);
    h_max_acumulada_alto_mw8 = h_max_acumulada_alto_mw8 + Mw8_alto{i,3};
    
    if Ms_acumulado_alto_mw8 >= Ms_9_5_alto
        indice_alto = i;
        alcanzado = true;
        break;
    end
end

% Mostrar resultados
if alcanzado
    fprintf('✅ Se necesitan %d eventos Mw 8 (alto slip) para igualar o superar el momento sísmico de Mw 9.5\n', indice_alto);
    fprintf('Momento acumulado: %.2e Nm\n', Ms_acumulado_alto_mw8);
else
    factor_alto = Ms_9_5_alto / Ms_acumulado_alto_mw8;
    fprintf('⚠️ No se alcanzó el momento sísmico de Mw 9.5 con todos los eventos disponibles.\n');
    fprintf('Se usaron los %d eventos disponibles (suma total: %.2e Nm)\n', length(Ms_8_alto), Ms_acumulado_alto_mw8);
    fprintf('El factor de escala necesario es: %.4f\n', factor_alto);
    h_max_acumulada_alto_mw8 = h_max_acumulada_alto_mw8 * factor_alto;
end

% Recorrer eventos y acumular
for i = 1:length(Ms_8_bajo)
    Ms_acumulado_bajo_mw8 = Ms_acumulado_bajo_mw8 + Ms_8_bajo(i);
    h_max_acumulada_bajo_mw8 = h_max_acumulada_bajo_mw8 + Mw8_bajo{i,3};
    
    if Ms_acumulado_bajo_mw8 >= Ms_9_5_bajo
        indice_bajo = i;
        alcanzado = true;
        break;
    end
end

% Mostrar resultados
if alcanzado
    fprintf('✅ Se necesitan %d eventos Mw 8 (bajo slip) para igualar o superar el momento sísmico de Mw 9.5\n', indice_bajo);
    fprintf('Momento acumulado: %.2e Nm\n', Ms_acumulado_bajo_mw8);
else
    factor_bajo = Ms_9_5_bajo / Ms_acumulado_bajo_mw8;
    fprintf('⚠️ No se alcanzó el momento sísmico de Mw 9.5 con todos los eventos disponibles.\n');
    fprintf('Se usaron los %d eventos disponibles (suma total: %.2e Nm)\n', length(Ms_8_bajo), Ms_acumulado_bajo_mw8);
    fprintf('El factor de escala necesario es: %.4f\n', factor_bajo);
    h_max_acumulada_bajo_mw8 = h_max_acumulada_bajo_mw8 * factor_bajo;
end

% --- Datos para el gráfico ---
alturas = [h_max_acumulada_alto_mw8, h_max_alto_mw9_5;
           h_max_acumulada_bajo_mw8, h_max_bajo_mw9_5];

labels = {'Alto Slip', 'Bajo Slip'};

% --- Gráfico de barras ---
figure;
bar(alturas)
set(gca, 'xticklabel', labels)
legend('Mw 8 acumulado', 'Mw 9.5 individual', 'Location', 'northwest')
ylabel('Altura máxima de tsunami (m)')
title('Comparación de alturas máximas de tsunami: Mw 8 acumulado vs Mw 9.5')
grid on

% --- Paso 1: Crear tabla combinada de eventos ---
eventos_alto = [Mw8_alto; alto];  % 7 eventos Mw 8 + 1 evento Mw 9.5 (alto slip)

figure(23)  % MISMA FIGURA QUE EL ORIGINAL

mags = unique(eventos_alto.Mw);  % Magnitudes únicas: debería ser [8.0; 9.5]

for i = 1:length(mags)
    idx = eventos_alto.Mw == mags(i);  % Filtro por magnitud actual
    scatter( ...
        ones(sum(idx),1) * mags(i) + eventos_alto.Dist_TC(idx) / (2 * max(eventos_alto.Dist_TC)), ...
        eventos_alto.Lat_S(idx), ...  % igual que el original, no usar lat media
        eventos_alto.Slip_max(idx) * 10, ...
        eventos_alto.H_max(idx), ...
        'filled', marker='s', MarkerEdgeColor='k');
    hold on
end

hold off
colorbar
title('Tamaño ∝ Slip × 10; Color ∝ H_{max} (altura máxima del tsunami)')
xlabel('Magnitud (con dispersión por distancia)')
ylabel('Latitud Sur')  % O la variable que estás graficando verticalmente
grid on

eventos_bajo = [Mw8_bajo; bajo];  % 7 eventos Mw 8 + 1 evento Mw 9.5 (bajo slip)

figure(24)  % Nueva figura (para que no se sobreescriba la 22 de alto)

mags_bajo = unique(eventos_bajo.Mw);  % Debería ser [8.0; 9.5]

for i = 1:length(mags_bajo)
    idx = eventos_bajo.Mw == mags_bajo(i);  % Filtro por magnitud
    scatter( ...
        ones(sum(idx),1) * mags_bajo(i) + eventos_bajo.Dist_TC(idx) / (2 * max(eventos_bajo.Dist_TC)), ...
        eventos_bajo.Lat_S(idx), ...
        eventos_bajo.Slip_max(idx) * 10, ...
        eventos_bajo.H_max(idx), ...
        'filled', marker='s', MarkerEdgeColor='k');
    hold on
end

hold off
colorbar
title('Tamaño ∝ Slip × 10; Color ∝ H_{max} (bajo slip)')
xlabel('Magnitud (con dispersión por distancia)')
ylabel('Latitud Sur')

