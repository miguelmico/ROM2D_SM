function [ReducedData] = craigBamptonReduction(MatrixData, ModelData, master_nodes, num_modes)
%CRAIGBAMPTONREDUCTION Realiza reducción Craig-Bampton del modelo FEM
%
% SINTAXIS:
%   [ReducedData] = craigBamptonReduction(MatrixData, ModelData, master_nodes, num_modes)
%
% ENTRADAS:
%   MatrixData    - Estructura con matrices K, M, DOF (salida de assembleFEMMatrices)
%   ModelData     - Estructura con datos del modelo (salida de createFEMModel)
%   master_nodes  - Vector con IDs de nodos maestros/interfaz [n1; n2; ...]
%   num_modes     - Número de modos normales a retener (opcional, por defecto auto)
%
% SALIDAS:
%   ReducedData   - Estructura con matrices reducidas y datos de reducción
%
% EJEMPLO:
%   master_nodes = [1; 7];  % Nodos maestros/interfaz
%   ReducedData = craigBamptonReduction(MatrixData, ModelData, master_nodes, 10);

% Validación de entradas
if nargin < 3
    error('Se requieren al menos 3 argumentos: MatrixData, ModelData, master_nodes');
end

if nargin < 4 || isempty(num_modes)
    num_modes = []; % Determinar automáticamente
end

fprintf('=== REDUCCIÓN CRAIG-BAMPTON ===\n');

%% Extraer datos
K = MatrixData.K;
M = MatrixData.M;
DOF = MatrixData.DOF;

fprintf('Matriz original: %dx%d\n', size(K, 1), size(K, 2));
fprintf('Nodos maestros: %s\n', mat2str(master_nodes'));

%% Determinar DOFs maestros y esclavos
fprintf('Determinando DOFs maestros y esclavos...\n');
[dof_m, dof_s, nodos_interfaz] = partitionDOFs(DOF, master_nodes, ModelData);

fprintf('DOFs maestros: %d\n', length(dof_m));
fprintf('DOFs esclavos: %d\n', length(dof_s));

%% Particionar matrices
fprintf('Particionando matrices...\n');
[Kmm, Kss, Kms, Ksm, Mmm, Mss, Mms, Msm] = partitionMatrices(K, M, dof_m, dof_s);

%% Condensación de Guyan (modos estáticos)
fprintf('Calculando modos estáticos de Guyan...\n');

% Verificar condicionamiento de Kss
if issparse(Kss)
    cond_Kss = 1 / condest(Kss);
else
    cond_Kss = rcond(Kss);
end

if cond_Kss < 1e-12
    warning('Matriz Kss mal condicionada (rcond ≈ %.2e). Usando pseudoinversa.', cond_Kss);
    if issparse(Kss)
        Kss_inv = pinv(full(Kss));
    else
        Kss_inv = pinv(Kss);
    end
else
    if issparse(Kss)
        Kss_inv = Kss \ eye(size(Kss, 1));  % Más eficiente para matrices dispersas
    else
        Kss_inv = inv(Kss);
    end
end

% Matriz de transformación de Guyan
TG = [eye(length(dof_m)); -Kss_inv * Ksm];

% Matrices condensadas de Guyan
KG = TG' * [Kmm, Kms; Ksm, Kss] * TG;
MG = TG' * [Mmm, Mms; Msm, Mss] * TG;

fprintf('✓ Matrices de Guyan calculadas (%dx%d)\n', size(KG, 1), size(KG, 2));

%% Cálculo de modos normales (Craig-Bampton)
fprintf('Calculando modos normales...\n');

% Problema de eigenvalores para DOFs esclavos con frontera fija
if issparse(Mss)
    cond_Mss = 1 / condest(Mss);
else
    cond_Mss = rcond(Mss);
end

if cond_Mss < 1e-12
    warning('Matriz Mss mal condicionada. Usando regularización.');
    if issparse(Mss)
        Mss_reg = Mss + 1e-10 * speye(size(Mss));
    else
        Mss_reg = Mss + 1e-10 * eye(size(Mss));
    end
else
    Mss_reg = Mss;
end

% Determinar número de modos automáticamente si no se especifica
if isempty(num_modes)
    num_modes = min(round(length(dof_s) * 0.1), 20); % 10% de DOFs esclavos, máximo 20
    num_modes = max(num_modes, 1); % Mínimo 1 modo
    fprintf('Número de modos determinado automáticamente: %d\n', num_modes);
end

% Asegurar que no pedimos más modos de los disponibles
num_modes = min(num_modes, length(dof_s));

try
    % Resolver problema de eigenvalores: Kss * phi = lambda * Mss * phi
    [phi, lambda] = eigs(Kss, Mss_reg, num_modes, 'smallestabs');
    
    % Ordenar por frecuencia
    [lambda_sorted, idx] = sort(diag(lambda));
    phi = phi(:, idx);
    lambda = lambda_sorted;
    
    % Frecuencias naturales
    frequencies = sqrt(abs(lambda)) / (2 * pi);
    
    fprintf('✓ %d modos calculados\n', num_modes);
    fprintf('Frecuencias (Hz): ');
    if length(frequencies) <= 5
        fprintf('%s\n', mat2str(frequencies', 3));
    else
        fprintf('%s ... %s\n', mat2str(frequencies(1:3)', 3), mat2str(frequencies(end)', 3));
    end
    
catch ME
    warning('Error en cálculo de eigenvalores: %s', ME.message);
    fprintf('Usando solo condensación de Guyan (sin modos normales)\n');
    phi = [];
    lambda = [];
    frequencies = [];
    num_modes = 0;
end

%% Construcción de matriz de transformación Craig-Bampton
fprintf('Construyendo matriz de transformación Craig-Bampton...\n');

if ~isempty(phi)
    % Matriz de transformación completa: [modos estáticos | modos normales]
    T_CB = [eye(length(dof_m)), zeros(length(dof_m), num_modes);
            -Kss_inv * Ksm, phi];
else
    % Solo modos estáticos (Guyan)
    T_CB = TG;
end

%% Matrices reducidas finales
if ~isempty(phi)
    % Matrices Craig-Bampton completas
    K_CB = T_CB' * [Kmm, Kms; Ksm, Kss] * T_CB;
    M_CB = T_CB' * [Mmm, Mms; Msm, Mss] * T_CB;
    
    fprintf('✓ Reducción Craig-Bampton completa: %dx%d → %dx%d\n', ...
            size(K, 1), size(K, 2), size(K_CB, 1), size(K_CB, 2));
else
    % Solo Guyan
    K_CB = KG;
    M_CB = MG;
    
    fprintf('✓ Reducción Guyan: %dx%d → %dx%d\n', ...
            size(K, 1), size(K, 2), size(K_CB, 1), size(K_CB, 2));
end

%% Organizar datos de salida
ReducedData = struct();

% Matrices reducidas
ReducedData.K_reduced = K_CB;
ReducedData.M_reduced = M_CB;
ReducedData.T_transformation = T_CB;

% Matrices de Guyan (para referencia)
ReducedData.K_Guyan = KG;
ReducedData.M_Guyan = MG;
ReducedData.T_Guyan = TG;

% Datos de partición
ReducedData.DOF_master = dof_m;
ReducedData.DOF_slave = dof_s;
ReducedData.nodes_interface = nodos_interfaz;
ReducedData.master_nodes = master_nodes;

% Datos modales
ReducedData.num_modes = num_modes;
ReducedData.frequencies = frequencies;
ReducedData.mode_shapes = phi;
ReducedData.eigenvalues = lambda;

% Información de reducción
ReducedData.reduction_ratio = size(K_CB, 1) / size(K, 1);
ReducedData.original_size = size(K, 1);
ReducedData.reduced_size = size(K_CB, 1);

% Metadatos
ReducedData.Info.Reduction_Date = datestr(now);
ReducedData.Info.Method = 'Craig-Bampton';
ReducedData.Info.Original_Matrix_Data = MatrixData.Info;

%% Resumen final
fprintf('\n=== RESUMEN DE REDUCCIÓN ===\n');
fprintf('Método: Craig-Bampton\n');
fprintf('Tamaño original: %d DOFs\n', ReducedData.original_size);
fprintf('Tamaño reducido: %d DOFs\n', ReducedData.reduced_size);
fprintf('Ratio de reducción: %.1f%%\n', (1 - ReducedData.reduction_ratio) * 100);
fprintf('Nodos interfaz: %d\n', length(master_nodes));
fprintf('Modos normales: %d\n', num_modes);
if ~isempty(frequencies)
    fprintf('Rango de frecuencias: %.2f - %.2f Hz\n', min(frequencies), max(frequencies));
end
fprintf('================================\n\n');

fprintf('✓ Reducción Craig-Bampton completada exitosamente.\n');

end

%% FUNCIONES AUXILIARES

function [dof_m, dof_s, nodos_interfaz] = partitionDOFs(DOF, master_nodes, ModelData)
%PARTITIONDOFS Particiona DOFs en maestros y esclavos

% Obtener coordenadas de nodos interfaz
Nodes = ModelData.Nodes;
nodos_interfaz = [];

for i = 1:length(master_nodes)
    node_id = master_nodes(i);
    node_idx = find(Nodes(:, 1) == node_id);
    if ~isempty(node_idx)
        nodos_interfaz = [nodos_interfaz; Nodes(node_idx, 2:4)];
    else
        warning('Nodo maestro %d no encontrado en la lista de nodos', node_id);
    end
end

% Determinar DOFs maestros
% Para análisis 2D con elementos BEAM: DOFs 01, 02, 06 por nodo
dof_m = [];
for i = 1:length(master_nodes)
    node_id = master_nodes(i);
    node_dofs = [node_id + 0.01; node_id + 0.02; node_id + 0.06];
    
    % Verificar que estos DOFs existen en el modelo
    for j = 1:length(node_dofs)
        dof_idx = find(abs(DOF - node_dofs(j)) < 1e-10);
        if ~isempty(dof_idx)
            dof_m = [dof_m; dof_idx];
        end
    end
end

% DOFs esclavos son todos los demás
all_dofs = (1:length(DOF))';
dof_s = all_dofs(~ismember(all_dofs, dof_m));

% Verificar que tenemos DOFs de ambos tipos
if isempty(dof_m)
    error('No se encontraron DOFs maestros válidos');
end

if isempty(dof_s)
    error('No se encontraron DOFs esclavos válidos');
end

end

function [Kmm, Kss, Kms, Ksm, Mmm, Mss, Mms, Msm] = partitionMatrices(K, M, dof_m, dof_s)
%PARTITIONMATRICES Particiona las matrices K y M

% Partición de matriz de rigidez
Kmm = K(dof_m, dof_m);
Kss = K(dof_s, dof_s);
Kms = K(dof_m, dof_s);
Ksm = K(dof_s, dof_m);

% Partición de matriz de masa
Mmm = M(dof_m, dof_m);
Mss = M(dof_s, dof_s);
Mms = M(dof_m, dof_s);
Msm = M(dof_s, dof_m);

% Verificar simetría
if norm(Kmm - Kmm', 'fro') > 1e-10 * norm(Kmm, 'fro')
    warning('Submatriz Kmm no es simétrica');
end

if norm(Kss - Kss', 'fro') > 1e-10 * norm(Kss, 'fro')
    warning('Submatriz Kss no es simétrica');
end

end