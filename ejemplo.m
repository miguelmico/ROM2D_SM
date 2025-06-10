clear; clc; close all;

%% 1. DEFINICIÓN DEL MODELO (solo elementos BEAM)
% Nodes=[NodID X  Y  Z] - Solo nodos únicos, sin duplicar
Nodes = [
1  0.0000         0.0000         0.0000;
2  0.0000         1.0000         0.0000;
3  0.7071         1.7071         0.0000;
4  1.4142         2.4142         0.0000;   % Este nodo tendrá restricción revolute
5  2.4142         2.4142         0.0000;
6  3.4142         2.4142         0.0000;
7  5.4142         2.4142         0.0000;
8  0.7778         1.6364         0.0000;   % Este nodo tendrá restricción revolute
9  2.4142         2.3142         0.0000;   % Este nodo tendrá restricción revolute
10 2 1 0
];

% Element types -> solo BEAM
Types=  {1        'beam'};

% Funciones auxiliares para secciones
function Iz = MomentoInerciaPerfilH(b,h,t_f,t_w)
    Iz = (((b*h^3)/12)-((b-t_w)*(h-2*t_f)^3/12));
end

function A = AreaPerfilH(b,h,t_f,t_w)
    A = 2 * (b * t_f) + (t_w * (h - 2 * t_f));
end

function Iz = MomentoInerciaPerfilCil(D,t)
    d = D - 2 * t;
    Iz = (pi / 64) * (D^4 - d^4);
end

function A = AreaPerfilCil(D,t)
    d = D - 2 * t;
    A = (pi / 4) * (D^2 - d^2);
end

% Sections=[SecID A      ky   kz   Ixx Iyy Izz       yt   yb   zt   zb]
Sections=  [1     AreaPerfilH(40,60,0.3,0.3)    Inf  Inf  0   0   MomentoInerciaPerfilH(40,60,0.3,0.3)  40/2  40/2  60/2  60/2;
            2     AreaPerfilH(40,40,0.3,0.3)    Inf  Inf  0   0   MomentoInerciaPerfilH(40,40,0.3,0.3)  40/2  40/2  40/2  40/2;
            3     AreaPerfilH(40,30,0.2,0.2)    Inf  Inf  0   0   MomentoInerciaPerfilH(40,30,0.2,0.2)  40/2  40/2  30/2  30/2;
            4     AreaPerfilH(40,20,0.2,0.2)    Inf  Inf  0   0   MomentoInerciaPerfilH(40,20,0.2,0.2)  40/2  40/2  20/2  20/2;
            5     3*10                          Inf  Inf  0   0   3*10^3/12                             3/2   3/2   10/2  10/2;
            6     AreaPerfilCil(10,0.4)         Inf  Inf  0   0   MomentoInerciaPerfilCil(10,0.4)       5     5     5     5];

Sections(:,2) = Sections(:,2) / 10^4;
Sections(:,7) = Sections(:,7) / 10^8;
Sections(:,8:11) = Sections(:,8:11) / 10^2;

% Materials=[MatID     E nu rho];
Materials=  [1      30e6 0.2 7850;               % concrete
             2     200e9 0.3 7850];              % steel

% Elements=[EltID TypID SecID MatID n1 n2 n3] - TODOS BEAM (tipo 1)
% NOTA: Usar nodos únicos, la función duplicará automáticamente para restricciones
Elements=  [1     1     1     2     1  2  10;              
            2     1     2     2     2  3  10;
            3     1     2     2     3  4  10;
            4     1     3     2     4  5  10;   % Este elemento usa nodo 4 (será duplicado)
            5     1     3     2     5  6  10;
            6     1     4     2     6  7  10;
            7     1     5     2     8  3  10;   % Este elemento usa nodo 8 (será duplicado)           
            8     1     5     2     9  5  10;   % Este elemento usa nodo 9 (será duplicado)
            9     1     6     2     8  9  10;   % Este elemento usa nodos 8 y 9 (serán duplicados)
          ];

% Restricciones revolute simplificadas
% Formato: [RestriccionID NodoID TipoRestriccion]
Constraints = {1, 4, 'revolute';    % Restricción revolute en nodo 4
               2, 8, 'revolute';    % Restricción revolute en nodo 8  
               3, 9, 'revolute'};   % Restricción revolute en nodo 9

%% 2. CREAR Y VALIDAR MODELO FEM
fprintf('=== PASO 1: CREANDO Y VALIDANDO MODELO FEM ===\n');
ModelData = createFEMModel(Nodes, Elements, Types, Sections, Materials, Constraints);

%% 3. ENSAMBLAJE SIMPLIFICADO PARA 2D
fprintf('=== PASO 2: ENSAMBLAJE PARA ANÁLISIS 2D ===\n');
MatrixData = assembleFEMMatrices(ModelData);

%% 4. VERIFICACIÓN DE PROCESAMIENTO AUTOMÁTICO
fprintf('=== VERIFICACIÓN DE PROCESAMIENTO ===\n');
fprintf('Datos originales vs procesados:\n');
fprintf('- Nodos originales: %d\n', size(ModelData.Original.Nodes, 1));
fprintf('- Nodos procesados: %d\n', size(ModelData.Nodes, 1));
fprintf('- Nodos duplicados: %d\n', length(ModelData.Info.DuplicateNodes));

fprintf('\nRestricciones revolute:\n');
for i = 1:size(ModelData.Original.Constraints, 1)
    constraint_id = ModelData.Original.Constraints{i, 1};
    node_id = ModelData.Original.Constraints{i, 2};
    constraint_type = ModelData.Original.Constraints{i, 3};
    fprintf('- Restricción %d: %s en nodo %d\n', constraint_id, constraint_type, node_id);
end

fprintf('\nRestricciones StaBIL generadas: %d\n', size(ModelData.Constraints, 1));
if ~isempty(ModelData.Constraints)
    fprintf('Primeras restricciones generadas:\n');
    for i = 1:min(6, size(ModelData.Constraints, 1))
        constr = ModelData.Constraints(i, :);
        fprintf('  [%g %g %g %g %g]\n', constr);
    end
    if size(ModelData.Constraints, 1) > 6
        fprintf('  ... y %d más\n', size(ModelData.Constraints, 1) - 6);
    end
end

fprintf('\nNodos duplicados generados:\n');
if ~isempty(ModelData.Info.DuplicateNodes)
    for i = 1:length(ModelData.Info.DuplicateNodes)
        dup_node = ModelData.Info.DuplicateNodes(i);
        fprintf('- Nodo duplicado: %d\n', dup_node);
    end
end

%% 5. VERIFICACIÓN DE MATRICES
fprintf('\n=== INFORMACIÓN DE MATRICES ===\n');
fprintf('Matriz ensamblada: %s\n', MatrixData.Info.Assembly_Date);
fprintf('- DOFs iniciales: %d\n', MatrixData.Properties.DOFs_Initial);
fprintf('- DOFs después de eliminar 2D: %d\n', MatrixData.Properties.DOFs_After_Removal);
fprintf('- DOFs activos: %d\n', MatrixData.Properties.DOFs_Active);
fprintf('- Tamaño matrices: %dx%d\n', size(MatrixData.K));

fprintf('\nPropiedades:\n');
fprintf('- K simétrica: %s\n', mat2str(MatrixData.Properties.Symmetric_K));
fprintf('- M simétrica: %s\n', mat2str(MatrixData.Properties.Symmetric_M));
if ~isnan(MatrixData.Properties.Condition_K)
    fprintf('- Condicionamiento K: %.2e\n', MatrixData.Properties.Condition_K);
    fprintf('- Condicionamiento M: %.2e\n', MatrixData.Properties.Condition_M);
end

%% 6. MOSTRAR DOFs ELIMINADOS
fprintf('\n=== DOFs ELIMINADOS PARA 2D ===\n');
fprintf('Se eliminaron automáticamente:\n');
for i = 1:length(MatrixData.DOF_removed)
    dof_num = round((MatrixData.DOF_removed(i) - floor(MatrixData.DOF_removed(i))) * 100);
    dof_names = {'', 'Ux', 'Uy', 'Uz', 'Rx', 'Ry', 'Rz'};
    if dof_num <= length(dof_names)
        fprintf('- DOF %02d (%s) de todos los nodos\n', dof_num, dof_names{dof_num+1});
    end
end

%% 7. REDUCCIÓN CRAIG-BAMPTON
fprintf('\n=== PASO 3: REDUCCIÓN CRAIG-BAMPTON ===\n');

% Definir nodos maestros/interfaz (como en tu ejemplo: nodos 1 y 7)
master_nodes = [1; 7];

% Realizar reducción Craig-Bampton
num_modes = 10; % Número de modos normales a retener
ReducedData = craigBamptonReduction(MatrixData, ModelData, master_nodes, num_modes);

%% 8. VERIFICACIÓN DE REDUCCIÓN
fprintf('\n=== VERIFICACIÓN DE REDUCCIÓN ===\n');
fprintf('Reducción completada: %s\n', ReducedData.Info.Reduction_Date);
fprintf('Matrices reducidas para Simscape:\n');
fprintf('- K_reduced: %dx%d\n', size(ReducedData.K_reduced));
fprintf('- M_reduced: %dx%d\n', size(ReducedData.M_reduced));
fprintf('- Ratio de reducción: %.1f%%\n', (1 - ReducedData.reduction_ratio) * 100);

if ~isempty(ReducedData.frequencies)
    fprintf('\nPrimeras frecuencias naturales (Hz):\n');
    for i = 1:min(5, length(ReducedData.frequencies))
        fprintf('  Modo %d: %.2f Hz\n', i, ReducedData.frequencies(i));
    end
end

fprintf('\nNodos interfaz (coordenadas):\n');
for i = 1:size(ReducedData.nodes_interface, 1)
    fprintf('  Nodo %d: [%.3f, %.3f, %.3f]\n', ReducedData.master_nodes(i), ...
            ReducedData.nodes_interface(i, :));
end

% Guardar datos completos
save('FEM_Model_2D_BEAM_Revolute_CB.mat', 'ModelData', 'MatrixData', 'ReducedData');
fprintf('\nDatos guardados en: FEM_Model_2D_BEAM_Revolute_CB.mat\n');

%% 8. VISUALIZACIÓN (si está disponible)
if exist('plotnodes', 'file') && exist('plotelem', 'file')
    figure('Name', 'Modelo 2D con restricciones revolute', 'Position', [100 100 1200 500]);
    
    % Subplot 1: Modelo original
    subplot(1,2,1);
    plotnodes(ModelData.Original.Nodes);
    hold on;
    plotelem(ModelData.Original.Nodes, ModelData.Original.Elements, Types);
    title('Modelo Original (nodos únicos)');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    grid on; axis equal;
    
    % Subplot 2: Modelo procesado
    subplot(1,2,2);
    plotnodes(ModelData.Nodes);
    hold on;
    plotelem(ModelData.Nodes, ModelData.Elements, Types);
    title('Modelo Procesado (con nodos duplicados)');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    grid on; axis equal;
end

fprintf('\n=== PROCESO COMPLETADO ===\n');
fprintf('Modelo 2D con restricciones revolute y reducción Craig-Bampton listo.\n');
fprintf('Matrices reducidas disponibles para Simscape Multibody ROM.\n');
fprintf('Siguiente paso: Export a formato Simscape\n');