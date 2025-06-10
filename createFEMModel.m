function [ModelData] = createFEMModel(Nodes, Elements, Types, Sections, Materials, Constraints)
%CREATEFEMMODEL Crea y valida completamente un modelo FEM con restricciones revolute
%
% SINTAXIS:
%   [ModelData] = createFEMModel(Nodes, Elements, Types, Sections, Materials, Constraints)
%
% ENTRADAS:
%   Nodes       - Matriz de nodos [NodID X Y Z] (obligatorio)
%   Elements    - Matriz de elementos [EltID TypID SecID MatID n1 n2 n3] (obligatorio)
%   Types       - Cell array de tipos {EltTypID EltName} (obligatorio)
%   Sections    - Matriz de secciones [SecID A ky kz Ixx Iyy Izz yt yb zt zb] (obligatorio)
%   Materials   - Matriz de materiales [MatID E nu rho] (obligatorio)
%   Constraints - Matriz de restricciones [RestriccionID NodoID TipoRestriccion] (opcional)
%                 TipoRestriccion: 'revolute' (otros tipos en futuras versiones)
%
% SALIDAS:
%   ModelData   - Estructura validada con todos los datos del modelo
%
% EJEMPLO:
%   % Restricción revolute en nodo 4
%   Constraints = [1 4 'revolute'; 2 8 'revolute'];
%   ModelData = createFEMModel(Nodes, Elements, Types, Sections, Materials, Constraints);

%% VALIDACIÓN DE ARGUMENTOS DE ENTRADA
if nargin < 5
    error('Se requieren al menos 5 argumentos: Nodes, Elements, Types, Sections, Materials');
end

if nargin < 6 || isempty(Constraints)
    Constraints = [];
end

fprintf('=== CREANDO Y VALIDANDO MODELO FEM ===\n');

%% VALIDACIÓN DE ENTRADAS BÁSICAS
fprintf('Validando entradas básicas...\n');
validateNodes(Nodes);
validateTypes(Types);
validateSections(Sections);
validateMaterials(Materials);
validateElements(Elements, Nodes, Types, Sections, Materials);

%% VALIDACIÓN DE RESTRICCIONES SIMPLIFICADAS
if ~isempty(Constraints)
    fprintf('Validando restricciones...\n');
    validateSimplifiedConstraints(Constraints, Nodes);
end

%% PROCESAMIENTO DE RESTRICCIONES REVOLUTE
[NodesProcessed, ElementsProcessed, ConstraintsStaBIL] = processRevoluteConstraints(Nodes, Elements, Constraints);

%% CREAR ESTRUCTURA DEL MODELO
fprintf('Creando estructura del modelo...\n');
ModelData = struct();

% Almacenar datos procesados (con nodos duplicados y restricciones StaBIL)
ModelData.Nodes = NodesProcessed;
ModelData.Elements = ElementsProcessed;
ModelData.Types = Types;
ModelData.Sections = Sections;
ModelData.Materials = Materials;
ModelData.Constraints = ConstraintsStaBIL;

% Almacenar datos originales para referencia
ModelData.Original.Nodes = Nodes;
ModelData.Original.Elements = Elements;
ModelData.Original.Constraints = Constraints;

% Información básica del modelo
ModelData.Info.NumNodes = size(NodesProcessed, 1);
ModelData.Info.NumElements = size(ElementsProcessed, 1);
ModelData.Info.NumTypes = size(Types, 1);
ModelData.Info.NumSections = size(Sections, 1);
ModelData.Info.NumMaterials = size(Materials, 1);
if ~isempty(ConstraintsStaBIL)
    ModelData.Info.NumConstraints = size(ConstraintsStaBIL, 1);
else
    ModelData.Info.NumConstraints = 0;
end

% Información sobre restricciones revolute
if ~isempty(Constraints)
    ModelData.Info.NumRevoluteConstraints = size(Constraints, 1);
    ModelData.Info.RevoluteNodes = Constraints(:, 2);
else
    ModelData.Info.NumRevoluteConstraints = 0;
    ModelData.Info.RevoluteNodes = [];
end

% Detectar nodos duplicados (ahora generados automáticamente)
node_ids = NodesProcessed(:, 1);
original_node_ids = Nodes(:, 1);
duplicated_nodes = setdiff(node_ids, original_node_ids);
ModelData.Info.DuplicateNodes = duplicated_nodes;

% Metadatos
ModelData.Info.CreationDate = datestr(now);
ModelData.Info.MATLABVersion = version;
ModelData.Info.Validated = true;

%% RESUMEN FINAL
fprintf('\n=== MODELO FEM CREADO Y VALIDADO ===\n');
fprintf('Nodos originales: %d\n', size(Nodes, 1));
fprintf('Nodos procesados: %d\n', ModelData.Info.NumNodes);
fprintf('Nodos duplicados: %d\n', length(ModelData.Info.DuplicateNodes));
fprintf('Elementos: %d\n', ModelData.Info.NumElements);
fprintf('Tipos de elementos: %d\n', ModelData.Info.NumTypes);
fprintf('Secciones: %d\n', ModelData.Info.NumSections);
fprintf('Materiales: %d\n', ModelData.Info.NumMaterials);
fprintf('Restricciones revolute: %d\n', ModelData.Info.NumRevoluteConstraints);
fprintf('Restricciones StaBIL: %d\n', ModelData.Info.NumConstraints);
fprintf('=============================================\n\n');

fprintf('✓ Modelo FEM creado y validado exitosamente.\n');

end

%% FUNCIONES DE VALIDACIÓN BÁSICAS (igual que antes)

function validateNodes(Nodes)
if ~isnumeric(Nodes)
    error('Nodes debe ser una matriz numérica');
end
if size(Nodes, 2) ~= 4
    error('Nodes debe tener exactamente 4 columnas: [NodID X Y Z]');
end
if size(Nodes, 1) < 1
    error('Debe haber al menos un nodo definido');
end
node_ids = Nodes(:, 1);
if any(node_ids <= 0) || any(mod(node_ids, 1) ~= 0)
    error('Los IDs de nodos deben ser enteros positivos');
end
coords = Nodes(:, 2:4);
if any(~isfinite(coords(:)))
    error('Las coordenadas de nodos deben ser números reales finitos');
end
fprintf('  ✓ Nodos: %d nodos validados\n', size(Nodes, 1));
end

function validateTypes(Types)
if ~iscell(Types)
    error('Types debe ser un cell array');
end
if size(Types, 2) ~= 2
    error('Types debe tener exactamente 2 columnas: {EltTypID EltName}');
end
if size(Types, 1) < 1
    error('Debe haber al menos un tipo de elemento definido');
end
for i = 1:size(Types, 1)
    if ~isnumeric(Types{i, 1}) || Types{i, 1} <= 0 || mod(Types{i, 1}, 1) ~= 0
        error('Types{%d,1} debe ser un entero positivo (ID del tipo)', i);
    end
    if ~ischar(Types{i, 2}) && ~isstring(Types{i, 2})
        error('Types{%d,2} debe ser un string (nombre del tipo)', i);
    end
    if isempty(Types{i, 2})
        error('Types{%d,2} no puede estar vacío (nombre del tipo)', i);
    end
end
type_ids = cell2mat(Types(:, 1));
if length(unique(type_ids)) ~= length(type_ids)
    error('Los IDs de tipos de elementos no pueden estar duplicados');
end
fprintf('  ✓ Tipos: %d tipos validados\n', size(Types, 1));
end

function validateSections(Sections)
if ~isnumeric(Sections)
    error('Sections debe ser una matriz numérica');
end
if size(Sections, 2) ~= 11
    error('Sections debe tener exactamente 11 columnas: [SecID A ky kz Ixx Iyy Izz yt yb zt zb]');
end
if size(Sections, 1) < 1
    error('Debe haber al menos una sección definida');
end
sec_ids = Sections(:, 1);
if any(sec_ids <= 0) || any(mod(sec_ids, 1) ~= 0)
    error('Los IDs de secciones deben ser enteros positivos');
end
if length(unique(sec_ids)) ~= length(sec_ids)
    error('Los IDs de secciones no pueden estar duplicados');
end
areas = Sections(:, 2);
if any(areas <= 0)
    error('Las áreas de secciones deben ser positivas');
end
moments = Sections(:, 7);
if any(moments < 0)
    error('Los momentos de inercia no pueden ser negativos');
end
% Verificar valores finitos (permitir Inf en ky, kz para análisis 2D)
sections_to_check = Sections;
sections_to_check(:, 3:4) = []; % Remover ky, kz de la verificación
if any(~isfinite(sections_to_check(:)))
    error('Las propiedades de secciones deben ser números finitos (excepto ky, kz que pueden ser Inf)');
end
fprintf('  ✓ Secciones: %d secciones validadas\n', size(Sections, 1));
end

function validateMaterials(Materials)
if ~isnumeric(Materials)
    error('Materials debe ser una matriz numérica');
end
if size(Materials, 2) ~= 4
    error('Materials debe tener exactamente 4 columnas: [MatID E nu rho]');
end
if size(Materials, 1) < 1
    error('Debe haber al menos un material definido');
end
mat_ids = Materials(:, 1);
if any(mat_ids <= 0) || any(mod(mat_ids, 1) ~= 0)
    error('Los IDs de materiales deben ser enteros positivos');
end
if length(unique(mat_ids)) ~= length(mat_ids)
    error('Los IDs de materiales no pueden estar duplicados');
end
E = Materials(:, 2);
nu = Materials(:, 3);
rho = Materials(:, 4);
if any(E <= 0)
    error('El módulo de Young debe ser positivo');
end
if any(nu <= -1) || any(nu >= 0.5)
    error('El coeficiente de Poisson debe estar en el rango (-1, 0.5)');
end
if any(rho <= 0)
    error('La densidad debe ser positiva');
end
if any(~isfinite(Materials(:)))
    error('Todas las propiedades de materiales deben ser números finitos');
end
fprintf('  ✓ Materiales: %d materiales validados\n', size(Materials, 1));
end

function validateElements(Elements, Nodes, Types, Sections, Materials)
if ~isnumeric(Elements)
    error('Elements debe ser una matriz numérica');
end
if size(Elements, 2) ~= 7
    error('Elements debe tener exactamente 7 columnas: [EltID TypID SecID MatID n1 n2 n3]');
end
if size(Elements, 1) < 1
    error('Debe haber al menos un elemento definido');
end

node_ids = Nodes(:, 1);
type_ids = cell2mat(Types(:, 1));
sec_ids = Sections(:, 1);
mat_ids = Materials(:, 1);

for i = 1:size(Elements, 1)
    elem = Elements(i, :);
    if elem(1) <= 0 || mod(elem(1), 1) ~= 0
        error('Element %d: El ID del elemento debe ser un entero positivo', i);
    end
    if ~ismember(elem(2), type_ids)
        error('Element %d: Tipo de elemento %d no está definido en Types', i, elem(2));
    end
    if ~ismember(elem(3), sec_ids)
        error('Element %d: Sección %d no está definida en Sections', i, elem(3));
    end
    if ~ismember(elem(4), mat_ids)
        error('Element %d: Material %d no está definido en Materials', i, elem(4));
    end
    element_nodes = elem(5:7);
    for j = 1:3
        if element_nodes(j) ~= 0 && ~ismember(element_nodes(j), node_ids)
            error('Element %d: Nodo %d no está definido en Nodes', i, element_nodes(j));
        end
    end
    valid_nodes = element_nodes(element_nodes ~= 0);
    if length(valid_nodes) < 2
        error('Element %d: Debe tener al menos 2 nodos válidos (no cero)', i);
    end
end

elem_ids = Elements(:, 1);
if length(unique(elem_ids)) ~= length(elem_ids)
    error('Los IDs de elementos no pueden estar duplicados');
end
fprintf('  ✓ Elementos: %d elementos validados\n', size(Elements, 1));
end

%% NUEVAS FUNCIONES PARA RESTRICCIONES REVOLUTE

function validateSimplifiedConstraints(Constraints, Nodes)
%VALIDATESIMPLIFIEDCONSTRAINTS Valida las restricciones en formato simplificado

if isempty(Constraints)
    return;
end

% Verificar formato básico
if size(Constraints, 2) ~= 3
    error('Constraints debe tener exactamente 3 columnas: [RestriccionID NodoID TipoRestriccion]');
end

node_ids = Nodes(:, 1);

for i = 1:size(Constraints, 1)
    % ID de restricción
    if ~isnumeric(Constraints{i, 1}) || Constraints{i, 1} <= 0 || mod(Constraints{i, 1}, 1) ~= 0
        error('Constraint %d: El ID de restricción debe ser un entero positivo', i);
    end
    
    % ID de nodo
    if ~isnumeric(Constraints{i, 2}) || ~ismember(Constraints{i, 2}, node_ids)
        error('Constraint %d: El nodo %d no está definido', i, Constraints{i, 2});
    end
    
    % Tipo de restricción
    if ~ischar(Constraints{i, 3}) && ~isstring(Constraints{i, 3})
        error('Constraint %d: El tipo de restricción debe ser un string', i);
    end
    
    constraint_type = lower(string(Constraints{i, 3}));
    if constraint_type ~= "revolute"
        error('Constraint %d: Solo se soporta tipo "revolute" por ahora', i);
    end
end

% Verificar que no hay IDs duplicados
constraint_ids = cell2mat(Constraints(:, 1));
if length(unique(constraint_ids)) ~= length(constraint_ids)
    error('Los IDs de restricciones no pueden estar duplicados');
end

fprintf('  ✓ Restricciones: %d restricciones revolute validadas\n', size(Constraints, 1));
end

function [NodesProcessed, ElementsProcessed, ConstraintsStaBIL] = processRevoluteConstraints(Nodes, Elements, Constraints)
%PROCESSREVOLUTECONSTRAINTS Procesa las restricciones revolute automáticamente

NodesProcessed = Nodes;
ElementsProcessed = Elements;
ConstraintsStaBIL = [];

if isempty(Constraints)
    return;
end

fprintf('Procesando restricciones revolute...\n');

% Para cada restricción revolute
for i = 1:size(Constraints, 1)
    constraint_id = Constraints{i, 1};
    node_id = Constraints{i, 2};
    constraint_type = lower(string(Constraints{i, 3}));
    
    if constraint_type == "revolute"
        fprintf('  Procesando restricción revolute en nodo %d...\n', node_id);
        
        % Encontrar elementos que usan este nodo
        elements_using_node = findElementsUsingNode(Elements, node_id);
        
        if length(elements_using_node) < 2
            warning('Nodo %d solo se usa en %d elemento(s). Restricción revolute requiere al menos 2 elementos.', ...
                    node_id, length(elements_using_node));
            continue;
        end
        
        % Crear nodos duplicados para cada elemento (excepto el primero)
        [NodesProcessed, ElementsProcessed, new_constraints] = createRevoluteConstraint(...
            NodesProcessed, ElementsProcessed, node_id, elements_using_node);
        
        % Agregar las nuevas restricciones StaBIL
        ConstraintsStaBIL = [ConstraintsStaBIL; new_constraints];
        
        fprintf('    ✓ Creados %d nodos duplicados y %d restricciones StaBIL\n', ...
                length(elements_using_node)-1, size(new_constraints, 1));
    end
end

fprintf('✓ Restricciones revolute procesadas correctamente.\n');
end

function elements_using_node = findElementsUsingNode(Elements, node_id)
%FINDELEMENTSUSINGNODE Encuentra todos los elementos que usan un nodo específico

elements_using_node = [];

for i = 1:size(Elements, 1)
    element_nodes = Elements(i, 5:7);
    if any(element_nodes == node_id)
        elements_using_node = [elements_using_node, Elements(i, 1)];
    end
end
end

function [NodesUpdated, ElementsUpdated, ConstraintsStaBIL] = createRevoluteConstraint(Nodes, Elements, node_id, elements_using_node)
%CREATEREVOLUTECONSTRAINT Crea la restricción revolute duplicando nodos

NodesUpdated = Nodes;
ElementsUpdated = Elements;
ConstraintsStaBIL = [];

% Encontrar las coordenadas del nodo original
node_idx = find(Nodes(:, 1) == node_id);
if isempty(node_idx)
    error('Nodo %d no encontrado', node_id);
end
original_coords = Nodes(node_idx, 2:4);

% El primer elemento mantiene el nodo original
% Los demás elementos usan nodos duplicados
duplicated_node_ids = [node_id]; % Incluir el nodo original para las restricciones

for i = 2:length(elements_using_node)
    element_id = elements_using_node(i);
    
    % Generar ID del nodo duplicado: node_id * 10 + element_id
    % Si hay conflicto, usar node_id * 100 + element_id
    new_node_id = node_id * 10 + element_id;
    
    % Verificar que el nuevo ID no existe
    while any(NodesUpdated(:, 1) == new_node_id)
        new_node_id = new_node_id + 100;
    end
    
    % Agregar el nodo duplicado con las mismas coordenadas
    NodesUpdated = [NodesUpdated; new_node_id, original_coords];
    duplicated_node_ids = [duplicated_node_ids, new_node_id];
    
    % Actualizar el elemento para usar el nodo duplicado
    elem_idx = find(Elements(:, 1) == element_id);
    if ~isempty(elem_idx)
        element_nodes = ElementsUpdated(elem_idx, 5:7);
        % Reemplazar todas las ocurrencias del nodo original
        element_nodes(element_nodes == node_id) = new_node_id;
        ElementsUpdated(elem_idx, 5:7) = element_nodes;
    end
end

% Crear restricciones StaBIL entre todos los nodos duplicados
% Para restricción revolute: igualar Ux y Uy entre todos los nodos
for i = 2:length(duplicated_node_ids)
    master_node = duplicated_node_ids(1);
    slave_node = duplicated_node_ids(i);
    
    % Restricción Ux: 0 = 1*master.01 - 1*slave.01
    constraint_ux = [0, 1, master_node + 0.01, -1, slave_node + 0.01];
    
    % Restricción Uy: 0 = 1*master.02 - 1*slave.02  
    constraint_uy = [0, 1, master_node + 0.02, -1, slave_node + 0.02];
    
    ConstraintsStaBIL = [ConstraintsStaBIL; constraint_ux; constraint_uy];
end

end