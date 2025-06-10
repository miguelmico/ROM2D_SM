function [MatrixData] = assembleFEMMatrices(ModelData)
%ASSEMBLEFEMMATRICES Ensambla matrices K y M usando funciones de StaBIL 
%
% SINTAXIS:
%   [MatrixData] = assembleFEMMatrices(ModelData)
%
% ENTRADAS:
%   ModelData     - Estructura del modelo (salida de defineFEMModel)
%
% SALIDAS:
%   MatrixData    - Estructura con matrices ensambladas y DOFs
%
% NOTA: Solo para análisis 2D con elementos BEAM
%       Elimina automáticamente DOFs 03, 04, 05 (Uz, Rx, Ry)
%       Mantiene DOFs 01, 02, 06 (Ux, Uy, Rz)
%
% EJEMPLO:
%   MatrixData = assembleFEMMatrices(ModelData);

fprintf('Ensamblando matrices para análisis 2D con elementos BEAM...\n');

%% Extraer datos del modelo
Nodes = ModelData.Nodes;
Elements = ModelData.Elements;
Types = ModelData.Types;
Sections = ModelData.Sections;
Materials = ModelData.Materials;
Constraints = ModelData.Constraints;

%% Usar funciones de StaBIL para obtener DOFs iniciales
fprintf('Obteniendo DOFs iniciales...\n');
DOF_initial = getdof(Elements, Types);

%% Eliminar DOFs para análisis 2D (elementos BEAM)
% Para BEAM en 2D: mantener solo DOFs 01, 02, 06 (Ux, Uy, Rz)
% Eliminar DOFs 03, 04, 05 (Uz, Rx, Ry) de todos los nodos
seldof_remove = [0.03; 0.04; 0.05];

fprintf('Eliminando DOFs para análisis 2D:\n');
fprintf('  - DOF 03 (Uz): desplazamiento en Z\n');
fprintf('  - DOF 04 (Rx): rotación en X\n');
fprintf('  - DOF 05 (Ry): rotación en Y\n');
fprintf('Manteniendo DOFs 01, 02, 06 (Ux, Uy, Rz)\n');

DOF = removedof(DOF_initial, seldof_remove);
fprintf('DOFs después de eliminación: %d (de %d originales)\n', length(DOF), length(DOF_initial));

%% Ensamblar matrices K y M con StaBIL
fprintf('Ensamblando matrices K y M...\n');
[K, M] = asmkm(Nodes, Elements, Types, Sections, Materials, DOF);

%% Aplicar restricciones multicuerpo si existen
if ~isempty(Constraints)
    fprintf('Aplicando %d restricciones multicuerpo...\n', size(Constraints, 1));
    % Crear vector de cargas vacío para usar addconstr
    P_dummy = zeros(length(DOF), 1);
    [K, ~, M] = addconstr(Constraints, DOF, K, P_dummy, M);
end

%% Organizar datos de salida
MatrixData = struct();
MatrixData.K = K;
MatrixData.M = M;
MatrixData.DOF = DOF;
MatrixData.DOF_initial = DOF_initial;
MatrixData.DOF_removed = seldof_remove;

%% Propiedades de las matrices
MatrixData.Properties.Size = size(K, 1);
MatrixData.Properties.DOFs_Initial = length(DOF_initial);
MatrixData.Properties.DOFs_After_Removal = length(DOF);
MatrixData.Properties.DOFs_Active = size(K, 1);
MatrixData.Properties.Symmetric_K = issymmetric(K);
MatrixData.Properties.Symmetric_M = issymmetric(M);

% Verificar definición positiva (solo para matrices no muy grandes)
if size(K, 1) < 500
    try
        eigs_K = eig(K);
        eigs_M = eig(M);
        MatrixData.Properties.Positive_Definite_K = all(eigs_K > 1e-10);
        MatrixData.Properties.Positive_Definite_M = all(eigs_M > 1e-10);
        MatrixData.Properties.Min_Eigenvalue_K = min(eigs_K);
        MatrixData.Properties.Min_Eigenvalue_M = min(eigs_M);
    catch
        MatrixData.Properties.Positive_Definite_K = NaN;
        MatrixData.Properties.Positive_Definite_M = NaN;
        MatrixData.Properties.Min_Eigenvalue_K = NaN;
        MatrixData.Properties.Min_Eigenvalue_M = NaN;
    end
else
    MatrixData.Properties.Positive_Definite_K = NaN;
    MatrixData.Properties.Positive_Definite_M = NaN;
    MatrixData.Properties.Min_Eigenvalue_K = NaN;
    MatrixData.Properties.Min_Eigenvalue_M = NaN;
end

% Condicionamiento (solo para matrices no muy grandes)
if size(K, 1) < 1000
    MatrixData.Properties.Condition_K = cond(K);
    MatrixData.Properties.Condition_M = cond(M);
else
    MatrixData.Properties.Condition_K = NaN;
    MatrixData.Properties.Condition_M = NaN;
end

%% Información de restricciones aplicadas
if ~isempty(Constraints)
    MatrixData.Constraints_Applied = Constraints;
    MatrixData.DOFs_Removed_Constraints = MatrixData.Properties.DOFs_After_Removal - MatrixData.Properties.DOFs_Active;
else
    MatrixData.Constraints_Applied = [];
    MatrixData.DOFs_Removed_Constraints = 0;
end

%% Metadatos
MatrixData.Info.Assembly_Date = datestr(now);
MatrixData.Info.Analysis_Type = '2D';
MatrixData.Info.Element_Type = 'BEAM';
MatrixData.Info.StaBIL_Functions = {'getdof', 'removedof', 'asmkm', 'addconstr'};
MatrixData.Info.Model_Reference = ModelData.Info;

%% Resumen
fprintf('\n=== MATRICES ENSAMBLADAS ===\n');
fprintf('Análisis: 2D\n');
fprintf('Elementos: BEAM\n');
fprintf('Tamaño matrices: %dx%d\n', size(K, 1), size(K, 2));
fprintf('DOFs iniciales: %d\n', MatrixData.Properties.DOFs_Initial);
fprintf('DOFs después de eliminación 2D: %d\n', MatrixData.Properties.DOFs_After_Removal);
fprintf('DOFs activos (después de restricciones): %d\n', MatrixData.Properties.DOFs_Active);
if ~isempty(Constraints)
    fprintf('DOFs eliminados por restricciones: %d\n', MatrixData.DOFs_Removed_Constraints);
end
fprintf('K simétrica: %s\n', mat2str(MatrixData.Properties.Symmetric_K));
fprintf('M simétrica: %s\n', mat2str(MatrixData.Properties.Symmetric_M));
if ~isnan(MatrixData.Properties.Condition_K)
    fprintf('Condicionamiento K: %.2e\n', MatrixData.Properties.Condition_K);
    fprintf('Condicionamiento M: %.2e\n', MatrixData.Properties.Condition_M);
end
fprintf('============================\n\n');

% Advertencias sobre definición positiva
if ~isnan(MatrixData.Properties.Positive_Definite_K) && ~MatrixData.Properties.Positive_Definite_K
    warning('La matriz de rigidez no es definida positiva (min eigenvalue: %.2e)', ...
            MatrixData.Properties.Min_Eigenvalue_K);
end

if ~isnan(MatrixData.Properties.Positive_Definite_M) && ~MatrixData.Properties.Positive_Definite_M
    warning('La matriz de masa no es definida positiva (min eigenvalue: %.2e)', ...
            MatrixData.Properties.Min_Eigenvalue_M);
end

fprintf('Matrices ensambladas correctamente para análisis 2D.\n');

end
