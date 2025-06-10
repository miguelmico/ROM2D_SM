# ROM2D_SM

FEM 2D to Simscape ROM Converter
A MATLAB toolbox for converting 2D beam finite element models into reduced-order matrices compatible with Simscape Multibody's Reduced Order Flexible Solid block.
Overview
This toolbox automates the process of creating reduced-order flexible body models from 2D beam structures for use in Simscape Multibody simulations. It handles the complete workflow from FEM model definition to Craig-Bampton reduction.
Features

2D Beam Analysis: Optimized for planar beam structures (Ux, Uy, Rz DOFs)
Automatic Revolute Joints: Simplified interface for defining joint connections
Craig-Bampton Reduction: Static condensation + normal modes for efficiency
StaBIL Integration: Uses StaBIL functions for matrix assembly
Input Validation: Comprehensive checks for model consistency
Simscape Ready: Outputs compatible with Reduced Order Flexible Solid block

Requirements

MATLAB R2019b+
StaBIL Toolbox (getdof, removedof, asmkm, addconstr)

Quick Start
matlab% 1. Define 2D beam model
Nodes = [1 0 0 0; 2 1 0 0; 3 2 0 0];  % [ID X Y Z]
Elements = [1 1 1 1 1 2 10];           % [ID Type Sec Mat n1 n2 ref]
Types = {1, 'beam'};
Sections = [1 0.01 Inf Inf 0 0 1e-6 0.1 0.1 0.1 0.1];
Materials = [1 200e9 0.3 7850];
Constraints = {1, 2, 'revolute'};      % Revolute joint at node 2

% 2. Process model
ModelData = createFEMModel(Nodes, Elements, Types, Sections, Materials, Constraints);
MatrixData = assembleFEMMatrices(ModelData);

% 3. Reduce for Simscape
master_nodes = [1; 3];  % Interface nodes
ReducedData = craigBamptonReduction(MatrixData, ModelData, master_nodes, 10);

% 4. Export matrices
K_reduced = ReducedData.K_reduced;
M_reduced = ReducedData.M_reduced;
File Structure

createFEMModel.m - Model creation and validation
assembleFEMMatrices.m - 2D matrix assembly using StaBIL
craigBamptonReduction.m - Craig-Bampton reduction implementation
example_usage.m - Complete usage example

Output
The toolbox generates matrices suitable for Simscape Multibody's Reduced Order Flexible Solid block:

K_reduced: Reduced stiffness matrix
M_reduced: Reduced mass matrix
Interface node coordinates for connection points
Modal frequencies and shapes for verification
