clc;clear;close all;

%####### ####### ####### ####### ####### ####### #######  
%####### images to meh save as 'xml' files can be read for FEnics
%####### ####### ####### ####### ####### ####### ####### 


%%  add necessary packages
addpath("Im2mesh_1_85/");
addpath("mesh2d/");

%% Design bianry images
[X1,Y1] = meshgrid(linspace(-1,1,200),linspace(-1,1,200));%% Images to mesh

% first images
ratio_b_a = 1.2;
a = 0.4;
b = ratio_b_a*a;
logic_geo= Cassini_oval(a,b,X1,Y1);
top1 = zeros(size(X1));
top1(logic_geo) = 1;

% second images
a2= 0.44;
ratio2_b2_a2 = 1.08; 
logic_geo2= Cassini_oval(a2,a2*ratio2_b2_a2,X1,Y1);
top2 = zeros(size(X1));
top2(logic_geo2) =1;

% combine
im = top1+top2;

figure; im = im+1; imagesc(im), axis equal  off, colormap("parula");







%% images to mesh
% parameters
hgrad = 1.5;                    % Mesh growth rate
tf_avoid_sharp_corner = false;  
tolerance = 1.;
hmax = 500;
 hmin = 0.01;                    % Target minimum mesh edge length
mesh_kind = 'delaunay';
grad_limit = +0.25;

select_phase_2 = [2 3];
[ vert,tria,tnum] = im2mesh( im, select_phase_2, tf_avoid_sharp_corner, tolerance, hmax, mesh_kind, grad_limit );
plotMeshes( vert, tria, tnum );

vertx_min = min(vert(:,1));
vertx_max = max(vert(:,1));
verty_min = min(vert(:,2));
verty_max = max(vert(:,2));
vert(:,1) = vert(:,1) -vertx_min;
vert(:,2) = vert(:,2) -verty_min;
xmlFileName = ['./','test','.xml'];
gnegate_xml_mesh(xmlFileName,vert,tria)
% generate_xml_mesh_with_materials(xmlFileName,vert,tria,tnum-1)
generate_materials_xml('./mesh.xml',tnum-1,length(tnum))

%% save as xml file



%% function to plot the geomatry
function [logicgeo] = Cassini_oval(a,b,X1,Y1)

%% 4 holes: central points (0.5,0.5); (0.5,-0.5); (-0.5,0.5); (-0.5,-0.5); 
r = 0.1;
% logic_holes = (((X1 - 0.5).^2 + (Y1 -0.5).^2)<r^2) +...
%               (((X1 - 0.5).^2 + (Y1 +0.5).^2)<r^2) +...
%               (((X1 + 0.5).^2 + (Y1 +0.5).^2)<r^2) +...
%               (((X1 + 0.5).^2 + (Y1 -0.5).^2)<r^2);
% logic_holes = logical(logic_holes);
% point (0.0)
logic_co = (X1.^2 +Y1.^2).^2 +  2* a^2 .*(Y1.^2 - X1.^2 ) + a^4 - b^4>0;
% point(1,1)
logic_co1 = ((X1-1).^2 +(Y1-1).^2).^2 +  2* a^2 .*((Y1-1).^2 -( X1-1).^2 ) + a^4 - b^4>0;
% point(-1,-1)
logic_co2 = ((X1+1).^2 +(Y1+1).^2).^2 +  2* a^2 .*((Y1+1).^2 -( X1+1).^2 ) + a^4 - b^4>0;
% point(1,-1)
logic_co3 = ((X1-1).^2 +(Y1+1).^2).^2 +  2* a^2 .*((Y1+1).^2 -( X1-1).^2 ) + a^4 - b^4>0;
% point(-1,1)
logic_co4 = ((X1+1).^2 +(Y1-1).^2).^2 +  2* a^2 .*((Y1-1).^2 -( X1+1).^2 ) + a^4 - b^4>0;
% point (0,-1) rota ang = 90;
logic_co5 = (X1.^2 +(Y1-1).^2).^2 +  2* a^2 .*(X1.^2 - (Y1-1).^2 ) + a^4 - b^4>0;
% point (0,1) rota ang = 90;
logic_co6 = (X1.^2 +(Y1+1).^2).^2 +  2* a^2 .*(X1.^2 - (Y1+1).^2 ) + a^4 - b^4>0;
% point (1,0) rota ang = 90;
logic_co7 =((X1-1).^2 +Y1.^2).^2 +  2* a^2 .*((X1-1).^2 - Y1.^2 ) + a^4 - b^4>0;
% point (-1,0) rota ang = 90;
logic_co8 =((X1+1).^2 +Y1.^2).^2 +  2* a^2 .*((X1+1).^2 - Y1.^2 ) + a^4 - b^4>0;


%% 
logicgeo=logic_co& logic_co1&logic_co2&logic_co3&logic_co4&logic_co5&logic_co6&logic_co7&logic_co8;
% logicgeo = logicgeo - logic_holes;
logicgeo = logical(logicgeo);
end


function gnegate_xml_mesh(xmlFileName,vertices,elements)
%% generate xml file

% Create the XML document object
docNode = com.mathworks.xml.XMLUtils.createDocument('dolfin');
rootNode = docNode.getDocumentElement();

% Create the mesh element and add it to the document
meshNode = docNode.createElement('mesh');
meshNode.setAttribute('celltype', 'triangle');
meshNode.setAttribute('dim', '2');
rootNode.appendChild(meshNode);

% Add the vertices to the mesh element
verticesNode = docNode.createElement('vertices');
verticesNode.setAttribute('size', num2str(size(vertices,1)));
meshNode.appendChild(verticesNode);
for i=1:size(vertices,1)
    vertexNode = docNode.createElement('vertex');
    vertexNode.setAttribute('index', num2str(i-1));
    vertexNode.setAttribute('x', num2str(vertices(i,1)));
    vertexNode.setAttribute('y', num2str(vertices(i,2)));
    vertexNode.setAttribute('z', '0.0');
    verticesNode.appendChild(vertexNode);
end

% Add the elements to the mesh element
elementsNode = docNode.createElement('cells');
elementsNode.setAttribute("size",num2str(size(elements)));

meshNode.appendChild(elementsNode);
for i=1:size(elements,1)
    elementNode = docNode.createElement('triangle');
    elementNode.setAttribute('index', num2str(i-1));
    elementNode.setAttribute('v0', num2str(elements(i,1)-1));
    elementNode.setAttribute('v1', num2str(elements(i,2)-1));
    elementNode.setAttribute('v2', num2str(elements(i,3)-1));
%     elementNode.setAttribute('mat', num2str(elemat(i)-1));
    elementsNode.appendChild(elementNode);
end

% Write the XML document to a file
xmlwrite(xmlFileName, docNode);

end

function generate_xml_mesh_with_materials(xmlFileName, vertices, elements, materials)
    %% Generate XML file with material IDs

    % Create the XML document object
    docNode = com.mathworks.xml.XMLUtils.createDocument('dolfin');
    rootNode = docNode.getDocumentElement();

    % Create the mesh element and add it to the document
    meshNode = docNode.createElement('mesh');
    meshNode.setAttribute('celltype', 'triangle');
    meshNode.setAttribute('dim', '2');
    rootNode.appendChild(meshNode);

    % Add the vertices to the mesh element
    verticesNode = docNode.createElement('vertices');
    verticesNode.setAttribute('size', num2str(size(vertices,1)));
    meshNode.appendChild(verticesNode);
    for i = 1:size(vertices, 1)
        vertexNode = docNode.createElement('vertex');
        vertexNode.setAttribute('index', num2str(i-1));
        vertexNode.setAttribute('x', num2str(vertices(i, 1)));
        vertexNode.setAttribute('y', num2str(vertices(i, 2)));
        vertexNode.setAttribute('z', '0.0');
        verticesNode.appendChild(vertexNode);
    end

    % Add the elements and material IDs to the mesh element
    elementsNode = docNode.createElement('cells');
    elementsNode.setAttribute('size', num2str(size(elements, 1)));
    meshNode.appendChild(elementsNode);
    for i = 1:size(elements, 1)
        elementNode = docNode.createElement('triangle');
        elementNode.setAttribute('index', num2str(i-1));
        elementNode.setAttribute('v0', num2str(elements(i, 1)-1));
        elementNode.setAttribute('v1', num2str(elements(i, 2)-1));
        elementNode.setAttribute('v2', num2str(elements(i, 3)-1));
        elementNode.setAttribute('material', num2str(materials(i))); % Add material ID
        elementsNode.appendChild(elementNode);
    end

    % Write the XML document to a file
    xmlwrite(xmlFileName, docNode);
end


function generate_materials_xml(xmlFileName, materials, numCells)
    %% Generate XML file for material IDs in FEniCS format

    % Create the XML document object
    docNode = com.mathworks.xml.XMLUtils.createDocument('dolfin');
    rootNode = docNode.getDocumentElement();

    % Create the mesh_function element and add it to the document
    meshFunctionNode = docNode.createElement('mesh_function');
    meshValueCollectionNode = docNode.createElement('mesh_value_collection');
    meshValueCollectionNode.setAttribute('type', 'uint');
    meshValueCollectionNode.setAttribute('dim', '2'); % Assuming 2D mesh
    meshValueCollectionNode.setAttribute('size', num2str(numCells));
    meshFunctionNode.appendChild(meshValueCollectionNode);
    rootNode.appendChild(meshFunctionNode);

    % Add material IDs to the mesh_value_collection
    for i = 1:numCells
        valueNode = docNode.createElement('value');
        valueNode.setAttribute('cell_index', num2str(i-1)); % Zero-based indexing for FEniCS
        valueNode.setAttribute('local_entity', '0');
        valueNode.setAttribute('value', num2str(materials(i)));
        meshValueCollectionNode.appendChild(valueNode);
    end

    % Write the XML document to a file
    xmlwrite(xmlFileName, docNode);
end
