
clear
load FILE
[elements,NNE]=size(elem);
[nodes,sp]=size(x_a);
df=sp;
sc=1;

x_a=x_0;

for cont=1:5:ste_p

    %Output file name
    filename=(['FRAC_4000_' num2str(cont) '.vtk']);
   

    %nr_of_elements=numel(x);
    fid = fopen(filename, 'w'); 

    %ASCII file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'vtk output\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid, ['POINTS ' num2str(nodes) ' float\n']);

    for i=1:nodes
        fprintf(fid,[num2str(x_a(i,1)+sc*d(i*df-1,cont)) ' ' num2str(x_a(i,2)+sc*d(i*df,cont)) ' 0\n']);
    end
    %grid	
    fprintf(fid, ['CELLS ' num2str(elements) ' ' num2str(elements*(NNE+1)) '\n']);
    for i=1:elements
        fprintf(fid,[num2str(NNE) ' ' num2str(elem(i,1)-1) ' ' num2str(elem(i,2)-1) ' ' num2str(elem(i,3)-1) '\n']);
    end
    fprintf(fid, ['CELL_TYPES ' num2str(elements) '\n']);
    for i=1:elements
        fprintf(fid,' 5 \n');
    end
    fprintf(fid, ['POINT_DATA ' num2str(nodes) '\n']);
    %d
    fprintf(fid, 'VECTORS d float \n');
    for i=1:nodes
        fprintf(fid,[num2str(d(i*df-1,cont)) ' ' num2str(d(i*df,cont)) ' 0\n']);
    end
    %V
    fprintf(fid, 'VECTORS v float \n');
    for i=1:nodes
        fprintf(fid,[num2str(v(i*df-1,cont)) ' ' num2str(v(i*df,cont)) ' 0\n']);
    end
    %A
    fprintf(fid, 'VECTORS a float \n');
    for i=1:nodes
        fprintf(fid,[num2str(a(i*df-1,cont)) ' ' num2str(a(i*df,cont)) ' 0\n']);
    end

    fprintf(fid, ['CELL_DATA ' num2str(elements) '\n']);
    %Pressure
    fprintf(fid,'SCALARS Pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ps(i,cont)) '\n']);
    end
%     %Strain Energy
     fprintf(fid,'SCALARS W float\n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:elements
         fprintf(fid,[num2str(W_elem(i,cont)) '\n']);
     end
     %Status
     fprintf(fid,'SCALARS Status float\n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:elements
         fprintf(fid,[num2str(Status(i,cont)) '\n']);
     end
    
    %E11
    fprintf(fid,'SCALARS E11 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-3,cont)) '\n']);
    end
    %E22
    fprintf(fid,'SCALARS E22 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-2,cont)) '\n']);
    end
    %E33
    fprintf(fid,'SCALARS E33 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-1,cont)) '\n']);
    end
    %G12
    fprintf(fid,'SCALARS G12 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4,cont)) '\n']);
    end
        
    %S11
    fprintf(fid,'SCALARS S11 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-3,cont)) '\n']);
    end
    %S22
    fprintf(fid,'SCALARS S22 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-2,cont)) '\n']);
    end
    %S33
    fprintf(fid,'SCALARS S33 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-1,cont)) '\n']);
    end
    %T12
    fprintf(fid,'SCALARS T12 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4,cont)) '\n']);
    end

    fclose(fid);
end
 
