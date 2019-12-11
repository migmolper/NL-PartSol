
clear
load 2-FILE_pastor_fbar
[elements,NNE]=size(elem);
[nodes,sp]=size(x_a);
df=sp;
sc=1;

x_a=x_0;

[ste_p,~]=size(tp);

for cont=1:10:ste_p

    %Output file name
    filename=(['VM_' num2str(cont) '.vtk']);
   

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
        fprintf(fid,[num2str(NNE) ' ' num2str(elem(i,1)-1) ' ' num2str(elem(i,2)-1) ' ' num2str(elem(i,3)-1) ' ' num2str(elem(i,4)-1) '\n']);
    end
    fprintf(fid, ['CELL_TYPES ' num2str(elements) '\n']);
    for i=1:elements
        fprintf(fid,' 9 \n');
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
    %Ep
    fprintf(fid, 'SCALARS gamma float \n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:nodes
        fprintf(fid,[num2str(Gamma_nds(i,cont)) '\n']);
    end
    %Aw
%     fprintf(fid, 'VECTORS a_w float \n');
%     for i=1:nodes
%         fprintf(fid,[num2str(v(i*df-1,cont)) ' ' num2str(a(i*df,cont)) ' 0\n']);
%     end
    
    fprintf(fid, ['CELL_DATA ' num2str(elements) '\n']);
    %Pressure
    fprintf(fid,'SCALARS Pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ps(i,cont)) '\n']);
    end
    %Sy
    fprintf(fid,'SCALARS Yield_str float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Sy_tot(i,cont)) '\n']);
    end
     %Plastic strain
     fprintf(fid,'SCALARS Ep float\n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:elements
         fprintf(fid,[num2str(Gamma_tot(i,cont)) '\n']);
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
    
    %E11_p
    fprintf(fid,'SCALARS E11_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-3,cont)) '\n']);
    end
    %E22_p
    fprintf(fid,'SCALARS E22_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-2,cont)) '\n']);
    end
    %E33_p
    fprintf(fid,'SCALARS E33_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-1,cont)) '\n']);
    end
    %G12_p
    fprintf(fid,'SCALARS G12_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4,cont)) '\n']);
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
 
