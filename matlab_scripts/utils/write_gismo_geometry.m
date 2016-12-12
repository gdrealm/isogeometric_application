%write the nurbs data to gismo file
function write_gismo_geometry(nurbs_list,fn)

%% write the header
fid = fopen(fn,'wt');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<xml>\n');

for i = 1:length(nurbs_list)
    %% prepare data for each geometry
    nurbs = nurbs_list{i};
    dim = length(nurbs.knots);
    geom_type = ['TensorNurbs',num2str(dim)];
    basis_type = ['TensorNurbsBasis',num2str(dim)];
    sub_basis_type = ['TensorBSplineBasis',num2str(dim)];
    sub2_basis_type = 'BSplineBasis';

    %% write data
    fprintf(fid,'\t<Geometry type="%s" id="%d">\n',geom_type,i-1);
    fprintf(fid,'\t\t<Basis type="%s">\n',basis_type);
    fprintf(fid,'\t\t\t<Basis type="%s">\n',sub_basis_type);
    for j = 1:dim
        fprintf(fid,'\t\t\t\t<Basis type="%s" index="%d">\n',sub2_basis_type,j-1);
        fprintf(fid,'\t\t\t\t\t<KnotVector degree="%d">\n',nurbs.order(j)-1);
        fprintf(fid,'\t\t\t\t\t\t');
        fprintf(fid,' %.16f',nurbs.knots{j});
        fprintf(fid,'\n\t\t\t\t\t</KnotVector>\n');
        fprintf(fid,'\t\t\t\t</Basis>\n');
    end
    fprintf(fid,'\t\t\t</Basis>\n');
    fprintf(fid,'\t\t\t<weights>\n');
    if dim == 2
        for i1 = 1:size(nurbs.coefs,3)
            for j1 = 1:size(nurbs.coefs,2)
                w = nurbs.coefs(4,j1,i1);
                fprintf(fid,'\t\t\t\t%.16f\n',w);
            end
        end
    elseif dim == 3
        for i1 = 1:size(nurbs.coefs,4)
            for j1 = 1:size(nurbs.coefs,3)
                for k1 = 1:size(nurbs.coefs,2)
                    w = nurbs.coefs(4,k1,j1,i1);
                    fprintf(fid,'\t\t\t\t%.16f\n',w);
                end
            end
        end
    end
    fprintf(fid,'\t\t\t</weights>\n');
    fprintf(fid,'\t\t</Basis>\n');
    fprintf(fid,'\t\t<coefs geoDim="3">\n');
    if dim == 2
        for i1 = 1:size(nurbs.coefs,3)
            for j1 = 1:size(nurbs.coefs,2)
                w = nurbs.coefs(4,j1,i1);
                x = nurbs.coefs(1,j1,i1) / w;
                y = nurbs.coefs(2,j1,i1) / w;
                z = 0.0;
                fprintf(fid,'\t\t\t%.16f %.16f %.16f\n',x,y,z);
            end
        end
    elseif dim == 3
        for i1 = 1:size(nurbs.coefs,4)
            for j1 = 1:size(nurbs.coefs,3)
                for k1 = 1:size(nurbs.coefs,2)
                    w = nurbs.coefs(4,k1,j1,i1);
                    x = nurbs.coefs(1,k1,j1,i1) / w;
                    y = nurbs.coefs(2,k1,j1,i1) / w;
                    z = nurbs.coefs(3,k1,j1,i1) / w;
                    fprintf(fid,'\t\t\t%.16f %.16f %.16f\n',x,y,z);
                end
            end
        end
    end
    fprintf(fid,'\t\t</coefs>\n');
    fprintf(fid,'\t</Geometry>\n');
end

%% finalize
fprintf(fid,'</xml>\n');
fclose(fid);
