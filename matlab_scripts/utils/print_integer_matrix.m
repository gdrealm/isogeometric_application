function print_integer_matrix(fid,M)
    fprintf(fid,'[');
    for i = 1:size(M,1)
        for j = 1:size(M,2)
            fprintf(fid,' %d,',M(i,j));
        end
    end
    fprintf(fid, ']\n');
end

