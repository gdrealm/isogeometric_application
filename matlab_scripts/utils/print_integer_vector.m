function print_integer_vector(fid,V)
    fprintf(fid,'[');
    for i = 1:length(V)
        fprintf(fid,' %d',V(i));
        if i < length(V)
            fprintf(fid,',');
        end
    end
    fprintf(fid, ']\n');
end
