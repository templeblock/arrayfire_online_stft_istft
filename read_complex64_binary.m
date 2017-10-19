function same_z = read_complex64_binary(filename, mat_size)

    %disp(mat_size);
    z = zeros(mat_size);
    
    f = fopen(filename);
    %disp(['read from ',filename]);
    %disp(size(z));
    same_real = fread(f, size(z(:)), 'double', 8);
    fseek(f, 8, 'bof');
    same_imag = fread(f, size(z(:)), 'double', 8);
    fclose(f);
    disp(['read ',filename, ' with size ', num2str(size(z(:)))]);
     
    same_z = complex(same_real, same_imag);
    same_z = reshape(same_z, mat_size);
    disp(['convert size ', num2str(size(z(:))), ' to ', num2str(mat_size)]);
    %disp(same_z);

    %same_z(:,1) % first col
    %same_z(1,:) % first row

end