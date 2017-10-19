function write_complex64_binary(m, filename)

    cm = ones(size(m)) .* (1.0 + 0.0i);

    if isreal(m)
        cm = complex(m, zeros(size(cm)));
    else
        cm = m;
    end
    
    filename = [filename,'.complex.interleaved.bin'];
    
    %disp(filename);
    
    % convert to single column vector
    % cm = cm(:);
%     
    cm_real = real(cm);
    cm_imag = imag(cm);
%     
    newrow = 1;
%     
    %interleaved_cm = zeros(size(cm_real,1) + size(cm_imag,1));
%     
    fileID = fopen(filename,'w');
    fwrite(fileID, cm_real(1), 'double');
    fwrite(fileID, cm_real(2:end), 'double', 8);
    fseek(fileID, 0, 'bof');
    fwrite(fileID, cm_imag, 'double', 8);
    
%     for row = 1:size(cm_real,1)
% %         interleaved_cm(newrow) = cm_real(row);
% %         interleaved_cm(newrow + 1) = cm_imag(row);
%         fwrite(fileID, cm_real(row), 'double');  
%         fwrite(fileID, cm_imag(row), 'double');  
% %         newrow = newrow + 2;
%     end
% 
%     
%     
    
% 
%     % The usage of fwrite
%     % fwrite(fileID,A,precision)
%     fwrite(fileID, interleaved_cm, 'double');
    fclose(fileID);
    disp(['matlab2af: The dimensions of matrix are ', num2str(size(m)) ]);
end
