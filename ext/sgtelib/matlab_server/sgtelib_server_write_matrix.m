function sgtelib_server_write_matrix(M,name,file)

% fid = fopen(file,'w');
% 
% newline = char(10);
% fwrite(fid,[name '=[' newline]);
% 
% for i=1:size(M,1)
%     Mi = num2str(M(i,:),12);
%     s = ['    ' Mi ' ; ' newline];
%     fwrite(fid,s);
% end
% fwrite(fid,'];');
% fclose(fid);

    
fid = fopen(file,'w');
fprintf(fid,'%s=[\n',name);

for i=1:size(M,1)
    fprintf(fid,'%12.12f  ',M(i,:));
    fprintf(fid,';\n');
end

fwrite(fid,'];');
fclose(fid);