function      st = cutStruc2Inds( st, indCutTo);
% cuts structure to indices
% CALL:   st = cutStruc2Inds( st, indCutTo)
%  INPUTS:
%    st -  a structure
%    indCutTo - list of indices
%  OUTPUTS: 
%     stCut
%  SEE ALSO: cutstruc 
%  WRITTEN & MODIFIED: 01/04/06 AK
%  � National ICT Australia 2006 (All Rights Reserved)
%  23/3/11 bwgoudey - added test for empty cases
%  24/3/11 bwgoudey - altered numeric case so that class of fields stays
%                     the same
%  6/9/11 bwgoudey - added case so we skip over chars.
%  18/10/12 AK - fixed cut for char arrays.
%  30/11/12 bwgoudey - skip elements with only one element
%  07/02/12 bwgoudey - case for cells. WAs causing errors. 



fields = fieldnames(st); 
for iF = 1 : length(fields)
    field = fields{iF};
    %HACK: This is only to remove some fields that accidentally got added
    if(strcmp(field,'process')) 
        continue
    end
    v = getfield( st, field);
    
    %bwgoudey - skip if empty or only one element
    if(isempty(v) || size(v, 1)==1)
        continue;
    end
    % st = setfield(st, field, v(indCutTo~=0));

    I = indCutTo ~=0;
    
    if ischar(v) || isstruct(v) || length(v)==1 || isa(v, 'cell')
      %  vv=v;
         vv = v(indCutTo(I),:);   % AK fixed this on 18 Oct, 2012
    %elseif isa(v, 'cell')
    %    vv = cell(size(I)); %Ben 
    %    vv(find(I),:) = v(indCutTo(I), :);
    %    %vv = v(indCutTo(I),:);
    else  %numeric
        %pre-allocation sometimes results in errors
        %vv = v(1:length(I));%NaN* zeros(size(I));        
        vv = v(indCutTo(I),:);
    end
    st = setfield(st, field, vv);
end
