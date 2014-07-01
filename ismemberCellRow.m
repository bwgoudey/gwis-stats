function [tf,I] = ismemberCellRow(set1, set2)
    set_all = union(set1,set2) ;
    [~,ai] = ismember(set1, set_all) ;
    [~,bi] = ismember(set2,set_all) ;
    [tf,I] = ismember(ai,bi,'rows');