function [ ] = run( scan_directory , feature_directory, label_directory, output_file, min_line)
%%
[components_preprocess, pins , netlist_mtx] = segment3(scan_directory, min_line);
[decision_2pin, decision_3pin] = KNN_sorted(fullfile(feature_directory,'feature_mtx.mat') ,fullfile(label_directory,'label_mtx.mat'), components_preprocess, pins);
%%
final_netlist = netlist_mtx;
index_2pin = 1;
index_3pin = 1;
for ii = 1:size(netlist_mtx,1)
    if (pins(ii) == 1)
        final_netlist(ii,1) = 4;
    elseif (pins(ii) == 2)
        final_netlist(ii,1) = decision_2pin(index_2pin);
        index_2pin = index_2pin + 1;
    elseif (pins(ii) == 3)
        final_netlist(ii,1) = 7;%decision_3pin(index_3pin);
        index_3pin = index_3pin +1;
    end       
end

for jj = 1:size(final_netlist,1)
    if(final_netlist(jj,1) == 1)
        final_netlist(jj,2) = 1;
        final_netlist(jj,5) = 1e3;
    elseif(final_netlist(jj,1) == 2)
        final_netlist(jj,2) = 15.9e-9;
    elseif(final_netlist(jj,1) == 3)
        final_netlist(jj,2) = 5;
    elseif(final_netlist(jj,1) == 10)
        final_netlist(jj,2) = 10e3;
    elseif(final_netlist(jj,1) == 7)
    end
end

output = netlist(final_netlist, output_file);
end