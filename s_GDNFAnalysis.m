clear; close all; clc;

vm_cell = load('Results-CellLocations.txt');
vm_GDNF = load('Results-GDNF.txt');
cSize = 100;
cT = 100;
cell_cells = cell(cT,1);
cell_GDNF = cell(cT,1);

for t = 0:cT-1
    t
    cell_cells{t+1} = vm_cell(1+t*cSize:(t+1)*cSize,:);
    cell_GDNF{t+1} = vm_GDNF(1+t*cSize:(t+1)*cSize,:);
end

for t = 1:t
    subplot(1,2,1),imagesc(cell_cells{t})
    subplot(1,2,2),imagesc(cell_GDNF{t})
    pause(0.1)
end
