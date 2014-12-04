clear; close all; clc;

vm_cell = load('Results-CellLocations.txt');
vm_GDNF = load('Results-GDNF.txt');
cSize = 100;
cT = 50;
cell_cells = cell(cT,1);
cell_GDNF = cell(cT,1);

for t = 0:cT-1
    t
    cell_cells{t+1} = vm_cell(1+t*cSize:(t+1)*cSize,:);
    cell_GDNF{t+1} = vm_GDNF(1+t*cSize:(t+1)*cSize,:);
end

for t = 1:t
    subplot(1,2,1),imagesc(cell_cells{t})
    m_GDNF = cell_GDNF{t};
    subplot(1,2,2),imagesc(m_GDNF)
    pause(0.01)
end
