% This function computes the dihedral angles of the
% protein structure
%
% Babak Alipanahi
% University of Waterloo
% May 10, 2010

function [phi psi] = ang_checker(X,Comp,plot_flag)

num_planes = size(Comp.planes,1);

num_residues = length(unique(Comp.residue));
phi = nan(num_residues,1);
psi = nan(num_residues,1);

psi_X = [X(:,1) X(:,Comp.planes(1,[1 2 4]))];
psi(1) = dicalc(psi_X);

for i = 1:num_planes-1
    if i < num_planes-1
        psi_atoms = [Comp.planes(i,[4 6]) Comp.planes(i+1,[2 4])];
        psi_X = X(:,psi_atoms);
        psi(i+1) = dicalc(psi_X);
    end
    
    phi_atoms = [Comp.planes(i,[2 4 6]) Comp.planes(i+1,2)];
    phi_X = X(:,phi_atoms);
    phi(i+1) = dicalc(phi_X);
    
    
end

psi = 180*psi/pi;
phi = 180*phi/pi;

if plot_flag
    regData = feval('ramachandranRegions');
    hp = zeros(1,numel(regData));
    for i = numel(regData):-1:1
        
        % print only contours
        %      hp(i) = patch(regData(i).Patch(1,:), regData(i).Patch(2,:),...
        %          regData(i).Color, 'EdgeColor', regData(i).Color,...
        %         'DisplayName', regData(i).Name, 'FaceColor', 'none');
        
        % print solid patches
        hp(i) = patch(regData(i).Patch(1,:), regData(i).Patch(2,:), regData(i).Color, ...
            'DisplayName', regData(i).Name, 'EdgeColor', 'none');
    end
    hold on
    hl = line([0 -180; 0 180], [-180 0; 180 0], 'Color', 'k', 'LineStyle', ':');
    set(hl, 'HitTest', 'off');
    h2 = line([-180 -180; -180 180], [-180 180; 180 180], 'Color', 'k', 'LineStyle', '-');
    set(h2,'HitTest', 'off');
    %scatter(phi,psi,'o')
    for i = 1:numel(phi)
        plot(phi(i),psi(i),'o','LineWidth',2,'MarkerSize',4);
        text(phi(i)+2,psi(i)+2,num2str(Comp.num_seq(i)),'FontSize',8);
    end
    axis square
    axis([-180 180 -180 180]);
    set(gca,'XTick',-180:60:180);
    set(gca,'YTick',-180:60:180);
    box off
    box on
    xlabel('\Phi (degrees)','FontSize',14)
    ylabel('\Psi (degrees)','FontSize',14)
    hold off;
end



