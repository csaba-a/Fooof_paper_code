% load surface 
surface_path = [analysis_location,'/atlas/brain_surf'];
[lhpialvert,lhpialface]=read_surf([surface_path '/lh.pial']);
lhpialface=lhpialface+1;
[rhpialvert,rhpialface]=read_surf([surface_path '/rh.pial']);
rhpialface=rhpialface+1;
lhpialvert(:,2)=lhpialvert(:,2)-18;
lhpialvert(:,3)=lhpialvert(:,3)+20;
rhpialvert(:,2)=rhpialvert(:,2)-18;
rhpialvert(:,3)=rhpialvert(:,3)+20;

%load color scheme