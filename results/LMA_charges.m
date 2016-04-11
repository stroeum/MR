% LMA charges %
for ii=1:Nb.Layers
    rectangle('Position',[-Layers.data(ii,5),(z_gnd+Layers.data(ii,4)-Layers.data(ii,6)/2),2*Layers.data(ii,5),Layers.data(ii,6)],'Curvature',[0,0],'LineWidth',Layers.Line.Width,'LineStyle',Layers.Line.Style,'EdgeColor',Layers.Edge.Color(ii,:));
end