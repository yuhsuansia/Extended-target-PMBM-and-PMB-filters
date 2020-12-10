function [handle_extent] = plot_extent_iw(x,X,line_style, color, line_width)

alpha = 0:pi/100:2*pi;
rotated = X*[cos(alpha);sin(alpha)];
xpoints = rotated(1,:) + x(1);
ypoints = rotated(2,:) + x(2);

handle_extent = plot(xpoints,ypoints,'LineStyle',line_style,'color',color,'LineWidth',line_width);

end