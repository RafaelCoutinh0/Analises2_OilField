
function ajustar_eixos(fig,margem)

if nargin < 1 || isempty(fig)
    fig = gcf;
end

if nargin < 2 || isempty(margem)
    margem = 0.15;
end

ax = findall(fig,'Type','axes');

for k = 1:numel(ax)

    curvas = findall(ax(k),'Type','line');
    valores = [];

    for j = 1:numel(curvas)
        y = curvas(j).YData;
        valores = [valores; y(:)];
    end

    valores = valores(isfinite(valores));

    if isempty(valores)
        continue;
    end

    ymin = min(valores);
    ymax = max(valores);

    amplitude = ymax-ymin;

    if amplitude == 0
        amplitude = max(abs(ymax)*0.1,1);
    end

    ylim(ax(k),[ymin-margem*amplitude, ...
        ymax+margem*amplitude]);

end
end
