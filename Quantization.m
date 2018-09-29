%funci�n que asocia un conjunto de datos a un conjunto de centroides y
%reemplaza el dato por el centroide de menor distancia euclidea

%datos es una matriz en la cual la fila es una observaci�n y las columnas
%son las variables.


function[datos_cuantizados]=Quantization(centroides,datos)


[Nd,~] = size(datos);         % N�mero de datos

Nc = size(centroides,1);      %N�mero de centroides

d = zeros(Nd,Nc);
for i=1:Nd
    for j = 1:Nc
        %for k = 1:dim
        d(i,j) = norm(centroides(j,:)-datos(i,:));
        %end
        %d(i,j) = sqrt(d(i,j));
    end
end

[~,datos_cuantizados]=min(d,[],2);