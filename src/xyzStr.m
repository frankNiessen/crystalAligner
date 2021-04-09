function str = xyzStr(vec)
%function str = xyzStr(vec)
%vec: vector3d [xvector | yvector | zvector]
    if vec.x == 1 && vec.y == 0 && vec.z == 0
       str = 'X-axis';                                                     %Assign axis string
    elseif vec.x == 0 && vec.y == 1 && vec.z == 0
       str = 'Y-axis';                                                     %Assign axis string
    elseif vec.x == 0 && vec.y == 0 && vec.z == 1
       str = 'Z-axis';                                                     %Assign axis string
    else
       str = num2str(vec.xyz);
       str = str(str ~= ' ');
    end
end