function iniExtLibs
%Initialize MTEX
try
    MTEXmenu;                                                              %Open m-tex menu
catch
    startup_mtex;                                                          %Startup m-tex
end
% Check for toolboxes
assert(logical(license('test','gads_toolbox')),'Global Optimization Toolbox not installed');
assert(logical(license('test','optimization_toolbox')),'Optimization Toolbox not installed');
end