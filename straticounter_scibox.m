function straticounter_scibox( settings_path, output_path )
%STRATICOUNTER_CSCIENCE Interface between the scibox software and
%straticounter function
%   Inputs:
%       settings_path (string): path to the .mat file from which the settings should be loaded
%       output_path (string): path to which the output of the straticounter() funciton should be written
%   Outputs:
%       NONE: This function does not directly produce any output or return any values. The outputs are written by the straticounter() function in the output_path
    load(settings_path);
    straticounter(Model, output_path)

end
