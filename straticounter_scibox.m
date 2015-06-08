function [ output_args ] = straticounter_scibox( settings_path, output_path )
%STRATICOUNTER_CSCIENCE Interface between the CScience software and
%straticounter fucntion
%   Detailed explanation goes here
    load(settings_path);
    straticounter(Model, output_path)

end

