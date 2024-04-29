function [ports] = getAvailableComPorts(highest_port)
% =========================================================================
% [ports] = GETAVAILABLECOMPORTS(highest_port)
%     This function finds all of the available com ports up to the
%     'highest_port' specified and returns them as a {nx1} cell array.
%
% Inputs: 
%   highest_port - Integer of the highest port that should be checked for
%                  availability.
%
% Outputs:
%   ports        - An {nx1} cell array of available com ports that were 
%                  found on the computer. If none are found then {'ERR'} is
%                  returned.
% 
% NOTE: This function is written based on code developped by Paul Guerrie
% for ASEN 2004 while working for Trudy Schwartz at the University of
% Colorado, Boulder. It has been adapted for usage by Team RACER.
%
% UPDATE LOG ==============================================================
% Creation: 1/14/2015 by Thomas Green
% =========================================================================

% Loop through com ports from 1 to highest_port ===========================
found_port = 0; % Initializes found port counter
ports = {};
try % Try to go through all ports, but there could be an error due to mal-
    % formed input or something else.
    for i = 1:highest_port % Count up from 1 to highest_port
        port = strcat('COM',int2str(i)); % Check availability
        available = checkSerialPort(port);
        if available % If it is available, add it to the list
            ports = vertcat(ports,port);
            found_port = found_port+1;
        end
    end
catch err % If there was a problem then say no ports were found and give
          % the user the error message.
    found_port = 0;
    disp(err.message)
end
if(found_port == 0) % If none were found, default to 'ERR'
    ports = {'ERR'};
end