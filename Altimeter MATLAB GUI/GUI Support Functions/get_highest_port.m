function [port] = get_highest_port()
%GET_HIGHEST_PORT Gets the highest available port on the computer.
%   PORT = GET_HIGHEST_PORT() gets the highest port and stores it in the
%   variable 'PORT'.

found_port = 0;%Initializes variable
max_port = 50;
wb = waitbar(0,'Detecting available COM ports...');
for i = max_port:-1:1%Start high and count backwards
    port = strcat('COM',int2str(i));%Check each com port
    found_port = check_port(port);
    if(found_port == 1)%If one works, exit loop
        break;
    end
    try
        waitbar((max_port-i)/max_port,wb,'Detecting available COM ports...');
    catch
        wb = waitbar((max_port-i)/max_port,'Detecting available COM ports...');
    end
end
close(wb)
if(found_port == 0)%if none were found, default to 'ERR'
    port = 'ERR';
end