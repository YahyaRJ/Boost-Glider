function [available] = check_port(port)
%CHECK_PORT Checks if port is available.
%   AVAILABLE = CHECK_PORT(PORT) checks to see if com port 'PORT' is
%   available.  Returns 1 if available, 0 if unavailable.

    available = 1;%Assume it is available
    s = serial(port);%Create serial port object
    try
        fopen(s);%Try to open serial port object for communications
    catch
        available = 0;%If there was an error, unavailable.
    end
    fclose(s);%Close object
    delete(s);%delete object
    clear s %clear variable for later use
end