function [ totalmem, freemem ] = linuxMemory()
% Usage: No input arguments needed.
% Returns total and free memory in Kilobytes.

    % probe system memory information
    [~,meminfo] = system('cat /proc/meminfo');

    % get total memory
    tokens = regexpi(meminfo,'^MemTotal:\s*(\d+)\s', 'tokens');
    totalmem = str2double(tokens{1}{1});

    % get available memory
    tokens = regexpi(meminfo,'^*MemFree:\s*(\d+)\s','tokens');
    freemem = str2double(tokens{1}{1});

end