function [pe_dot] = pe_dot(pe, dt) % 18/03/2020 10:01

    pe_dot = (pe(2:end) - pe(1:end-1)) / dt;

end