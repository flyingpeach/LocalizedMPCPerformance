function connected = check_connected(sys)
% Check if A matrix is fully connected

if ~isempty(sys.AComm)
    commsAdj = sys.AComm; % Use the specified communication structure
else
    commsAdj = sys.A ~= 0; % Use comm structure implied by A
end

supp      = commsAdj^(sys.Nx) > 0;
connected = isempty(find(supp==0, 1));

end