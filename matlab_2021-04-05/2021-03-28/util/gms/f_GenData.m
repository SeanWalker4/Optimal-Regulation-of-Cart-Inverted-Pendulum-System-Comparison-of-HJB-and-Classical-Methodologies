function [DATArz, DATAdz, DATAniz] = f_GenData(P_d, weights, n_xi)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% GENERATE PROBLEM DATA FOR GMS OPTIMIZATION
%
% Karan Puttannaiah
% Brent Wallace  
%
% 12/14/2020
%
% Extract data from problem setup for GMS optimization.
%
% ***** CALL SYNTAX:
%
% [DATArz, DATAdz, DATAniz] = f_GenData(P_d, weights, n_xi)
%
% ***** INPUTS:
%
%   P_d                 : Design plant (state space object).
%   weights             : Data structure containing weightings (constructed
%                           in gms_main.m).
%   n_xi                : Number of states being fed back in inner loop
%                           (integer). Zero if inner-outer loop structure
%                           not used.
%
% ***** OUTPUTS:
%
%   DATArz, DATAdz, DATAniz     : Each a data structure with the following
%                                   fields:
%       ObjVec          : Vector which points to the row indices
%                           corresponding to the weighted sensitivity
%                           portion of the optimization (at the respective
%                           loop breaking point).
%       ConVec          : Cell array of vector(s) which point to the row
%                           indices corresponding to the constraint
%                           functional(s) associated with the respective
%                           loop breaking point.
%       ConNam          : Cell array containing the type of constraint of
%                           the respective constraint functional (i.e., the
%                           .Fun entry of the constraint functional's
%                           struct array, see gms_main.m for further
%                           documumentation).
%       ConVal          : Cell array containing the value of the respective
%                           constraint functional (i.e., the .Val entry of
%                           the constraint functional's struct array, see
%                           gms_main.m for further documumentation).
%       ConNum          : Total number of constraint functionals associated
%                           with the respective loop breaking point.
%
% *************************************************************************
% *************************************************************************
% *************************************************************************


% *************************************************************************
% *************************************************************************
%
% INITIALIZATION
%
% *************************************************************************
% *************************************************************************

% ***********************
%
% PLANT DIMENSIONS
%

% Number of inputs, outputs to plant.
[n_op, n_u] = size(P_d);

% Number of measured signals in outer loop (relevant for inner-outer only).
n_yp = n_op - n_xi;

% ***********************
%
% UNPACK WEIGHTING FUNCTIONS
%
W1 = weights.W1; 
W2 = weights.W2; 
W3 = weights.W3; 
Wd1 = weights.Wd1; 
Wd2 = weights.Wd2; 
Wd3 = weights.Wd3; 
Wni1 = weights.Wni1; 
Wni2 = weights.Wni2; 
Wni3 = weights.Wni3; 

% ***********************
%
% UNPACK CONSTRAINT FUNCTIONALS
%
W1c = weights.W1c; 
W2c = weights.W2c; 
W3c = weights.W3c; 
Wd1c = weights.Wd1c; 
Wd2c = weights.Wd2c; 
Wd3c = weights.Wd3c;
Wni1c = weights.Wni1c; 
Wni2c = weights.Wni2c; 
Wni3c = weights.Wni3c; 


%%
% *************************************************************************
% *************************************************************************
%
% PROBLEM DATA: r -> z
%
% *************************************************************************
% *************************************************************************

% Initialize row dimension (r -> z).
nObj = 0;

% *************************************************************************
%
% WEIGHTING FUNCTIONS
%

% Weighting function W1 (S_e).
if weights.exist_W1
    [noutput, ninput] = size(W1);
    if noutput ~= ninput
        error('Error: W1 is not square')
    end
    if ninput ~= n_yp
        error('Error: Dimension mismatch in W1')
    end
    nObj = nObj + n_yp;
end

% Weighting function W2 (K*S_e).
if weights.exist_W2
    [noutput, ninput] = size(W2);
    if noutput ~= ninput
        error('Error: W2 is not square')
    end  
    if ninput ~= n_u
        error('Error: Dimension mismatch in W2')
    end
    nObj = nObj + n_u;
end

% Weighting function W3 (T_e).
if weights.exist_W3
    [noutput, ninput] = size(W3);
    if noutput ~= ninput
        error('Error: W3 is not square')
    end 
    if ninput ~= n_yp
        error('Error: Dimension mismatch in W3')
    end
    nObj = nObj + n_yp;
end

% Row pointer for weighted closed loop maps (r -> z).
DATArz.ObjVec = 1:nObj;

% Update total row dimension counter (r -> z).
TotalRows = nObj;

% *************************************************************************
%
% CONSTRAINT FUNCTIONALS
%

% Initialize number of constraint functionals in optimization (r -> z).
ConstraintCounter = 0;

% Constraint functional W1c (S_e).
[~, nCol]=size(W1c);
for i=1:nCol
    if ~isempty(W1c{i}.tfm)
        [noutput, ninput] = size(W1c{i}.tfm);
        if noutput ~= ninput
            error(['Error: W1c{' num2str(i) '} is not square'])
        end
        if ninput ~= n_yp
            error(['Error: Dimension mismatch in W1c{' num2str(i) '}'])
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATArz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_yp;
        DATArz.ConNam{ConstraintCounter} = W1c{i}.Fun;
        DATArz.ConVal{ConstraintCounter} = W1c{i}.Val;
        TotalRows = TotalRows + n_yp;
    end
end

% Constraint functional W2c (K*S_e).
[~, nCol]=size(W2c);
for i=1:nCol
    if ~isempty(W2c{i}.tfm)
        [noutput, ninput] = size(W2c{i}.tfm);
        if noutput ~= ninput
            error(['Error: W2c{' num2str(i) '} is not square'])
        end
        if ninput ~= n_u
            error(['Error: Dimension mismatch in W2c{' num2str(i) '}'])
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATArz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_u;
        DATArz.ConNam{ConstraintCounter} = W2c{i}.Fun;
        DATArz.ConVal{ConstraintCounter} = W2c{i}.Val;
        TotalRows = TotalRows + n_u;
    end
end

% Constraint functional W3c (T_e).
[~, nCol]=size(W3c);
for i=1:nCol
    if ~isempty(W3c{i}.tfm)
        [noutput, ninput] = size(W3c{i}.tfm);
        if noutput ~= ninput
            error(['Error: W3c{' num2str(i) '} is not square'])
        end
        if ninput ~= n_yp
            error(['Error: Dimension mismatch in W3c{' num2str(i) '}'])
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATArz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_yp;
        DATArz.ConNam{ConstraintCounter} = W3c{i}.Fun;
        DATArz.ConVal{ConstraintCounter} = W3c{i}.Val;
        TotalRows = TotalRows + n_yp;
    end
end

% Total number of constraints (r -> z).
DATArz.ConNum = ConstraintCounter;



%%
% *************************************************************************
% *************************************************************************
%
% PROBLEM DATA: d -> z
%
% *************************************************************************
% *************************************************************************

% Initialize row dimension (d -> z).
nObj = 0;

% *************************************************************************
%
% WEIGHTING FUNCTIONS
%

% Weighting function Wd1 (S_c).
if weights.exist_Wd1
    [noutput, ninput] = size(Wd1);
    if noutput ~= ninput
        error('Error: Wd1 is not square')
    end
    if ninput ~= n_u
        error('Error: Dimension mismatch in Wd1')
    end
    nObj = nObj + n_u;
end

% Weighting function Wd2 (P*S_c).
if weights.exist_Wd2
    [noutput, ninput] = size(Wd2);
    if noutput ~= ninput
        error('Error: Wd2 is not square')
    end  
    if ninput ~= n_yp
        error('Error: Dimension mismatch in Wd2')
    end
    nObj = nObj + n_yp;
end

% Weighting function Wd3 (T_c).
if weights.exist_Wd3
    [noutput, ninput] = size(Wd3);
    if noutput ~= ninput
        error('Error: Wd3 is not square')
    end 
    if ninput ~= n_u
        error('Error: Dimension mismatch in Wd3')
    end
    nObj = nObj + n_u;
end

% Row pointer for weighted closed loop maps (d -> z).
DATAdz.ObjVec = 1:nObj;

% Update total row dimension counter (d -> z).
TotalRows = nObj;

% *************************************************************************
%
% CONSTRAINT FUNCTIONALS
%

% Initialize number of constraint functionals in optimization (d -> z).
ConstraintCounter = 0;

% Constraint functional Wd1c (S_c).
[~, nCol]=size(Wd1c);
for i=1:nCol
    if ~isempty(Wd1c{i}.tfm)
        [noutput, ninput] = size(Wd1c{i}.tfm);
        if noutput ~= ninput
            error(['Error: Wd1c{' num2str(i) '} is not square'])
        end
        if ninput ~= n_u
            error(['Error: Dimension mismatch in Wd1c{' num2str(i) '}'])
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATAdz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_u;
        DATAdz.ConNam{ConstraintCounter} = Wd1c{i}.Fun;
        DATAdz.ConVal{ConstraintCounter} = Wd1c{i}.Val;
        TotalRows = TotalRows + n_u;
    end
end

% Constraint functional Wd2c (P*S_c).
[~, nCol]=size(Wd2c);
for i=1:nCol
    if ~isempty(Wd2c{i}.tfm)
        [noutput, ninput] = size(Wd2c{i}.tfm);
        if noutput ~= ninput
            error(['Error: Wd2c{' num2str(i) '} is not square'])
        end
        if ninput ~= n_yp
            error(['Error: Dimension mismatch in Wd2c{' num2str(i) '}'])
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATAdz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_yp;
        DATAdz.ConNam{ConstraintCounter} = Wd2c{i}.Fun;
        DATAdz.ConVal{ConstraintCounter} = Wd2c{i}.Val;
        TotalRows = TotalRows + n_yp;
    end
end

% Constraint functional Wd3c (T_c).
[~, nCol]=size(Wd3c);
for i=1:nCol
    if ~isempty(Wd3c{i}.tfm)
        [noutput, ninput] = size(Wd3c{i}.tfm);
        if noutput ~= ninput
            error(['Error: Wd3c{' num2str(i) '} is not square'])
        end
        if ninput ~= n_u
            error(['Error: Dimension mismatch in Wd3c{' num2str(i) '}'])
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATAdz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_u;
        DATAdz.ConNam{ConstraintCounter} = Wd3c{i}.Fun;
        DATAdz.ConVal{ConstraintCounter} = Wd3c{i}.Val;
        TotalRows = TotalRows + n_u;
    end
end

% Total number of constraints (d -> z).
DATAdz.ConNum = ConstraintCounter;


%%
% *************************************************************************
% *************************************************************************
%
% PROBLEM DATA: ni -> z
%
% *************************************************************************
% *************************************************************************

% Initialize row dimension (ni -> z).
nObj = 0;

% *************************************************************************
%
% WEIGHTING FUNCTIONS
%

% Weighting function Wni1 (S_ni).
if weights.exist_Wni1
    [noutput, ninput] = size(Wni1);
    if noutput ~= ninput
        error('Error: Wni1 is not square')
    end
    if ninput ~= n_xi
        error('Error: Dimension mismatch in Wni1')
    end
    nObj = nObj + n_xi;
end

% Weighting function Wni2 (K*S_ni).
if weights.exist_Wni2
    [noutput, ninput] = size(Wni2);
    if noutput ~= ninput
        error('Error: Wni2 is not square')
    end  
    if ninput ~= n_u
        error('Error: Dimension mismatch in Wni2')
    end
    nObj = nObj + n_u;
end

% Weighting function Wni3 (T_ni).
if weights.exist_Wni3
    [noutput, ninput] = size(Wni3);
    if noutput ~= ninput
        error('Error: Wni3 is not square')
    end 
    if ninput ~= n_yp
        error('Error: Dimension mismatch in Wni3')
    end
    nObj = nObj + n_yp;
end

% Row pointer for weighted closed loop maps (ni -> z).
DATAniz.ObjVec = 1:nObj;

% Update total row dimension counter (ni -> z).
TotalRows = nObj;

% *************************************************************************
%
% CONSTRAINT FUNCTIONALS
%

% Initialize number of constraint functionals in optimization (ni -> z).
ConstraintCounter = 0;

% Constraint functional Wni1c (S_ni).
[~, nCol]=size(Wni1c);
for i=1:nCol
    if ~isempty(Wni1c{i}.tfm)
        [noutput, ninput] = size(Wni1c{i}.tfm);
        if noutput ~= ninput
            error(['Error: Wni1c{' num2str(i) '} is not square'])
        end
        if ninput ~= n_xi
            error(['Error: Dimension mismatch in Wni1c{' num2str(i) '}'])
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATAniz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_xi;
        DATAniz.ConNam{ConstraintCounter} = Wni1c{i}.Fun;
        DATAniz.ConVal{ConstraintCounter} = Wni1c{i}.Val;
        TotalRows = TotalRows + n_xi;
    end
end

% Constraint functional Wni2c (K*S_ni).
[~, nCol]=size(Wni2c);
for i=1:nCol
    if ~isempty(Wni2c{i}.tfm)
        [noutput, ninput] = size(Wni2c{i}.tfm);
        if noutput ~= ninput
            error(['Error: Wni2c{' num2str(i) '} is not square'])
        end
        if ninput ~= n_u
            error(['Error: Dimension mismatch in Wni2c{' num2str(i) '}'])
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATAniz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_u;
        DATAniz.ConNam{ConstraintCounter} = Wni2c{i}.Fun;
        DATAniz.ConVal{ConstraintCounter} = Wni2c{i}.Val;
        TotalRows = TotalRows + n_u;
    end
end

% Constraint functional Wni3c (T_ni).
[~, nCol]=size(Wni3c);
for i=1:nCol
    if ~isempty(Wni3c{i}.tfm)
        [noutput, ninput] = size(Wni3c{i}.tfm);
        if noutput ~= ninput
            error(['Error: Wni3c{' num2str(i) '} is not square'])
        end
        if ninput ~= n_yp
            error(['Error: Dimension mismatch in Wni3c{' num2str(i) '}'])
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATAniz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_yp;
        DATAniz.ConNam{ConstraintCounter} = Wni3c{i}.Fun;
        DATAniz.ConVal{ConstraintCounter} = Wni3c{i}.Val;
        TotalRows = TotalRows + n_yp;
    end
end

% Total number of constraints (ni -> z).
DATAniz.ConNum = ConstraintCounter;



