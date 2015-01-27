classdef invertA
% INVERTA   Inverse of a matrix A.
%   You still need to write better documentation, implement with vectors
%   instead of permutation matrices, and implement for dense matrices, also
%   for symmetric matrices and symmetric positive definite matrices!

properties
    amISparse;      % boolean to determine sparsity of A
    tranpose;       % boolean to indicate whether we store tranpose of A
    L; U; P; Q;     % quantities to factor A
    Lt; Ut; Pt; Qt; % quantities to factor A'
end

methods
    function this=invertA(A,transpose)
        if nargin<2, transpose=0; end
        this.amISparse=issparse(A);
        if this.amISparse
            [this.L,this.U,this.P,this.Q]=lu(A);
        else
            [this.L,this.U,this.P,this.Q]=lu(A);
        end
        if transpose
            if this.amISparse
                this.Lt=this.L';
                this.Ut=this.U';
                this.Pt=this.P';
                this.Qt=this.Q';
            else
            end
        end
    end
    
    function y=apply(this,x)
        if this.amISparse
            y=this.Q*(this.U\(this.L\(this.P*x)));
        else
        end
    end
    
    function y=applyt(this,x)
        if this.amISparse
            y=this.Pt*(this.Lt\(this.Ut\(this.Qt*x)));
        else
        end
    end
end
end