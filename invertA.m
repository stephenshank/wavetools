classdef invertA
%INVERTA   Inverse of a matrix A.
%   AINV = INVERTA(A) computes and stores a sparse-LUPQ factorization of
%   the the sparse matrix A. Linear solves (i.e., vectors X such that
%   A*X = B for a given vector B) can be performed by calling the apply
%   method, i.e., X = AINV.APPLY(B).
%
%   AINV = INVERT(A,1) allows solves with the transpose of A to also be
%   performed, and stores transposes of the LUPQ factors. Linear solves
%   with A^T (i.e., vectors X such that A^T*X = B for a given vector B) can
%   be performed by calling the apply method, i.e., X = AINV.APPLYT(B).

properties
	transpose;		% boolean to indicate whether we store tranpose of A
	L; U; P; Q;		% quantities to factor A
	Lt; Ut; Pt; Qt;	% quantities to factor A'
end

methods
	function this=invertA(A,transpose)
		if nargin<2, transpose=0; end
		this.transpose = transpose;
		% Compute L, U, P, and Q so that P*A*Q = L*U
		[this.L,this.U,this.P,this.Q]=lu(A);
		% Store transposes of these (if user specifies so)
		if transpose
			this.Lt=this.L';
			this.Ut=this.U';
			this.Pt=this.P';
			this.Qt=this.Q';
		end
	end

	function y=apply(this,x)
		% Do a linear solve with A via, LUPQ factors
		y=this.Q*(this.U\(this.L\(this.P*x)));
	end

	function y=applyt(this,x)
		% Check whether transposed factors were stored
		if ~this.transpose
			error('Transposed factor were not stored. Use AINV = INVERTA(A,1).')
		end
		% Do a linear solve with A^T via transposes of LUPQ factors
		y=this.Pt*(this.Lt\(this.Ut\(this.Qt*x)));
	end
end

end