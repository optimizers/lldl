classdef opLLDL < opSpot
%OPLLDL  Operator representing the preconditioner
%
%          L |D| L'
%
%        where L and D are incomplete LDL factors of the SQD matrix K.
%
%   opLLDL(K,p) creates an operator for multiplication by preconditioned
%   operator, where L and D are limited-memory factors with factor p.
%   Note that K is an explicit matrix.
%
%   See also ldl.
%
%   D. Orban, 2013.
%
%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Properties
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   properties( SetAccess = private )
      L             % "Inverse" lower triangular incomplete factor
      D             % "Inverse" absolute value of diagonal incomplete factor
      Dinv          % "Inverse" diagonal incomplete factor (for normest)
      nnz           % Number of nonzeros in L
      p             % Limited-memory factor
      shift         % The shift necessary to complete the factorization.
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - Public
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % opLDL. Constructor
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function op = opLLDL(K, p)
         if nargin ~= 2
            error('Invalid number of arguments.');
         end

         % Get size of input matrix
         m = size(K,1);
         n = size(K,2);
         if m ~= n
            error('Input matrix must be square.');
         end

         % Construct operator
         op = op@opSpot('LLDL', n, n);
         op.cflag         = false;
         op.p             = max(p,0);
         [L, D, op.shift] = lldl(sparse(tril(K,-1)), full(diag(K)), op.p);
         op.nnz           = nnz(L);
         op.L             = inv(opMatrix(L + speye(size(L))));
         op.D             = opDiag(1./abs(D));
         op.Dinv          = opDiag(1./D);           % For op.normest().
         op.sweepflag     = true;
      end % function opLLDL

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % transpose
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function opOut = transpose(op)
         opOut = op;
      end

      function n = normest(op)
         % This is typically not accurate.
         n = normest(op.L' * op.Dinv * op.L);
      end

   end % methods - public

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - protected
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods( Access = protected )

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % multiply
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function y = multiply(op, x, mode)
         y = op.L' * op.D * op.L * x;
      end % function multiply

   end % methods - protected

end % classdef
