classdef opLLDL < opSpot
%OPLLDL  Operator representing the preconditioner
%
%          L D L'   or   L |D| L'
%
%        where L and D are incomplete LDL factors of the SQD matrix K.
%
%   opLLDL(K,p) creates an operator for multiplication by preconditioned
%   operator, where L and D are limited-memory factors with factor p.
%
%   opLLDL(K,p,shift) is equivalent to opLLDL(K + shift * E, p), where
%   E is diagonal with +1s and -1s on the diagonal, with the same
%   definiteness pattern as K. The default value of shift is zero.
%
%   Options may be set by passing a structure as last argument:
%   opLLDL(K,p,shift,opts) where
%     opts.indef = whether the operator is allowed to be indefinite (true)
%                  or not (false). If set to `true`, the diagonal D
%                  resulting from the factorization is used as is when
%                  applying the preconditioner. If set to `false`, then
%                  |D| is used.
%     opts.fill_factor = whether the argument `p` should be interpreted
%                        as an additive (nnz + p*n) or multiplicative
%                        (p * nnz) memory specification. In the latter
%                        case (`opts.fill_factor = true`), an equivalent
%                        additive value of p is computed to perform the
%                        factorization.
%
%   Note that K is an explicit matrix. The public properties are:
%
%      L         : "Inverse" lower triangular incomplete factor
%      D         : "Inverse" absolute value of diagonal incomplete factor
%      Dinv      : "Inverse" diagonal incomplete factor (for normest)
%      nnz       : Number of nonzeros in L
%      p         : Limited-memory factor
%      shift     : The shift necessary to complete the factorization
%      growth    : Growth factor max(max |L|, max |D|) / max |K|
%      minpivot  : Smallest pivot in absolute value.
%
%   Public methods:
%
%      normest : estimate the norm of inv(L*D*L') (note that D appears,
%                not its absolute value). This may be useful to
%                estimate the condition number of the input operator,
%                though this estimate is rarely accurate.
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
      D             % Currently used variant of the diagonal (signed or unsigned)
      Dabs          % "Inverse" absolute value of diagonal incomplete factor
      Dinv          % "Inverse" diagonal incomplete factor (for normest)
      nnz           % Number of nonzeros in L
      p             % Limited-memory factor
      shift         % The shift necessary to complete the factorization
      growth        % Growth factor max(max |L|, max |D|) / max |K|
      minpivot      % Smallest pivot in absolute value
      indef         % Whether D should be replaced with |D|.
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - Public
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % opLDL. Constructor
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function op = opLLDL(K, p, shift, opts)
         if nargin < 2 || nargin > 4
            error('Invalid number of arguments.');
         end
         if ~exist('shift', 'var') || isempty(shift)
            shift = 0;
         end
         indef = false;
         fill_factor = false;
         if exist('opts', 'var')
            if isfield(opts, 'indef')
               if islogical(opts.indef)
                  indef = opts.indef;
               end
            end
            if isfield(opts, 'fill_factor')
               if islogical(opts.fill_factor)
                  fill_factor = opts.fill_factor;
               end
            end
         end

         % Get size of input matrix
         m = size(K,1);
         n = size(K,2);
         if m ~= n
            error('Input matrix must be square.');
         end

         % Construct operator
         op = op@opSpot('LLDL', n, n);
         op.cflag   = false;
         op.indef   = indef;

         % If opts.fill_factor, p specifies the fill factor instead of
         % the absolute additional memory.
         % The total number of nonzers allowed in the factors is as in iLDL:
         % fill/n * nnz(K). This gives
         % fill/n * nnz(K) = nnz(K) + pn, i.e., p = (fill/n - 1) * nnz(K) / n
         T = tril(K, -1);
         if fill_factor
           nnzK = nnz(T);
           op.p = max(round((p * (n + 2 * nnzK) - nnzK) / n), 0);
         else
           op.p = max(p, 0);
         end
         [L, D, op.shift] = lldl(sparse(T), full(diag(K)), op.p, shift);
         absD             = abs(D);
         op.nnz           = nnz(L);

         % Growth factor: max(|L * sqrt(|D|)) / max(|K|).
         L                = L + speye(n);
         LD               = L * spdiags(sqrt(abs(D)), 0, n, n);
         op.growth        = full(max(max(abs(LD))) / max(max(abs(tril(K)))));

         op.minpivot      = min(absD);
         op.L             = inv(opMatrix(L));
         op.Dabs          = opDiag(1./absD);
         op.Dinv          = opDiag(1./D);
         if indef
           op.D           = op.Dinv;
         else
           op.D           = op.Dabs;
         end
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
