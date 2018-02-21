% Iplementation : Antoine Scherrer
% antoine.scherrer@ens-lyon.fr
% Apply clustering after :
% "Fast unfolding of community hierarchies in large networks"
% Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte,
% Etienne Lefebvre
% http://arxiv.org/abs/0803.0476
%
% NON ORIENTED VERSION USING SYMETRIC MATRIX A = M + M^t INSTEAD OF 
% (POSSIBLY NON SYMETRIC) INPUT MATRIX M
%
% FULL MATLAB VERSION (SLOWER)
%
% Inputs : 
% M : weight matrix (the matrix is symetrized with
% the sum of weights in both directions)
% s : 1 = Recursive computation
%   : 0 = Just one level computation
% self : 1 = Use self weights
%        0 = Do not use self weights
% debug   : 1 = outputs some debug messages
% verbose : 1 = outputs some messages
%
% Output :
% COMTY, structure with the following information
% for each level i :
%   COMTY.COM{i} : vector of community IDs (sorted by community sizes)
%   COMTY.SIZE{i} : vector of community sizes
%   COMTY.MOD(i) : modularity of clustering
%   COMTY.Niter(i) : Number of iteration before convergence
%
function [COMTY ending] = cluster_jl(M,s,self,debug,verbose)

if nargin < 1
  error('not enough argument');
end

if nargin < 2
  s = 1;
end

if nargin < 3
  self = 1;
end

if nargin < 4
  debug = 0;
end

if nargin < 5
  verbose = 0;
end

S = size(M);
N = S(1);

ddebug = 0;
ending = 0;

tic
[COM Niter]=ClusterMex(M,logical(debug),logical(self));
toc

[COM COMSIZE] = reindex_com(COM);
COMTY.COM{1} = COM;
COMTY.SIZE{1} = COMSIZE;
COMTY.MOD(1) = compute_modularity(COM,M);
COMTY.Niter(1) = Niter;

% Perform part 2
if (s == 1)
  
  Mnew = M;
  Mold = Mnew;
  COMcur = COM;
  COMfull = COM;
  k = 2;

  if (debug)
    Nco2 = length(COMSIZE(COMSIZE>1));
    fprintf('Pass number 1 - %d com (%d iterations)\n',Nco2,Niter);
  end
  while 1
    Mold = Mnew;
    S2 = size(Mold);
    Nnode = S2(1);
    
    COMu = unique(COMcur);
    Ncom = length(COMu);
    ind_com = zeros(Ncom,Nnode);
    ind_com_full = zeros(Ncom,N);

    for p=1:Ncom
      ind = find(COMcur==p);
      ind_com(p,1:length(ind)) = ind;
    end
    for p=1:Ncom
      ind = find(COMfull==p);
      ind_com_full(p,1:length(ind)) = ind;
    end
    
    Mnew = zeros(Ncom,Ncom);
    for m=1:Ncom
      for n=m:Ncom
        ind1 = ind_com(m,:);
        ind2 = ind_com(n,:);
        Mnew(m,n) = sum(sum(Mold(ind1(ind1>0),ind2(ind2>0))));
        Mnew(n,m) = sum(sum(Mold(ind1(ind1>0),ind2(ind2>0))));
      end
    end
    
    [COMt e] = cluster_jl(Mnew,0,self,debug,verbose);
    if (e ~= 1)
      COMfull = zeros(1,N);
      COMcur = COMt.COM{1};
      for p=1:Ncom
        ind1 = ind_com_full(p,:);
        COMfull(ind1(ind1>0)) = COMcur(p);
      end
      [COMfull COMSIZE] = reindex_com(COMfull);
      COMTY.COM{k} = COMfull;
      COMTY.SIZE{k} = COMSIZE;
      COMTY.MOD(k) = compute_modularity(COMfull,M);
      COMTY.Niter(k) = COMt.Niter;
      Nco2 = length(COMSIZE(COMSIZE>1));
      if (debug)
        fprintf('Pass number %d - %d com\n',k,Nco2);
      end
      Ind = (COMfull == COMTY.COM{k-1});
      if (sum(Ind) == length(Ind))
        if (debug)
          fprintf('Identical segmentation => End\n');
        end
        return;
      end
    else
      if (debug)
        fprintf('Empty matrix => End\n');
      end
      return;
    end
    k = k + 1;
  end
end

end

% Re-index community IDs
function [C Ss] = reindex_com(COMold)

C = zeros(1,length(COMold));
COMu = unique(COMold);
S = zeros(1,length(COMu));
for l=1:length(COMu)
    S(l) = length(COMold(COMold==COMu(l)));
end
[Ss INDs] = sort(S,'descend');

for l=1:length(COMu)
    C(COMold==COMu(INDs(l))) = l;
end

end

%Compute modulartiy
function MOD = compute_modularity(C,Mat)

m = sum(sum(Mat));
MOD = 0;
COMu = unique(C);
for j=1:length(COMu)
    Cj = find(C==COMu(j));
    Ec = sum(sum(Mat(Cj,Cj)));
    Et = sum(sum(Mat(Cj,:)));
    if Et>0
        MOD = MOD + Ec/m-(Et/m)^2;
    end
end

end

