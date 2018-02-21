function [COM,Niter]=ClusterMex(M,debug,self);
COM=0;
S = size(M);
N = S(1);
ddebug = 0;
ending = 0;

M = M + M';
if (self == 0)
  M((N+1).*[0:N-1]+1) = 0;
end
M2 = M;
M2((N+1).*[0:N-1]+1) = 0;

m = sum(sum(M));
Niter = 1;

if m==0 | N == 1
  fprintf('No more possible decomposition\n');
  ending = 1;
  COMTY = 0;
  return;
end

% Main loop
K = sum(M); % Sum of wieght incident to node i
SumTot = sum(M);
SumIn = diag(M); % Sum of weight inside community i
COM = 1:S(1); % Community of node i
Neighbor=cell(1,N);
for k=1:N
  Neighbor{k} = find(M2(k,:));
end


%coder.varsize('Cnew'); coder.varsize('Cnew_t');
Cnew=0; Cnew_t=0; %Ck=[];
sCost = 10;
gain = 1;
while (gain == 1)
  Cost = zeros(1,N);
  gain = 0;
  for i=1:N
    Ci = COM(i);
    NB = Neighbor{i};
    G = zeros(1,N); % Gain vector
    best_increase = -1;
    Cnew = Ci;
    COM(i) = -1;
    SumTot(Ci) = SumTot(Ci) - K(i);
    CNj1 = find(COM==Ci);
    SumIn(Ci) = SumIn(Ci) - 2*sum(M(i,CNj1)) - M(i,i);
    for j=1:length(NB)
      Cj = COM(NB(j));
      if (G(Cj) == 0)
        CNj = find(COM==Cj);
        Ki_in = 2*sum(M(i,CNj));
        G(Cj) = Ki_in/m - 2*K(i)*SumTot(Cj)/(m*m);
        if (ddebug)
          fprintf('Gaim for comm %d => %g\n',Cj-1,G(Cj));
        end
        if G(Cj) > best_increase;
          best_increase = G(Cj);
          Cnew_t = Cj;
        end
      end
    end
    if best_increase > 0
      Cnew = Cnew_t;
      if (debug)
        %fprintf('Move %d => %d\n',i-1,Cnew-1);
      end
      Cost(i) = best_increase;
    end
%      sel=COM==Cnew;
%      indices=[1:N];
%      Ck=indices(sel);
    Ck = find(COM==Cnew);
%     for ii=1:size(COM,2)
%         if COM(ii)==Cnew
%             Ck=ii;
%         end
%     end
    
    SumIn(Cnew) = SumIn(Cnew) + 2*sum(M(i,Ck));
    SumTot(Cnew) = SumTot(Cnew) + K(i);
    COM(i) = Cnew;
    if (Cnew ~= Ci)
      gain = 1;
    end
    
  end
  sCost = sum(Cost);
  [C2 S2] = reindex_com(COM);
  Nco = length(unique(COM));
  Nco2 = length(S2(S2>1));
  mod = compute_modularity(COM,M);
  if (debug)
    %fprintf('It %d - Mod=%f %d com (%d non isolated)\n',Niter,mod,Nco,Nco2);
  end
  Niter = Niter + 1;
end



end

function MOD = compute_modularity(C,Mat)

m = sum(sum(Mat));
MOD = 0;
COMu = unique(C);
for j=1:length(COMu)
    Cj = find(C==COMu(j));
    Ec = sum(sum(Mat(Cj,Cj),1),2);
    Et = sum(sum(Mat(Cj,:),1),2);
    if Et>0
        MOD = MOD + Ec/m-(Et/m)^2;
    end
end

end

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