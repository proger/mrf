% MATLAB example of usage the MEX function maxsum.cpp
% (c) Tomas Werner, Oct 2005, Center for Machine Perception, Prague

% Create input image
I = double(imread('russian_p.png'))/255;
I = I + randn(size(I))/3;

% Define the max-sum problem
q = -3*[0 0 0 0 0 1 1];
L = [0 0 1 1 1 1 1];
G = [1 1 1 1 1 1 1 1 1 1 1 1 1   2 2 2 2 2 2 2 2 2 2 2
     1 1 1 1 2 2 3 4 4 5 5 6 7   1 1 1 2 2 3 3 4 5 6 7
     1 3 6 7 1 2 2 1 4 1 5 4 5   1 4 6 2 5 3 7 2 1 3 1
     0 0 0 0 0 0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 0 0 0 0];

K = max(max(G(2:3,:)));
for k = 1:K
  Q(k,:,:) = q(k) - (L(k)-I).^2;
end


% Solve the max-sum problem

cmap = [ linspace(0,1,126)'*[1 1 1]; [1 0 0] ];
colormap(cmap);
subplot(211); image( 126*(I-min(min(I)))/(max(I(:))-min(I(:))) ); drawnow

E = grid_graph([size(Q,2) size(Q,3)]);
G(1:3,:) = G(1:3,:) - 1;
FACT= 1e4;
G(4,:)= FACT*G(4,:);
f = int32(zeros(K,2,size(E,2)));
for epsilon = 5:-3:-1 % initial thresholds to test for maximality of nodes and edges
  disp(epsilon)
  Ik= reshape(maxsum(E,repmat(uint16(K),[1 prod(size(I))]),int32(G),int32(FACT*Q),f,uint32(2^epsilon)),size(Q));
  UNDEF= 10000;  AMBIG= 10001;
  J= repmat(UNDEF,size(I));
  for k= 1:K
    ik= Ik(k,:,:); ik= find(ik(:));
    J(ik( J(ik)~=UNDEF & J(ik)~=L(k) ))= AMBIG;
    J(ik( J(ik)==UNDEF & J(ik)~=AMBIG ))= L(k);
  end
  subplot(212); image( 126*J.*(J~=AMBIG) + 127*(J==AMBIG) ); drawnow;
end

%imwrite(1-I,'in.png');
%imwrite( 126*(1-J).*(J~=AMBIG) + 127*(J==AMBIG), cmap, 'out.png');
return
