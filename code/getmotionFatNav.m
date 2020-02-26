function [ mot_mean ] = getmotionFatNav( map_corr, brd_crds, MPos, voxRes)
%Calculates motion in ROIs as displacement

G = gradientvector(map_corr);
G_matx = im2mat(G{1});
G_maty = im2mat(G{2});
G_matz = im2mat(G{3});
grad = [G_matx(brd_crds)',G_maty(brd_crds)',G_matz(brd_crds)'];

[x, y, z] = ind2sub(size(map_corr),brd_crds);

spacing=mean(voxRes); % approximation
offset=[-(size(map_corr,1)/2) -(size(map_corr,2)/2) -(size(map_corr,3)/2)];

c = spacing.*[x(:)+offset(1) y(:)+offset(2) z(:)+offset(3)]; %manual changing of rotation axes
c(:,4) = 1;

% apply transform
crall = zeros(numel(x),3,size(MPos.mats,3));
for ii=1:size(MPos.mats,3)
  m=MPos.mats(:,:,ii);
  cr=(m*c')';
  crall(:,:,ii)=cr(:,1:3);
end

nf=size(crall,3);
mid=round(nf/2);
idx=(mid-5):(mid+5);

% take center 11 time points as reference
crallrel=bsxfun(@minus,crall,mean(crall(:,:,idx),3));

% weighting vector, stronger weight in center
w=exp(-linspace(-1.8,1.8,nf).^2);
w=w/sum(w);

score=squeeze(mean(sqrt(sum(crallrel.^2,2))));
wscore=w*score;

% output scalar
mot_mean=wscore;

end

