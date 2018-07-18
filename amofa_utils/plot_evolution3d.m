function plot_evolution3d(trainingHistory,data,rr,ax)
% Function to plot 2D and 3D mixtures
% trainingHistory: obtained from amofa
% data: n x d dataset
% rr: n x 1 labels 
% If the labels are provided it will use them to color the classes
% -----------------------------------------------------------------------
% Copyleft (2014): Heysem Kaya and Albert Ali Salah
%
% This software is distributed under the terms
% of the GNU General Public License Version 3
% 
% Permission to use, copy, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
models=trainingHistory.models;
k_l=ceil(numel(models)/3);
[N,d]=size(data);
if (d>3 || d<=1)
    disp('Dimensionality must be 2 or 3!');
    return;
end
if nargin<4
    ax=[-10 10 -10 10];
end
figure,
for k=1:numel(models)
   
   count = 0;

   model=models{k};
   model.numSamples=N;
   clear mixture
   K= size(model.Mu,2);
   mixture=cell(K,1);
   for i = 1:K
        mixture{i}.Lambda = model.Lambda(:,count+1:count+model.numFactors(i));
        mixture{i}.Mu= model.Mu(:,i);
        mixture{i}.Pi =  model.Pi(i);
        mixture{i}.Psi =  model.Psi(:,i);
        mixture{i}.numFactors =  model.numFactors(i);
        count = count + model.numFactors(i);
    end

   %subplot(k_l,3,k)
   loglik=loglike3(data,model.Lambda,model.Psi,model.Mu,model.Pi,model.numFactors);
   descLen=getDescLen(model,loglik,1);
   
   
   axes1=subplot(k_l,3,k);
   
   if (d==3)
       view(axes1,[44 22]);
       if nargin>=3
         plot_mixture3d(mixture,data,rr)
       else
         plot_mixture3d(mixture,data);
       end
   else
       if nargin>=3
         plot_mixture(mixture,data,rr)
       else
         plot_mixture(mixture,data);
       end
   end
   axis(ax)
   title(['K=' num2str(K) ', DL=' num2str(descLen,'%.2f')])

end