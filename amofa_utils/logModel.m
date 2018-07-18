function [models]=logModel(models,Mu,Lambda,Pi,Psi,numFactors,logsTra,dl,action,tLog,tLogN)

i=numel(models)+1;
models{i}.Mu=Mu;
models{i}.Lambda=Lambda;
models{i}.Pi=Pi;
models{i}.Psi=Psi;
models{i}.numFactors=numFactors;
models{i}.logsTra = logsTra;
models{i}.dl=dl;
models{i}.action=action;
models{i}.tLog=tLog;
models{i}.tLogN=tLogN;

end