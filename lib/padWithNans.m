function data=padWithNans(varargin)
sizesOfVectors=zeros(1,numel(varargin));
for i=1:length(sizesOfVectors)
    sizesOfVectors(i)=numel(varargin{i});
    if size(varargin{i},1) ~=1 & size(varargin{i},2)==1
        varargin{i}=varargin{i}';
    elseif min(size(varargin{i}))~=1
        error('Error: Input data should be a vector');
    end    
        
        
end

data=nan(numel(sizesOfVectors),max(sizesOfVectors));
for i=1:length(sizesOfVectors)
   data(i,1:sizesOfVectors(i))=varargin{i}; 
end

data=data';

end