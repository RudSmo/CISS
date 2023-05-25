function [out] = pagemtimes(A,B)
  
  A_layers = size(A,3);
  B_layers = size(B,3);
  
  out = zeros(size(A,1),size(B,2),A_layers);
  
 if A_layers == 1 
  
    for n = 1:B_layers
    
    out(:,:,n) = A(:,:) * B(:,:,n);
  
     end
  
 elseif B_layers == 1
  
      for n = 1:A_layers
    
      out(:,:,n) = A(:,:,n) * B(:,:);
  
      end
else
      
     for n = 1:B_layers
    
    out(:,:,n) = A(:,:,n) * B(:,:,n);
  
       end

       
 end
 
 end