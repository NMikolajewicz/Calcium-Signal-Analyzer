function y = AT(w) 
y=( sum(w) * ones( length( w ), 1 ) - [ 0; cumsum( w( 1 : end - 1 ) ) ] ); 
end
