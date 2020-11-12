# Function to calculate I squared statistic
# Source: https://github.com/MRCIEU/Health-and-Wellbeing-MR/blob/master/Two-sample%20MR%20Base%20script.R

Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}