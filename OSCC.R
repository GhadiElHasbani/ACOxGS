calculate_OSCC <- function(x, y, full_res = FALSE){
  require(dplyr)
  # X and Y should have the same length
  if(is.vector(x) & is.vector(y) & length(x) == length(y)) {
    N <- length(x)
    if(N > 1){
      input <- data.frame(x = x, y = y, id = 1:N)
      
      sorted_x <- input %>% arrange(x, id) %>% dplyr::select(-id)#id for stable sort to break ties consistently
      sorted_y <- input %>% arrange(y, id) %>% dplyr::select(y)
      
      numerator <- 0
      denominator <- 0
      for(i in 1:N){
        common_term <- sorted_x$x[i] - sorted_x$x[c(N-i+1)]
        #if(!is.finite(common_term)){
        #  print("===================================================================================CT IS NULL:")
        #  print(paste(common_term, "<-", sorted_x$x[i], "-",  sorted_x$x[c(N-i+1)]))
        #}
        numerator <- numerator + common_term * sorted_x$y[i]
        denominator <- denominator + common_term * sorted_y$y[i]
      }
      
      if(denominator == 0){
        oscc <- 0
      } else {
        oscc <- numerator/denominator
      }
      if(full_res){
        return(list(oscc = oscc, numerator = numerator, denominator = denominator, common_term = common_term))
      }else{
        return(oscc)
      }
    } else {
      print("Warning: OSCC calculated with vectors of length 1! Input vectors should have length N > 1. Default output: 1")
      if(full_res){
        return(list(oscc = 1, numerator = 1, denominator = 1, common_term = 1))
      }else{
        return(1)
      }
    }
    
  } else {
    e <- simpleError("Input Error... usage: calculate_OSCC(x,y) where x and y are vectors of the same length N > 1")
    stop(e)
  }
}
