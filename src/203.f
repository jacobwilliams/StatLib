      subroutine mixture(a, k, m, c, x, n, alpha, mean, sd, tol, nobvs,         
     +    logl, counter, putout, ifault)                                        
c                                                                               
c     Translation from Algol 60 by Alan Miller of AS 203.                       
c                                                                               
c     ALGORITHM AS 203 APPL. STATIST. (1984) VOL.33, NO. 3                      
c                                                                               
c     This algorithm calculates the maximum likelihood estimates of the         
c     parameters of a mixture of normal or exponential or Poisson or            
c     binomial distributions.   These parameters are the mixing                 
c     proportions, the means and (in the normal distribution case)              
c     standard deviations.   It also calculates the log-likelihood              
c     function and the number of iterations taken to satisfy the                
c     tolerance value.                                                          
c                                                                               
c     Users must provide their own routine `putout' which prints the            
c     estimated parameter values after each iteration.   The form of            
c     this routine is:                                                          
c        subroutine putout(k, alpha, mean, sd, logl)                            
c        integer k                                                              
c        real alpha(k), mean(k), sd(k), logl                                    
c                                                                               
c     Notes:                                                                    
c     The original (Algol) variable names have been retained.   This            
c     means that some variable names have more than 6 characters, such          
c     as `counter', `newalpha', etc.                                            
c     The authors treat the normal distribution as if it were a discrete        
c     distribution.                                                             
c     In Fortran 77 it is necessary to give explicit dimensions to some         
c     of the temporary arrays which are not passed as arguments.                
c     Maximum values have been set for k and m (kmax and mmax) in the           
c     parameter statement under `Local variables' below.                        
c                                                                               
      integer a, k, m, c, n(m), nobvs, counter, ifault                          
      real x(m), alpha(k), mean(k), sd(k), tol, logl                            
      external putout                                                           
c                                                                               
c     Local variables                                                           
c                                                                               
      integer kmax, mmax, i, j                                                  
      parameter (kmax = 20, mmax = 100)                                         
      logical test                                                              
      real sumalpha, part, oldlogl, newalpha(kmax), newmean(kmax),              
     +    newsd(kmax), dt(kmax), nt(kmax), vt(kmax), g(mmax),                   
     +    f(mmax, kmax), zero, one, half                                        
      data zero /0.0/, one /1.0/, half /0.5/                                    
c                                                                               
      if (a .lt. 1 .or. a .gt. 4) then                                          
        ifault = 1                                                              
        return                                                                  
      else                                                                      
        ifault = 0                                                              
      end if                                                                    
                                                                                
      do 10 i = 1, m-1                                                          
        if (x(i) .gt. x(i+1)) then                                              
          ifault = 5                                                            
          return                                                                
        end if                                                                  
   10 continue                                                                  
                                                                                
      if (nobvs .lt. 2 * m) then                                                
        ifault = 6                                                              
        return                                                                  
      end if                                                                    
                                                                                
      do 20 i = 1, m                                                            
        if (n(i) .lt. 0) then                                                   
          ifault = 6                                                            
          return                                                                
        end if                                                                  
        if (a .ne. 1) then                                                      
          if (x(i) .lt. zero) then                                              
            ifault = 7                                                          
            return                                                              
          end if                                                                
        end if                                                                  
   20 continue                                                                  
                                                                                
      oldlogl = zero                                                            
      counter = 0                                                               
c                                                                               
c     Start of iterative cycle                                                  
c                                                                               
   30 counter = counter + 1                                                     
      do 40 j = 1, k                                                            
        if (alpha(j) .gt. one .or. alpha(j) .lt. zero) then                     
          ifault = 2                                                            
          return                                                                
        end if                                                                  
        if (mean(j) .ge. x(m) .or. mean(j) .le. x(1)) then                      
          ifault = 3                                                            
          return                                                                
        end if                                                                  
        if (a .eq. 1) then                                                      
          if (sd(j) .le. zero) then                                             
            ifault = 4                                                          
            return                                                              
          end if                                                                
        end if                                                                  
   40 continue                                                                  
                                                                                
      do 60 i = 1, k-1                                                          
        do 50 j = i+1, k                                                        
          if (mean(i) .eq. mean(j)) then                                        
            if (a .eq. 1) then                                                  
              if (sd(i) .eq. sd(j)) then                                        
                ifault = 9                                                      
                return                                                          
              end if                                                            
            else                                                                
              ifault = 8                                                        
              return                                                            
            end if                                                              
          end if                                                                
   50   continue                                                                
   60 continue                                                                  
                                                                                
      logl = zero                                                               
      do 80 i = 1, m                                                            
        g(i) = zero                                                             
        do 70 j = 1, k                                                          
c                                                                               
c     a = 1 denotes normal mixture                                              
c     a = 2 denotes exponential mixture                                         
c     a = 3 denotes Poisson mixture                                             
c     a = 4 denotes binomial mixture                                            
c                                                                               
          if (a .eq. 1) then                                                    
            f(i,j) = exp(-half*((x(i) - mean(j))/sd(j))**2) / sd(j)             
          else if (a .eq. 2) then                                               
            f(i,j) = exp(-x(i)/mean(j)) / mean(j)                               
          else if (a .eq. 3) then                                               
            if (x(i) .eq. x(1)) then                                            
              f(i,j) = exp(-mean(j)) * mean(j)**x(i)                            
            else                                                                
              f(i,j) = f(i-1,j) * mean(j)                                       
            end if                                                              
          else                                                                  
            if (x(i) .eq. x(1)) then                                            
              f(i,j) = (one - mean(j) / x(m))**x(m) * (mean(j) /                
     +                  (x(m) - mean(j)))**x(i)                                 
            else                                                                
              f(i,j) = f(i-1,j) * (mean(j) / (x(m) - mean(j)))                  
            end if                                                              
          end if                                                                
          g(i) = g(i) + alpha(j) * f(i,j)                                       
   70   continue                                                                
        logl = logl + n(i) * log(g(i))                                          
   80 continue                                                                  
c                                                                               
c     Calculation of the probability densities of the sub-populations           
c     which form the mixture, and the log-likelihood function.                  
c                                                                               
      test = .false.                                                            
      sumalpha = zero                                                           
      do 100 j = 1, k                                                           
        nt(j) = zero                                                            
        dt(j) = zero                                                            
        vt(j) = zero                                                            
        do 90 i = 1, m                                                          
          part = f(i,j) * n(i) / g(i)                                           
          dt(j) = dt(j) + part                                                  
          nt(j) = nt(j) + part * x(i)                                           
          if (a .eq. 1) vt(j) = vt(j) + part * (x(i) - mean(j))**2              
   90   continue                                                                
c                                                                               
c     Calculation of denominators and numerators of new estimates.              
c                                                                               
        newmean(j) = nt(j) / dt(j)                                              
        if (j .ne. k) then                                                      
          newalpha(j) = alpha(j) * dt(j) / nobvs                                
          sumalpha = sumalpha + newalpha(j)                                     
        else                                                                    
          newalpha(k) = one - sumalpha                                          
        end if                                                                  
        if (a .eq. 1) newsd(j) = sqrt(vt(j) / dt(j))                            
c                                                                               
c     Convergence test.                                                         
c                                                                               
        if (abs(oldlogl - logl) .gt. tol) test = .true.                         
        oldlogl = logl                                                          
        alpha(j) = newalpha(j)                                                  
        mean(j) = newmean(j)                                                    
        if (a .eq. 1) sd(j) = newsd(j)                                          
  100 continue                                                                  
c                                                                               
      if (c .gt. 0) then                                                        
        if ((counter/c)*c .eq. counter) then                                    
          call putout(k, alpha, mean, sd, logl)                                 
        end if                                                                  
      end if                                                                    
                                                                                
      if (test) go to 30                                                        
      return                                                                    
      end                                                                       
