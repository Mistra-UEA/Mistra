C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C                                                                  
C Sparse Data Header File                                          
C                                                                  
C Generated by KPP-2.2.4 for Mistra symbolic chemistry Kinetics PreProcessor
C       (http://www.cs.vt.edu/~asandu/Software/Kpp)                
C KPP is distributed under GPL, the general public licence         
C       (http://www.gnu.org/copyleft/gpl.html)                     
C (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           
C (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            
C     With important contributions from:                           
C        M. Damian, Villanova University, USA                      
C        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
C                                                                  
C File                 : gas_Sparse.h                              
C Time                 : Wed Jul 14 18:26:57 2021                  
C Working directory    : /local/josue/Mistra_2019/src/mech         
C Equation file        : gas.k                                     
C Output root filename : gas                                       
C                                                                  
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




C  ----------> Sparse Jacobian Data                                

C LU_ICOL_g - Column indexes of the LU Jacobian of variables         
      INTEGER LU_ICOL_g(1110)
      COMMON /SDATA_g/ LU_ICOL_g
C LU_CROW_g - Compressed row indexes of the LU Jacobian of variables 
      INTEGER LU_CROW_g(103)
      COMMON /SDATA_g/ LU_CROW_g
C LU_DIAG_g - Diagonal indexes of the LU Jacobian of variables       
      INTEGER LU_DIAG_g(103)
      COMMON /SDATA_g/ LU_DIAG_g

