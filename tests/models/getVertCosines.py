from math import cos,sin,sqrt,pi
import opensees as ops
def getVertCosines(eleTag  vecxz): # CONVERT-COMPLETE
       nodeTags = [eleNodes eleTag]
       n1Crds =       [nodeCoord '[lindex', nodeTags 0]]
       n2Crds =       [nodeCoord '[lindex', nodeTags 1]]
       eleLocX = (0)      , ([lindex n2Crds 0]-[lindex n1Crds 0])
       eleLocX = (1)      , ([lindex n2Crds 1]-[lindex n1Crds 1])
       eleLocX = (2)      , ([lindex n2Crds 2]-[lindex n1Crds 2])
       eleLocXMag =, (sqrt(eleLocX(0)**2+eleLocX(1)**2+eleLocX(2)**2))
       eleLocXNorm = (0)      , (eleLocX(0)/eleLocXMag)
       eleLocXNorm = (1)      , (eleLocX(1)/eleLocXMag)
       eleLocXNorm = (2)      , (eleLocX(2)/eleLocXMag)
       
       vecxzMag =, (sqrt([lindex vecxz 0]**2+[lindex vecxz 1]**2+[lindex vecxz 2]**2))

       array vecxzNorm = "
              0       [expr '[lindex', vecxz 0]/vecxzMag]
              1       [expr '[lindex', vecxz 1]/vecxzMag]
              2       [expr '[lindex', vecxz 2]/vecxzMag]
       "
       array eleLocYNorm = "
              0      , (      vecxzNorm(1)*eleLocXNorm(2)-vecxzNorm(2)*eleLocXNorm(1)  )
              1      , (-1 * (vecxzNorm(0)*eleLocXNorm(2)-vecxzNorm(2)*eleLocXNorm(0)) )
              2      , (      vecxzNorm(0)*eleLocXNorm(1)-vecxzNorm(1)*eleLocXNorm(0)  )
       "
       array eleLocZNorm = "
              0      , (-1 * (eleLocYNorm(1)*eleLocXNorm(2)-eleLocYNorm(2)*eleLocXNorm(1)) )
              1      , ( 1 * (eleLocYNorm(0)*eleLocXNorm(2)-eleLocYNorm(2)*eleLocXNorm(0)) )
              2      , (-1 * (eleLocYNorm(0)*eleLocXNorm(1)-eleLocYNorm(1)*eleLocXNorm(0)) )
       "
#       print("eleLocXNorm(0) eleLocXNorm(1) eleLocXNorm(2)")
#       print("eleLocYNorm(0) eleLocYNorm(1) eleLocYNorm(2)")
#       print("eleLocZNorm(0) eleLocZNorm(1) eleLocZNorm(2)")
       return "eleLocYNorm(2) eleLocZNorm(2) eleLocXNorm(2) eleTag eleLocXMag"

