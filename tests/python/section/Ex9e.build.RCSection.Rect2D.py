from math import cos,sin,sqrt,pi
import opensees as ops
# --------------------------------------------------------------------------------------------------
# build a section
#                Silvia Mazzoni & Frank McKenna, 2006
#

# SET UP ----------------------------------------------------------------------------
from opensees.units.english import * # define units

# MATERIAL parameters -------------------------------------------------------------------
IDconcCore = 1                          # material, ID, tag -- confined, core, concrete
IDconcCover = 2                         # material, ID, tag -- unconfined, cover, concrete
IDreinf = 3                             # material, ID, tag -- reinforcement
# nominal concrete compressive strength
fc = -4.0*ksi                   # CONCRETE, Compressive, Strength, ksi   (+Tension, -Compression)
Ec = 57*ksi*sqrt(-fc/psi)       # Concrete, Elastic, Modulus
# confined concrete
Kfc = 1.3                    # ratio of confined to unconfined concrete strength
fc1C = Kfc*fc                # CONFINED concrete (mander model) maximum stress
eps1C = 2.*fc1C/Ec           # strain at maximum stress 
fc2C = 0.2*fc1C              # ultimate stress
eps2C = 5*eps1C              # strain at ultimate stress 
# unconfined concrete
fc1U = fc                    # UNCONFINED concrete (todeschini parabolic model) maximum stress
eps1U = -0.003               # strain at maximum strength of unconfined concrete
fc2U = 0.2*fc1U              # ultimate stress
eps2U = -0.01                # strain at ultimate stress
_lambda = 0.1                # ratio between unloading slope at eps2 and initial slope Ec
# tensile-strength properties
ftC = -0.14*fc1C             # tensile, strength +tension
ftU = -0.14*fc1U             # tensile, strength +tension
Ets = ftU/0.002              # tension, softening, stiffness
# -----------
Fy = 66.8*ksi                # STEEL, yield, stress
Es = 29000.*ksi              # modulus, of, steel
Bs = 0.01                    # strain-hardening, ratio 
R0 = 18                      # control, the, transition, from, elastic, to, plastic, branches
cR1 = 0.925                  # control, the, transition, from, elastic, to, plastic, branches
cR2 = 0.15                   # control, the, transition, from, elastic, to, plastic, branches
uniaxialMaterial.Concrete02(IDconcCore,  fc1C, eps1C, fc2C, eps2C, _lambda, ftC, Ets)      # build, 'core', concrete (confined)
uniaxialMaterial.Concrete02(IDconcCover, fc1U, eps1U, fc2U, eps2U, _lambda, ftU, Ets)      # build, 'cover', concrete (unconfined)
uniaxialMaterial.Steel02(IDreinf, Fy, Es, Bs, R0, cR1, cR2)                           # build, 'reinforcement', 'material'

# section GEOMETRY -------------------------------------------------------------
HSec = 5.*ft                # Column Depth
BSec = 3.*ft                # Column Width
coverH = 5.*in              # Column cover to reinforcing steel NA parallel to H
coverB = 3.5*in             # Column cover to reinforcing steel NA parallel to B
numBarsTop = 16             # number of longitudinal-reinforcement bars in steel layer. -- top
numBarsBot = 16             # number of longitudinal-reinforcement bars in steel layer. -- bot
numBarsIntTot = 6           # number of longitudinal-reinforcement bars in steel layer. -- total intermediate skin reinforcement symm about y-axis
barAreaTop = 2.25*in2       # area of longitudinal-reinforcement bars -- top
barAreaBot = 2.25*in2       # area of longitudinal-reinforcement bars -- bot
barAreaInt = 2.25*in2       # area of longitudinal-reinforcement bars -- intermediate skin reinf
SecTag = 1                  # tag = for symmetric section

# FIBER SECTION properties -------------------------------------------------------------
#
#                        y
#                        ^
#                        |     
#             ---------------------   --   --
#             |  o      o      o  |   |    -- coverH
#             |                   |   |
#             |  o             o  |   |
  #    z <--- |         +         |  Hsec
#             |  o             o  |   |
#             |                   |   |
#             |  o o o o o o o o  |   |    -- coverH
#             ---------------------   --   --
#             |--------Bsec-------|
#             |--|    coverB   |--|
#
#                       y
#                       ^
#                       |    
#             ----------------------
#             |\      cover       /|
#             | \------Top-------/ |
#             |c|                |c|
#             |o|                |o|
#     z <-----|v|      core      |v|  Hsec
#             |e|                |e|
#             |r|                |r|
#             | /-------Bot------\ |
#             |/       cover      \|
#             ----------------------
#                       Bsec
#
# Notes
#    The core concrete ends at the NA of the reinforcement
#    The center of the section is at (0,0) in the local axis system

coverY = HSec/2.0              # The, distance, from, the, section, z-axis, to, the, edge, of, the, cover, concrete -- outer, edge, of, cover, concrete
coverZ = BSec/2.0              # The, distance, from, the, section, y-axis, to, the, edge, of, the, cover, concrete -- outer, edge, of, cover, concrete
coreY = coverY-coverH          # The, distance, from, the, section, z-axis, to, the, edge, of, the, core, concrete --  edge, of, the, core, concrete/inner, edge, of, cover, concrete
coreZ = coverZ-coverB          # The, distance, from, the, section, y-axis, to, the, edge, of, the, core, concrete --  edge, of, the, core, concrete/inner, edge, of, cover, concretenfY = 16                     # number, of, fibers, for, concrete, in, y-direction
nfY = 16                       # number, of, fibers, for, concrete, in, y-direction
nfZ = 4                            # number, of, fibers, for, concrete, in, z-direction
numBarsInt = numBarsIntTot/2       # number, of, intermediate, bars, per, side

section.FiberSection(SecTag, None, {
   patch.quad( core, [nfZ,nfY],[ -coreY coreZ -coreY -coreZ coreY -coreZ coreY coreZ),        # Define 'the', core 'patch'
   patch.quad(cover, [1,nfY], [-coverY coverZ -coreY coreZ coreY coreZ coverY coverZ),       # Define 'the', four 'cover', 'patches'
   patch.quad(cover, [1,nfY], [coreY=True, coreZ=coverY, =True, coverZ=coverY,  coverZ=coreY,  coreZ=True),
   patch.quad(cover, [nfZ,1], [coverY=coverZ,  coverY=True, coverZ=coreY, =True, coreZ=, -coreY coreZ),
   patch.quad(cover, [nfZ,1], [coreY coreZ coreY coreZ=coverY,  coverZ=coverY,  coverZ),
   layer.line(reinf, numBarsInt, barAreaInt, [[-coreY, coreZ],[ coreY, coreZ]]),       # intermediate 'skin', reinf. +z
   layer.line(reinf, numBarsInt, barAreaInt, [[-coreY,-coreZ],[ coreY,-coreZ]]),       # intermediate 'skin', reinf. -z
   layer.line(reinf, numBarsTop, barAreaTop, [[ coreY, coreZ],[ coreY,-coreZ]]),       # top 'layer', 'reinfocement'
   layer.line(reinf, numBarsBot, barAreaBot, [[-coreY, coreZ],[-coreY,-coreZ]]),       # bottom 'layer', 'reinforcement'
})

