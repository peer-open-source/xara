from math import cos,sin,sqrt,pi
import opensees as ops

# RAYLEIGH DAMPING
loadConst -time, 0.0; # maintain, 'constant', gravity, 'loads', and, 'reset', time, 'to', 'zero'
wipeAnalysis
ana.rayleigh(){*}[brace2::rayleigh_alpha {1, 0.03} {2, 0.03}]

# DYNAMIC EQ ANALYSIS
wipeAnalysis
dtfact = 1;
SF = 1;
IT = NM;
SOS = 1;
file, 'mkdir', dataDir/DynResponse;
Tol =                            1.0e-8;
maxNumIter =              100;
printFlag =              0;
TestType =              EnergyIncr;
NewmarkGamma =       0.50;
NewmarkBeta =              0.25;
ana.constraints(Transformation);
ana.numberer(RCM);
ana.test(TestType, Tol, maxNumIter, printFlag);
if IT == "NM":
  # algorithmType = NewtonLineSearch;
  algorithmType =   Newton;
  # system SparseGeneral -piv;
  ana.system(BandGeneral);
  ana.integrator('Newmark', NewmarkGamma NewmarkBeta);

inFilelong = "./Records/GM_long_global2.txt";
inFiletrans = "./Records/GM_trans_global2.txt";
inFilevert = "./Records/GM_vert_global2.txt";
iGMfile = "inFilelong inFiletrans inFilevert";
iGMdirection = "1 2 3";
iGMfact = "SF SF SF";
IDloadTag =               1;
iloop = "1 2 3";
dt =                            , 0.005*sec
NumPts =                12600;
ana.algorithm(algorithmType);
ana.analysis(Transient);
foreach, 'Sloop', iloop, 'GMdirection', iGMdirection, 'GMfile', iGMfile, 'GMfact', iGMfact {
       incr IDloadTag;
       inFile = "GMfile.txt";
       ana.timeSeries('Path', Sloop dt=dt,  filePath=GMfile,  factor='[expr', GMfact)*g];
       pattern 'UniformExcitation',  IDloadTag  GMdirection -accel  Sloop;              # create 'Uniform', 'excitation'
       if Sloop == 1:
              print("Horizontal component 1 checked")

       if {Sloop == 2} {
              print("Horizontal component 2 checked")

       if {Sloop == 3} {
              print("Vertical component checked")


DtAnalysis =              , dt/dtfact
TmaxAnalysis =       , dt*NumPts
Nsteps =                     , int(TmaxAnalysis/DtAnalysis)
print("Ground Motion: dt= DtAnalysis, NumPts=, (NumPts*dtfact), TmaxAnalysis= TmaxAnalysis");
print("Running dynamic ground motion analysis...")
t0 = [clock, clicks -millisec];  # Time, the, analysis
# DYNAMIC OUTPUT RECORDERS
# py -m fiberRecorders dataDir/modelDetails.json dataDir/dsr1 -e 4010 -d Dcol -s 0,-1 -l dsr1;
# py -m fiberRecorders dataDir/modelDetails.json dataDir/dsr2 -e 4010 -d Dcol -s 0,-1 -l dsr2;
# py -m fiberRecorders dataDir/modelDetails.json dataDir/dsr3 -e 4010 -d Dcol -s 0,-1 -l dsr3;
# py -m fiberRecorders dataDir/modelDetails.json dataDir/dsr4 -e 4010 -d Dcol -s 0,-1 -l dsr4;
# py -m fiberRecorders dataDir/modelDetails.json dataDir/dsr5 -e 4010 -d Dcol -s 0,-1 -l dsr5;
# py -m fiberRecorders dataDir/modelDetails.json dataDir/dsr6 -e 4010 -d Dcol -s 0,-1 -l dsr6;
dynRec = ""
dynRec.append([ana.recorder('Node', -file, dataDir)/DynResponse/Col4010TopDisp.txt, time=True, -model.node(402 -dof, 1 2, disp]))
dynRec.append([ana.recorder('Element', -xml, dataDir)/DynResponse/eleDef1.txt -ele, 2010, 2020, 3010, 3020, 4010, 4020, 5010, 5020, 6010, 6020, 7010, 7020, 8010, 8020, 9010, 9020, 10010, 10020, 11010, 11020, 12010, 12020, 12030, 13010, 13020, 13040, 14010, 14020, 14030, 14040, section, 1 deformation])
dynRec.append([ana.recorder('Element', -xml, dataDir)/DynResponse/eleDef4.txt -ele, 2010, 2020, 3010, 3020, 4010, 4020, 5010, 5020, 6010, 6020, 7010, 7020, 8010, 8020, 9010, 9020, 10010, 10020, 11010, 11020, 12010, 12020, 12030, 13010, 13020, 13040, 14010, 14020, 14030, 14040, section, 4 deformation])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch02&03_X.out -timeSeries 1 -time -node 1031 -dof 1 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch02&03_Y.out -timeSeries 2 -time -node 1031 -dof 2 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch12&13_X.out -timeSeries 1 -time -node 1030 -dof 1 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch12&13_Y.out -timeSeries 2 -time -node 1030 -dof 2 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch06&07_X.out -timeSeries 1 -time -node 307 -dof 1 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch06&07_Y.out -timeSeries 2 -time -node 307 -dof 2 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch14&15_X.out -timeSeries 1 -time -node 304 -dof 1 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch14&15_Y.out -timeSeries 2 -time -node 304 -dof 2 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch17&18_X.out -timeSeries 1 -time -node 401 -dof 1 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch17&18_Y.out -timeSeries 2 -time -node 401 -dof 2 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch19&20_X.out -timeSeries 1 -time -node 402 -dof 1 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch19&20_Y.out -timeSeries 2 -time -node 402 -dof 2 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch24&25_X.out -timeSeries 1 -time -node 407 -dof 1 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch24&25_Y.out -timeSeries 2 -time -node 407 -dof 2 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch22&23_X.out -timeSeries 1 -time -node 405 -dof 1 accel])
# dynRec.append([recorder Node -file dataDir/DynResponse/Ch22&23_Y.out -timeSeries 2 -time -node 405 -dof 2 accel])
dynRec.append([ana.recorder('Node', -file, dataDir)/DynResponse/AX4.out -ana.timeSeries(1, time=True, -model.node)(402, 405 -dof, 1 accel]))
dynRec.append([ana.recorder('Node', -file, dataDir)/DynResponse/AY4.out -ana.timeSeries(2, time=True, -model.node)(402, 405 -dof, 2 accel]))
dynRec.append([ana.recorder('Node', -file, dataDir)/DynResponse/AZ4.out -ana.timeSeries(3, time=True, -model.node)(402, 405 -dof, 3 accel]))
dynRec.append([ana.recorder('Node', -file, dataDir)/DynResponse/AZ3.out -ana.timeSeries(3, time=True, -model.node)(30003 -dof, 3 accel]))
# dynRec.append([recorder Node -binary dataDir/DynResponse/allNodes.bin              -precision 5       -time       -dof 1 2 3 4 5 6 disp])
print, JSON=True, file=dataDir, /DynResponse/modelDetails.json

for {ik = 1} {ik <= Nsteps} {ik += 1} {
   ok =  ana.analyze(1, DtAnalysis);
       # Convergence
   if ok != 0:
     print("Trying Bisection ...");
     ana.algorithm('NewtonLineSearch', type=Bisection, );
     ok = [ana.analyze(1 DtAnalysis])

   if ok != 0:
          print("Trying Secant ...");
          ana.algorithm('NewtonLineSearch', type=Secant, );
          ok = [ana.analyze(1 DtAnalysis])

   if ok != 0:
          print("Trying RegulaFalsi ...");
          ana.algorithm('NewtonLineSearch', type=RegulaFalsi, );
          ok = [ana.analyze(1 DtAnalysis])

   if ok != 0:
      print("Trying KrylovNewton ...");
      ana.algorithm(KrylovNewton);
      ok = [ana.analyze(1 DtAnalysis]);
      ana.algorithm(algorithmType);

   if ok != 0:
      print("Trying Newton...");
      ana.algorithm(Newton);
      ok = [ana.analyze(1 DtAnalysis]);
      ana.algorithm(algorithmType);

   if ok != 0:
      print("Trying Broyden ...");
      ana.algorithm(Broyden);
      ok = [ana.analyze(1 DtAnalysis]);
      ana.algorithm(algorithmType);

   if ok != 0:
      print("Trying BFGS ...");
      ana.algorithm(BFGS);
      ok = [ana.analyze(1 DtAnalysis]);
      ana.algorithm(algorithmType);

   if SOS == 1:
       if ok != 0:
            print("Trying OS ...");
            ana.integrator('AlphaOS', 1.00);
            ana.algorithm(Linear);
            ok = [ana.analyze(1 DtAnalysis]);
            ana.integrator('Newmark', NewmarkGamma NewmarkBeta);
            ana.algorithm(algorithmType);


       if ok != 0:
              Nstepsmax =, ik-1
              break;


if [expr ik-1] == Nsteps:
       AnalysisA =, (1)
       print("Analysis complete")
else:
       AnalysisA =, (0) };
puts stderr "The analysis time was [expr {([clock clicks -millisec]-t0)/1000.}] seconds";
