#pragma once

enum class Field {
  Unit,
  UnitY,
  UnitZ,
  UnitXX,
  UnitZZ,
  UnitYY,
  UnitYZ,
  UnitCentroidYY,
  UnitCentroidZZ,

  Density,
  DensityXX,
  DensityYY,
  DensityZZ,
  DensityXY,
  DensityXZ,

  PolarInertia,      // NOTE: not J0; this includes density

  // Energy

  StrainEnergy,
  KineticEnergy,

};

