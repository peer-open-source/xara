#include <cmath>
#include <vector>
#include <Parsing.h>
#include <Logging.h>
#include <ArgumentTracker.h>
#include <ModelRegistry.h>
#include <uniaxial/ASDConcrete1DMaterial.h>
#include <damage/ASDConcrete3DMaterial.h>

// anonymous namespace for utilities
namespace {

  double bezier3(double xi,
    double x0, double x1, double x2,
    double y0, double y1, double y2)
  {
    double A = x0 - 2.0 * x1 + x2;
    double B = 2.0 * (x1 - x0);
    double C = x0 - xi;
    if (fabs(A) < 1.0e-12) {
      x1 = x1 + 1.0E-6 * (x2 - x0);
      A = x0 - 2.0 * x1 + x2;
      B = 2.0 * (x1 - x0);
      C = x0 - xi;
    }
    if (A == 0.0)
      return 0.0;

    double D = B * B - 4.0 * A * C;
    double t = (sqrt(D) - B) / (2.0 * A);

    return (y0 - 2.0 * y1 + y2) * t * t + 2.0 * (y1 - y0) * t + y0;
  }

struct ConcreteBackbone {
  std::vector<double> Te;
  std::vector<double> Ts;
  std::vector<double> Td;
  std::vector<double> Ce;
  std::vector<double> Cs;
  std::vector<double> Cd;
};

}



static int
CreateConcreteBackbone(
  ConcreteBackbone& backbone,
  double E,
  double fc,
  double ft,
  double Gc,
  double Gt,
  bool have_lch_ref,
  double& lch_ref)
{
  double ec = 2.0*fc/E;

  if (!have_lch_ref) {
    //
    // _get_lch_ref from ASDConcrete[13]D_MakeLaws.py
    //

    // min lch for tension
    double et_el = ft/E;
    double Gt_min = 0.5*ft*et_el;
    double hmin_t = 0.01*Gt/Gt_min;

    // min lch for compression
    double ec1 = fc/E;
    double ec_pl = (ec - ec1)*0.4 + ec1;
    double Gc_min = 0.5*fc*(ec - ec_pl);
    double hmin_c = 0.01*Gc/Gc_min;

    lch_ref = std::min(hmin_c, hmin_t);
  }

  //
  // _make_tension from ASDConcrete[13]D_MakeLaws.py
  //

  Gt = Gt/lch_ref;

  double f0 = 0.9*ft;
  double f1 = ft;
  double e0 = f0/E;
  double e1 = 1.5*f1/E;
  double ep = e1 - f1/E;
  double f2 = 0.2*ft;
  double f3 = 1.0e-3*ft;
  double w2 = Gt/ft;
  double w3 = 5.0*w2;
  double e2 = w2 + f2/E + ep;
  if (e2 <= e1)
    e2 = 1.001*e1;
  double e3 = w3 + f3/E + ep;
  if (e3 <= e2)
    e3 = 1.001*e2;
  double e4 = 10.0*e3;

  backbone.Te = {0.0, e0, e1, e2, e3, e4};
  backbone.Ts = {0.0, f0, f1, f2, f3, f3};
  backbone.Td = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double Tpl[6] = {0.0, 0.0, ep, 0.9*e2, 0.8*e3, 0.8*e3};

  for (int i = 2; i < 6; i++) {
    double xi = backbone.Te[i];
    double si = backbone.Ts[i];
    double xipl = Tpl[i];
    double xipl_max = xi - si/E;
    xipl = std::min(xipl, xipl_max);
    double qi = (xi - xipl)*E;
    backbone.Td[i] = 1.0 - si/qi;
  }

  //
  // _make_compression from ASDConcrete[13]D_MakeLaws.py
  //

  Gc = Gc/lch_ref;

  double fc0 = 0.5*fc;
  double ec0 = fc0/E;
  double ec1 = fc/E;
  double fcr = 0.1*fc;
  double ec_pl = (ec - ec1)*0.4 + ec1;
  double Gc1 = 0.5*fc*(ec - ec_pl);
  double Gc2 = std::max(0.01*Gc1, Gc - Gc1);
  double ecr = ec + 2.0*Gc2/(fc + fcr);

  const int nc = 10;

  backbone.Ce.resize(nc + 3);
  backbone.Cs.resize(nc + 3);

  backbone.Ce[0] = 0.0;
  backbone.Ce[1] = ec0;

  backbone.Cs[0] = 0.0;
  backbone.Cs[1] = fc0;

  double Cpl[nc + 3];
  Cpl[0] = 0.0;
  Cpl[1] = 0.0;

  double dec = (ec - ec0)/(nc - 1);

  for (int i = 0; i < nc - 1; i++) {
    double iec = ec0 + (i + 1)*dec;
    backbone.Ce[i + 2] = iec;
    backbone.Cs[i + 2] = bezier3(iec, ec0, ec1, ec, fc0, fc, fc);
    Cpl[i + 2] = Cpl[i + 1] + 0.7*(iec - Cpl[i + 1]);
  }

  backbone.Ce[nc + 1] = ecr;
  backbone.Cs[nc + 1] = fcr;
  Cpl[nc + 1] = Cpl[nc] + 0.7*(ecr - Cpl[nc]);

  backbone.Ce[nc + 2] = ecr + ec0;
  backbone.Cs[nc + 2] = fcr;
  Cpl[nc + 2] = Cpl[nc + 1];

  backbone.Cd.resize(nc + 3);
  backbone.Cd[0] = 0.0;
  backbone.Cd[1] = 0.0;

  for (int i = 2; i < nc + 3; i++) {
    double xi = backbone.Ce[i];
    double si = backbone.Cs[i];
    double xipl = Cpl[i];
    double xipl_max = xi - si/E;
    xipl = std::min(xipl, xipl_max);
    double qi = (xi - xipl)*E;
    backbone.Cd[i] = 1.0 - si/qi;
  }

  return TCL_OK;
}


int
TclCommand_addASDConcrete1D(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  // Accept:
  //   nDMaterial ASDConcrete1D tag? E?
  //     <-fc fc?> <-ft ft?> <-Gc Gc?> <-Gt Gt?>
  //     <-Te Te?> <-Ts Ts?> <-Td Td?>
  //     <-Ce Ce?> <-Cs Cs?> <-Cd Cd?>
  //     <-implex>
  //     <-implexControl implexErrorTolerance? implexTimeReductionLimit?>
  //     <-implexAlpha alpha?>
  //     <-eta eta?> <-tangent> <-autoRegularization lch_ref?>

  static bool first_done = false;
  if (!first_done) {
    opslog << "Using ASDConcrete1D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
    first_done = true;
  }

  const char* info =
    "nDMaterial ASDConcrete1D tag? E? "
    "<-fc fc?> <-ft ft?> <-Gc Gc?> <-Gt Gt?> "
    "<-Te Te?> <-Ts Ts?> <-Td Td?> "
    "<-Ce Ce?> <-Cs Cs?> <-Cd Cd?> "
    "<-implex> "
    "<-implexControl implexErrorTolerance? implexTimeReductionLimit?> "
    "<-implexAlpha alpha?> "
    "<-eta eta?> <-tangent> <-autoRegularization lch_ref?>";

  if (argc < 4) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: " << info << "\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid 'tag'.\n";
    return TCL_ERROR;
  }

  double E;
  if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid E " << argv[3] << "\n";
    return TCL_ERROR;
  }


  bool implex = false;
  bool implex_control = false;
  double implex_error_tolerance = 0.05;
  double implex_time_redution_limit = 0.01;
  double implex_alpha = 1.0;
  double eta = 0.0;
  bool tangent = false;
  bool auto_regularization = false;
  double lch_ref = 1.0;

  double fc = 0.0;
  double ft = 0.0;
  double Gc = 0.0;
  double Gt = 0.0;

  bool have_fc = false;
  bool have_ft = false;
  bool have_Gc = false;
  bool have_Gt = false;
  bool have_lch_ref = false;

  ConcreteBackbone backbone;

  int i = 4;
  while (i < argc) {
    const char* value = argv[i++];

    if ((strcmp(value, "-fc") == 0) || (strcmp(value, "-Fc") == 0)) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-fc' given without the next argument fc.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &fc) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'fc'.\n";
        return TCL_ERROR;
      }
      have_fc = true;
      i++;
    }
    else if ((strcmp(value, "-ft") == 0) || (strcmp(value, "-Ft") == 0)) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-ft' given without the next argument ft.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &ft) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'ft'.\n";
        return TCL_ERROR;
      }
      have_ft = true;
      i++;
    }
    else if (strcmp(value, "-Gc") == 0) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-Gc' given without the next argument Gc.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Gc) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'Gc'.\n";
        return TCL_ERROR;
      }
      have_Gc = true;
      i++;
    }
    else if (strcmp(value, "-Gt") == 0) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-Gt' given without the next argument $Gt.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Gt) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'Gt'.\n";
        return TCL_ERROR;
      }
      have_Gt = true;
      i++;
    }
    else if (strcmp(value, "-implex") == 0) {
      implex = true;
    }
    else if (strcmp(value, "-implexControl") == 0) {
      implex_control = true;

      if (argc - i < 2) {
        opserr << OpenSees::PromptValueError << "'-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i], &implex_error_tolerance) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'implexErrorTolerance'.\n";
        return TCL_ERROR;
      }
      i++;

      if (Tcl_GetDouble(interp, argv[i], &implex_time_redution_limit) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'implexTimeReductionLimit'.\n";
        return TCL_ERROR;
      }
      i++;
    }
    else if (strcmp(value, "-implexAlpha") == 0) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-implexAlpha' given without the next argument $alpha.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &implex_alpha) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'alpha'.\n";
        return TCL_ERROR;
      }
      i++;
    }
    else if (strcmp(value, "-eta") == 0) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-eta' given without the next argument $eta.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &eta) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'eta'.\n";
        return TCL_ERROR;
      }
      i++;
    }
    else if (strcmp(value, "-tangent") == 0) {
      tangent = true;
    }
    else if (strcmp(value, "-autoRegularization") == 0) {
      auto_regularization = true;

      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-autoRegularization' given without the next 1 argument $lch_ref.\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i], &lch_ref) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'lch_ref'.\n";
        return TCL_ERROR;
      }

      have_lch_ref = true;
      i++;
    }
    else if (
      strcmp(value, "-Te") == 0 ||
      strcmp(value, "-Ts") == 0 ||
      strcmp(value, "-Td") == 0 ||
      strcmp(value, "-Ce") == 0 ||
      strcmp(value, "-Cs") == 0 ||
      strcmp(value, "-Cd") == 0)
    {
      std::vector<double>* list = nullptr;
      const char* list_name = value + 1;

      if (strcmp(value, "-Te") == 0)
        list = &backbone.Te;
      else if (strcmp(value, "-Ts") == 0)
        list = &backbone.Ts;
      else if (strcmp(value, "-Td") == 0)
        list = &backbone.Td;
      else if (strcmp(value, "-Ce") == 0)
        list = &backbone.Ce;
      else if (strcmp(value, "-Cs") == 0)
        list = &backbone.Cs;
      else
        list = &backbone.Cd;

      list->clear();

      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'" << list_name << "' requested but not provided.\n";
        return TCL_ERROR;
      }

      double item;
      if (Tcl_GetDouble(interp, argv[i], &item) == TCL_OK) {
        while (i < argc) {
          if (Tcl_GetDouble(interp, argv[i], &item) != TCL_OK) {
            Tcl_ResetResult(interp);
            break;
          }

          list->push_back(item);
          i++;
        }
      }
      else {
        Tcl_ResetResult(interp);

        TCL_Char** listArgv = nullptr;
        Tcl_Size listArgc = 0;

        if (Tcl_SplitList(interp, argv[i], &listArgc, &listArgv) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "cannot parse the '" << list_name << "' list.\n";
          return TCL_ERROR;
        }

        for (Tcl_Size j = 0; j < listArgc; j++) {
          if (Tcl_GetDouble(interp, listArgv[j], &item) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "cannot parse the '" << list_name << "' list.\n";
            Tcl_Free((char*)listArgv);
            return TCL_ERROR;
          }

          list->push_back(item);
        }

        Tcl_Free((char*)listArgv);
        i++;
      }
    }
    else {
      opserr << OpenSees::PromptValueError << "unknown option '" << value << "'.\n";
      return TCL_ERROR;
    }
  }

  if (E <= 0.0) {
    opserr << OpenSees::PromptValueError << "invalid value for 'E' (" << E << "). It should be strictly positive.\n";
    return TCL_ERROR;
  }

  if (!have_fc && (have_Gc || have_Gt)) {
    opserr << OpenSees::PromptValueError << "'-Gc' and '-Gt' require '-fc'.\n";
    return TCL_ERROR;
  }

  if (have_fc && !have_ft)
    ft = 0.1*fc;

  if (have_fc) {
    if (!have_Gt)
      Gt = 0.073*pow(fc, 0.18);

    if (!have_Gc)
      Gc = 2.0*Gt*(fc*fc)/(ft*ft);

    if (CreateConcreteBackbone(backbone, E, fc, ft, Gc, Gt, have_lch_ref, lch_ref) != TCL_OK)
      return TCL_ERROR;
  }

  if (backbone.Te.size() < 1) {
    opserr << OpenSees::PromptValueError << "'Te' list is empty. At least 1 non-zero value should be provided.\n";
    return TCL_ERROR;
  }

  if (backbone.Ts.size() != backbone.Te.size()) {
    opserr << OpenSees::PromptValueError << "'Te' (size = " <<
      static_cast<int>(backbone.Te.size()) << ") and 'Ts' (size = " <<
      static_cast<int>(backbone.Ts.size()) << ") lists should have the same size.\n";
    return TCL_ERROR;
  }

  if (backbone.Td.size() == 0) {
    backbone.Td.resize(backbone.Te.size(), 0.0);
  }
  else if (backbone.Td.size() != backbone.Te.size()) {
    opserr << OpenSees::PromptValueError << "'Te' (size = " <<
      static_cast<int>(backbone.Te.size()) << ") and 'Td' (size = " <<
      static_cast<int>(backbone.Td.size()) << ") lists should have the same size.\n";
    return TCL_ERROR;
  }

  if (backbone.Ce.size() < 1) {
    opserr << OpenSees::PromptValueError 
           << "'Ce' list is empty. At least 1 non-zero value should be provided.\n";
    return TCL_ERROR;
  }

  if (backbone.Cs.size() != backbone.Ce.size()) {
    opserr << OpenSees::PromptValueError << "'Ce' (size = " <<
      static_cast<int>(backbone.Ce.size()) << ") and 'Cs' (size = " <<
      static_cast<int>(backbone.Cs.size()) << ") lists should have the same size.\n";
    return TCL_ERROR;
  }

  if (backbone.Cd.size() == 0) {
    backbone.Cd.resize(backbone.Ce.size(), 0.0);
  }
  else if (backbone.Cd.size() != backbone.Ce.size()) {
    opserr << OpenSees::PromptValueError << "'Ce' (size = " <<
      static_cast<int>(backbone.Ce.size()) << ") and 'Cd' (size = " <<
      static_cast<int>(backbone.Cd.size()) << ") lists should have the same size.\n";
    return TCL_ERROR;
  }

  ASDConcrete1DMaterial::HardeningLaw HT(
    tag,
    ASDConcrete1DMaterial::HardeningLawType::Tension,
    E,
    backbone.Te,
    backbone.Ts,
    backbone.Td);

  if (!HT.isValid()) {
    opserr << OpenSees::PromptValueError << "Tensile hardening law is not valid.\n";
    return TCL_ERROR;
  }

  ASDConcrete1DMaterial::HardeningLaw HC(
    tag,
    ASDConcrete1DMaterial::HardeningLawType::Compression,
    E,
    backbone.Ce,
    backbone.Cs,
    backbone.Cd);

  if (!HC.isValid()) {
    opserr << OpenSees::PromptValueError << "Compressive hardening law is not valid.\n";
    return TCL_ERROR;
  }

  UniaxialMaterial* instance = new ASDConcrete1DMaterial(
    tag,
    E, eta,
    implex, implex_control, implex_error_tolerance, implex_time_redution_limit, implex_alpha,
    tangent, auto_regularization, lch_ref,
    HT, HC);

  if (builder->addTaggedObject<UniaxialMaterial>(*instance) != TCL_OK) {
    delete instance;
    opserr << OpenSees::PromptValueError << "could not add material with tag " << tag << ".\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}


int
TclCommand_addASDConcrete3D(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  // Accept:
  //   nDMaterial ASDConcrete3D tag? E? v?
  //     <-fc fc?> <-ft ft?> <-Gc Gc?> <-Gt Gt?>
  //     <-Te Te?> <-Ts Ts?> <-Td Td?>
  //     <-Ce Ce?> <-Cs Cs?> <-Cd Cd?>
  //     <-rho rho?> <-Kc Kc?>
  //     <-implex>
  //     <-implexControl implexErrorTolerance? implexTimeReductionLimit?>
  //     <-implexAlpha alpha?>
  //     <-cdf cdf?>
  //     <-crackPlanes nct? ncc? smoothingAngle?>
  //     <-eta eta?> <-tangent> <-autoRegularization lch_ref?>

  static bool first_done = false;
  if (!first_done) {
    opslog << "Using ASDConcrete3D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
    first_done = true;
  }

  const char* info =
    "nDMaterial ASDConcrete3D tag? E? v? "
    "<-fc fc?> <-ft ft?> <-Gc Gc?> <-Gt Gt?> "
    "<-Te Te?> <-Ts Ts?> <-Td Td?> "
    "<-Ce Ce?> <-Cs Cs?> <-Cd Cd?> "
    "<-rho rho?> <-Kc Kc?> "
    "<-implex> "
    "<-implexControl implexErrorTolerance? implexTimeReductionLimit?> "
    "<-implexAlpha alpha?> "
    "<-cdf cdf?> "
    "<-crackPlanes nct? ncc? smoothingAngle?> "
    "<-eta eta?> <-tangent> <-autoRegularization lch_ref?>";

  if (argc < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: " << info << "\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid 'tag'.\n";
    return TCL_ERROR;
  }

  double E;
  if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid 'E'.\n";
    return TCL_ERROR;
  }

  if (E <= 0.0) {
    opserr << OpenSees::PromptValueError << "invalid value for 'E' (" << E << "). It should be strictly positive.\n";
    return TCL_ERROR;
  }

  double v;
  if (Tcl_GetDouble(interp, argv[4], &v) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid 'v'.\n";
    return TCL_ERROR;
  }

  double rho = 0.0;
  bool implex = false;
  bool implex_control = false;
  double implex_error_tolerance = 0.05;
  double implex_time_redution_limit = 0.01;
  double implex_alpha = 1.0;
  double eta = 0.0;
  double Kc = 2.0/3.0;
  bool tangent = false;
  bool auto_regularization = false;
  double lch_ref = 1.0;
  double cdf = 0.0;
  int nct = 0;
  int ncc = 0;
  double smoothing_angle = 45.0;

  double fc = 0.0;
  double ft = 0.0;
  double Gc = 0.0;
  double Gt = 0.0;

  bool have_fc = false;
  bool have_ft = false;
  bool have_Gc = false;
  bool have_Gt = false;
  bool have_lch_ref = false;

  ConcreteBackbone backbone;

  int i = 5;
  while (i < argc) {
    const char* value = argv[i++];

    if (strcmp(value, "-rho") == 0) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-rho' given without the next argument $rho.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &rho) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'rho'.\n";
        return TCL_ERROR;
      }
      i++;
    }
    else if ((strcmp(value, "-fc") == 0) || (strcmp(value, "-Fc") == 0)) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-fc' given without the next argument $fc.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &fc) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'fc'.\n";
        return TCL_ERROR;
      }
      have_fc = true;
      i++;
    }
    else if ((strcmp(value, "-ft") == 0) || (strcmp(value, "-Ft") == 0)) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-ft' given without the next argument $ft.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &ft) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'ft'.\n";
        return TCL_ERROR;
      }
      have_ft = true;
      i++;
    }
    else if (strcmp(value, "-Gc") == 0) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-Gc' given without the next argument $Gc.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Gc) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'Gc'.\n";
        return TCL_ERROR;
      }
      have_Gc = true;
      i++;
    }
    else if (strcmp(value, "-Gt") == 0) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-Gt' given without the next argument $Gt.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Gt) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'Gt'.\n";
        return TCL_ERROR;
      }
      have_Gt = true;
      i++;
    }
    else if (strcmp(value, "-Kc") == 0) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-Kc' given without the next argument $Kc.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Kc) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'Kc'.\n";
        return TCL_ERROR;
      }
      if (Kc < 2.0/3.0 || Kc > 1.0) {
        opserr << OpenSees::PromptValueError << "'Kc' (" << Kc << ") should be >= 2/3 and <= 1.\n";
        return TCL_ERROR;
      }
      i++;
    }
    else if (strcmp(value, "-implex") == 0) {
      implex = true;
    }
    else if (strcmp(value, "-implexControl") == 0) {
      implex_control = true;

      if (argc - i < 2) {
        opserr << OpenSees::PromptValueError << "'-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i], &implex_error_tolerance) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'implexErrorTolerance'.\n";
        return TCL_ERROR;
      }
      i++;

      if (Tcl_GetDouble(interp, argv[i], &implex_time_redution_limit) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'implexTimeReductionLimit'.\n";
        return TCL_ERROR;
      }
      i++;
    }
    else if (strcmp(value, "-implexAlpha") == 0) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-implexAlpha' given without the next argument $alpha.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &implex_alpha) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'alpha'.\n";
        return TCL_ERROR;
      }
      i++;
    }
    else if (strcmp(value, "-eta") == 0) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-eta' given without the next argument $eta.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &eta) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'eta'.\n";
        return TCL_ERROR;
      }
      i++;
    }
    else if (strcmp(value, "-tangent") == 0) {
      tangent = true;
    }
    else if (strcmp(value, "-autoRegularization") == 0) {
      auto_regularization = true;

      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-autoRegularization' given without the next 1 argument $lch_ref.\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i], &lch_ref) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'lch_ref'.\n";
        return TCL_ERROR;
      }

      have_lch_ref = true;
      i++;
    }
    else if (strcmp(value, "-cdf") == 0) {
      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'-cdf' given without the next 1 argument $cdf.\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i], &cdf) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'cdf'.\n";
        return TCL_ERROR;
      }

      i++;
    }
    else if (strcmp(value, "-crackPlanes") == 0) {
      if (argc - i < 3) {
        opserr << OpenSees::PromptValueError << "'-crackPlanes' given without the next 3 arguments $nct $ncc and $smoothingAngle.\n";
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[i], &nct) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'nct'.\n";
        return TCL_ERROR;
      }
      i++;

      if (Tcl_GetInt(interp, argv[i], &ncc) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'ncc'.\n";
        return TCL_ERROR;
      }
      i++;

      if (Tcl_GetDouble(interp, argv[i], &smoothing_angle) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to get 'smoothingAngle'.\n";
        return TCL_ERROR;
      }
      i++;
    }
    else if (
      strcmp(value, "-Te") == 0 ||
      strcmp(value, "-Ts") == 0 ||
      strcmp(value, "-Td") == 0 ||
      strcmp(value, "-Ce") == 0 ||
      strcmp(value, "-Cs") == 0 ||
      strcmp(value, "-Cd") == 0)
    {
      std::vector<double>* list = nullptr;
      const char* list_name = value + 1;

      if (strcmp(value, "-Te") == 0)
        list = &backbone.Te;
      else if (strcmp(value, "-Ts") == 0)
        list = &backbone.Ts;
      else if (strcmp(value, "-Td") == 0)
        list = &backbone.Td;
      else if (strcmp(value, "-Ce") == 0)
        list = &backbone.Ce;
      else if (strcmp(value, "-Cs") == 0)
        list = &backbone.Cs;
      else
        list = &backbone.Cd;

      list->clear();

      if (i >= argc) {
        opserr << OpenSees::PromptValueError << "'" << list_name << "' requested but not provided.\n";
        return TCL_ERROR;
      }

      double item;
      if (Tcl_GetDouble(interp, argv[i], &item) == TCL_OK) {
        while (i < argc) {
          if (Tcl_GetDouble(interp, argv[i], &item) != TCL_OK) {
            Tcl_ResetResult(interp);
            break;
          }

          list->push_back(item);
          i++;
        }
      }
      else {
        Tcl_ResetResult(interp);

        TCL_Char** listArgv = nullptr;
        Tcl_Size listArgc = 0;

        if (Tcl_SplitList(interp, argv[i], &listArgc, &listArgv) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "cannot parse the '" << list_name << "' list.\n";
          return TCL_ERROR;
        }

        for (Tcl_Size j = 0; j < listArgc; j++) {
          if (Tcl_GetDouble(interp, listArgv[j], &item) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "cannot parse the '" << list_name << "' list.\n";
            Tcl_Free((char*)listArgv);
            return TCL_ERROR;
          }

          list->push_back(item);
        }

        Tcl_Free((char*)listArgv);
        i++;
      }
    }
    else {
      opserr << OpenSees::PromptValueError << "unknown option '" << value << "'.\n";
      return TCL_ERROR;
    }
  }

  if (!have_fc && (have_Gc || have_Gt)) {
    opserr << OpenSees::PromptValueError << "'-Gc' and '-Gt' require '-fc'.\n";
    return TCL_ERROR;
  }

  if (have_fc && !have_ft)
    ft = 0.1*fc;

  if (have_fc) {
    if (!have_Gt)
      Gt = 0.073*pow(fc, 0.18);

    if (!have_Gc)
      Gc = 2.0*Gt*(fc*fc)/(ft*ft);

    if (CreateConcreteBackbone(backbone, E, fc, ft, Gc, Gt, have_lch_ref, lch_ref) != TCL_OK)
      return TCL_ERROR;
  }

  if (backbone.Te.size() < 1) {
    opserr << OpenSees::PromptValueError << "'Te' list is empty. At least 1 non-zero value should be provided.\n";
    return TCL_ERROR;
  }

  if (backbone.Ts.size() != backbone.Te.size()) {
    opserr << OpenSees::PromptValueError << "'Te' (size = " <<
      static_cast<int>(backbone.Te.size()) << ") and 'Ts' (size = " <<
      static_cast<int>(backbone.Ts.size()) << ") lists should have the same size.\n";
    return TCL_ERROR;
  }

  if (backbone.Td.size() == 0) {
    backbone.Td.resize(backbone.Te.size(), 0.0);
  }
  else if (backbone.Td.size() != backbone.Te.size()) {
    opserr << OpenSees::PromptValueError << "'Te' (size = " <<
      static_cast<int>(backbone.Te.size()) << ") and 'Td' (size = " <<
      static_cast<int>(backbone.Td.size()) << ") lists should have the same size.\n";
    return TCL_ERROR;
  }

  if (backbone.Ce.size() < 1) {
    opserr << OpenSees::PromptValueError << "'Ce' list is empty. At least 1 non-zero value should be provided.\n";
    return TCL_ERROR;
  }

  if (backbone.Cs.size() != backbone.Ce.size()) {
    opserr << OpenSees::PromptValueError << "'Ce' (size = " <<
      static_cast<int>(backbone.Ce.size()) << ") and 'Cs' (size = " <<
      static_cast<int>(backbone.Cs.size()) << ") lists should have the same size.\n";
    return TCL_ERROR;
  }

  if (backbone.Cd.size() == 0) {
    backbone.Cd.resize(backbone.Ce.size(), 0.0);
  }
  else if (backbone.Cd.size() != backbone.Ce.size()) {
    opserr << OpenSees::PromptValueError << "'Ce' (size = " 
           << static_cast<int>(backbone.Ce.size()) << ") and 'Cd' (size = " 
           << static_cast<int>(backbone.Cd.size()) << ") lists should have the same size.\n";
    return TCL_ERROR;
  }

  ASDConcrete3DMaterial::HardeningLaw HT(
    tag,
    ASDConcrete3DMaterial::HardeningLawType::Tension,
    E,
    backbone.Te,
    backbone.Ts,
    backbone.Td);

  if (!HT.isValid()) {
    opserr << OpenSees::PromptValueError << "Tensile hardening law is not valid.\n";
    return TCL_ERROR;
  }

  ASDConcrete3DMaterial::HardeningLaw HC(
    tag,
    ASDConcrete3DMaterial::HardeningLawType::Compression,
    E,
    backbone.Ce,
    backbone.Cs,
    backbone.Cd);

  if (!HC.isValid()) {
    opserr << OpenSees::PromptValueError << "Compressive hardening law is not valid.\n";
    return TCL_ERROR;
  }

  NDMaterial* instance = new ASDConcrete3DMaterial(
    tag,
    E, v, rho, eta, Kc,
    implex, implex_control, implex_error_tolerance, implex_time_redution_limit, implex_alpha,
    tangent, auto_regularization, lch_ref,
    HT, HC,
    cdf, nct, ncc, smoothing_angle);

  if (builder->addTaggedObject<NDMaterial>(*instance) != TCL_OK) {
    delete instance;
    opserr << OpenSees::PromptValueError << "could not add material with tag " << tag << ".\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}