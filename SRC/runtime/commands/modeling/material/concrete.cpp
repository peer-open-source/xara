//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
// April 2025
//
#include <tcl.h>
#include <set>
#include <string.h>
#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include <ModelRegistry.h>
#include "isotropy.h"

#include <damage/FariaPlasticDamage3d.h>

int
TclCommand_newConcreteMaterial(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);
  enum class Position {
    Tag, E, Nu, PeakTension, PeakCompression, EndRequired,
    Beta, Ap, An, Bn,
    Density,
    G, K, Lambda,
    End,
  };
  ArgumentTracker<Position> tracker;
  std::set<int> positional;


  // Count the number of required isotropic parameters.
  // This is needed to accommodate uniaxial materials, which generally
  // only require one isotropic parameter (ie, E).
  int niso = (
    (Position::E      < Position::EndRequired) +
    (Position::G      < Position::EndRequired) +
    (Position::Nu     < Position::EndRequired) +
    (Position::K      < Position::EndRequired) +
    (Position::Lambda < Position::EndRequired)
  );

  // 
  // Values we're parsing for
  //
  int tag;
  double density = 0.0;
  // Isotropy
  IsotropicConstants consts {};
  // Plasticity
  double Fc, Ft=0;
  // Hardening
  double beta = 0.6;
  double Ap = 0.5,
         An = 2.0,
         Bn = 0.75;
  

  //
  // 1. Keyword arguments
  //

  // Isotropy
  IsotropicParse iso {consts, niso};
  if (TclCommand_setIsotropicParameters((ClientData)&iso, interp, argc, argv) == TCL_OK) {
    tracker.consume(Position::E);
    tracker.consume(Position::G);
    tracker.consume(Position::Nu);
    tracker.consume(Position::K);
    tracker.consume(Position::Lambda);
  }

  // Other arguments
  for (int i=2; i<argc; i++) {
    if (iso.positions.find(i) != iso.positions.end()) {
      continue;
    }

    else if (strcmp(argv[i], "-rho") == 0 || strcmp(argv[i], "-density") == 0) {
      if (++i >= argc) {
          opserr << OpenSees::PromptValueError
                 << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &density) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "Invalid density value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
    }

    // Compression
    else if (strcmp(argv[i], "-Fc") == 0 || 
             strcmp(argv[i], "-fc") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
               << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Fc) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "Invalid " << &argv[i-1][1] << " value " << argv[i] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::PeakCompression);
    }
    // Tension
    else if ((strcmp(argv[i], "-Ft") == 0)) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
               << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Ft) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "Invalid " << &argv[i-1][1] << " value " << argv[i] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::PeakTension);
    }
    //
    // Hardening
    //
    else if (strcmp(argv[i], "-beta") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
               << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &beta) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Beta);
    }
    //
    //
    else if (strcmp(argv[i], "-Ap") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
               << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Ap) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Ap);
    }
    else if (strcmp(argv[i], "-An") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
               << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &An) != TCL_OK) {
        opserr << OpenSees::PromptValueError
                << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::An);
    }
    else if (strcmp(argv[i], "-Bn") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
                << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Bn) != TCL_OK) {
        opserr << OpenSees::PromptValueError
                << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Bn);
    }
    else
      positional.insert(i);
  }

  //
  // 2) Positional arguments
  //
  for (int i : positional) {
  
    if (tracker.current() == Position::EndRequired)
      tracker.increment();

    switch (tracker.current()) {
      // General
      case Position::Tag :
        if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid section Elastic tag.\n";
            return TCL_ERROR;           
        } else {
          tracker.increment();
          break;
        }
      case Position::Density:
        if (Tcl_GetDouble (interp, argv[i], &density) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid density.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      // Isotropy
      case Position::E:
        if (Tcl_GetDouble (interp, argv[i], &consts.E) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid E.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      case Position::G:
        if (Tcl_GetDouble (interp, argv[i], &consts.G) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid G.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      case Position::K:
        if (Tcl_GetDouble (interp, argv[i], &consts.K) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid K.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      case Position::Nu:
        if (Tcl_GetDouble (interp, argv[i], &consts.nu) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid nu.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      
      case Position::Lambda:
        if (Tcl_GetDouble (interp, argv[i], &consts.lambda) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Lame lambda.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
    
      // Yielding
      case Position::PeakCompression:
        if (Tcl_GetDouble(interp, argv[i], &Fc) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Fc.\n";
            return TCL_ERROR;           
        } else {
          tracker.increment();
          break;
        }

      case Position::PeakTension:
        if (Tcl_GetDouble (interp, argv[i], &Ft) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Ft.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Beta:
        if (Tcl_GetDouble (interp, argv[i], &beta) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid beta.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Ap:
        if (Tcl_GetDouble (interp, argv[i], &Ap) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Ap.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::An:
        if (Tcl_GetDouble (interp, argv[i], &An) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid An.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Bn:
        if (Tcl_GetDouble (interp, argv[i], &Bn) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Bn.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
    

      case Position::EndRequired:
        // This will not be reached
        break;

      case Position::End:
        opserr << OpenSees::PromptParseError << "unexpected argument " << argv[i] << ".\n";
        return TCL_ERROR;
    }
  }


  //
  // 3. Check for required arguments
  //
  if (tracker.current() < Position::EndRequired) {
    opserr << OpenSees::PromptParseError
            << "missing required arguments: ";
    while (tracker.current() != Position::EndRequired) {
      switch (tracker.current()) {
        case Position::Tag :
          opserr << "tag ";
          break;
        // Isotropy
        case Position::E:
          opserr << "E ";
          break;
        case Position::G:
          opserr << "G ";
          break;
        case Position::K:
          opserr << "K ";
          break;
        case Position::Nu:
          opserr << "nu ";
          break;
        case Position::Lambda:
          opserr << "lambda ";
          break;
        // Yielding
        case Position::PeakCompression:
          opserr << "Fc ";
          break;
        case Position::PeakTension:
          opserr << "Ft ";
          break;
        case Position::Ap:
          opserr << "Ap ";
          break;
        // Hardening
        case Position::An:
          opserr << "An ";
          break;
        case Position::Bn:
          opserr << "Bn ";
          break;
        case Position::Beta:
          opserr << "beta ";
          break;

        case Position::EndRequired:
        case Position::End:
        default:
          break;
      }

      if (tracker.current() == Position::EndRequired)
        break;

      tracker.consume(tracker.current());
    }

    opserr << OpenSees::SignalMessageEnd;

    return TCL_ERROR;
  }

  //
  // Create the material (TODO)
  //
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);
  if ((strcmp(argv[1], "PlasticDamageConcrete3d") == 0) ||
      (strcasecmp(argv[1], "PlasticDamageConcrete") == 0) ||
      (strcmp(argv[1], "FariaPlasticDamage") == 0)) {

    NDMaterial* theMaterial = new FariaPlasticDamage3d(tag, consts.E, consts.nu, Ft, Fc, 
                                                        beta, Ap, An, Bn, density);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  return TCL_ERROR;
}


#include <elementAPI.h>
#include <vector>
#include <uniaxial/ASDConcrete1DMaterial.h>
#include <damage/ASDConcrete3DMaterial.h>

// anonymous namespace for utilities
namespace {

	/**
	Converts a string into a vector of doubles using whitespace as delimiter
	*/
	bool string_to_double(const std::string& text, double& num) {
		num = 0.0;
		try {
			num = std::stod(text);
			return true;
		}
		catch (...) {
			return false;
		}
	}
	bool string_to_list_of_doubles(const std::string& text, char sep, std::vector<double>& out) {
		if (out.size() > 0) out.clear();
		std::size_t start = 0, end = 0;
		double value;
		while (true) {
			end = text.find(sep, start);
			if (end == std::string::npos) {
				if (start < text.size()) {
					if (!string_to_double(text.substr(start), value))
						return false;
					out.push_back(value);
				}
				break;
			}
			std::string subs = text.substr(start, end - start);
			if (subs.size() > 0) {
				if (!string_to_double(subs, value))
					return false;
				out.push_back(value);
			}
			start = end + 1;
		}
		return true;
	}


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

}




void* OPS_ADD_RUNTIME_VPV(OPS_ASDConcrete1DMaterial)
{
	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opslog << "Using ASDConcrete1D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}

	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 2) {
		opserr <<
			"nDMaterial ASDConcrete1D Error: Few arguments (< 2).\n"
			"nDMaterial ASDConcrete1D $tag $E "
			"<-fc $fc> <-ft $ft> "
			"<-Te $Te -Ts $Ts <-Td $Td>> <-Ce $Ce -Cs $Cs <-Cd $Cd>> "
			"<-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> <-implexAlpha $alpha>"
			"<-eta $eta> <-tangent> <-autoRegularization $lch_ref>\n";
		return nullptr;
	}

	// numData
	int numData = 1;

	// data
	int tag;
	double E;
	bool implex = false;
	bool implex_control = false;
	double implex_error_tolerance = 0.05;
	double implex_time_redution_limit = 0.01;
	double implex_alpha = 1.0;
	double eta = 0.0;
	bool tangent = false;
	bool auto_regularization = false;
	double lch_ref = 1.0;
	std::vector<double> Te, Ts, Td, Ce, Cs, Cd;

	// get tag
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "nDMaterial ASDConcrete1D Error: invalid 'tag'.\n";
		return nullptr;
	}

	// get Elasticity arguments
	if (OPS_GetDouble(&numData, &E) != 0) {
		opserr << "nDMaterial ASDConcrete1D Error: invalid 'E'.\n";
		return nullptr;
	}
	if (E <= 0.0) {
		opserr << "nDMaterial ASDConcrete1D Error: invalid value for 'E' (" << E << "). It should be strictly positive.\n";
		return nullptr;
	}

  //
	// utilities
  //
	auto lam_optional_double = [&numData](const char* variable, double& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetDouble(&numData, &value) < 0) {
				opserr << "nDMaterial ASDConcrete1D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDConcrete1D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};

	auto lam_optional_list = [&numData](const char* variable, std::vector<double>& value) -> bool {
		// first try expanded list like {*}$the_list,
		// also used in python like *the_list
		value.clear();
		while (OPS_GetNumRemainingInputArgs() > 0) {
			double item;
			auto old_num_rem = OPS_GetNumRemainingInputArgs();
			if (OPS_GetDoubleInput(&numData, &item) < 0) {
				auto new_num_rem = OPS_GetNumRemainingInputArgs();
				if (new_num_rem < old_num_rem)
					OPS_ResetCurrentInputArg(-1);
				break;
			}
			value.push_back(item);
		}
		// try Tcl list (it's a string after all...)
		if (value.size() == 0 && OPS_GetNumRemainingInputArgs() > 0) {
			std::string list_string = OPS_GetString();
			if (!string_to_list_of_doubles(list_string, ' ', value)) {
				opserr << "nDMaterial ASDConcrete1D Error: cannot parse the '" << variable << "' list.\n";
				return false;
			}
		}
		return true;
	};

	double fc;
	double ft;
	bool have_fc = false;
	bool have_ft = false;
	bool have_lch_ref = false;

	// optional parameters
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-fc") == 0) {
			if (!lam_optional_double("fc", fc))
				return nullptr;
			have_fc = true;
		}
		else if (strcmp(value, "-ft") == 0) {
			if (!lam_optional_double("ft", ft))
				return nullptr;
			have_ft = true;
		}
		else if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		else if (strcmp(value, "-implexControl") == 0) {
			implex_control = true;
			if (OPS_GetNumRemainingInputArgs() < 2) {
				opserr << "nDMaterial ASDConcrete1D Error: '-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
				return nullptr;
			}
			if (!lam_optional_double("implexErrorTolerance", implex_error_tolerance))
				return nullptr;
			if (!lam_optional_double("implexTimeReductionLimit", implex_time_redution_limit))
				return nullptr;
		}
		else if (strcmp(value, "-implexAlpha") == 0) {
			if (!lam_optional_double("alpha", implex_alpha))
				return nullptr;
		}
		else if (strcmp(value, "-eta") == 0) {
			if (!lam_optional_double("eta", eta))
				return nullptr;
		}
		else if (strcmp(value, "-tangent") == 0) {
			tangent = true;
		}
		else if (strcmp(value, "-autoRegularization") == 0) {
			auto_regularization = true;
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "nDMaterial ASDConcrete1D Error: '-autoRegularization' given without the next 1 argument $lch_ref.\n";
				return nullptr;
			}
			if (!lam_optional_double("lch_ref", lch_ref))
				return nullptr;
			have_lch_ref = true;
		}
		else if (strcmp(value, "-Te") == 0) {
			if (!lam_optional_list("Te", Te))
				return nullptr;
		}
		else if (strcmp(value, "-Ts") == 0) {
			if (!lam_optional_list("Ts", Ts))
				return nullptr;
		}
		else if (strcmp(value, "-Td") == 0) {
			if (!lam_optional_list("Td", Td))
				return nullptr;
		}
		else if (strcmp(value, "-Ce") == 0) {
			if (!lam_optional_list("Ce", Ce))
				return nullptr;
		}
		else if (strcmp(value, "-Cs") == 0) {
			if (!lam_optional_list("Cs", Cs))
				return nullptr;
		}
		else if (strcmp(value, "-Cd") == 0) {
			if (!lam_optional_list("Cd", Cd))
				return nullptr;
		}
	}

	// Set a default value of tension strength if none specified
	if (have_fc && !have_ft)
		ft = 0.1 * fc;

	if (have_fc) {
		double ec = 2 * fc / E;
		double Gt = 0.073 * pow(fc, 0.18);
		double Gc = 2 * Gt * (fc * fc) / (ft * ft);


		if (!have_lch_ref) {
			//
			// _get_lch_ref from ASDConcrete1D_MakeLaws.py
			//

			// min lch for tension
			double et_el = ft / E;
			double Gt_min = 0.5 * ft * et_el;
			double hmin_t = 0.01 * Gt / Gt_min;

			// min lch for compression
			double ec1 = fc / E;
			double ec_pl = (ec - ec1) * 0.4 + ec1;
			double Gc_min = 0.5 * fc * (ec - ec_pl);
			double hmin_c = 0.01 * Gc / Gc_min;

			lch_ref = std::min(hmin_c, hmin_t);
		}

		//
		// _make_tension from ASDConcrete1D_MakeLaws.py
		//

		Gt = Gt / lch_ref;

		double f0 = 0.9 * ft;
		double f1 = ft;
		double e0 = f0 / E;
		double e1 = 1.5 * f1 / E;
		double ep = e1 - f1 / E;
		double f2 = 0.2 * ft;
		double f3 = 1.0e-3 * ft;
		double w2 = Gt / ft;
		double w3 = 5.0 * w2;
		double e2 = w2 + f2 / E + ep;
		if (e2 <= e1)
			e2 = 1.001 * e1;
		double e3 = w3 + f3 / E + ep;
		if (e3 <= e2)
			e3 = 1.001 * e2;
		double e4 = 10.0 * e3;
		Te.resize(6); Te = { 0.0, e0, e1, e2, e3, e4 };
		Ts.resize(6); Ts = { 0.0, f0, f1, f2, f3, f3 };
		Td.resize(6); Td = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		double Tpl[6] = { 0.0, 0.0, ep, 0.9 * e2, 0.8 * e3, 0.8 * e3 };

		for (int i = 2; i < 6; i++) {
			double xi = Te[i];
			double si = Ts[i];
			double xipl = Tpl[i];
			double xipl_max = xi - si / E;
			xipl = std::min(xipl, xipl_max);
			double qi = (xi - xipl) * E;
			Td[i] = 1.0 - si / qi;
		}

		//
		// _make_compression from ASDConcrete1D_MakeLaws.py
		//

		Gc = Gc / lch_ref;

		double fc0 = 0.5 * fc;
		double ec0 = fc0 / E;
		double ec1 = fc / E;
		double fcr = 0.1 * fc;
		double ec_pl = (ec - ec1) * 0.4 + ec1;
		double Gc1 = 0.5 * fc * (ec - ec_pl);
		double Gc2 = std::max(0.01 * Gc1, Gc - Gc1);
		double ecr = ec + 2.0 * Gc2 / (fc + fcr);
		const int nc = 10;
		Ce.resize(nc + 3); Ce[0] = 0.0; Ce[1] = ec0;
		Cs.resize(nc + 3); Cs[0] = 0.0; Cs[1] = fc0;
		double Cpl[nc + 3]; Cpl[0] = 0.0; Cpl[1] = 0.0;
		double dec = (ec - ec0) / (nc - 1);
		for (int i = 0; i < nc - 1; i++) {
			double iec = ec0 + (i + 1) * dec;
			Ce[i + 2] = iec;
			Cs[i + 2] = bezier3(iec, ec0, ec1, ec, fc0, fc, fc);
			Cpl[i + 2] = Cpl[i + 1] + 0.7 * (iec - Cpl[i + 1]);
		}
		Ce[nc + 1] = ecr;
		Cs[nc + 1] = fcr;
		Cpl[nc + 1] = Cpl[nc] + 0.7 * (ecr - Cpl[nc]);
		Ce[nc + 2] = ecr + ec0;
		Cs[nc + 2] = fcr;
		Cpl[nc + 2] = Cpl[nc + 1];
		Cd.resize(nc + 3); Cd[0] = 0.0; Cd[1] = 0.0;
		for (int i = 2; i < nc + 3; i++) {
			double xi = Ce[i];
			double si = Cs[i];
			double xipl = Cpl[i];
			double xipl_max = xi - si / E;
			xipl = std::min(xipl, xipl_max);
			double qi = (xi - xipl) * E;
			Cd[i] = 1.0 - si / qi;
		}
	}

	// check lists
	if (Te.size() < 1) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Te' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Ts.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Ts' (size = " <<
			static_cast<int>(Ts.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Td.size() == 0) {
		Td.resize(Te.size(), 0.0);
	}
	else if (Td.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Td' (size = " <<
			static_cast<int>(Td.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Ce.size() < 1) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Tc' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Cs.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cs' (size = " <<
			static_cast<int>(Cs.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Cd.size() == 0) {
		Cd.resize(Ce.size(), 0.0);
	}
	else if (Cd.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cd' (size = " <<
			static_cast<int>(Cd.size()) << ") lists should have the same size.\n";
		return nullptr;
	}

	// build the hardening laws
	ASDConcrete1DMaterial::HardeningLaw HT(tag, ASDConcrete1DMaterial::HardeningLawType::Tension, E, Te, Ts, Td);
	if (!HT.isValid()) {
		opserr << "nDMaterial ASDConcrete1D Error: Tensile hardening law is not valid.\n";
		return nullptr;
	}
	ASDConcrete1DMaterial::HardeningLaw HC(tag, ASDConcrete1DMaterial::HardeningLawType::Compression, E, Ce, Cs, Cd);
	if (!HC.isValid()) {
		opserr << "nDMaterial ASDConcrete1D Error: Compressive hardening law is not valid.\n";
		return nullptr;
	}

	// create the material
	UniaxialMaterial* instance = new ASDConcrete1DMaterial(
		tag,
		E, eta,
		implex, implex_control, implex_error_tolerance, implex_time_redution_limit, implex_alpha,
		tangent, auto_regularization, lch_ref,
		HT, HC);
	return instance;
}


void *OPS_ADD_RUNTIME_VPV(OPS_ASDConcrete3DMaterial_V2)
{
	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opslog << "Using ASDConcrete3D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}

	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 3) {
		opserr << 
			"nDMaterial ASDConcrete3D Error: Few arguments (< 3).\n"
			"nDMaterial ASDConcrete3D $tag $E $v "
			"-Te $Te -Ts $Ts <-Td $Td> -Ce $Ce -Cs $Cs <-Cd $Cd> "
			"<-rho $rho> <-Kc $Kc>"
			"<-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> <-implexAlpha $alpha>"
			"<-cdf $cdf>"
			"<-crackPlanes $nct $ncc $smoothingAngle>"
			"<-eta $eta> <-tangent> <-autoRegularization $lch_ref>\n";
		return nullptr;
	}

	// numData
	int numData = 1;

	// data
	int tag;
	double E;
	double v;
	double rho = 0.0;
	bool implex = false;
	bool implex_control = false;
	double implex_error_tolerance = 0.05;
	double implex_time_redution_limit = 0.01;
	double implex_alpha = 1.0;
	double eta = 0.0;
	double Kc = 2.0 / 3.0; // default suggested by Lubliner et al.
	bool tangent = false;
	bool auto_regularization = false;
	double lch_ref = 1.0;
	std::vector<double> Te, Ts, Td, Ce, Cs, Cd;
	double cdf = 0.0;
	int nct = 0;
	int ncc = 0;
	double smoothing_angle = 45.0;
	
	// get tag
	if (OPS_GetInt(&numData, &tag) != 0)  {
		opserr << "nDMaterial ASDConcrete3D Error: invalid 'tag'.\n";
		return nullptr;
	}

	// get Elasticity arguments
	if (OPS_GetDouble(&numData, &E) != 0) {
		opserr << "nDMaterial ASDConcrete3D Error: invalid 'E'.\n";
		return nullptr;
	}
	if (E <= 0.0) {
		opserr << "nDMaterial ASDConcrete3D Error: invalid value for 'E' (" << E << "). It should be strictly positive.\n";
		return nullptr;
	}
	if (OPS_GetDouble(&numData, &v) != 0) {
		opserr << "nDMaterial ASDConcrete3D Error: invalid 'v'.\n";
		return nullptr;
	}

	// utilities (code re-use)
	auto lam_optional_int = [&numData](const char* variable, int& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetInt(&numData, &value) < 0) {
				opserr << "nDMaterial ASDConcrete3D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDConcrete3D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	auto lam_optional_double = [&numData](const char* variable, double& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetDouble(&numData, &value) < 0) {
				opserr << "nDMaterial ASDConcrete3D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDConcrete3D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	auto lam_optional_list = [&numData](const char* variable, std::vector<double>& value) -> bool {
		// first try expanded list like {*}$the_list,
		// also used in python like *the_list
		value.clear();
		while (OPS_GetNumRemainingInputArgs() > 0) {
			double item;
			auto old_num_rem = OPS_GetNumRemainingInputArgs();
			if (OPS_GetDoubleInput(&numData, &item) < 0) {
				auto new_num_rem = OPS_GetNumRemainingInputArgs();
				if (new_num_rem < old_num_rem)
					OPS_ResetCurrentInputArg(-1);
				break;
			}
			value.push_back(item);
		}
		// try Tcl list (it's a string after all...)
		if (value.size() == 0 && OPS_GetNumRemainingInputArgs() > 0) {
			std::string list_string = OPS_GetString();
			if (!string_to_list_of_doubles(list_string, ' ', value)) {
				opserr << "nDMaterial ASDConcrete3D Error: cannot parse the '" << variable << "' list.\n";
				return false;
			}
		}
		return true;
	};

	double fc;
	double ft;
  double density = 0.0;
	bool have_fc = false;
	bool have_ft = false;	
	bool have_lch_ref = false;
	
	// optional parameters
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-rho") == 0) {
			if (!lam_optional_double("rho", rho))
				return nullptr;
		}
		else if (strcmp(value, "-fc") == 0) {
			if (!lam_optional_double("fc", fc))
				return nullptr;
			have_fc = true;
		}
		else if (strcmp(value, "-ft") == 0) {
			if (!lam_optional_double("ft", ft))
				return nullptr;
			have_ft = true;
		}		
		else if (strcmp(value, "-Kc") == 0) {
			if (!lam_optional_double("Kc", Kc))
				return nullptr;
			if (Kc < 2.0 / 3.0 || Kc > 1.0) {
				opserr << "nDMaterial ASDConcrete3D Error: 'Kc' (" << Kc << ") double be >= 2/3 and <= 1.\n";
				return nullptr;
			}
		}
		else if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		else if (strcmp(value, "-implexControl") == 0) {
			implex_control = true;
			if (OPS_GetNumRemainingInputArgs() < 2) {
				opserr << "nDMaterial ASDConcrete3D Error: '-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
				return nullptr;
			}
			if (!lam_optional_double("implexErrorTolerance", implex_error_tolerance))
				return nullptr;
			if (!lam_optional_double("implexTimeReductionLimit", implex_time_redution_limit))
				return nullptr;
		}
		else if (strcmp(value, "-implexAlpha") == 0) {
			if (!lam_optional_double("alpha", implex_alpha))
				return nullptr;
		}
		else if (strcmp(value, "-eta") == 0) {
			if (!lam_optional_double("eta", eta))
				return nullptr;
		}
		else if (strcmp(value, "-tangent") == 0) {
			tangent = true;
		}
		else if (strcmp(value, "-autoRegularization") == 0) {
			auto_regularization = true;
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "nDMaterial ASDConcrete3D Error: '-autoRegularization' given without the next 1 argument $lch_ref.\n";
				return nullptr;
			}
			if (!lam_optional_double("lch_ref", lch_ref))
				return nullptr;
			have_lch_ref = true;
		}
		else if (strcmp(value, "-Te") == 0) {
			if (!lam_optional_list("Te", Te))
				return nullptr;
		}
		else if (strcmp(value, "-Ts") == 0) {
			if (!lam_optional_list("Ts", Ts))
				return nullptr;
		}
		else if (strcmp(value, "-Td") == 0) {
			if (!lam_optional_list("Td", Td))
				return nullptr;
		}
		else if (strcmp(value, "-Ce") == 0) {
			if (!lam_optional_list("Ce", Ce))
				return nullptr;
		}
		else if (strcmp(value, "-Cs") == 0) {
			if (!lam_optional_list("Cs", Cs))
				return nullptr;
		}
		else if (strcmp(value, "-Cd") == 0) {
			if (!lam_optional_list("Cd", Cd))
				return nullptr;
		}
		else if (strcmp(value, "-crackPlanes") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 3) {
				opserr << "nDMaterial ASDConcrete3D Error: '-crackPlanes' given without the next 3 arguments $nct $ncc and $smoothingAngle.\n";
				return nullptr;
			}
			if (!lam_optional_int("nct", nct))
				return nullptr;
			if (!lam_optional_int("ncc", ncc))
				return nullptr;
			if (!lam_optional_double("smoothingAngle", smoothing_angle))
				return nullptr;
		}
		else if (strcmp(value, "-cdf") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "nDMaterial ASDConcrete3D Error: '-cdf' given without the next 1 argument $cdf.\n";
				return nullptr;
			}
			if (!lam_optional_double("cdf", cdf))
				return nullptr;
		}
	}

	// Set a default value of tension strength if none specified
	if (have_fc && !have_ft)
	  ft = 0.1*fc;

	if (have_fc) {
	  double ec = 2*fc/E;
	  double Gt = 0.073*pow(fc,0.18);
	  double Gc = 2*Gt*(fc*fc)/(ft*ft);
	  
	  
	  if (!have_lch_ref) {
	    //
	    // _get_lch_ref from ASDConcrete3D_MakeLaws.py
	    //
	    
	    // min lch for tension
	    double et_el = ft/E;
	    double Gt_min = 0.5*ft*et_el;
	    double hmin_t = 0.01*Gt/Gt_min;
	    
	    // min lch for compression
	    double ec1 = fc/E;
	    double ec_pl = (ec-ec1)*0.4 + ec1;
	    double Gc_min = 0.5*fc*(ec-ec_pl);
	    double hmin_c = 0.01*Gc/Gc_min;
	    
	    lch_ref = std::min(hmin_c,hmin_t);
	  }
	  
	  //
	  // _make_tension from ASDConcrete3D_MakeLaws.py
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
	  Te.resize(6); Te = {0.0, e0, e1, e2, e3, e4};
	  Ts.resize(6); Ts = {0.0, f0, f1, f2, f3, f3};
	  Td.resize(6); Td = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	  double Tpl[6] = {0.0, 0.0, ep, 0.9*e2, 0.8*e3, 0.8*e3};
	  
	  for (int i = 2; i < 6; i++) {
	    double xi = Te[i];
	    double si = Ts[i];
	    double xipl = Tpl[i];
	    double xipl_max = xi-si/E;
	    xipl = std::min(xipl, xipl_max);
	    double qi = (xi-xipl)*E;
	    Td[i] = 1.0 - si/qi;
	  }

	  
	  //
	  // _make_compression from ASDConcrete3D_MakeLaws.py
	  //

	  Gc = Gc/lch_ref;
	  
	  double fc0 = 0.5*fc;
	  double ec0 = fc0/E;
	  double ec1 = fc/E;
	  double fcr = 0.1*fc;
	  double ec_pl = (ec-ec1)*0.4 + ec1;
	  double Gc1 = 0.5*fc*(ec-ec_pl);
	  double Gc2 = std::max(0.01*Gc1,Gc-Gc1);
	  double ecr = ec + 2.0*Gc2/(fc+fcr);
	  const int nc = 10;
	  Ce.resize(nc+3); Ce[0] = 0.0; Ce[1] = ec0;
	  Cs.resize(nc+3); Cs[0] = 0.0; Cs[1] = fc0;
	  double Cpl[nc+3]; Cpl[0] = 0.0; Cpl[1] = 0.0;
	  double dec = (ec-ec0)/(nc-1);
	  for (int i = 0; i < nc-1; i++) {
	    double iec = ec0 + (i+1)*dec;
	    Ce[i+2] = iec;
	    Cs[i+2] = bezier3(iec,  ec0, ec1, ec,  fc0, fc, fc);
	    Cpl[i+2] = Cpl[i+1] + 0.7*(iec-Cpl[i+1]);
	  }
	  Ce[nc+1] = ecr;
	  Cs[nc+1] = fcr;
	  Cpl[nc+1] = Cpl[nc] + 0.7*(ecr-Cpl[nc]);
	  Ce[nc+2] = ecr + ec0;
	  Cs[nc+2] = fcr;
	  Cpl[nc+2] = Cpl[nc+1];
	  Cd.resize(nc+3); Cd[0] = 0.0; Cd[1] = 0.0;
	  for (int i = 2; i < nc+3; i++) {
	    double xi = Ce[i];
	    double si = Cs[i];
	    double xipl = Cpl[i];
	    double xipl_max = xi-si/E;
	    xipl = std::min(xipl, xipl_max);
	    double qi = (xi-xipl)*E;
	    Cd[i] = 1.0 - si/qi;
	  }
	}
	
	// check lists
	if (Te.size() < 1) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Te' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Ts.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Ts' (size = " <<
			static_cast<int>(Ts.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Td.size() == 0) {
		Td.resize(Te.size(), 0.0);
	}
	else if (Td.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Td' (size = " <<
			static_cast<int>(Td.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Ce.size() < 1) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Tc' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Cs.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cs' (size = " <<
			static_cast<int>(Cs.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Cd.size() == 0) {
		Cd.resize(Ce.size(), 0.0);
	}
	else if (Cd.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cd' (size = " <<
			static_cast<int>(Cd.size()) << ") lists should have the same size.\n";
		return nullptr;
	}

	// build the hardening laws
	ASDConcrete3DMaterial::HardeningLaw HT(tag, ASDConcrete3DMaterial::HardeningLawType::Tension, E, Te, Ts, Td);
	if (!HT.isValid()) {
		opserr << "nDMaterial ASDConcrete3D Error: Tensile hardening law is not valid.\n";
		return nullptr;
	}
	ASDConcrete3DMaterial::HardeningLaw HC(tag, ASDConcrete3DMaterial::HardeningLawType::Compression, E, Ce, Cs, Cd);
	if (!HC.isValid()) {
		opserr << "nDMaterial ASDConcrete3D Error: Compressive hardening law is not valid.\n";
		return nullptr;
	}

	// create the material
	return new ASDConcrete3DMaterial(
		tag, 
		E, v, rho, eta, Kc,
		implex, implex_control, implex_error_tolerance, implex_time_redution_limit, implex_alpha,
		tangent, auto_regularization, lch_ref,
		HT, HC,
		cdf, nct, ncc, smoothing_angle);
}
