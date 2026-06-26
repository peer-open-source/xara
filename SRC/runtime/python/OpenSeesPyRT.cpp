//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
// Author: cmp
//
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <G3_Runtime.h>
#include <elementAPI.h> // G3_getRuntime/SafeBuilder
#include <runtime/runtime/ModelRegistry.h>

#include <Domain.h>
#include <Vector.h>
#include <Node.h>
#include <NodeData.h>
#include <Element.h>
#include <SectionForceDeformation.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>


#define ARRAY_FLAGS py::array::c_style|py::array::forcecast


#if 0
std::unique_ptr<G3_Runtime, py::nodelete> 
getRuntime(py::object interpaddr) {
      void *interp_addr;
      interp_addr = (void*)PyLong_AsVoidPtr(interpaddr.ptr());
      return std::unique_ptr<G3_Runtime, py::nodelete>(G3_getRuntime((Tcl_Interp*)interp_addr));
} // , py::return_value_policy::reference
#endif


std::unique_ptr<ModelRegistry, py::nodelete> 
get_builder(py::object interpaddr) {
    void *interp_addr;
    interp_addr = PyLong_AsVoidPtr(interpaddr.ptr());

    Tcl_Interp* interp = (Tcl_Interp*)interp_addr;
    Tcl_InitStubs(interp, "8.6", 0);

    Tcl_CmdInfo info;
    void *builder_addr = nullptr;
    if (Tcl_GetCommandInfo(interp, "uniaxialMaterial", &info)==1) {
      builder_addr = info.clientData;
    } else {
      return nullptr;
    }

//    void *builder_addr = G3_getSafeBuilder(G3_getRuntime((Tcl_Interp*)interp_addr));
    return std::unique_ptr<ModelRegistry, py::nodelete>((ModelRegistry*)builder_addr);
} // , py::return_value_policy::reference

class Channel;
class FEM_ObjectBroker;

class PyUniaxialMaterial : public UniaxialMaterial {
public:
  ~PyUniaxialMaterial()
  {

  }
  PyUniaxialMaterial(const UniaxialMaterial& m) :
    // m_object(py::cast(this)),
    UniaxialMaterial(m.getTag(), 11) 
  {

  }
  PyUniaxialMaterial(const py::object &obj, int tag) :
    // m_object(obj), 
    UniaxialMaterial(tag,10)
  { 

  }

  void Print(OPS_Stream &s, int flag) {
    // TODO:
    s << "" << "\n";
  }

  /* Trampoline (need one for each virtual function) */
  double getStress() override {
      PYBIND11_OVERRIDE_PURE(
          double, 
          UniaxialMaterial, 
          getStress,
      );
  }

  double getStrain() override {
    PYBIND11_OVERRIDE_PURE(double, UniaxialMaterial, getStrain);
  }
  double getTangent() override {
    PYBIND11_OVERRIDE_PURE(double, UniaxialMaterial, getTangent);
  }
  double getInitialTangent() override {
    PYBIND11_OVERRIDE_PURE(double, UniaxialMaterial, getInitialTangent);
  }
  int commitState() override {
      PYBIND11_OVERRIDE_PURE(
          int,                      // Return type
          UniaxialMaterial,         // Parent class
          commitState               // Name of function in C++ (must match Python name)
      );
  }
  int revertToLastCommit() override {
      PYBIND11_OVERRIDE_PURE(
          int,                      /* Return type */
          UniaxialMaterial,         /* Parent class */
          revertToLastCommit        /* Name of function in C++ (must match Python name) */
      );
  }
  int revertToStart() override {
      PYBIND11_OVERRIDE_PURE(int, UniaxialMaterial, revertToStart);
  }
  UniaxialMaterial *getCopy() override {
    py::gil_scoped_acquire acquire;

    auto self = py::cast(this);
    auto cloned = self.attr("getCopy")();

    auto keep_python_state_alive = std::make_shared<py::object>(cloned);
    // auto ptr = cloned.cast<UniaxialMaterial*>();
    auto ptr = py::cast<PyUniaxialMaterial*>(cloned.release());

    py::gil_scoped_release release;
    return ptr;
    // return std::shared_ptr<UniaxialMaterial>(keep_python_state_alive, ptr);
    // PYBIND11_OVERRIDE_PURE(UniaxialMaterial*, UniaxialMaterial, getCopy);
    // return this->clone();
  }
    
  int setTrialStrain(double strain, double strainRate=0) override {
      PYBIND11_OVERRIDE_PURE(
          double,                 /* Return type */
          UniaxialMaterial,       /* Parent class */
          setTrialStrain,         /* Name of function in C++ (must match Python name) */
          strain,                 /* Argument(s) */
          strainRate
      );
  }

// private:
//   const py::object m_object;
};



static py::array_t<double>
copy_vector(Vector vector)
{
  py::array_t<double> array(vector.Size());
  double *ptr = static_cast<double*>(array.request().ptr);
  for (int i=0; i<vector.Size(); i++)
    ptr[i] = vector(i);
  return array;
}

static Vector *
new_vector(py::array_t<double> array)
{
  py::buffer_info info = array.request();
  return new Vector(static_cast<double*>(info.ptr),(int)info.shape[0]);
}

py::array_t<double>
copy_matrix(Matrix matrix)
{
  int nr = matrix.noRows();
  int nc = matrix.noCols();
  py::array_t<double> array({nr,nc},{nr*nc*sizeof(double), nc*sizeof(double)});
  double *ptr = static_cast<double*>(array.request().ptr);

  for (int i=0; i<matrix.noRows(); i++)
    for (int j=0; j<matrix.noCols(); j++)
      ptr[i*nc+j] = matrix(i,j);
  return array;
}


void
init_obj_module(py::module &m)
{
  py::class_<Vector, std::unique_ptr<Vector, py::nodelete>> PyVector(m, "Vector", py::buffer_protocol());

  PyVector.def (py::init([](
        py::array_t<double, py::array::c_style|py::array::forcecast> array
      ) -> Vector {
        const bool verbose = false;
        py::buffer_info info = array.request();
        if (verbose){
          py::print("ptr\t",      info.ptr);
          py::print("itemsize\t", info.itemsize);
          py::print("format\t",   info.format);
          py::print("ndim\t",     info.ndim);
          py::print("shape\t",    py::cast(info.shape));
          py::print("strides\t",  py::cast(info.strides));
          py::print("array\t",    array);
          printf("%lf\n",         *((double*)info.ptr));
        }
        return Vector(static_cast<double*>(info.ptr),(int)info.shape[0]);
      }))
      /* Allow reference by numpy array; requires access to Vector.theData */
      // .def_buffer([](Vector& v) -> py::buffer_info{
      //       return py::buffer_info(
      //       );
      // })

      .def (py::init([](
        py::array_t<double, py::array::c_style|py::array::forcecast> array,
        int assert_size
      ) -> Vector {
        py::buffer_info info = array.request();
        if (info.shape[0] != assert_size)
          throw std::runtime_error("Incompatible buffer dimension.");
        return Vector(static_cast<double*>(info.ptr), static_cast<int>(info.shape[0]));
      }))
  ;

  py::class_<Matrix, std::unique_ptr<Matrix, py::nodelete>>(m, "Matrix", py::buffer_protocol())
    .def (py::init([](
        py::array_t<double, ARRAY_FLAGS> array,
        int assert_size
    ) -> Matrix {
      bool verbose = true;
      py::buffer_info info = array.request();
      if (verbose){
        py::print("ptr\t",info.ptr);
        py::print("itemsize\t", info.itemsize);
        py::print("format\t", info.format);
        py::print("ndim\t", info.ndim);
        py::print("shape\t", py::cast(info.shape));
        py::print("strides\t", py::cast(info.strides));
      }
      if (info.shape[0] != assert_size)
        throw std::runtime_error("Incompatible buffer dimension.");
      return Matrix(
        static_cast<double*>(info.ptr), 
        static_cast<int>(info.shape[0]),
        static_cast<int>(info.shape[1])
      );
    }))
  ; 
  py::class_<Node,    std::unique_ptr<Node,py::nodelete>>(m, "_Node")
  ;
  py::class_<Element, std::unique_ptr<Element,py::nodelete>>(m, "_Element")
    .def ("commitState",           &Element::commitState)
    .def ("revertToStart",         &Element::revertToStart)
    .def ("revertToLastCommit",    &Element::revertToLastCommit)
  ;

  py::class_<SectionForceDeformation, std::unique_ptr<SectionForceDeformation, py::nodelete> > (m, "_SectionForceDeformation")
    .def ("getSectionTangent",          [](SectionForceDeformation& section) {
      return copy_matrix(section.getSectionTangent());
    })
    .def ("getInitialTangent",          [](SectionForceDeformation& section) {
      return copy_matrix(section.getInitialTangent());
    })
    .def ("getSectionFlexibility",      [](SectionForceDeformation& section) {
      return copy_matrix(section.getSectionFlexibility());
    })
    .def ("getInitialFlexibility",      [](SectionForceDeformation& section) {
      return copy_matrix(section.getInitialFlexibility());
    })
    .def ("setTrialSectionDeformation", [](SectionForceDeformation& section,  
        py::array_t<double, py::array::c_style|py::array::forcecast> deformation) {
      return section.setTrialSectionDeformation(*new_vector(deformation));
    }) 
    .def ("setTrialSectionDeformation", [](SectionForceDeformation& section, Vector &deformation) {
        return section.setTrialSectionDeformation(deformation);
    }) 
    .def ("getSectionDeformation", &SectionForceDeformation::getSectionDeformation)

    .def ("getStressResultant",    [](SectionForceDeformation &section, py::array_t<double> deformation, bool commit=false) {
        section.setTrialSectionDeformation(*new_vector(deformation));
        if (commit) section.commitState();
        return copy_vector(section.getStressResultant());
    })
    .def ("getStressResultant",    [](SectionForceDeformation &section) {
        return copy_vector(section.getStressResultant());
    })
    .def ("commitState",           &SectionForceDeformation::commitState)
    .def ("revertToStart",         &SectionForceDeformation::revertToStart)
    .def ("revertToLastCommit",    &SectionForceDeformation::revertToLastCommit)
  ;

  py::class_<UniaxialMaterial, PyUniaxialMaterial, std::shared_ptr<UniaxialMaterial>>(m, "_UniaxialMaterial", py::multiple_inheritance())
    .def (py::init<py::object &, int>())
    .def (py::init<const UniaxialMaterial &>())
    .def (py::init_alias<const UniaxialMaterial &>())
    .def ("setTrialStrain",        [](UniaxialMaterial &material, double strain) {
        return material.setTrialStrain(strain);
      }
    )
    .def ("getStress",             &UniaxialMaterial::getStress)
    .def ("getStress",            [](UniaxialMaterial &material, double strain, bool commit=false){
        material.setTrialStrain(strain);
        if (commit)
          material.commitState();
        return material.getStress();
      }, 
      py::arg("strain"), py::arg("commit")
    )
    .def ("getTangent",            &UniaxialMaterial::getTangent)
    .def ("getDampTangent",        &UniaxialMaterial::getDampTangent)
    .def ("getStrainRate",         &UniaxialMaterial::getStrainRate)
    .def ("commitState",           &UniaxialMaterial::commitState)
    .def ("revertToStart",         &UniaxialMaterial::revertToStart)
    .def ("revertToLastCommit",    &UniaxialMaterial::revertToLastCommit)
  ;

  // py::class_<HystereticBackbone,   PyHystereticBackbone>(m, "HystereticBackbone")
  //   .def("getStress", &HystereticBackbone::getStress);
  // ;

  // py::class_<ManderBackbone, HystereticBackbone>(m, "PopovicsBackbone")
  //   .def(py::init<int, double, double, double>(),
  //        py::arg("tag"), py::arg("f"), py::arg("e"), py::arg("E")
  //   )
  //   .def("getStress", &ManderBackbone::getStress)
  // ;


  py::class_<ModelRegistry, std::unique_ptr<ModelRegistry, py::nodelete> >(m, "TclModelRegistry")
    .def (py::init([](py::object interpaddr)->std::unique_ptr<ModelRegistry, py::nodelete>{
        void *interp_addr;
        interp_addr = (void*)PyLong_AsVoidPtr(interpaddr.ptr());
        void *builder_addr = Tcl_GetAssocData((Tcl_Interp*)interp_addr, "OPS::theModelRegistry", NULL);
        return std::unique_ptr<ModelRegistry, py::nodelete>((ModelRegistry*)builder_addr);
      }) // , py::return_value_policy::reference
    )
    .def ("getSection", [](ModelRegistry& builder, int id){
        return builder.getTypedObject<SectionForceDeformation>(id);
    })

    .def ("getNDMaterial", [](ModelRegistry& builder, int tag){
        return builder.getTypedObject<NDMaterial>(tag);
    })
    /*
    .def ("getUniaxialMaterial", [](ModelRegistry& builder, py::str tag){
        return builder.getTypedObject<UniaxialMaterial>(tag);
    })
    */
    .def ("getUniaxialMaterial", [](ModelRegistry& builder, int tag){
        return builder.getTypedObject<UniaxialMaterial>(tag);
    })
    .def ("addPythonObject", [](ModelRegistry& builder, int tag, PyUniaxialMaterial& material){
        return builder.addTypedObject<UniaxialMaterial>(tag, &material);
    })
    // .def ("getHystereticBackbone", [](ModelRegistry& builder, int tag){
    //     return std::unique_ptr<HystereticBackbone, py::nodelete>(builder.getTypedObject<HystereticBackbone>(tag));
    // })
  ;

  py::class_<Domain>(m, "_Domain")
    // .def ("getElementResponse", &Domain::getElementResponse)
    .def ("getNodeResponse", [](Domain& domain, int node, std::string type) {
        NodeData typ;
        if (type == "displ") typ = NodeData::Disp;
        if (type == "accel") typ = NodeData::Accel;
        if (type == "veloc") typ = NodeData::Vel;
        if (type == "react") typ = NodeData::Reaction;
      return copy_vector(*domain.getNodeResponse(node, typ));
    })
    .def ("getTime", &Domain::getCurrentTime)
  ;
  
  py::class_<G3_Runtime>(m, "_Runtime")
  ;

  //
  // Module-Level Functions
  //
  m.def ("get_builder", &get_builder);
  m.def ("get_domain", [](G3_Runtime *rt)->std::unique_ptr<Domain, py::nodelete>{
      Domain *domain_addr = rt->m_domain;
      return std::unique_ptr<Domain, py::nodelete>((Domain*)domain_addr);
    } // , py::return_value_policy::reference
  );

}

PYBIND11_MODULE(OpenSeesPyRT, m) {
  init_obj_module(m);
}

