#include <comp.hpp>
#include <python_comp.hpp>

using namespace ngcomp;

class HuLinShiSpace;
class HLSFiniteElement;



class HLSFiniteElement : public FiniteElement
{
  const HuLinShiSpace * space;
  Array<int> edges;
  Matrix<> basistrafo{14,14};
  
public:
  HLSFiniteElement(const HuLinShiSpace * _space, Array<int> _edges)
    : FiniteElement(14,1), space{_space}, edges{_edges}
  {
    CalcBasisTrafo();
  }
  
  ELEMENT_TYPE ElementType() const override { return ET_TET; }

  void CalcShapes (const MappedIntegrationPoint<3,3> & mip,
                   FlatMatrix<double> shapes) const    // 14 x 9 matrix
  {
    Matrix<> shapes1(14,9);
    CalcShapes1 (mip.GetPoint(), shapes1);
    shapes = basistrafo * shapes1;   
  }

  void CalcSymCurlShapes (const MappedIntegrationPoint<3,3> & mip,
                          FlatMatrix<double> symcurl_shapes) const    // 14 x 9 matrix
  {
    symcurl_shapes = 0;
  }


private:
  // calc some basis
  void CalcShapes1 (Vec<3> p, FlatMatrix<double> shapes) const
  {
    double x=p(0), y=p(1), z=p(2);
    shapes.Row(0) = Vector ( { 1, 0, 0,  0, 0, 0,  0, 0, 0 } );
    shapes.Row(8) = Vector ( { 0, -z, y,  0, 0, 0,  0, 0, 0 } );
  }

  void CalcBasisTrafo (); // implemented below
};





class HuLinShiSpace : public FESpace
{
  Array<Vec<3>> edge_t, edge_n1, edge_n2;
  
public:
  HuLinShiSpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    cout << "Created a Hu-Lin-Shi finite element space for HCurlSym" << endl;
  }
  
  string GetClassName () const override { return "HuLinShiFESpace"; }

  static DocInfo GetDocu()
  {
    // just copied from example, WIP
    auto docu = FESpace::GetDocu();
    docu.short_docu = "My own FESpace.";
    docu.long_docu =
      R"raw_string(My own FESpace provides first and second order triangular elements.
)raw_string";      
      
    docu.Arg("secondorder") = "bool = False\n"
      "  Use second order basis functions";
    return docu;
  };
  
  // organzize the FESpace, called after every mesh update
  void Update() override
  {
    // compute edge tangential and normal vectors:
    edge_t.SetSize(ma->GetNEdges());

    SetNDof (2*ma->GetNEdges() + 2*ma->GetNE());
  }
    
  // dof-numbers for element-id ei
  void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override
  {
    dnums.SetSize(0);
    if (ei.VB() != VOL) return;
    
    // first dofs are vertex numbers:
    for (auto e : ma->GetElement(ei).Edges())
      {
        dnums.Append (2*e);
        dnums.Append (2*e+1);
      }
    dnums.Append (2*ma->GetNEdges()+2*ei.Nr());
    dnums.Append (2*ma->GetNEdges()+2*ei.Nr()+1);
  }
    
  // generate FiniteElement for element-id ei
  FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override
  {
    switch (ma->GetElement(ei).GetType())
      {
        case ET_TET:
          return * new (alloc) HLSFiniteElement(this, Array<int>{ma->GetElement(ei).Edges()});
      default:
        throw Exception("element type nto implemented");
      }
  }

  tuple<Vec<3>, Vec<3>, Vec<3>> Get_t_n1_n2 (int edgenr) const
  {
    return { edge_t[edgenr], edge_n1[edgenr], edge_n2[edgenr] };
  }

};








void HLSFiniteElement :: CalcBasisTrafo ()
{
  for (int e = 0; e < 6; e++)
    {
      auto [t,n1,n2] = space->Get_t_n1_n2(edges[e]);
    }
}




   
extern "C" void HLSmodule(py::object & res)
{
  cout << "imported Hu-Lin-Shi FESpace" << endl;
  auto ngs = py::module::import("ngsolve");    

  static py::module::module_def def;    
  py::module m = py::module::create_extension_module("", "", &def);    
    
  ExportFESpace<HuLinShiSpace>(m, "HuLinShiSpace", true);
  res = m;    
}    

