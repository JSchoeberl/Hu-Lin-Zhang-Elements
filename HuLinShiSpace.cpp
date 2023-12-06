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
    Matrix<> sc_shapes1(14,9);
    CalcSymCurlShapes1 (mip.GetPoint(), sc_shapes1);
    symcurl_shapes = basistrafo * sc_shapes1;   
  }


private:
  // calc some basis, 14x9 matrix
  void CalcShapes1 (Vec<3> p, FlatMatrix<double> shapes) const
  {
    double x=p(0), y=p(1), z=p(2);
    // 8 components of T
    shapes.Row(0) = Vector ( { 1, 0, 0,  0, -1, 0,  0, 0, 0 } );
    shapes.Row(1) = Vector ( { 0, 0, 0,  0, -1, 0,  0, 0, 1 } );
    shapes.Row(2) = Vector ( { 0, 1, 0,  0, 0, 0,  0, 0, 0 } );    
    shapes.Row(3) = Vector ( { 0, 0, 1,  0, 0, 0,  0, 0, 0 } );    
    shapes.Row(4) = Vector ( { 0, 0, 0,  1, 0, 0,  0, 0, 0 } );    
    shapes.Row(5) = Vector ( { 0, 0, 0,  0, 0, 1,  0, 0, 0 } );    
    shapes.Row(6) = Vector ( { 0, 0, 0,  0, 0, 0,  1, 0, 0 } );    
    shapes.Row(7) = Vector ( { 0, 0, 0,  0, 0, 0,  0, 1, 0 } );    

    // 6 components of S \cross x
    shapes.Row(8)  = Vector ( {  0, -z,  y,    0,  0,  0,    0,  0,  0 } );
    shapes.Row(9)  = Vector ( {  z,  0, -x,    0, -z,  y,    0,  0,  0 } );
    shapes.Row(10) = Vector ( { -y,  x,  0,    0,  0,  0,    0, -z,  y } );
    shapes.Row(11) = Vector ( {  0,  0,  0,   -y,  x,  0,    z,  0, -x } );
    shapes.Row(12) = Vector ( {  0,  0,  0,    z,  0, -x,    0,  0,  0 } );
    shapes.Row(13) = Vector ( {  0,  0,  0,    0,  0,  0,   -y,  x,  0 } );
    // todo
  }

  void CalcSymCurlShapes1 (Vec<3> /* p */, FlatMatrix<double> symcurl_shapes) const
  {
    // double x=p(0), y=p(1), z=p(2);
    
    symcurl_shapes = 0;
    symcurl_shapes.Row(8)  = Vector ( { 2, 0, 0,  0, 0, 0,  0, 0, 0 } );
    symcurl_shapes.Row(9)  = Vector ( { 0, 2, 0,  2, 0, 0,  0, 0, 0 } );
    symcurl_shapes.Row(10) = Vector ( { 0, 0, 2,  0, 0, 0,  2, 0, 0 } );
    symcurl_shapes.Row(11) = Vector ( { 0, 0, 0,  0, 0, 2,  0, 2, 0 } );                
    symcurl_shapes.Row(12) = Vector ( { 0, 0, 0,  0, 2, 0,  0, 0, 0 } );                
    symcurl_shapes.Row(13) = Vector ( { 0, 0, 0,  0, 0, 0,  0, 0, 2 } );                
  }

  
  void CalcBasisTrafo (); // implemented below
};






/// Identity operator, covariant transformation
class DiffOpIdHLS : public DiffOp<DiffOpIdHLS>
{
public:
  enum { DIM = 1 };
  enum { DIM_SPACE = 3 };
  enum { DIM_ELEMENT = 3 };
  enum { DIM_DMAT = 9 };
  enum { DIFFORDER = 0 };
  
  static auto & Cast (const FiniteElement & fel) 
  { return static_cast<const HLSFiniteElement&> (fel); }
  
  template <typename MIP, typename MAT>
  static void GenerateMatrix (const FiniteElement & fel, 
                              const MIP & mip,
                              MAT && mat, LocalHeap & lh)  // mat is Col-major
  {
    Cast(fel).CalcShapes (mip, Trans(mat));
  }
};




class HuLinShiSpace : public FESpace
{
  Array<Vec<3>> edge_t, edge_n1, edge_n2;
  
public:
  HuLinShiSpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    cout << "Created a Hu-Lin-Shi finite element space for HCurlSym" << endl;

    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHLS>>();
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
    edge_n1.SetSize(ma->GetNEdges());
    edge_n2.SetSize(ma->GetNEdges());

    for (size_t i = 0; i < ma->GetNEdges(); i++)
      {
        auto [v1,v2] = ma->GetEdgePNums(i);
        Vec<3> p1 = ma->GetPoint<3> (v1);
        Vec<3> p2 = ma->GetPoint<3> (v2);
        Vec<3> t = p2-p1;
        t /= L2Norm(t);
        Vec<3> n1;
        if (fabs(t(0)) > 0.5 || fabs(t(1)) > 0.5)
          n1 = Vec<3> (t(1),-t(0), 0);
        else
          n1 = Vec<3> (0, t(2), -t(1));
        Vec<3> n2 = Cross(t, n1);
        
        edge_t[i] = t;
        edge_n1[i] = n1;
        edge_n2[i] = n2;
      }

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
  Matrix<> shapes1(14, 9);
  const shared_ptr<MeshAccess> & ma = space->GetMeshAccess();
  for (int e = 0; e < 6; e++)
    {
      auto [t,n1,n2] = space->Get_t_n1_n2(edges[e]);

      auto [v1,v2] = ma->GetEdgePNums(edges[e]);
      Vec<3> p = 0.5 * (ma->GetPoint<3> (v1) + ma->GetPoint<3> (v2));
      
      Vec<9> tn1, tn2;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          {
            tn1(3*i+j) = t(i)*n1(j);   // maybe transpose ? 
            tn2(3*i+j) = t(i)*n2(j);
          }
      CalcShapes1 (p, shapes1);
      basistrafo.Row(2*e)   = shapes1 * tn1;
      basistrafo.Row(2*e+1) = shapes1 * tn2;
    }
  
  // what are dofs 13 and 14 ?
  for (int i = 0; i < 14; i++)
    {
      basistrafo(12, i) = sin(i);
      basistrafo(13, i) = cos(i);
    }

  // cout << "basistrafo = " << endl << basistrafo << endl;
  CalcInverse (basistrafo);
  // cout << "basistrafo = " << endl << basistrafo << endl;
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

