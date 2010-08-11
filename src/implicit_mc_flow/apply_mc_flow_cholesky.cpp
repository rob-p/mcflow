// OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>

/** Standard Includes */
#include <map>
#include <set>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <numeric>

/** Boost Includes */
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

/** Matrix Library **/
#include <gmm/gmm.h>

extern "C" {
#include "taucs.h"
}
// taucs pisses me off
#undef max
#undef min

/** Local Includes */
#include <arg_parser.hpp>

/** Namespace uses */
using std::string;
using std::ifstream;
using std::map;
using std::ofstream;
using std::ios;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::exception;
using std::abort;
using boost::progress_display;
using boost::shared_ptr;
using namespace OpenMesh;

/** Namespace declarations */
namespace po = boost::program_options;

/** Typedefs **/
struct MyTraits : public OpenMesh::DefaultTraits 
{
  VertexAttributes( OpenMesh::Attributes::Normal |
                    OpenMesh::Attributes::Color );

  FaceAttributes( OpenMesh::Attributes::Normal );
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> cmesh_t;
typedef OpenMesh::VPropHandleT<double> vfloat_handle_t;
typedef OpenMesh::EPropHandleT<double> efloat_handle_t;

typedef gmm::row_matrix< gmm::wsvector<double> > compressed_matrix_t;
typedef gmm::col_matrix< std::vector<double> > dense_matrix_t;


template <typename mesh_t>
double vfarea(mesh_t& mesh, typename mesh_t::VertexHandle vh) {
  typedef typename mesh_t::Point p_t;
  typedef typename mesh_t::VertexHandle vh_t;
  typedef typename mesh_t::HalfedgeHandle hh_t;
  typedef typename mesh_t::Scalar scalar_t;

  typename mesh_t::VertexFaceIter vf_it;
  typename mesh_t::FaceVertexIter fv_it;

  double area = 0.0;
  for (vf_it = mesh.vf_iter(vh); vf_it; ++vf_it) {
    fv_it = mesh.fv_iter(vf_it);
			
    const p_t& P = mesh.point(fv_it);  ++fv_it;
    const p_t& Q = mesh.point(fv_it);  ++fv_it;
    const p_t& R = mesh.point(fv_it);
			
    area += ((Q-P)%(R-P)).norm() * 0.5f;
  }
  return fabs(area) > std::numeric_limits<double>::min() ? area : 1e-8;
}

template <typename mesh_t>
double varea(mesh_t& mesh, typename mesh_t::VertexHandle vh) {
  typedef typename mesh_t::Point p_t;
  typedef typename mesh_t::VertexHandle vh_t;
  typedef typename mesh_t::HalfedgeHandle hh_t;
  typedef typename mesh_t::Scalar scalar_t;

  typename mesh_t::VertexFaceIter vf_it;
  typename mesh_t::FaceVertexIter fv_it;

  double area = 0.0;
  for (vf_it = mesh.vf_iter(vh); vf_it; ++vf_it) {
    fv_it = mesh.fv_iter(vf_it);
			
    const p_t& P = mesh.point(fv_it);  ++fv_it;
    const p_t& Q = mesh.point(fv_it);  ++fv_it;
    const p_t& R = mesh.point(fv_it);
			
    area += ((Q-P)%(R-P)).norm() * 0.5f * 0.3333f;
  }
  return fabs(area) > std::numeric_limits<double>::min() ? area : 1e-8;
}

/*
template <typename mesh_t>
m4_mesh_t::Point normalized_mean_curvature_flow(m4_mesh_t& mesh, m4_mesh_t::VertexHandle vh) {
  typedef m4_mesh_t::Point p_t;
  typedef m4_mesh_t::VertexHandle vh_t;
  typedef m4_mesh_t::HalfedgeHandle hh_t;
  typedef m4_mesh_t::Scalar scalar_t;

  p_t p0, p1, p2, d0, d1;
  vh_t v0, v1;
  hh_t h0, h1, h2;
  scalar_t w, area0, area1, b(0.99);

  size_t valence = 0;

  p0 = mesh.point(vh);
  p_t cog(0.0, 0.0, 0.0);
  v0 = vh;
  area0 = varea(mesh, vh);
  double total_weight = 0.0;
  m4_mesh_t::VertexOHalfedgeIter voh_it;
  for (voh_it = mesh.voh_iter(vh); voh_it; ++voh_it) {
    w = 0.0;

    v1 = mesh.to_vertex_handle(voh_it);
    p1 = mesh.point(v1);

    h0 = voh_it.handle();
    h1 = mesh.opposite_halfedge_handle(voh_it.handle());


    if (!mesh.is_boundary(h0)) {
      h2 = mesh.next_halfedge_handle(h0);
      p2 = mesh.point(mesh.to_vertex_handle(h2));
      d0 = (p0 - p2).normalize();
      d1 = (p1 - p2).normalize();
      w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0|d1)))));
    }

    if (!mesh.is_boundary(h1)) {
      h2 = mesh.next_halfedge_handle(h1);
      p2 = mesh.point(mesh.to_vertex_handle(h2));
      d0 = (p0 - p2).normalize();
      d1 = (p1 - p2).normalize();
      w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0|d1)))));
    }
    w = std::max(static_cast<double>(w), 0.0);
    total_weight += w;
    //area1 = varea(mesh, v1);
    cog += (w) * (p0 - p1);
    ++valence;
  }
  std::cerr << "cog = " << cog << "\n";
  //    return (p0 - cog) * (1.0 / 2.0 * area0);
  return cog * (1.0 / total_weight);
}
*/

template <typename mesh_t>
double mesh_volume ( mesh_t& mesh ) {
  typedef typename mesh_t::HalfedgeHandle hh_t;
  typedef typename mesh_t::VertexHandle vh_t;
  typedef typename mesh_t::Point p_t;

  double volume = 0.0;
  double three_inv = 1.0 / 3.0;
  typename mesh_t::FaceIter f_it;
  for (f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
    hh_t hh = mesh.halfedge_handle(f_it);
    hh_t hhn = mesh.next_halfedge_handle(hh);

    vh_t v0(mesh.from_vertex_handle(hh));
    vh_t v1(mesh.to_vertex_handle(hh));
    vh_t v2(mesh.to_vertex_handle(hhn));

    p_t p0(mesh.point(v0));
    p_t p1(mesh.point(v1));
    p_t p2(mesh.point(v2));
      
    p_t g = (p0 + p1 + p2) * three_inv; 
    p_t n = (p1 - p0) % (p2 - p0);
    volume += (g|n);
  }

  volume *= 1.0 / 6.0;
  return volume;
}

typedef struct {
  int     flags;

  char    uplo;     /* 'u' for upper, 'l' for lower, ' ' don't know; prefer lower. */
  int     n;        /* size of matrix */
  int     n_sn;     /* number of supernodes */

  int* parent;      /* supernodal elimination tree */
  int* first_child; 
  int* next_child;

  int* sn_size;     /* size of supernodes (diagonal block) */
  int* sn_up_size;  /* size of subdiagonal update blocks   */
  int** sn_struct;  /* row structure of supernodes         */

  int* sn_blocks_ld;  /* lda of supernode blocks */
  taucs_datatype** sn_blocks; /* supernode blocks        */
    
  int* up_blocks_ld;  /* lda of update blocks    */
  taucs_datatype** up_blocks; /* update blocks           */
} supernodal_factor_matrix;

template <typename mesh_t>
void implicit_fairing_curvature_flow( mesh_t& mesh, const float dt,
                                      void*& F, int*& perm, int*& iperm, std::ofstream& of) {
  typedef typename mesh_t::Point p_t;
  typedef typename mesh_t::VertexHandle vh_t;
  typedef typename mesh_t::HalfedgeHandle hh_t;
  typedef typename mesh_t::Scalar scalar_t;

  p_t p0, p1, p2, d0, d1;
  vh_t v0, v1;
  hh_t h0, h1, h2;
  scalar_t w, area0, area1, b(0.99);

  size_t valence = 0;
  size_t nnz = mesh.n_vertices();

  std::vector<double> B0(mesh.n_vertices());
  std::vector<double> B1(mesh.n_vertices());
  std::vector<double> B2(mesh.n_vertices());

  std::vector<double> pxvec(mesh.n_vertices());
  std::vector<double> pbvec(mesh.n_vertices());

  std::vector<double> X0(mesh.n_vertices());
  std::vector<double> X1(mesh.n_vertices());
  std::vector<double> X2(mesh.n_vertices());

  gmm::row_matrix< gmm::wsvector<double> > K(mesh.n_vertices(), mesh.n_vertices());

  double vol_init = mesh_volume(mesh);
  gmm::clean(K, std::numeric_limits<double>::min() );
  //    gmm::clean(A, 1e-12);
  
    
  typename mesh_t::VertexIter v_it;
  for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
    size_t idx = v_it.handle().idx();
    double area = varea(mesh, v_it);
    //A(idx, idx) = vfarea(mesh, v_it);
    //Ai(idx, idx) = 1.0 / (vfarea(mesh, v_it));
    K(idx, idx) = area;//A(idx,idx);
    B0[idx] = area * static_cast<double>(mesh.point(v_it)[0]);
    B1[idx] = area * static_cast<double>(mesh.point(v_it)[1]);
    B2[idx] = area * static_cast<double>(mesh.point(v_it)[2]);

    X0[idx] = 0.0;
    X1[idx] = 0.0;
    X2[idx] = 0.0;
  }

  typename mesh_t::EdgeIter e_it;
  for (e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
    w = 0;
    h0 = mesh.halfedge_handle(e_it, 0);
    h1 = mesh.halfedge_handle(e_it, 1);
    v0 = mesh.from_vertex_handle(h0);
    v1 = mesh.to_vertex_handle(h0);
    p0 = mesh.point(v0);
    p1 = mesh.point(v1);

    if (!mesh.is_boundary(h0)) {
      h2 = mesh.next_halfedge_handle(h0);
      p2 = mesh.point(mesh.to_vertex_handle(h2));
      d0 = (p0 - p2).normalize();
      d1 = (p1 - p2).normalize();
      w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0|d1)))));
    }

    if (!mesh.is_boundary(h1)) {
      h2 = mesh.next_halfedge_handle(h1);
      p2 = mesh.point(mesh.to_vertex_handle(h2));
      d0 = (p0 - p2).normalize();
      d1 = (p1 - p2).normalize();
      w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0|d1)))));
    }
    w = std::max( static_cast<double>(w), 1e-10);

     

    w *= 0.25 * dt;
    ++nnz;
    if ( v0.idx() > v1.idx() ) {
      K(v0.idx(), v1.idx()) = -w;
    } else {
      K(v1.idx(), v0.idx()) = -w;
    }
    K(v0.idx(), v0.idx()) += w;
    K(v1.idx(), v1.idx()) += w;
  }

  gmm::csc_matrix< double > CK;
  gmm::copy(K, CK);

  taucs_ccs_matrix* tK = NULL;
  taucs_ccs_matrix _tK;// = ccs_create( K.ncols(), K.nrows(), nnz, TAUCS_DOUBLE | TAUCS_HERMITIAN | TAUCS_LOWER );

  _tK.n = K.ncols();
  _tK.flags = TAUCS_DOUBLE | TAUCS_HERMITIAN | TAUCS_LOWER;
  _tK.colptr = reinterpret_cast<int*>(&CK.jc[0]);
  _tK.rowind = reinterpret_cast<int*>(&CK.ir[0]);
  _tK.values.d = reinterpret_cast<taucs_double*>(&CK.pr[0]);

  /*
  if ( F == NULL ) {
    // Compute ordering
    string ordering("metis");
    taucs_ccs_order( &_tK, &perm, &iperm, const_cast<char*>(ordering.c_str()) );
      
    tK = taucs_ccs_permute_symmetrically( &_tK, perm, iperm );

    char* opts[] = {"taucs.factor.LLT=true", "taucs.factor.mf=true", "taucs.factor.symbolic=true", 
                    "taucs.factor.numeric=false", "taucs.maxdepth=1000", "taucs.factor.ordering=identity", NULL};
    int r = taucs_linsolve(tK, &F, 0, NULL, NULL, opts, NULL );
    // Perform symbolic factorization
    // F = taucs_ccs_factor_llt_symbolic( tK );
    if ( F == NULL ) {
      cerr << "failed symbolic factorization\n";
      exit(-1);
    } else {
      cerr << "performed symbolic factor\n";
    }

  } else {
    tK = taucs_ccs_permute_symmetrically( &_tK, perm, iperm );
  }

  char* opts[] = {"taucs.factor.LLT=true", "taucs.factor.mf=true", "taucs.factor.symbolic=false", 
                "taucs.factor.numeric=true", "taucs.maxdepth=1000", "taucs.factor.ordering=identity", NULL};
  // Solve X
  // Permute RHS
  taucs_vec_permute( pbvec.size(), TAUCS_DOUBLE, &B0[0], &pbvec[0], perm );
  int r = taucs_linsolve(tK, &F, 1, &pxvec[0], &pbvec[0], opts, NULL );
  // Inverse Permute result
  taucs_vec_permute( pxvec.size(), TAUCS_DOUBLE, &pxvec[0], &X0[0], iperm );

  // Solve Y
  // Permute RHS
  taucs_vec_permute( pbvec.size(), TAUCS_DOUBLE, &B1[0], &pbvec[0], perm );
  r = taucs_linsolve(tK, &F, 1, &pxvec[0], &pbvec[0], opts, NULL );
  // Inverse Permute result
  taucs_vec_permute( pxvec.size(), TAUCS_DOUBLE, &pxvec[0], &X1[0], iperm );

  // Solve Z
  // Permute RHS
  taucs_vec_permute( pbvec.size(), TAUCS_DOUBLE, &B2[0], &pbvec[0], perm );
  r = taucs_linsolve(tK, &F, 1, &pxvec[0], &pbvec[0], opts,  NULL );
  // Inverse Permute result
  taucs_vec_permute( pxvec.size(), TAUCS_DOUBLE, &pxvec[0], &X2[0], iperm );

  taucs_ccs_free( tK );
  */

  // If we need to compute a permutation and factorization; do it
  if ( F == NULL ) {
    // Compute ordering
    string ordering("metis");
    taucs_ccs_order( &_tK, &perm, &iperm, const_cast<char*>(ordering.c_str()) );
      
    tK = taucs_ccs_permute_symmetrically( &_tK, perm, iperm );

    // Perform symbolic factorization
    F = taucs_ccs_factor_llt_symbolic( tK );
    if ( F == NULL ) {
      cerr << "failed symbolic factorization\n";
      exit(-1);
    }

  } else {
    tK = taucs_ccs_permute_symmetrically( &_tK, perm, iperm );  
  }

  _tK.values.d = reinterpret_cast<taucs_double*>(&CK.pr[0]);

  // Perform numeric factorization
  int num_fact = taucs_ccs_factor_llt_numeric( tK, F );
  taucs_ccs_free( tK );
  if ( num_fact == TAUCS_ERROR ||
       num_fact == TAUCS_ERROR_NOMEM ||  
       num_fact == TAUCS_ERROR_BADARGS ||  
       num_fact == TAUCS_ERROR_MAXDEPTH ||  
       num_fact == TAUCS_ERROR_INDEFINITE ) {  
    cerr << "failed numeric factorization\n";
    exit(-1);
  } else
    if (num_fact != TAUCS_SUCCESS) {
      cerr << "WTF?\n";
      exit(-1);
    }

  // Solve X
  // Permute RHS
  taucs_vec_permute( pbvec.size(), TAUCS_DOUBLE, &B0[0], &pbvec[0], perm );
  int r = taucs_supernodal_solve_llt( F, &pxvec[0], &pbvec[0] );
  // Inverse Permute result
  taucs_vec_permute( pxvec.size(), TAUCS_DOUBLE, &pxvec[0], &X0[0], iperm );

  // Solve Y
  // Permute RHS
  taucs_vec_permute( pbvec.size(), TAUCS_DOUBLE, &B1[0], &pbvec[0], perm );
  r = taucs_supernodal_solve_llt( F, &pxvec[0], &pbvec[0] );
  // Inverse Permute result
  taucs_vec_permute( pxvec.size(), TAUCS_DOUBLE, &pxvec[0], &X1[0], iperm );

  // Solve Z
  // Permute RHS
  taucs_vec_permute( pbvec.size(), TAUCS_DOUBLE, &B2[0], &pbvec[0], perm );
  r = taucs_supernodal_solve_llt( F, &pxvec[0], &pbvec[0] );
  // Inverse Permute result
  taucs_vec_permute( pxvec.size(), TAUCS_DOUBLE, &pxvec[0], &X2[0], iperm );

  // Release the numeric factorization
   taucs_supernodal_factor_free_numeric( F );
  //taucs_supernodal_factor_free( F );

  // Copy over the new coordinates
   for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
     size_t idx = v_it.handle().idx();
     typename mesh_t::Point p(X0[idx], X1[idx], X2[idx]);
     mesh.set_point(v_it, p);     
   }

  double vol_final = mesh_volume(mesh);

  double beta = pow((vol_init / vol_final), (1.0 / 3.0));
  // std::cerr << "init volume = " << vol_init 
  //           << ", final volume = " << vol_final 
  //           << ", scaling factor = " << beta << "\n";
  for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
    typename mesh_t::Point p(mesh.point(v_it));
    mesh.set_point(v_it, p*beta);
  }
}



template <typename mesh_t>
void write_gaussian_curvatures( mesh_t& mesh, std::ostream& of ) {

  typedef typename mesh_t::Point p_t;
  typedef typename mesh_t::VertexHandle vh_t;
  typedef typename mesh_t::HalfedgeHandle hh_t;
  typedef typename mesh_t::Scalar scalar_t;
  typedef typename mesh_t::VertexIter vit_t;

  p_t p0, p1, p2, d0, d1;
  vh_t v0, v1;
  hh_t h0, h1, h2;
  scalar_t w, area0, area1, b(0.99);
  vit_t v_it;

  mesh.update_normals();

  for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
    typename mesh_t::VertexIHalfedgeIter he_it;

    double area = varea( mesh, v_it );
    double angle = 2.0 * M_PI;

    for (he_it = mesh.vih_iter(v_it); he_it; ++he_it) {
      angle -= mesh.calc_sector_angle(he_it);
    }
    of.write( reinterpret_cast<char*>(&angle), sizeof(double) ); 
  }
}


template <typename mesh_t>
void write_mean_curvatures( mesh_t& mesh, std::ostream& of ) {

  typedef typename mesh_t::Point p_t;
  typedef typename mesh_t::VertexHandle vh_t;
  typedef typename mesh_t::HalfedgeHandle hh_t;
  typedef typename mesh_t::Scalar scalar_t;
  typedef typename mesh_t::VertexIter vit_t;

  p_t p0, p1, p2, d0, d1;
  vh_t v0, v1;
  hh_t h0, h1, h2;
  scalar_t w, area0, area1, b(0.99);
  vit_t v_it;

  mesh.update_normals();

  for (v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
    typename mesh_t::VertexEdgeIter e_it;
    double w = 0.0;
    double tw = 0.0;
    p_t mcn(0.0, 0.0, 0.0);
    p_t po;
    double area = varea( mesh, v_it );
    for (e_it = mesh.ve_iter(v_it); e_it; ++e_it) {

      h0 = mesh.halfedge_handle(e_it, 0);
      h1 = mesh.halfedge_handle(e_it, 1);
      v0 = mesh.from_vertex_handle(h0);
      v1 = mesh.to_vertex_handle(h0);
      p0 = mesh.point(v0);
      p1 = mesh.point(v1);

      if (!mesh.is_boundary(h0)) {
        h2 = mesh.next_halfedge_handle(h0);
        p2 = mesh.point(mesh.to_vertex_handle(h2));
        d0 = (p0 - p2).normalize();
        d1 = (p1 - p2).normalize();
        w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0|d1)))));
      }

      if (!mesh.is_boundary(h1)) {
        h2 = mesh.next_halfedge_handle(h1);
        p2 = mesh.point(mesh.to_vertex_handle(h2));
        d0 = (p0 - p2).normalize();
        d1 = (p1 - p2).normalize();
        w += 1.0 / tan(acos(std::max(-b, std::min(b, (d0|d1)))));
      }

      po = ( v_it.handle().idx() == v0.idx() ) ? p1 : p0;
      mcn += (mesh.point(v_it) - po) * w;
    }
    mcn *= (1.0 / (4.0*area));
    tw = mcn.norm() * 0.5;

    mcn.normalize();
    tw *= (mcn|mesh.normal(v_it));

    of.write( reinterpret_cast<char*>(&tw), sizeof(double) ); 
  }
}

int main( int argc, char **argv ){
  po::options_description desc("Allowed Options");

  string curve_type;

  desc.add_options()
    ("help,h", "produce this message")
    ("mesh_input_file,i", po::value< string >(), "mesh input file")
    ("mesh_output_file,o", po::value< string >(), "mesh output file")
    ("signature_output_file,s", po::value< string >(), "signature output file")
    ("curvature_type,c", po::value< string >(&curve_type)->default_value("mean"), "curvature type")
    ("timestep,t", po::value< float >(), "delta t")
    ("smoothing_steps,n", po::value< size_t >(), "# of smoothing steps")
    ;

  try {
    ArgParser ap( desc );
    ap.parse_args( argc, argv );

    string mesh_input_filename = ap["mesh_input_file"].as< string >();
    string mesh_output_filename = ap["mesh_output_file"].as< string >();
    float dt = ap["timestep"].as< float >();
    size_t n = ap["smoothing_steps"].as< size_t >();

    // The type of curvature we will write out
    curve_type = ap["curvature_type"].as< string >();
    enum { mean, gaussian, unknown } curvature_t;

    if (curve_type == "mean" or curve_type == "m") {
      curvature_t = mean;
    } else 
      if (curve_type == "gaussian" or curve_type == "g") {
        curvature_t = gaussian;
      } else {
        curvature_t = unknown;
        throw runtime_error("unknown curvature type requested for signature");
      }
    
    bool write_sig = false;
    string signature_output_filename;

    po::variable_value sfn = ap["signature_output_file"];
    if (!sfn.empty()) {
     signature_output_filename = ap["signature_output_file"].as< string >();
     write_sig = true;
    }

    cerr << "reading model: [" << mesh_input_filename << "] ... " << flush;
 
    cmesh_t mesh; 
    OpenMesh::IO::Options opt;

    opt += OpenMesh::IO::Options::VertexColor;

    if( !OpenMesh::IO::read_mesh(mesh, mesh_input_filename, opt) ) {
      cerr << "Cannot read mesh: " << mesh_input_filename << endl;
      abort();
    }

    /// Computing Bounding Box
    cmesh_t::ConstVertexIter vIt(mesh.vertices_begin()), vEnd(mesh.vertices_end());
    cmesh_t::Point bbMin, bbMax;

    bbMin = bbMax = mesh.point(vIt);

    for (; vIt!=vEnd; ++vIt) { bbMin.minimize(mesh.point(vIt)); bbMax.maximize(mesh.point(vIt)); } 

    float ubbox = ( bbMax - bbMin ).norm() * 0.01f;

    dt *= ubbox;

    cerr << "Applying " << n 
         << " smoothing steps with dt = " << dt << " ("
	 << "1% bbox diagonal = " << ubbox << ")\n";

    ofstream fout;
    if (write_sig) {
      fout.open(signature_output_filename.c_str(), ios::binary);
      size_t nverts = mesh.n_vertices();
      size_t ntimesteps = n;
      fout.write( (const char*) &nverts, sizeof(size_t) );
      fout.write( (const char*) &ntimesteps, sizeof(size_t) );
    }

    void* F = NULL; // Symbolic factorization of matrix; persistent between integration steps
    int* perm = NULL; // Matrix column permutation and it's inverse
    int* invperm = NULL; 

    progress_display show_progress( n );
    for (size_t i = 0; i < n; ++i)  {
      implicit_fairing_curvature_flow( mesh, dt, F, perm, invperm, fout );

      if (write_sig) {
        switch (curvature_t) {
        case gaussian:
          write_gaussian_curvatures( mesh, fout );
          break;
        case mean:
          write_mean_curvatures( mesh, fout );
          break;
        default:
          break;
        }
      }

      ++show_progress;
    }
    
    if (write_sig) {
      fout.close();
    }

    if( !OpenMesh::IO::write_mesh(mesh, mesh_output_filename, opt) ) {
      cerr << "Cannot write mesh: " << mesh_output_filename << endl;
      abort();
    }

    cerr << "writing model: [" << mesh_input_filename << "] ... " << flush;

  } catch(exception& e) {

    cerr << "Caught Exception: [" << e.what() << "]\n";
    abort();

  }

  return 0;
}
