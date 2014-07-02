
#include "icp.h"

bool register_pc_point_plane(const Mat& cloud_source, const Mat& cloud_target, Transformation& T)
{
  // Check the input
  // n < n_min already checked in the icp main loop
  const size_t n = cloud_source->size ();
  if (cloud_target->size () != n)
  {
    std::cerr << "ERROR in icp.cpp: Input must have the same size!\n";
    return (false);
  }

  // For numerical stability
  // - Low: Linear Least-Squares Optimization for Point-to-Plane ICP Surface Registration (2004), in the discussion: "To improve the numerical stability of the computation, it is important to use a unit of distance that is comparable in magnitude with the rotation angles. The simplest way is to rescale and move the two input surfaces so that they are bounded within a unit sphere or cube centered at the origin."
  // - Gelfand et al.: Geometrically Stable Sampling for the ICP Algorithm (2003), in sec 3.1: "As is common with PCA methods, we will shift the center of mass of the points	to the origin." ... "Therefore, af- ter shifting the center of mass, we will scale the point set so that the average distance of points	from the origin is 1."
  // - Hartley, Zisserman: - Multiple View Geometry (2004), page 109: They normalize to sqrt(2)
  // TODO: Check the resulting C matrix for the conditioning.

  // Subtract the centroid and calculate the scaling factor
  Eigen::Vector4f c_s (0.f, 0.f, 0.f, 1.f);
  Eigen::Vector4f c_t (0.f, 0.f, 0.f, 1.f);
  pcl::compute3DCentroid (*cloud_source, c_s); c_s.w () = 1.f;
  pcl::compute3DCentroid (*cloud_target, c_t); c_t.w () = 1.f;

  // - There is no need to carry the rgb information along
  // - The normals are only needed for the target
  typedef std::vector <Eigen::Vector4f, Eigen::aligned_allocator <Eigen::Vector4f> > Vec4Xf;

  Vec4Xf xyz_s, xyz_t, nor_t;
  xyz_s.reserve (n);
  xyz_t.reserve (n);
  nor_t.reserve (n);

  CloudProcessed::const_iterator it_s = cloud_source->begin ();
  CloudProcessed::const_iterator it_t = cloud_target->begin ();

  float accum = 0.f;
  for (; it_s!=cloud_source->end (); ++it_s, ++it_t)
  {
    // Subtract the centroid
    const Eigen::Vector4f pt_s = it_s->getVector4fMap () - c_s;
    const Eigen::Vector4f pt_t = it_t->getVector4fMap () - c_t;

    xyz_s.push_back (pt_s);
    xyz_t.push_back (pt_t);
    nor_t.push_back (it_t->getNormalVector4fMap ());

    // Calculate the radius (L2 norm) of the bounding sphere through both shapes and accumulate the average
    // TODO: Change to squared norm and adapt the rest accordingly
    accum += pt_s.head <3> ().norm () + pt_t.head <3> ().norm ();
  }

  // Inverse factor (do a multiplication instead of division later)
  const float factor         = 2.f * static_cast <float> (n) / accum;
  const float factor_squared = factor*factor;

  // Covariance matrix C
  Eigen::Matrix <float, 6, 6> C;

  // Right hand side vector b
  Eigen::Matrix <float, 6, 1> b;

  // For Eigen vectorization: use 4x4 submatrixes instead of 3x3 submatrixes
  // -> top left 3x3 matrix will form the final C
  // Same for b
  Eigen::Matrix4f C_tl    = Eigen::Matrix4f::Zero(); // top left corner
  Eigen::Matrix4f C_tr_bl = Eigen::Matrix4f::Zero(); // top right / bottom left
  Eigen::Matrix4f C_br    = Eigen::Matrix4f::Zero(); // bottom right

  Eigen::Vector4f b_t     = Eigen::Vector4f::Zero(); // top
  Eigen::Vector4f b_b     = Eigen::Vector4f::Zero(); // bottom

  Vec4Xf::const_iterator it_xyz_s = xyz_s.begin ();
  Vec4Xf::const_iterator it_xyz_t = xyz_t.begin ();
  Vec4Xf::const_iterator it_nor_t = nor_t.begin ();

  for (; it_xyz_s!=xyz_s.end (); ++it_xyz_s, ++it_xyz_t, ++it_nor_t)
  {
    const Eigen::Vector4f cross = it_xyz_s->cross3 (*it_nor_t);

    C_tl           += cross     * cross.    transpose ();
    C_tr_bl        += cross     * it_nor_t->transpose ();
    C_br           += *it_nor_t * it_nor_t->transpose ();

    const float dot = (*it_xyz_t-*it_xyz_s).dot (*it_nor_t);

    b_t            += cross     * dot;
    b_b            += *it_nor_t * dot;
  }

  // Scale with the factor and copy the 3x3 submatrixes into C and b
  C_tl    *= factor_squared;
  C_tr_bl *= factor;

  C << C_tl.  topLeftCorner <3, 3> ()            , C_tr_bl.topLeftCorner <3, 3> (),
      C_tr_bl.topLeftCorner <3, 3> ().transpose(), C_br.   topLeftCorner <3, 3> ();

  b << b_t.head <3> () * factor_squared,
      b_b. head <3> () * factor;

  // Solve C * x = b with a Cholesky factorization with pivoting
  // x = [alpha; beta; gamma; trans_x; trans_y; trans_z]
  Eigen::Matrix <float, 6, 1> x = C.selfadjointView <Eigen::Lower> ().ldlt ().solve (b);

  // The calculated transformation in the scaled coordinate system
  const float
      sa = std::sin (x (0)),
      ca = std::cos (x (0)),
      sb = std::sin (x (1)),
      cb = std::cos (x (1)),
      sg = std::sin (x (2)),
      cg = std::cos (x (2)),
      tx = x (3),
      ty = x (4),
      tz = x (5);

  Eigen::Matrix4f TT;
  TT << cg*cb, -sg*ca+cg*sb*sa,  sg*sa+cg*sb*ca, tx,
      sg*cb  ,  cg*ca+sg*sb*sa, -cg*sa+sg*sb*ca, ty,
      -sb    ,  cb*sa         ,  cb*ca         , tz,
      0.f    ,  0.f           ,  0.f           , 1.f;

  // Transformation matrixes into the local coordinate systems of model/data
  Eigen::Matrix4f T_s, T_t;

  T_s << factor, 0.f   , 0.f   , -c_s.x () * factor,
      0.f      , factor, 0.f   , -c_s.y () * factor,
      0.f      , 0.f   , factor, -c_s.z () * factor,
      0.f      , 0.f   , 0.f   ,  1.f;

  T_t << factor, 0.f   , 0.f   , -c_t.x () * factor,
      0.f      , factor, 0.f   , -c_t.y () * factor,
      0.f      , 0.f   , factor, -c_t.z () * factor,
      0.f      , 0.f   , 0.f   ,  1.f;

  // Output transformation T
  T = T_t.inverse () * TT * T_s;

  return (true);
}
