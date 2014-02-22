#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "curvature.h"
using namespace OpenMesh;
using namespace Eigen;

#define esp 0.0001

void FindKTFromMatrix(Matrix3d m, CurvatureInfo *info) {

        EigenSolver<Matrix3d> solver;
        solver.compute(m, true);

        Vector3d e1, e2, e3;
        e1 = solver.eigenvectors().col(1);
        e2 = solver.eigenvectors().col(2);
        e3 = solver.eigenvectors().col(3);

        MatrixType::Scalar v1, v2, v3;
        v1 = solver.eigenvalues()(0);
        v2 = solver.eigenvalues()(1);
        v3 = solver.eigenvalues()(2);

//T1, T2, k1, k2 corresponding to algorithm in paper
        Vector3d T1, T2;
        double k1, k2;

//find the eigen value == 0
//thus removing the normal vector
        if(v1 < esp && v1 > -esp) {
            T1 = e2;
            T2 = e3;
        
            k1 = v2;
            k2 = v3;
        } else if (v2 < esp && v2 > -esp) {
            T1 = e1;
            T2 = e3;
            k1 = v1;
            k2 = v3;
        } else if (v3 < esp && v3 > -esp) {
            T1 = e1;
            T2 = e2;
            k1 = v1;
            k2 = v2;
        }
        

		// In the end you need to fill in this struct
		

		info->curvatures[0] = k1;
		info->curvatures[1] = k2;
		info->directions[0] = Vec3f(T1[0], T1[1], T1[2]);
		info->directions[1] = Vec3f(T2[0], T2[1], T2[2]);    
}


void computeCurvature(Mesh &mesh, OpenMesh::VPropHandleT<CurvatureInfo> &curvature) {

	for (Mesh::VertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
		// WRITE CODE HERE TO COMPUTE THE CURVATURE AT THE CURRENT VERTEX ----------------------------------------------
		Vec3f normal = mesh.normal(it.handle());
		Vector3d N(normal[0],normal[1],normal[2]); // example of converting to Eigen's vector class for easier math

        Vec3f meshVi = mesh.point(it.handle());
        Vector3d vi(meshVi[0], meshVi[1], meshVi[2]);

        double TotalArea = 0;

        Mesh::VertexOHalfedgeIter area_iter;
        area_iter = mesh.voh_iter(it.handle());

        while(area_iter) {
            TotalArea += mesh.calc_sector_area(area_iter.handle());
        }

        Matrix3d m = Matrix3d::Zero();
    
        Mesh::VertexOHalfedgeIter he_iter;
        he_iter = mesh.voh_iter(it.handle());

        while(he_iter) {
            Vec3f meshVj = mesh.point(he_iter.value_type());
            Vector3d vj(meshVj[0], meshVj[1], meshVj[2]); 

            Vector3d Tij = (Matrix3d::Identity(3, 3) - N * N.transpose()) * (vi - vj);
            Tij.normalize();

            Vector3d vjminvi = vj - vi;
            Matrix3d Kij = 2 * N.transpose * (vj - vi)/(vjminvi.norm() * vjminvi.norm());
            double A1 = mesh.calc_sector_area(he_iter.handle());
            double A2 = mesh.calc_sector_area(he_iter.handle().flip());

            double Wij = (A1 + A2) / TotalArea;

            m += Wij * Kij * Tij * Tij.transpose();
        }

		CurvatureInfo info;
        FindKTFromMatrix(m, &info);

		mesh.property(curvature,it) = info;
		// -------------------------------------------------------------------------------------------------------------
	}
}

