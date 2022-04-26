#pragma once

#include <Eigen/Geometry>
#include <Eigen/Core>
#include "../../qpOASES/qpOASES.hpp"
#include <iostream>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_v_aux_ext_dep.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <blasfeo_i_aux_ext_dep.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include <blasfeo_d_aux_ext_dep.h>
#include <hpipm_d_dense_qp_ipm.h>
#include <hpipm_d_dense_qp_dim.h>
#include <hpipm_d_dense_qp.h>
#include <hpipm_d_dense_qp_sol.h>
#include <hpipm_timing.h>

inline bool isEven(int n) {
    return not (n % 2);
}

inline bool isOdd(int n) {
    return (n % 2);
}

inline double wrapToPi(double angle){
  double ret=angle;
  while(ret>M_PI){
    ret-=2*M_PI;
  }

  while(ret<=-M_PI){
    ret+=2*M_PI;
  }

  return ret;
}

inline double angDiff(double thetaD, double theta){
  double alpha=0;
  Eigen::Vector2d nD,n;

  nD<<cos(thetaD), sin(thetaD);
  n<<cos(theta), sin(theta);

  double alphaAbs = acos(nD.transpose()*n);

  Eigen::Vector3d n3,nD3;

  n3<<n(0),n(1),0;
  nD3<<nD(0),nD(1),0;

  Eigen::Vector3d nC3;

  nC3=n3.cross(nD3);

  if(nC3(2)>0){
    alpha=alphaAbs;
  }
  else{
    alpha=-alphaAbs;
  }

  return alpha;

}

inline Eigen::MatrixXd matrixPower(Eigen::MatrixXd &A, int exp){

	Eigen::MatrixXd result = Eigen::MatrixXd::Identity(A.rows(),A.cols());

	for (int i=0; i<exp;++i)
        	result *= A;

	return result;
}

inline double sign(double x){
	if(x>0) return +1;
	if(x<0) return -1;
	return -1;
}

inline Eigen::VectorXd solveQP(Eigen::MatrixXd costFunctionH, Eigen::VectorXd costFunctionF, 
		Eigen::MatrixXd AConstraint, Eigen::VectorXd bConstraintMin, Eigen::VectorXd bConstraintMax) {

	int nVariables = costFunctionH.rows();
	int nConstraints = AConstraint.rows();

	qpOASES::QProblem qp;

	qpOASES::real_t H[nVariables*nVariables];
	qpOASES::real_t g[nVariables];

	qpOASES::real_t A[nConstraints*nVariables];
	qpOASES::real_t lb[nConstraints];
	qpOASES::real_t ub[nConstraints];

	for(int i=0;i<nVariables;++i){
		for(int j=0;j<nVariables;++j){
			H[i*nVariables+j] = costFunctionH(i,j);
		}
		g[i] = costFunctionF(i);
	}

	for(int i=0;i<nConstraints;++i){
		for(int j=0;j<nVariables;++j){
			A[i*nVariables+j] = AConstraint(i,j);
		}
		lb[i] = bConstraintMin(i);
		ub[i] = bConstraintMax(i);
	}

	qpOASES::real_t xOpt[nVariables];

	qpOASES::Options options;
	options.setToMPC();
	options.printLevel=qpOASES::PL_NONE;
	qpOASES::int_t nWSR = 300;

	qp = qpOASES::QProblem(nVariables, nConstraints);
	qp.setOptions(options);
	qp.init(H,g,A,0,0,lb,ub,nWSR,NULL,NULL,NULL,NULL,NULL,NULL);

	qp.getPrimalSolution(xOpt);

	Eigen::VectorXd decisionVariables(nVariables);

	for(int i=0;i<nVariables;++i){
		decisionVariables(i) = xOpt[i];
	}

        //if (qp.isInfeasible()) std::cout << "INFEASIBLE " << std::endl;
	return decisionVariables;
}

inline Eigen::VectorXd solveQP_hpipm(Eigen::MatrixXd costFunctionH, Eigen::VectorXd costFunctionF, 
		Eigen::MatrixXd AConstraint, Eigen::VectorXd bConstraintMin, Eigen::VectorXd bConstraintMax, int n_eq) {

    int n_variables = costFunctionH.cols();
    int n_constr = AConstraint.rows() - n_eq;

    int nv = n_variables;
    int ne = n_eq;  //number of equality constraints
    //if (itr == 20 ) ne = 0;

    int nb = 0; 
    int ng = n_constr;
    int ns = 0; 
    int nsb = 0; 
    int nsg = 0;
    int idxb[n_variables] = {};
    double H[n_variables*n_variables] = {};
    double g[n_variables] = {};
    double d_lb[n_constr] = {};
    double d_ub[n_constr] = {};
    double C[n_constr*n_variables] = {};

    double A[ne*n_variables] = {};
    double b[ne] = {};

    for (int i = 0; i<n_variables; i++) {
        //idxb[i] = i;
        g[i] = costFunctionF(i);

        for (int j = 0; j<n_variables; j++) {
            H[j*n_variables+i] = costFunctionH(i,j);
        }
    }

    for (int k = 0; k<n_constr; k++) {
        d_lb[k] = bConstraintMin(k+ne); //bConstraintMin(k);
        d_ub[k] = bConstraintMax(k+ne); //bConstraintMax(k);

        for (int j = 0; j<n_variables; j++) {
            C[j*n_constr+k] = AConstraint(k+ne,j); //AConstraint(k,j);
        }
    }

    for (int k = 0; k<ne; k++) {
        b[k] = bConstraintMax(k); //beq(k);

        for (int j = 0; j<n_variables; j++) {
            A[j*ne+k] = AConstraint(k,j); //Aeq(k,j);
        }
    }
    // allocate memory for QP struct

    int dim_size = d_dense_qp_dim_memsize();
    void *dim_mem = malloc(dim_size);
    struct d_dense_qp_dim dim;
    d_dense_qp_dim_create(&dim, dim_mem);

    d_dense_qp_dim_set_all(nv, ne, nb, ng, nsb, nsg, &dim);

    int qp_size = d_dense_qp_memsize(&dim);
    void *qp_mem = malloc(qp_size);
    struct d_dense_qp qp;
    d_dense_qp_create(&dim, &qp, qp_mem);

    d_dense_qp_set_H(H, &qp);
    d_dense_qp_set_g(g, &qp);

    //if (itr == 20 ) { // do nothing
    //} else {
    if (ne == 0) {
        d_dense_qp_set_A(A, &qp);
        d_dense_qp_set_b(b, &qp);
    }
    //}

    d_dense_qp_set_C(C, &qp);
    d_dense_qp_set_lg(d_lb, &qp);
    d_dense_qp_set_ug(d_ub, &qp);

    // allocate memory for the solution

    int qp_sol_size = d_dense_qp_sol_memsize(&dim);
    void *qp_sol_mem = malloc(qp_sol_size);
    struct d_dense_qp_sol qp_sol;
    d_dense_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

    // allocate memory for ipm solver and its workspace

    int ipm_arg_size = d_dense_qp_ipm_arg_memsize(&dim);
    void *ipm_arg_mem = malloc(ipm_arg_size);
    struct d_dense_qp_ipm_arg arg;
    d_dense_qp_ipm_arg_create(&dim, &arg, ipm_arg_mem);
    enum hpipm_mode mode = ROBUST;//ROBUST;   // set mode ROBUST, SPEED, BALANCE, SPEED_ABS
    d_dense_qp_ipm_arg_set_default(mode, &arg);

    int ipm_size = d_dense_qp_ipm_ws_memsize(&dim, &arg);
    void *ipm_mem = malloc(ipm_size);
    struct d_dense_qp_ipm_ws workspace;
    d_dense_qp_ipm_ws_create(&dim, &arg, &workspace, ipm_mem);

    // solve QP

    d_dense_qp_ipm_solve(&qp, &qp_sol, &arg, &workspace);

    // store solution in an Eigen Vector

    int nu_max = n_variables;   // see here what is missing
    Eigen::VectorXd u_store = Eigen::VectorXd::Zero(nu_max);
    double *u = (double *) malloc(nu_max*sizeof(double));
    d_dense_qp_sol_get_v(&qp_sol, u);
    for (int i=0; i<n_variables; i++) u_store(i) = u[i];

    // free memory

    free(u);
    free(dim_mem);
    free(qp_mem);
    free(qp_sol_mem);
    free(ipm_arg_mem);
    free(ipm_mem);

    return u_store;
/**/


/*
    int n_variables = costFunctionH.cols();
    int n_constr = AConstraint.rows();

    int nv = n_variables;
    int ne = 0;  //number of equality constraints
    //if (itr == 20 ) ne = 0;

    int nb = 0; 
    int ng = n_constr;
    int ns = 0; 
    int nsb = 0; 
    int nsg = 0;
    int idxb[n_variables] = {};
    double H[n_variables*n_variables] = {};
    double g[n_variables] = {};
    double d_lb[n_constr] = {};
    double d_ub[n_constr] = {};
    double C[n_constr*n_variables] = {};

    double A[ne*n_variables] = {};
    double b[ne] = {};

    for (int i = 0; i<n_variables; i++) {
        //idxb[i] = i;
        g[i] = costFunctionF(i);

        for (int j = 0; j<n_variables; j++) {
            H[j*n_variables+i] = costFunctionH(i,j);
        }
    }

    for (int k = 0; k<n_constr; k++) {
        d_lb[k] = bConstraintMin(k+ne); //bConstraintMin(k);
        d_ub[k] = bConstraintMax(k+ne); //bConstraintMax(k);

        for (int j = 0; j<n_variables; j++) {
            C[j*n_constr+k] = AConstraint(k+ne,j); //AConstraint(k,j);
        }
    }


    // allocate memory for QP struct

    int dim_size = d_dense_qp_dim_memsize();
    void *dim_mem = malloc(dim_size);
    struct d_dense_qp_dim dim;
    d_dense_qp_dim_create(&dim, dim_mem);

    d_dense_qp_dim_set_all(nv, ne, nb, ng, nsb, nsg, &dim);

    int qp_size = d_dense_qp_memsize(&dim);
    void *qp_mem = malloc(qp_size);
    struct d_dense_qp qp;
    d_dense_qp_create(&dim, &qp, qp_mem);

    d_dense_qp_set_H(H, &qp);
    d_dense_qp_set_g(g, &qp);

    //if (itr == 20 ) { // do nothing
    //} else {
    //    d_dense_qp_set_A(A, &qp);
    //    d_dense_qp_set_b(b, &qp);
    //}

    d_dense_qp_set_C(C, &qp);
    d_dense_qp_set_lg(d_lb, &qp);
    d_dense_qp_set_ug(d_ub, &qp);

    // allocate memory for the solution

    int qp_sol_size = d_dense_qp_sol_memsize(&dim);
    void *qp_sol_mem = malloc(qp_sol_size);
    struct d_dense_qp_sol qp_sol;
    d_dense_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

    // allocate memory for ipm solver and its workspace

    int ipm_arg_size = d_dense_qp_ipm_arg_memsize(&dim);
    void *ipm_arg_mem = malloc(ipm_arg_size);
    struct d_dense_qp_ipm_arg arg;
    d_dense_qp_ipm_arg_create(&dim, &arg, ipm_arg_mem);
    enum hpipm_mode mode = ROBUST;//ROBUST;   // set mode ROBUST, SPEED, BALANCE, SPEED_ABS
    d_dense_qp_ipm_arg_set_default(mode, &arg);

    int ipm_size = d_dense_qp_ipm_ws_memsize(&dim, &arg);
    void *ipm_mem = malloc(ipm_size);
    struct d_dense_qp_ipm_ws workspace;
    d_dense_qp_ipm_ws_create(&dim, &arg, &workspace, ipm_mem);

    // solve QP

    d_dense_qp_ipm_solve(&qp, &qp_sol, &arg, &workspace);

    // store solution in an Eigen Vector

    int nu_max = n_variables;   // see here what is missing
    Eigen::VectorXd u_store = Eigen::VectorXd::Zero(nu_max);
    double *u = (double *) malloc(nu_max*sizeof(double));
    d_dense_qp_sol_get_v(&qp_sol, u);
    for (int i=0; i<n_variables; i++) u_store(i) = u[i];

    // free memory

    free(u);
    free(dim_mem);
    free(qp_mem);
    free(qp_sol_mem);
    free(ipm_arg_mem);
    free(ipm_mem);

    return u_store;
/**/
}


inline Eigen::VectorXd solveQP_hpipm_z(Eigen::MatrixXd costFunctionH, Eigen::VectorXd costFunctionF, 
		Eigen::MatrixXd AConstraint, Eigen::VectorXd bConstraintMin, Eigen::VectorXd bConstraintMax, Eigen::MatrixXd AeqZ, Eigen::VectorXd beqZ) {



    int n_variables = costFunctionH.cols();
    int n_constr = AConstraint.rows() ;

    int nv = n_variables;
    int ne = AeqZ.rows();

    int nb = 0; 
    int ng = n_constr;
    int ns = 0; 
    int nsb = 0; 
    int nsg = 0;
    int idxb[n_variables] = {};
    double H[n_variables*n_variables] = {};
    double g[n_variables] = {};
    double d_lb[n_constr] = {};
    double d_ub[n_constr] = {};
    double _lb[n_constr] = {};
    double _ub[n_constr] = {};
    double C[n_constr*n_variables] = {};
    
    double A[ne*n_variables] = {};
    double b[ne] = {};

    for (int i = 0; i<n_variables; i++) {
        //idxb[i] = i;
        g[i] = costFunctionF(i);

        for (int j = 0; j<n_variables; j++) {
            H[j*n_variables+i] = costFunctionH(i,j);
        }
    }


    for (int k = 0; k<n_constr; k++) {
        d_lb[k] = bConstraintMin(k); //bConstraintMin(k);
        d_ub[k] = bConstraintMax(k); //bConstraintMax(k);

        for (int j = 0; j<n_variables; j++) {
            C[j*n_constr+k] = AConstraint(k,j); //AConstraint(k,j);
        }
    }

    for (int k = 0; k<ne; k++) {
        b[k] = beqZ(k); 
        for (int j = 0; j<n_variables; j++) {
            A[j*ne+k] = AeqZ(k,j); //AConstraint(k,j);
        }
    }

    int dim_size = d_dense_qp_dim_memsize();
    void *dim_mem = malloc(dim_size);
    struct d_dense_qp_dim dim;
    d_dense_qp_dim_create(&dim, dim_mem);

    d_dense_qp_dim_set_all(nv, ne, nb, ng, nsb, nsg, &dim);

    int qp_size = d_dense_qp_memsize(&dim);
    void *qp_mem = malloc(qp_size);
    struct d_dense_qp qp;
    d_dense_qp_create(&dim, &qp, qp_mem);

    d_dense_qp_set_H(H, &qp);
    d_dense_qp_set_g(g, &qp);

    d_dense_qp_set_C(C, &qp);
    d_dense_qp_set_lg(d_lb, &qp);
    d_dense_qp_set_ug(d_ub, &qp);

    d_dense_qp_set_A(A, &qp);
    d_dense_qp_set_b(b, &qp);

    // allocate memory for the solution

    int qp_sol_size = d_dense_qp_sol_memsize(&dim);
    void *qp_sol_mem = malloc(qp_sol_size);
    struct d_dense_qp_sol qp_sol;
    d_dense_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

    // allocate memory for ipm solver and its workspace

    int ipm_arg_size = d_dense_qp_ipm_arg_memsize(&dim);
    void *ipm_arg_mem = malloc(ipm_arg_size);
    struct d_dense_qp_ipm_arg arg;
    d_dense_qp_ipm_arg_create(&dim, &arg, ipm_arg_mem);
    enum hpipm_mode mode = SPEED;//ROBUST;   // set mode ROBUST, SPEED, BALANCE, SPEED_ABS
    d_dense_qp_ipm_arg_set_default(mode, &arg);

    int ipm_size = d_dense_qp_ipm_ws_memsize(&dim, &arg);
    void *ipm_mem = malloc(ipm_size);
    struct d_dense_qp_ipm_ws workspace;
    d_dense_qp_ipm_ws_create(&dim, &arg, &workspace, ipm_mem);

    // solve QP

    d_dense_qp_ipm_solve(&qp, &qp_sol, &arg, &workspace);

    // store solution in an Eigen Vector

    int nu_max = n_variables;   // see here what is missing
    Eigen::VectorXd u_store = Eigen::VectorXd::Zero(nu_max);
    double *u = (double *) malloc(nu_max*sizeof(double));
    d_dense_qp_sol_get_v(&qp_sol, u);
    for (int i=0; i<n_variables; i++) u_store(i) = u[i];

    // free memory

    free(u);
    free(dim_mem);
    free(qp_mem);
    free(qp_sol_mem);
    free(ipm_arg_mem);
    free(ipm_mem);

    return u_store;
}

inline Eigen::VectorXd solveQP_hpipm_xy_piecewiseconstantZMP(Eigen::MatrixXd costFunctionH, Eigen::VectorXd costFunctionF, Eigen::MatrixXd Aeq, Eigen::VectorXd beq, Eigen::MatrixXd AConstraint, Eigen::VectorXd bZMPMin, Eigen::VectorXd bZMPMax) {

    int n_variables = costFunctionH.cols();
    int n_constr = bZMPMin.rows();


    int nv = n_variables;
    int ne = 1 ;  //number of equality constraints
    //if (itr == 20 ) ne = 0;

    int nb = 0*n_constr; 
    int ng = n_constr;
    int ns = 0; 
    int nsb = 0; 
    int nsg = 0;
    int idxb[n_variables] = {};
    double H[n_variables*n_variables] = {};
    double g[n_variables] = {};
    double d_lb[n_constr] = {};
    double d_ub[n_constr] = {};
    double A[ne*n_variables] = {};
    double b[ne] = {};
    double C[n_constr*n_variables] = {};

    for (int i = 0; i<n_variables; i++) {
        //idxb[i] = i;
        g[i] = costFunctionF(i);

        for (int j = 0; j<n_variables; j++) {
            H[j*n_variables+i] = costFunctionH(i,j);
        }
    }

    for (int k = 0; k<n_constr; k++) {
        d_lb[k] = bZMPMin(k); //bConstraintMin(k);
        d_ub[k] = bZMPMax(k); //bConstraintMax(k);
        for (int j = 0; j<n_variables; j++) {
            C[j*n_constr+k] = AConstraint(k,j); //AConstraint(k,j);
        }
    }

/*
    b[0] = beq(0); //beq(k);
    for (int k = 0; k<n_variables; k++) {
        A[k] = Aeq(0,k);
    }
/**/
  
    for (int k = 0; k<ne; k++) {
        b[k] = beq(k); 
        for (int j = 0; j<n_variables; j++) {
            A[j*ne+k] = Aeq(k,j); //AConstraint(k,j);
        }
    }


    // allocate memory for QP struct

    int dim_size = d_dense_qp_dim_memsize();
    void *dim_mem = malloc(dim_size);
    struct d_dense_qp_dim dim;
    d_dense_qp_dim_create(&dim, dim_mem);

    d_dense_qp_dim_set_all(nv, ne, nb, ng, nsb, nsg, &dim);

    int qp_size = d_dense_qp_memsize(&dim);
    void *qp_mem = malloc(qp_size);
    struct d_dense_qp qp;
    d_dense_qp_create(&dim, &qp, qp_mem);

    d_dense_qp_set_H(H, &qp);
    d_dense_qp_set_g(g, &qp);

    //if (itr == 20 ) { // do nothing
    //} else {
        d_dense_qp_set_A(A, &qp);
        d_dense_qp_set_b(b, &qp);
    //}

    d_dense_qp_set_C(C, &qp);
    d_dense_qp_set_lg(d_lb, &qp);
    d_dense_qp_set_ug(d_ub, &qp);

    // allocate memory for the solution

    int qp_sol_size = d_dense_qp_sol_memsize(&dim);
    void *qp_sol_mem = malloc(qp_sol_size);
    struct d_dense_qp_sol qp_sol;
    d_dense_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

    // allocate memory for ipm solver and its workspace

    int ipm_arg_size = d_dense_qp_ipm_arg_memsize(&dim);
    void *ipm_arg_mem = malloc(ipm_arg_size);
    struct d_dense_qp_ipm_arg arg;
    d_dense_qp_ipm_arg_create(&dim, &arg, ipm_arg_mem);
    enum hpipm_mode mode = BALANCE;//ROBUST;   // set mode ROBUST, SPEED, BALANCE, SPEED_ABS
    d_dense_qp_ipm_arg_set_default(mode, &arg);

    int ipm_size = d_dense_qp_ipm_ws_memsize(&dim, &arg);
    void *ipm_mem = malloc(ipm_size);
    struct d_dense_qp_ipm_ws workspace;
    d_dense_qp_ipm_ws_create(&dim, &arg, &workspace, ipm_mem);

    // solve QP

    d_dense_qp_ipm_solve(&qp, &qp_sol, &arg, &workspace);

    // store solution in an Eigen Vector

    int nu_max = n_variables;   // see here what is missing
    Eigen::VectorXd u_store = Eigen::VectorXd::Zero(nu_max);
    double *u = (double *) malloc(nu_max*sizeof(double));
    d_dense_qp_sol_get_v(&qp_sol, u);
    for (int i=0; i<n_variables; i++) u_store(i) = u[i];

    // free memory

    free(u);
    free(dim_mem);
    free(qp_mem);
    free(qp_sol_mem);
    free(ipm_arg_mem);
    free(ipm_mem);

    return u_store;
}
inline Eigen::Matrix3d rotx(double ax){
	
	Eigen::Matrix3d rx;
	rx.setZero();

	rx(0,0) = 1;
	rx(0,1) = 0;
	rx(0,2) = 0;
	rx(1,0) = 0;
	rx(1,1) = cos(ax);
	rx(1,2) = -sin(ax);
	rx(2,0) = 0;
	rx(2,1) = sin(ax);
	rx(2,2) = cos(ax);
	
	return rx;
}


inline Eigen::Matrix3d roty(double ay){
	
	Eigen::Matrix3d ry;
	ry.setZero();

	ry(0,0) = cos(ay);
	ry(0,1) = 0;
	ry(0,2) = sin(ay);
	ry(1,0) = 0;
	ry(1,1) = 1;
	ry(1,2) = 0;
	ry(2,0) = -sin(ay);
	ry(2,1) = 0;
	ry(2,2) = cos(ay);
	
	return ry;
}


inline Eigen::Matrix3d rotz(double az){
	
	Eigen::Matrix3d rz;
	rz.setZero();
	
	rz(0,0) = cos(az);
	rz(0,1) = -sin(az);
	rz(0,2) = 0;
	rz(1,0) = sin(az);
	rz(1,1) = cos(az);
	rz(1,2) = 0;
	rz(2,0) = 0;
	rz(2,1) = 0;
	rz(2,2) = 1;
	
	return rz;
}

inline Eigen::Matrix3d rot(Eigen::Vector3d eul){
	
	Eigen::Matrix3d r;

	Eigen::Matrix3d rx = rotx(eul(0));
	Eigen::Matrix3d ry = roty(eul(1));
	Eigen::Matrix3d rz = rotz(eul(2));

	r = rx*ry*rz;

	return r;		
}

inline Eigen::Matrix4d v2t(Eigen::VectorXd v){
	
	Eigen::Matrix4d m = Eigen::Matrix4d::Identity();

	Eigen::Vector3d eul = v.head(3);
	Eigen::Matrix3d r = rot(eul);

	m.block<3,3>(0,0) = r;
	m.block<3,1>(0,3) = v.tail(3);
	
	return m;	
}

inline Eigen::VectorXd t2v(Eigen::Matrix4d m) {
	Eigen::VectorXd v(6);

	float beta = atan2( m(0,2), sqrt(pow(m(0,0), 2) + pow(m(0,1), 2)) );
	float alpha = atan2( -m(1,2)/cos(beta), m(2,2)/cos(beta) );
	float gamma = atan2( -m(0,1)/cos(beta), m(0,0)/cos(beta) );
	
	v(0) = alpha;
	v(1) = beta;
	v(2) = gamma;  
	v(3) = m(0,3);
	v(4) = m(1,3);
	v(5) = m(2,3);  

	return v;
}

// Express v2 in the frame of v1
inline Eigen::VectorXd vvRel(Eigen::VectorXd v2, Eigen::VectorXd v1) {
	return t2v(v2t(v1).inverse()*v2t(v2));
}
